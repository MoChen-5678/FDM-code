!*************************************************************************************************!
! Module: FDM_Solver
! Logic Source: Detgff.f90 (Shooting Method Structure)
! Core Solver: 7-Point Finite Difference Matrix
! Key Features:
!   1. Mimics the loop structure of Detgff (it -> ib -> i0).
!   2. Calculates Bulk Properties (r2, rho, gg, ff) exactly like Detgff.
!   3. Uses Matrix Diagonalization to find eigenvalues instead of Shooting.
!*************************************************************************************************!

Module FDM_Solver
    Use Define
    Use BASE, only: NBD
    Use RHFlib, only: well, lev, wav, chunk, simps
    Use PotelHF, only: epotl
    Use Eigenvalues         ! Your Matrix Solver (Cdiag)
    Use omp_lib
    Implicit None

Contains

    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    ! Subroutine: Dirac_FDM
    ! Note: Named Dirac_FDM to match your calling interface, but logic mimics Detgff.
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    Subroutine Dirac_FDM(lpr)
        Logical, Intent(in) :: lpr

        Integer :: it, ib, i0, n0, nmax, n, N_Grid, N_Dim, info
        Integer :: i, j
        Double Precision :: h, r_val, norm_fac, Mass_Eff, inv_hbc
        Double Precision :: t1, t2, eig_MeV
        Double Precision :: val_vms_raw, val_vms_physics
        Double Precision :: rel ! For radius calculation

        ! Matrix Arrays
        Double Precision, Allocatable :: Ham(:,:)         
        Double Precision, Allocatable :: D_Plus(:,:), D_Minus(:,:) 
        Double Precision, Allocatable :: EigVal(:)
        Double Precision, Allocatable :: EigVec(:,:) 
        
        ! Temporary arrays for integration
        Double Precision, Dimension(MSD) :: fun, dens_integ

        t1 = omp_get_wtime()

        ! 1. Basic Constants
        N_Grid = well%npt
        h      = well%h
        ! Only Mass needs conversion (MeV -> fm^-1). Potentials assumed fm^-1.
        inv_hbc = 1.0d0 / 197.327d0 

        If (N_Grid > MSD) Stop 'Error[FDM]: N_Grid > MSD'

        ! 2. Memory Allocation
        N_Dim = 2 * N_Grid
        Allocate(Ham(N_Dim, N_Dim))
        Allocate(D_Plus(N_Grid, N_Grid), D_Minus(N_Grid, N_Grid))
        Allocate(EigVal(N_Dim), EigVec(N_Dim, N_Dim))

        If (lpr) Write(6, '(A)') ' [FDM] Solving Dirac Eq (Matrix Method, Detgff Logic)...'

        !======================================================================
        ! Loop over Particle Types (Proton/Neutron) [Matches Detgff source: 141]
        !======================================================================
        Do it = 1, 2
            
            Mass_Eff = pset%amu(it) * inv_hbc 

            !==================================================================
            ! Loop over Blocks (ib) [Matches Detgff source: 142]
            !==================================================================
            Do ib = 1, chunk(it)%nb
                
                n0   = chunk(it)%ia(ib)
                nmax = chunk(it)%id(ib) ! Number of states in this block (1s, 2s...)

                ! Initialize Matrix
                Ham = 0.0d0; D_Plus = 0.0d0; D_Minus = 0.0d0

                ! 1. Build Kinetic Matrices (7-point Stencil)
                Call Build_Kinetic_Matrices(chunk(it)%ka(ib), N_Grid, h, D_Plus, D_Minus)

                ! 2. Fill Hamiltonian (Kinetic + Potential + Mass)
                !    Note: We build the matrix ONCE per block, which contains all radial excitations.
                Call Build_Hamiltonian(it, ib, N_Grid, h, Mass_Eff, Ham, D_Plus, D_Minus)

                ! 3. Diagonalize (Solve All States in Block)
                Call Cdiag(N_Dim, N_Dim, Ham, EigVal, EigVec, info)
                
                If (info /= N_Dim) Then
                    Write(6,*) 'Error: Cdiag failed at it=', it, ' ib=', ib
                    Stop
                End If

                !==================================================================
                ! Loop over States (n) [Matches Detgff source: 143]
                ! Extraction of Eigenvalues and Wavefunctions
                !==================================================================
                Do n = 1, nmax
                    i0 = n0 + n
                    
                    ! In Matrix FDM, states are ordered by energy.
                    ! The n-th state corresponds to index (N_Grid + n) in EigVal.
                    ! (Skipping N_Grid negative energy states).
                    
                    ! 1. Extract Energy (fm^-1 -> MeV)
                    eig_MeV = (EigVal(N_Grid + n) * 197.327d0) - pset%amu(it)
                    lev(it)%ee(i0) = eig_MeV

                    ! 2. Extract Wavefunctions
                    wav(i0,it)%xg(:) = 0.0d0
                    wav(i0,it)%xf(:) = 0.0d0
                    wav(i0,it)%xg(1:N_Grid) = EigVec(1:N_Grid, N_Grid + n)
                    wav(i0,it)%xf(1:N_Grid) = EigVec(N_Grid+1:N_Dim, N_Grid + n)

                    ! 3. Normalize [Matches Detgff source: 121]
                    fun(1:N_Grid) = wav(i0,it)%xg(1:N_Grid)**2 + wav(i0,it)%xf(1:N_Grid)**2
                    Call simps(fun, N_Grid, h, norm_fac)
                    
                    If (norm_fac > 1.0d-20) Then
                        norm_fac = 1.0d0 / dsqrt(norm_fac)
                        wav(i0,it)%xg = wav(i0,it)%xg * norm_fac
                        wav(i0,it)%xf = wav(i0,it)%xf * norm_fac
                    End If

                    ! Log output similar to Detgff
                    If (lpr) Write(6, '(2x,A, I4, I2, A, F12.6)') &
                        lev(it)%tb(i0), i0, it, ' E=', eig_MeV

                End Do ! End State Loop (n)

            End Do ! End Block Loop (ib)

            !==================================================================
            ! Calculate Bulk Properties [Matches Detgff source: 153-156]
            ! This section is copied logic from Detgff to ensure consistency.
            !==================================================================
            Do i0 = 1, lev(it)%nt
                
                ! RMS Radius Calculation 
                fun = ( wav(i0,it)%xg**2 + wav(i0,it)%xf**2 ) * well%xr**2
                Call simps(fun, N_Grid, well%h, rel)
                wav(i0,it)%r2 = dsqrt(rel)
                
                ! Density Calculation 
                wav(i0,it)%rho = 0.0d0
                Where(well%xr > 1.0d-10)
                    wav(i0,it)%rho = (wav(i0,it)%xg**2 + wav(i0,it)%xf**2) / (4.0d0 * pi * well%xr**2)
                End Where

                ! Component Integration (gg/ff) 
                fun = wav(i0,it)%xg**2
                Call simps(fun, N_Grid, well%h, wav(i0,it)%gg)
                
                fun = wav(i0,it)%xf**2
                Call simps(fun, N_Grid, well%h, wav(i0,it)%ff)

            End Do

        End Do ! End Particle Loop (it)

        ! Cleanup
        If (Allocated(Ham))     Deallocate(Ham)
        If (Allocated(D_Plus))  Deallocate(D_Plus)
        If (Allocated(D_Minus)) Deallocate(D_Minus)
        If (Allocated(EigVal))  Deallocate(EigVal)
        If (Allocated(EigVec))  Deallocate(EigVec)

        t2 = omp_get_wtime()
        If (lpr) Write(6, '(A, F10.4, A)') ' [FDM] Time used: ', t2-t1, 's'

    End Subroutine Dirac_FDM


    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    ! Helper: Build Hamiltonian Matrix
    ! Encapsulates the 7-point stencil + Potential Logic
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    Subroutine Build_Hamiltonian(it, ib, N, h, mass, Ham, Dp, Dm)
        Integer, Intent(in) :: it, ib, N
        Double Precision, Intent(in) :: h, mass
        Double Precision, Intent(inout) :: Ham(2*N, 2*N)
        Double Precision, Intent(in) :: Dp(N,N), Dm(N,N)
        
        Integer :: i, j
        Double Precision :: r_val, val_vms_raw, val_vms_phys
        
        ! [Step 1] Kinetic Energy
        Do j = 1, N
            Do i = 1, N
                Ham(i, j + N) = -Dm(i, j)      ! Top-Right
                Ham(i + N, j) = Dp(i, j)       ! Bottom-Left
            End Do
        End Do

        ! [Step 2] Local Potentials & Mass
        Do i = 1, N
            r_val = dble(i) * h
            
            ! Upper Component (G): V + S + M
            Ham(i, i) = Ham(i, i) + dpotl(it)%vps(i) + mass

            ! Lower Component (F): V - S - M
            ! Auto-correct V-S sign
            val_vms_raw = dpotl(it)%vms(i)
            If (val_vms_raw < 0.0d0) Then
                val_vms_phys = -val_vms_raw ! Convert (S-V) to (V-S)
            Else
                val_vms_phys = val_vms_raw  ! Already (V-S)
            End If
            
            Ham(i + N, i + N) = Ham(i + N, i + N) + val_vms_phys - mass

            ! Off-Diagonal (Kappa/r + Tensor)
            Ham(i, i + N) = Ham(i, i + N) + (dble(chunk(it)%ka(ib))/r_val) + dpotl(it)%vtt(i)
            Ham(i + N, i) = Ham(i + N, i) + (dble(chunk(it)%ka(ib))/r_val) + dpotl(it)%vtt(i)
        End Do

        ! [Step 3] Exchange Potentials (RHF)
        If (pset%IE == 2) Then
            Do j = 1, N
                Do i = 1, N
                    Ham(i, j) = Ham(i, j) + epotl(ib,it)%YG(i, j)
                    Ham(i+N, j+N) = Ham(i+N, j+N) + epotl(ib,it)%XF(i, j)
                    Ham(i, j+N) = Ham(i, j+N) + epotl(ib,it)%YF(i, j)
                    Ham(i+N, j) = Ham(i+N, j) + epotl(ib,it)%XG(i, j)
                End Do
            End Do
        End If

        ! [Step 4] Symmetrization
        Do j = 1, 2*N
            Do i = j + 1, 2*N
                Ham(i, j) = 0.5d0 * (Ham(i, j) + Ham(j, i))
                Ham(j, i) = Ham(i, j)
            End Do
        End Do

    End Subroutine Build_Hamiltonian


    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    ! 7-Point Stencil Matrices (Keeping User's Logic)
    !+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
    Subroutine Build_Kinetic_Matrices(kappa, N, h, M_Plus, M_Minus)
        Integer, Intent(in) :: kappa, N
        Double Precision, Intent(in) :: h
        Double Precision, Intent(out) :: M_Plus(N,N), M_Minus(N,N)
        
        If (kappa < 0) Then
            Call Build_D_Plus_Kernel(M_Plus, N, h)   
            Call Build_D_Minus_Kernel(M_Minus, N, h) 
        Else
            Call Build_D_Minus_Kernel(M_Plus, N, h)  
            Call Build_D_Plus_Kernel(M_Minus, N, h)  
        End If
    End Subroutine Build_Kinetic_Matrices

    Subroutine Build_D_Plus_Kernel(M, N, h)
        Integer, Intent(in) :: N
        Double Precision, Intent(in) :: h
        Double Precision, Intent(out) :: M(N,N)
        Double Precision :: c(0:6), inv_60h
        Integer :: i, k
        M = 0.0d0
        inv_60h = 1.0d0 / (60.0d0 * h)
        c(0) = -147.0d0; c(1) = 360.0d0;  c(2) = -450.0d0
        c(3) = 400.0d0;  c(4) = -225.0d0; c(5) = 72.0d0; c(6) = -10.0d0
        Do i = 1, N
            If (i <= N - 6) Then
                Do k = 0, 6
                    M(i, i+k) = c(k) * inv_60h
                End Do
            Else
                Do k = 0, 6
                    M(i, i-k) = -c(k) * inv_60h 
                End Do
            End If
        End Do
    End Subroutine Build_D_Plus_Kernel

    Subroutine Build_D_Minus_Kernel(M, N, h)
        Integer, Intent(in) :: N
        Double Precision, Intent(in) :: h
        Double Precision, Intent(out) :: M(N,N)
        Double Precision :: c(0:6), inv_60h
        Integer :: i, k
        M = 0.0d0
        inv_60h = 1.0d0 / (60.0d0 * h)
        c(0) = -147.0d0; c(1) = 360.0d0;  c(2) = -450.0d0
        c(3) = 400.0d0;  c(4) = -225.0d0; c(5) = 72.0d0; c(6) = -10.0d0
        Do i = 1, N
            If (i >= 7) Then
                Do k = 0, 6
                    M(i, i-k) = -c(k) * inv_60h
                End Do
            Else
                Do k = 0, 6
                    M(i, i+k) = c(k) * inv_60h
                End Do
            End If
        End Do
    End Subroutine Build_D_Minus_Kernel

End Module FDM_Solver