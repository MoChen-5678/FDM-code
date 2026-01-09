!*************************************************************************************************!
!                                                                                                 !
Module DiracB                                                                                     !
!                                                                                                 !
!*************************************************************************************************!
Use Define
Use BASE, only: BASIS, NBD
Use RHFlib
Use Eigenvalues, only: Cdiag
Use PotelHF, only: STM, Stiff
Use omp_lib
Implicit None

Contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine Dirac(lpr)                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Calculate the eigenvalues and wave functions
!
!-------------------------------------------------------------------------------------------------!
!
    logical, intent(in) :: lpr

    integer :: it, ib, i0, n, ka, ibn, nn, iw, kan, nps, kap, i1, i2, i12, NA, ND, ise
    double precision :: rel, t1, t2
    double precision, dimension(MSD) :: fun 
    double precision, allocatable, dimension(:,:) :: xa, AA
    double precision, allocatable, dimension(  :) :: eig,eia

!---    Calculate the Stiff matrix
    t1  = omp_get_wtime ( )
    Call Stiff
    t2  = omp_get_wtime ( )

    write( 6, 200) 'Time used in calculating the stiff matrix:         ', t2-t1, 's'
    write(l6, 200) 'Time used in calculating the stiff matrix:         ', t2-t1, 's'
200 format(2x, a, f12.6, a, f12.6, a)
        
!--- Find the solutions with negative and positive energies 
    do it = 1, 2

    !--- Loop over blocks
        ise = 1
    !$OMP PARALLEL Do private(ib,ka,nn, ND,NA,AA,xa,eig,eia, i1,i2,i12, i0,iw)
        do ib = 1, chunk(it)%nb

            ka  = chunk(it)%ka(ib)
            nn  = chunk(it)%id(ib)
            if( nn .gt. BASIS(ka,it)%NF ) Stop 'Not enough eigen states in Dirac!'

        !--- Load the stiff matrix
            ND  = BASIS(ka,it)%ND
            NA  = BASIS(ka,it)%NDF
            if( .not.allocated( AA  ) ) allocate( AA(NA,NA) )
            if( .not.allocated( xa  ) ) allocate( xa(NA,NA) )
            if( .not.allocated( eig ) ) allocate( eig(  NA) )
           ! if( .not.allocated( eia ) ) allocate( eia(  NA) )
            i12 = 0
            do i1 = 1, NA
                do i2 = i1, NA
                    i12 = i12 + 1
                    AA(i1,i2)   = STM(ib,it)%AH(i12)
                !    AA(i2,i1)   = AA(i1,i2)
                end do
            end do
            
        !--- Sovle the eigen system: eig is the set of eigenvalues and xa is the set of eigenvectors
            Call Cdiag(NA, NA, AA, eig, xa, ise)
            !Call sdiag(NA, NA, AA, eig, xa, eia, ise)
            
        !--- Eigenstates from NBD + 1 to NBD + NMX: Positive energy
            do nn = 1, chunk(it)%id(ib)
            !--- Index of single-particle orbit
                i0  = chunk(it)%ia(ib) + nn

            !--- For positive energy states: single particle energy and wave functions
                lev(it)%ee(i0)  = eig(ND+nn)*hbc
                
                wav(i0,it)%xg   = zero;             wav(i0,it)%xf   = zero
                do iw = 1, ND
                    i1  = iw 
                    wav(i0,it)%xg(:)    = wav(i0,it)%xg(:) + BASIS(ka,it)%XG(:,i1)*xa(iw,ND+nn)
                    wav(i0,it)%xf(:)    = wav(i0,it)%xf(:) + BASIS(ka,it)%XF(:,i1)*xa(iw,ND+nn)
                end do
                do iw = ND+1, NA
                    i1  = iw - ND + NBD
                    wav(i0,it)%xg(:)    = wav(i0,it)%xg(:) + BASIS(ka,it)%XG(:,i1)*xa(iw,ND+nn)
                    wav(i0,it)%xf(:)    = wav(i0,it)%xf(:) + BASIS(ka,it)%XF(:,i1)*xa(iw,ND+nn)
                end do
            end do

            if( allocated( AA  ) ) deallocate( AA  )
            if( allocated( xa  ) ) deallocate( xa  )
            if( allocated( eig ) ) deallocate( eig )
           ! if( allocated( eia ) ) deallocate( eia )
        end do
    !$OMP END PARALLEL Do
        
    !--- Calculate rms radii and probability densities of each orbits
        do i0 = 1, lev(it)%nt
            
            fun = ( wav(i0,it)%xg**2 + wav(i0,it)%xf**2 )*well%xr**2
            Call simps(fun, well%npt, well%h, rel)
            wav(i0,it)%r2  = dsqrt(rel)
            
            wav(i0,it)%rho  = (wav(i0,it)%xg**2 + wav(i0,it)%xf**2)/(4.0*pi*well%xr**2)
        end do

    end do
    
    return

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine Dirac                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!*************************************************************************************************!
!                                                                                                 !
End Module DiracB                                                                                 !
!                                                                                                 !
!*************************************************************************************************!
