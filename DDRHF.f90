!*************************************************************************************************!
!                                                                                                 !
Program DDRHF                                                                                     !
!                                                                                                 !
!*************************************************************************************************!
!
Use INOUTS
Use DiracR
Use Expectation
Use configuration
Use Density
Use Define
Use Base
Use Combination
Use PotelHF
Use Gogny
Use DiracB
Use FDM_Solver          ! [USE CORRECT MODULE]
Use omp_lib
implicit none
    
    integer ::  ip, ik, ipa, idn, isolve
    double precision :: xmx0, si0
    double precision :: time1, time2, xtb, xte, tpot, tdir

    type (Lagrangian), dimension(0:8) :: PARA

!
    Call cpu_time(xtb)
    xtb     = omp_get_wtime()
!--- Effective interaction and pairing force
    PARA(0) = PKA1
    PARA(1) = PKO1;     PARA(2) = PKO2;         PARA(3) = PKO3
    PARA(4) = DDME1;    PARA(5) = DDME2
    PARA(6) = PKDD;     PARA(7) = TW99;         PARA(8) = DDLZ1

    write(*,'(a)') 'Which Model: 0 --> PKA1, 1-3 --> PKOx, 4-5 --> DD-MEx, 6 --> PKDD, 7 --> TW99, 8 --> DD-LZ1'
    read(*,*) ipa 
    pset    = PARA(ipa)

    write(*,'(a)') 'Which Solver: 0 --> Basis Expansion (DiracB), 1 --> Finite Difference (FDM)'
    read(*,*) isolve
    if (isolve /= 1) isolve = 0 

!--- Coupling constants for rho-tensor couplings
    pset%grtn   = pset%grho*pset%grtn/(pset%amu(1)/hbc + pset%amu(2)/hbc)

!--- Input from rhf.dat
    Call Reader 

!--- Determine the cannoical basis of the calculated nuclei
    Call Config(.false.)

!--- Calculate the multiexpansions for the propagators and 3j symbols
    time1   = omp_get_wtime()
    Call PreMedia
    time2   = omp_get_wtime()
    write(6, 200) 'Time used in calculating propagator and j symbols: ', time2-time1, 's'

!--- Calculate the propagator terms in spherical coordinate space
    time1   = omp_get_wtime()
    Call PrePotels
    time2   = omp_get_wtime()
    write(6, 200) 'Time used in preparing the interacting matrix:     ', time2-time1, 's'

!--- Calculate the Gogny force
    if( pair%is(1).eq.4 .or. pair%is(2).eq.4 ) Call MediaG

!--- Loop over isototic or isotonic chain: IDT = 1 or 2
   
    xmx0    = xmix
    
    write(*,'(a, 2a1, 3i4)') 'The calculated nucleus:', nuc%nucnam, nuc%npr

!--- Output the the input information from rhf.dat
    Call output
    
!--- If inin.eq.0, input the wave functions and the occupation of each orbit from the files
    if(inin.eq.0) Call inout(1)
    
!--- Initialize the single-particle configurations with a Dirac Woods-Saxon Basis and prepare the DWS BASIS
    Call DBASE

!--- Initialize the variables of the self-consistent iteration
    si      = 1.0d0;    ii      = 1
    tpot    = zero;     tdir    = zero    
    
    do while(ii.lt.maxi .and. si.gt.epsi) 
        si0     = si

    !--- For the case using the existing single-particle solutions
        if(inin.eq.0 .and. ii.eq.1) xmix = 1.0d0

    !--- Print out the results of the current step
        write( 6,100) ii, si, ene%etot(0)/nuc%amas, fermi%ala, ene%rms(0), ene%rch, pair%del, xmix
        write(l6,100) ii, si, ene%etot(0)/nuc%amas, fermi%ala, ene%rms(0), ene%rch, pair%del, xmix
100     format(i4,' si =',f16.10,'  E/A =',f9.4, ' Ef =', 2f9.3, ' R =', f8.4, '  Rc = ', f8.4,' Del = ', 2f8.4, '  mix =',f5.2) 

    !--- Determine the configuration with BCS pairing: BCS pairing
        Call occup(.false.)

    !--- Calculate Density: local densities, non-local ones, and coupling constants
        Call Densit 
        
    !--- Calculate the self-energies: local and non-local potentials in Dirac equation
        time1   = omp_get_wtime()
        Call Potel(.false.)
        time2   = omp_get_wtime()

        tpot    = tpot + time2 - time1
        write( 6, 200) 'Time used in calculating the potentials:           ', time2-time1, 's'
        write(l6, 200) 'Time used in calculating the potentials:           ', time2-time1, 's'

    !--- Calculate the bulk quantities: radii, energy functional and the detailed contributions
        time1   = omp_get_wtime()
        Call Expect( .true.)
        time2   = omp_get_wtime()
        write( 6, 200) 'Time used in calculating the bulk quanities:       ', time2-time1, 's'
        write(l6, 200) 'Time used in calculating the bulk quanities:       ', time2-time1, 's'

        If (isolve == 1) Then
            ! 使用有限差分矩阵法 (FDM) 求解 Dirac 方程
            Call Dirac_FDM(.false.)
        Else
            ! 使用基矢展开法 (Basis Expansion) 求解 Dirac 方程
            Call Dirac(.false.)
        End If
        
        time2   = omp_get_wtime()

        tdir    = tdir + time2 - time1
        write( 6, 200) 'Time used in solving the Dirac-HF equations:       ', time2-time1, 's'
        write(l6, 200) 'Time used in solving the Dirac-HF equations:       ', time2-time1, 's'
        write( 6, *)
        write(l6, *)
        
    !--- Adjust the mixture of the neighboring steps
        if(dabs(si).lt.dabs(si0) .and. xmix .lt.1.0d0) then
            xmix    = xmix*1.0618
            if(xmix.gt.1.0d0) xmix  = 1.0d0
        else if(dabs(si).gt.dabs(si0)) then
            xmix    = xmx0
        end if
    
    !--- Iterative steps
        ii  = ii + 1
    end do

 
    if(ii.ge.maxi .and. si.gt.epsi) write(*,*) 'The calculation is not converged to required accuracy!'
    tpot    = tpot + time2 - time1
    write( 6, 200) 'Time used in calculating the potentials:           ', time2-time1, 's'
    write(l6, 200) 'Time used in calculating the potentials:           ', time2-time1, 's'
        
!--- OUTPUT the single-particle configuration for convenience.
    Call inout(2)

    xte = omp_get_wtime()
    
    write( 6, 200) 'The total time used in solve the Dirac-HF equation:', tdir, 's', tdir/60, 'm'
    write( 6, 200) 'The total time used in calculating the potentials: ', tpot, 's', tpot/60, 'm'
    write( 6, 200) 'The total time used in the whole calculation:      ', xte-xtb, 's', (xte-xtb)/60, 'm'
    write(l6, 200) 'The total time used in solve the Dirac-HF equation:', tdir, 's', tdir/60, 'm'
    write(l6, 200) 'The total time used in calculating the potentials: ', tpot, 's', tpot/60, 'm'
    write(l6, 200) 'The total time used in the whole calculation:      ', xte-xtb, 's', (xte-xtb)/60, 'm'
200 format(2x, a, f12.6, a, f12.6, a)
    Close(l6)

    Stop 'The final Stop!'
!*************************************************************************************************!
!                                                                                                 !
End Program DDRHF                                                                                 !
!                                                                                                 !
!*************************************************************************************************!