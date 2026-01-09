!*************************************************************************************************!
!                                                                                                 !
Module Expectation                                                                                !
!                                                                                                 !
!*************************************************************************************************!
!
!--- This module contains the subroutine to solve the radial Dirac equation
!
Use Define
Use rhflib
Use Meanfield
Use Configuration, only: pair, fermi
Use Center
implicit none

!--- Bulk quantities of calculate nuclei
    type Energy
        
    !--- Charge radius: rch and rcc <==> w/ and w/o CoM; rcd <==> charge density (w/ CoM)
        double precision :: rch, rcc, rcd
        
    !--- Matter radius and particle numbers: total, neutron and proton
        double precision, dimension(0:IBX) :: rms, rmn, xn                      !--- rmn: w/o CoM
        
    !--- Energy: the total and contributions from various channels
        double precision, dimension(0:IBX) :: etot, ecom                        !--- ecom: CoM correction
        double precision, dimension(0:IBX) :: ekin, dsig, dome, drho, dcou      !--- kinetic and direct terms
        double precision, dimension(0:IBX) :: epio, esig, eome, erho, ecou      !--- Exchange terms
        double precision, dimension(0:IBX) :: drtn, drvt
        double precision, dimension(0:IBX) :: ertn, ervt, epai                  !--- epai: pairing energy 

    end type Energy

    type (Energy) :: ene

    double precision, dimension(MSD), private :: fun, df, dg, r2, r
    double precision, dimension(MSD), private :: sig, ome, rho, rtn, pio, rvt, cou, fvt, ftv
    double precision, dimension(MSD, MSD), private :: gg, gf, fg, ff

Contains

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine Expect(lpr)                                                                            !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- This subroutine is used to calculate the bulk properties of nuclei.
!
    logical, intent(in) :: lpr
    
    character(len=6) :: nucnam
    integer :: it, i, npt, i0, ia, iap, ka, ib, n
    double precision :: h, enl, b2, b3, ea, sigs, sig0, sigt, vvs, rch, vv, xmu
    double precision, dimension(0:IBX) :: epart, ekt1, era
    double precision, dimension(IBX) :: drms

    if(lpr) write(l6, *) '***************************** Begin Expect ***********************************'
    
    npt     = well%npt
    h       = well%h
    r2      = well%xr**2
    r       = well%xr

!--- Particle Number
    ene%xn  = zero
    do it = 1, IBX
        fun(1:npt)  = dens(it)%rv(1:npt)*r2(1:npt)
        call simps(fun, npt, h, ene%xn(it))
        ene%xn(it)  = 4.*pi*ene%xn(it)
        ene%xn(0)   = ene%xn(it) + ene%xn(0)
    end do

!--- Center-of-Mass correction with microscopic treatment
    Call CoM(1,.false.);            ene%ecom = ECM

!--- calculate the charge-distribtion
    Call distrib(PCM(0))

!--- RMS radii with microscopic center-of-mass correction
    ene%rms = zero
    do it = 1, IBX
        fun(1:npt)  = dens(it)%rv(1:npt)*r2(1:npt)**2
        call simps(fun, npt, h, ene%rms(it))
        ene%rms(it) = ene%rms(it)*4.*pi
        
        ene%rms(0)  = ene%rms(it)/ene%xn(0) + ene%rms(0)
    end do

!--- Consider the corrections from the center-of-mass motion
    do it = 1, IBX
        ene%rms(it) = ene%rms(it)/ene%xn(it)
        drms(it)    = (ene%rms(it)*two - ene%rms(0)) + (RCM(it)*two - RCM(0))
        ene%rmn(it) = dsqrt( ene%rms(it) )
        ene%rms(it) = dsqrt( ene%rms(it) - drms(it)/ene%xn(0) )        
    end do
    
    ene%rms(0)  = zero
    ene%rmn(0)  = zero
    do it = 1, IBX
        ene%rms(0)  = ene%rms(it)**2*ene%xn(it)/ene%xn(0) + ene%rms(0)
        ene%rmn(0)  = ene%rmn(it)**2*ene%xn(it)/ene%xn(0) + ene%rmn(0)
    end do
    ene%rms(0)  = dsqrt( ene%rms(0) )
    ene%rmn(0)  = dsqrt( ene%rmn(0) )
    ene%rch     = dsqrt(ene%rms(2)**2 + 0.862**2 - 0.336**2*ene%xn(1)/ene%xn(2))
    ene%rcc     = dsqrt(ene%rmn(2)**2 + 0.862**2 - 0.336**2*ene%xn(1)/ene%xn(2))


!--- Calculate the charge radius from the charge density distributions
    do i = 1, npt
        fun(i)  = den%rc(i)*well%xr(i)**4
    end do
    Call simps(fun, npt, h, rch)
    ene%rcd = dsqrt(rch*4.*pi/ene%xn(2))

!--- Pairing-energy   !（BCS中的对能公式）
    ene%epai    = zero
    do it = 1, IBX
        do i = 1, lev(it)%nt
            xmu             = lev(it)%mu(i);            if( lev(it)%lpb(i) ) xmu = xmu - 2
            ene%epai(it)    = -lev(it)%del(i)*lev(it)%spk(i)*xmu + ene%epai(it)
        end do
        
        ene%epai(0) = ene%epai(0) + ene%epai(it)
    enddo

!--- Single-Particle Energy
    epart   = zero
    do it = 1, IBX
        do i = 1, lev(it)%nt
            epart(it)   = epart(it) + lev(it)%vv(i)*lev(it)%ee(i)*lev(it)%mu(i)
        end do
        
        epart(0)    = epart(0) + epart(it)
    end do

!--- Potential Energy of the Direct part
    ene%dsig    = zero;         ene%dome    = zero;         ene%drho    = zero
    ene%dcou    = zero;         ene%drtn    = zero;         ene%drvt    = zero
    do it = 1, IBX
        sig(1:npt)  = dens(it)%rs(1:npt)*dself%sig(1:npt)*r2(1:npt)*cct%gsig(1:npt)
        ome(1:npt)  = dens(it)%rv(1:npt)*dself%ome(1:npt)*r2(1:npt)*cct%gome(1:npt)
        rho(1:npt)  = dens(it)%rv(1:npt)*dself%rho(1:npt)*r2(1:npt)*cct%grho(1:npt)
        cou(1:npt)  = dens(it)%rv(1:npt)*dself%cou(1:npt)*r2(1:npt)
        rtn(1:npt)  = dens(it)%rt(1:npt)*dself%rtn(1:npt)*r2(1:npt)*cct%grtn(1:npt)
        rvt(1:npt)  = dens(it)%rt(1:npt)*dself%rtv(1:npt)*r2(1:npt)*cct%grtn(1:npt) + &
                    & dens(it)%rv(1:npt)*dself%rvt(1:npt)*r2(1:npt)*cct%grho(1:npt)

        call simps(sig, npt, h, ene%dsig(it));          ene%dsig(it)    = ene%dsig(it)*hbc*two*pi
        call simps(ome, npt, h, ene%dome(it));          ene%dome(it)    = ene%dome(it)*hbc*two*pi
        call simps(rho, npt, h, ene%drho(it));          ene%drho(it)    = ene%drho(it)*hbc*two*pi*tauz(it)
        call simps(cou, npt, h, ene%dcou(it));          ene%dcou(it)    = ene%dcou(it)*two*pi*hbc*tauc(it)
        call simps(rtn, npt, h, ene%drtn(it));          ene%drtn(it)    = ene%drtn(it)*hbc*two*pi*tauz(it)
        call simps(rvt, npt, h, ene%drvt(it));          ene%drvt(it)    = ene%drvt(it)*hbc*two*pi*tauz(it)        

        ene%dsig(0) = ene%dsig(0) + ene%dsig(it)
        ene%dome(0) = ene%dome(0) + ene%dome(it)
        ene%drho(0) = ene%drho(0) + ene%drho(it)
        ene%dcou(0) = ene%dcou(0) + ene%dcou(it)
        ene%drtn(0) = ene%drtn(0) + ene%drtn(it)
        ene%drvt(0) = ene%drvt(0) + ene%drvt(it)        
    end do

!--- Potential energies in exchange channels
    ene%esig    = zero;     ene%eome    = zero;     ene%erho    = zero
    ene%ecou    = zero;     ene%epio    = zero;     ene%ertn    = zero
    ene%ervt    = zero
    do it = 1, IBX
        
        if(pset%IE.eq.1) cycle

        sig = zero;         ome = zero;         rho = zero
        cou = zero;         pio = zero;         rtn = zero
        fvt = zero;         ftv = zero
        do ib = 1, chunk(it)%nb
            gg  = denf(ib,it)%gg;           gf  = denf(ib,it)%gf
            fg  = denf(ib,it)%fg;           ff  = denf(ib,it)%ff
        !$OMP PARALLEL DO private(i) SCHEDULE(STATIC)
            do n = 1, npt
                do i = 1, npt
                    sig(n)  = sig(n) + gg(n,i)*FSIG(ib,it)%XPP(n,i) + gf(n,i)*FSIG(ib,it)%XPM(n,i) + &
                            &          fg(n,i)*FSIG(ib,it)%XMP(n,i) + ff(n,i)*FSIG(ib,it)%XMM(n,i)
                    ome(n)  = ome(n) + gg(n,i)*FOME(ib,it)%XPP(n,i) + gf(n,i)*FOME(ib,it)%XPM(n,i) + &
                            &          fg(n,i)*FOME(ib,it)%XMP(n,i) + ff(n,i)*FOME(ib,it)%XMM(n,i)
                    rho(n)  = rho(n) + gg(n,i)*FRHO(ib,it)%XPP(n,i) + gf(n,i)*FRHO(ib,it)%XPM(n,i) + &
                            &          fg(n,i)*FRHO(ib,it)%XMP(n,i) + ff(n,i)*FRHO(ib,it)%XMM(n,i)
                    cou(n)  = cou(n) + gg(n,i)*FCOU(ib,it)%XPP(n,i) + gf(n,i)*FCOU(ib,it)%XPM(n,i) + &
                            &          fg(n,i)*FCOU(ib,it)%XMP(n,i) + ff(n,i)*FCOU(ib,it)%XMM(n,i)
                    pio(n)  = pio(n) + gg(n,i)*FPIO(ib,it)%XPP(n,i) + gf(n,i)*FPIO(ib,it)%XPM(n,i) + &
                            &          fg(n,i)*FPIO(ib,it)%XMP(n,i) + ff(n,i)*FPIO(ib,it)%XMM(n,i)
                    rtn(n)  = rtn(n) + gg(n,i)*FRTN(ib,it)%XPP(n,i) + gf(n,i)*FRTN(ib,it)%XPM(n,i) + &
                            &          fg(n,i)*FRTN(ib,it)%XMP(n,i) + ff(n,i)*FRTN(ib,it)%XMM(n,i)
                    fvt(n)  = fvt(n) + gg(n,i)*FRVT(ib,it)%XPP(n,i) + gf(n,i)*FRVT(ib,it)%XPM(n,i) + &
                            &          fg(n,i)*FRVT(ib,it)%XMP(n,i) + ff(n,i)*FRVT(ib,it)%XMM(n,i)
                    ftv(n)  = ftv(n) + gg(n,i)*FRTV(ib,it)%XPP(n,i) + gf(n,i)*FRTV(ib,it)%XPM(n,i) + &
                            &          fg(n,i)*FRTV(ib,it)%XMP(n,i) + ff(n,i)*FRTV(ib,it)%XMM(n,i)

                end do
            end do
        !$OMP END PARALLEL DO
        end do

        sig(1:npt) = sig(1:npt)*cct%gsig(1:npt);        call simps( sig, npt, h,  ene%esig(it) )          !--- sigma field
        ome(1:npt) = ome(1:npt)*cct%gome(1:npt);        call simps( ome, npt, h,  ene%eome(it) )          !--- omega field
        rho(1:npt) = rho(1:npt)*cct%grho(1:npt);        call simps( rho, npt, h,  ene%erho(it) )          !--- rho-vector coupling
                                                        call simps( cou, npt, h,  ene%ecou(it) )          !--- Coulomb field
        pio(1:npt) = pio(1:npt)*cct%fpio(1:npt);        call simps( pio, npt, h,  ene%epio(it) )          !--- pion pseudo-vector coupling
        rtn(1:npt) = rtn(1:npt)*cct%grtn(1:npt);        call simps( rtn, npt, h,  ene%ertn(it) )          !--- rho-tensor coupling 
        
        rvt(1:npt) = fvt(1:npt)*cct%grho(1:npt) + ftv(1:npt)*cct%grtn(1:npt)
        call simps( rvt, npt, h,  ene%ervt(it) )                                                          !--- rho vector-tensor coupling

        ene%esig(it)    =  ene%esig(it)*hbc*half
        ene%eome(it)    =  ene%eome(it)*hbc*half
        ene%erho(it)    =  ene%erho(it)*hbc*half 
        ene%ecou(it)    =  ene%ecou(it)*hbc*half*tauc(it)
        ene%epio(it)    =  ene%epio(it)*hbc*half
        ene%ertn(it)    =  ene%ertn(it)*hbc*half
        ene%ervt(it)    =  ene%ervt(it)*hbc*half
        
        ene%esig(0) = ene%esig(it) + ene%esig(0)
        ene%eome(0) = ene%eome(it) + ene%eome(0)
        ene%erho(0) = ene%erho(it) + ene%erho(0)
        ene%ecou(0) = ene%ecou(it) + ene%ecou(0)
        ene%epio(0) = ene%epio(it) + ene%epio(0)
        ene%ertn(0) = ene%ertn(it) + ene%ertn(0)
        ene%ervt(0) = ene%ervt(it) + ene%ervt(0)

    end do

!--- Potential energy corresponding with the rearrangement term
    era = zero
    do it = 1, IBX
        fun(1:npt)  = dens(it)%rv(1:npt)*SigR(1:npt)*r2(1:npt)
        Call simps(fun, npt, well%h, era(it))
        era(it) = era(it)*hbc*pi*two
        
        era(0)  = era(0) + era(it)
    end do

!--- Kinetic energy: following relation can be derived from the Dirac equation
    ene%ekin    = epart - two*(ene%dsig + ene%dome + ene%drho + ene%dcou + ene%drtn + ene%drvt + era + &
                &              ene%esig + ene%eome + ene%erho + ene%ecou + ene%ertn + ene%ervt + ene%epio)

!--- Total binding energy
    ene%etot    = ene%ekin + ene%dsig + ene%dome + ene%drho + ene%dcou + ene%drtn + ene%drvt + ene%ecom + &
                & ene%epai + ene%esig + ene%eome + ene%erho + ene%ecou + ene%ertn + ene%ervt + ene%epio

!--- Binding energy per nucleon   
    ea          = ene%etot(0)/nuc%amas

!--- Determine the total binding energy in another way, which can be used to check the convergence
    ekt1        = zero
    do it = 1, IBX
        fun     = zero
        do i0 = 1, lev(it)%nt
            ka  = lev(it)%nk(i0)
            Call tderiv(wav(i0,it)%xg(2), dg(2), npt-1, h)
            Call tderiv(wav(i0,it)%xf(2), df(2), npt-1, h)
            vv          = lev(it)%vv(i0)*lev(it)%mu(i0)
            fun(2:npt)  = fun(2:npt) + vv*( wav(i0, it)%xg(2:npt)*(-df(2:npt) + ka/r(2:npt)*wav(i0, it)%xf(2:npt)) &
                        &                 + wav(i0, it)%xf(2:npt)*( dg(2:npt) + ka/r(2:npt)*wav(i0, it)%xg(2:npt)) &
                        &                 + pset%amu(it)/hbc*(wav(i0,it)%xg(2:npt)**2 - wav(i0,it)%xf(2:npt)**2) )
        end do
        Call simps(fun(2), npt-1, h, ekt1(it))
        ekt1(it)    = ekt1(it)*hbc - pset%amu(it)*ene%xn(it)
        
        ekt1(0)     = ekt1(it) + ekt1(0)
    end do

    enl     =  half*(ekt1(0) + epart(0)) - era(0) + ene%ecom(0) + ene%epai(0) 

!---- Output the corresponding expectations
    if(lpr) then
        write(l6,'(32x, a, 12x,a, 12x, a)') 'neutron','proton','total'

    !---  particle number
        write(l6,'(a, 3f18.6)')             ' particle number       ', ene%xn (1:IBX), ene%xn (0)

    !---  rms-Radius    
        write(l6,'(/,a, 3f18.6)')           ' rms-Radius w/o CoM    ', ene%rmn(1:IBX), ene%rmn(0)

    !---  charge-Radius    
        write(l6,'(a, 18x, f18.6)')         ' charge-Radius w/o CoM ', ene%rcc

    !---  rms-Radius    
        write(l6,'(/,a, 3f18.6)')           ' rms-Radius w/ CoM     ', ene%rms(1:IBX), ene%rms(0)

    !---  charge-Radius    
        write(l6,'(a, 18x, f18.6)')         ' charge-Radius w/ CoM  ', ene%rch

    !--- Fermi energy
        write(l6,'(/,a, 3f18.6)')           ' Fermi Energy          ', fermi%ala

    !--- Effective pairing gap
        write(l6,'(a, 3f18.6)')             ' effective pairing     ', pair%del

    !--- single-particle energy
        write(l6,'(/,a, 3f18.6)')           ' Particle Energy       ', epart(1:IBX), epart(0)

    !--- kinetic energy
        write(l6,'(/,a, 3f18.6)')           ' Kinetic Energy        ', ene%ekin(1:IBX), ene%ekin(0)
        write(l6,'(/,a, 3f18.6)')           ' Kinetic Energy        ', ekt1(1:IBX), ekt1(0)
 
    !--- sigma energy 
        write(l6,'(/,a, 3f18.6)')           ' Direct E-sigma        ', ene%dsig(1:IBX), ene%dsig(0)
        write(l6,'(a, 3f18.6)')             ' Exchange E-sigma      ', ene%esig(1:IBX), ene%esig(0)  
         
    !--- omega energy  
        write(l6,'(/, a, 3f18.6)')          ' Direct E-omega        ', ene%dome(1:IBX), ene%dome(0)
        write(l6,'(a, 3f18.6)')             ' Exchange E-omega      ', ene%eome(1:IBX), ene%eome(0) 

    !--- rho-energy(Vector) 
        write(l6,'(/, a, 3f18.6)')          ' Direct E-rho(V)       ', ene%drho(1:IBX), ene%drho(0) 
        write(l6,'(a, 3f18.6)')             ' Exchange E-rho(V)     ', ene%erho(1:IBX), ene%erho(0) 

    !--- rho-energy(Tensor)
        write(l6,'(/, a, 3f18.6)')          ' Direct E-rho(T)       ', ene%drtn(1:IBX), ene%drtn(0) 
        write(l6,'(a, 3f18.6)')             ' Exchange E-rho(T)     ', ene%ertn(1:IBX), ene%ertn(0) 

    !--- rho-energy(Vector-Tensor)
        write(l6,'(/, a, 3f18.6)')          ' Direct E-rho(VT)      ', ene%drvt(1:IBX), ene%drvt(0)
        write(l6,'(a, 3f18.6)')             ' Exchange E-rho(VT)    ', ene%ervt(1:IBX), ene%ervt(0)

    !--- pion-energy 
        write(l6,'(/, a, 3f18.6)')          ' Direct E-pion         ', zero, zero, zero
        write(l6,'(a, 3f18.6)')             ' Exchange E-pion       ', ene%epio(1:IBX), ene%epio(0) 

    !--- The rearrangement terms
        write(l6,'(/, a, 3f18.6)')          ' Rearrangement         ', era(1:IBX), era(0)

    !--- Coulomb energy (direct part)
        write(l6,'(/, a, 3f18.6)')          ' Direct Coulomb        ', ene%dcou(1:IBX), ene%dcou(0)
        write(l6,'(a, 3f18.6)')             ' Exchange Coulomb      ', ene%ecou(1:IBX), ene%ecou(0)

    !--- Pairing Energy
        write(l6,'(/, a, 3f18.6)')          ' Pairing Energy        ', ene%epai(1:IBX), ene%epai(0)

    !--- center of mass correction
        write(l6,'(/, a, 3f18.6)')          ' E-cm                  ', ene%ecom(1:IBX), ene%ecom(0)

    !--- total energy
        write(l6,'(/, a, 3f18.6)')          ' Total Energy          ', ene%etot(1:IBX), ene%etot(0)
        write(l6,'(a, 36x, f18.6)')         ' Total Energy          ', enl 

    !--- energy per particle
        write(l6,'(a, 36x,f18.6)')          ' Energy per Particle   ', ea

        write(l6, *) '***************************** End Expect *************************************'
    end if

    return

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine Expect                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!*************************************************************************************************!
!                                                                                                 !
End module Expectation                                                                            !
!                                                                                                 !
!*************************************************************************************************!
