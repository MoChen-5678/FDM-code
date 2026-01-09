!*************************************************************************************************!
!
Module Meanfield
!
!*************************************************************************************************!
!
!--- This module provides the subroutine to calculate the Hartree and Fock mean field, also includes
!    the rearrangement terms induced by the density dependent coupling strength
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
Use Define
Use Combination
Use density
Use RHFlib
Implicit None

!--- Direct part of the mean fields 
    type dselfs
        double precision, dimension(MSD) :: sig, ome, rho, cou, rtn, rvt, rtv 
    end type dselfs
    type (dselfs) :: dself

!--- Rearrangement term
    double precision, dimension(MSD) :: SigR 

    TYPE (NLMF_T), DIMENSION(NBX, IBX) :: FSIG, FOME, FRHO, FCOU, FPIO, FRTN, FRTV, FRVT

Contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine Hartree                                                                                !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- This subroutine is used to calculate the direct part of fields with
!--- Green function technical
!
!--- The direct part of meson and phonon fields
    double precision, dimension(MSD) ::r2
    double precision :: h, h6
    integer :: npt, i, n

    h           = well%h;                   npt         = well%npt;             r2          = well%xr**2
    h6          = h/6.0d0

!--- Performing the integration over r', and n for the index of r
!$OMP PARALLEL DO PRIVATE(n, i) SCHEDULE(STATIC)
    do n = 1, npt
        
        dself%sig(n)    = zero;         dself%ome(n)    = zero;         dself%rho(n)    = zero
        dself%rtn(n)    = zero;         dself%rvt(n)    = zero;         dself%rtv(n)    = zero
        dself%cou(n)    = zero 
        do i = 1, npt   
            dself%sig(n)    = dself%sig(n) - XBIS(0)%wsig(n,i)*den%rs (i)*cct%gsig(i)*r2(i)*h6
            dself%ome(n)    = dself%ome(n) + XBIS(0)%wome(n,i)*den%rv (i)*cct%gome(i)*r2(i)*h6
            dself%rho(n)    = dself%rho(n) + XBIS(0)%wrho(n,i)*den%rv3(i)*cct%grho(i)*r2(i)*h6
            dself%cou(n)    = dself%cou(n) + XBIS(0)%wcou(n,i)*dens(2)%rv(i)*r2(i)/alphi*4.0*pi*h6
        
            if(pset%grtn.eq.zero) cycle
            dself%rtn(n)    = dself%rtn(n) - XBIT(0,4)%wrtn(n,i)*den%rt3(i)*cct%grtn(i)*r2(i)*h6
            dself%rvt(n)    = dself%rvt(n) + XBIV(0,2)%wrvt(n,i)*den%rt3(i)*cct%grtn(i)*r2(i)*h6
            dself%rtv(n)    = dself%rtv(n) - XBIV(0,2)%wrtv(n,i)*den%rv3(i)*cct%grho(i)*r2(i)*h6

        end do
        dself%rtn(n)   = dself%rtn(n)*(pset%amrho/hbc)**2 + den%rt3(n)*cct%grtn(n)*two*third
        dself%rvt(n)   = dself%rvt(n)*(pset%amrho/hbc)
        dself%rtv(n)   = dself%rtv(n)*(pset%amrho/hbc)
    end do
!$OMP END PARALLEL DO
    return

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine Hartree                                                                            !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine Fock                                                                                   !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Calculate the non-local mean field contributed by the Fock terms
!
!-------------------------------------------------------------------------------------------------!

    integer :: it, ia, ib, ka, kb, i

    do it = 1, IBX

    !$OMP PARALLEL Do private( ia, ib, ka, kb, i )
        do ia = 1, chunk(it)%nb
            ka  = chunk(it)%ka(ia) 
        
        !--- Initialize
            FSIG(ia,it)%XPP = zero;     FOME(ia,it)%XPP = zero;     FRHO(ia,it)%XPP = zero;     FCOU(ia,it)%XPP = zero
            FSIG(ia,it)%XPM = zero;     FOME(ia,it)%XPM = zero;     FRHO(ia,it)%XPM = zero;     FCOU(ia,it)%XPM = zero
            FSIG(ia,it)%XMP = zero;     FOME(ia,it)%XMP = zero;     FRHO(ia,it)%XMP = zero;     FCOU(ia,it)%XMP = zero
            FSIG(ia,it)%XMM = zero;     FOME(ia,it)%XMM = zero;     FRHO(ia,it)%XMM = zero;     FCOU(ia,it)%XMM = zero

            FPIO(ia,it)%XPP = zero;     FRTN(ia,it)%XPP = zero;     FRVT(ia,it)%XPP = zero;     FRTV(ia,it)%XPP = zero
            FPIO(ia,it)%XPM = zero;     FRTN(ia,it)%XPM = zero;     FRVT(ia,it)%XPM = zero;     FRTV(ia,it)%XPM = zero
            FPIO(ia,it)%XMP = zero;     FRTN(ia,it)%XMP = zero;     FRVT(ia,it)%XMP = zero;     FRTV(ia,it)%XMP = zero
            FPIO(ia,it)%XMM = zero;     FRTN(ia,it)%XMM = zero;     FRVT(ia,it)%XMM = zero;     FRTV(ia,it)%XMM = zero


        !--- Loop over occupied states
            do ib = 1, chunk(it)%nb
                kb  = chunk(it)%ka(ib)
                
            !--- Calculate the time-component terms for sigma(scalar), omega(vector) and rho(vector), and Coulomb field
                FSIG(ia,it)%XPP = FSIG(ia,it)%XPP + denf(ib,it)%gg*CSIG(ka,kb)%XPP
                FSIG(ia,it)%XPM = FSIG(ia,it)%XPM - denf(ib,it)%gf*CSIG(ka,kb)%XPP
                FSIG(ia,it)%XMP = FSIG(ia,it)%XMP - denf(ib,it)%fg*CSIG(ka,kb)%XPP
                FSIG(ia,it)%XMM = FSIG(ia,it)%XMM + denf(ib,it)%ff*CSIG(ka,kb)%XPP

                FOME(ia,it)%XPP = FOME(ia,it)%XPP + denf(ib,it)%gg*COVT(ka,kb)%XPP + denf(ib,it)%ff*COVS(ka,kb)%XPP
                FOME(ia,it)%XPM = FOME(ia,it)%XPM + denf(ib,it)%gf*COVT(ka,kb)%XPP + denf(ib,it)%fg*COVS(ka,kb)%XPM
                FOME(ia,it)%XMP = FOME(ia,it)%XMP + denf(ib,it)%fg*COVT(ka,kb)%XPP + denf(ib,it)%gf*COVS(ka,kb)%XPM
                FOME(ia,it)%XMM = FOME(ia,it)%XMM + denf(ib,it)%ff*COVT(ka,kb)%XPP + denf(ib,it)%gg*COVS(ka,kb)%XMM

                if(it.eq.2) then
                    FCOU(ia,it)%XPP = FCOU(ia,it)%XPP + denf(ib,it)%gg*CAVT(ka,kb)%XPP + denf(ib,it)%ff*CAVS(ka,kb)%XPP
                    FCOU(ia,it)%XPM = FCOU(ia,it)%XPM + denf(ib,it)%gf*CAVT(ka,kb)%XPP + denf(ib,it)%fg*CAVS(ka,kb)%XPM
                    FCOU(ia,it)%XMP = FCOU(ia,it)%XMP + denf(ib,it)%fg*CAVT(ka,kb)%XPP + denf(ib,it)%gf*CAVS(ka,kb)%XPM
                    FCOU(ia,it)%XMM = FCOU(ia,it)%XMM + denf(ib,it)%ff*CAVT(ka,kb)%XPP + denf(ib,it)%gg*CAVS(ka,kb)%XMM
                end if

                FRHO(ia,it)%XPP = FRHO(ia,it)%XPP + denv(ib,it)%gg*CRVT(ka,kb)%XPP + denv(ib,it)%ff*CRVS(ka,kb)%XPP
                FRHO(ia,it)%XPM = FRHO(ia,it)%XPM + denv(ib,it)%gf*CRVT(ka,kb)%XPP + denv(ib,it)%fg*CRVS(ka,kb)%XPM
                FRHO(ia,it)%XMP = FRHO(ia,it)%XMP + denv(ib,it)%fg*CRVT(ka,kb)%XPP + denv(ib,it)%gf*CRVS(ka,kb)%XPM
                FRHO(ia,it)%XMM = FRHO(ia,it)%XMM + denv(ib,it)%ff*CRVT(ka,kb)%XPP + denv(ib,it)%gg*CRVS(ka,kb)%XMM

        !--- Pion (pseudo-vector couplings)
                if(pset%fpio.ne.zero) then
                    FPIO(ia,it)%XPP = FPIO(ia,it)%XPP + denv(ib,it)%gg*CPIO(ka,kb)%XPP
                    FPIO(ia,it)%XPM = FPIO(ia,it)%XPM + denv(ib,it)%gf*CPIO(ka,kb)%XPM
                    FPIO(ia,it)%XMP = FPIO(ia,it)%XMP + denv(ib,it)%fg*CPIO(ka,kb)%XMP
                    FPIO(ia,it)%XMM = FPIO(ia,it)%XMM + denv(ib,it)%ff*CPIO(ka,kb)%XMM
                end if

        !--- Calculate the 0 component of rho-Tensor coupling (T1) and (VT1); Rho tensor couplings (T2) and (VT2)
                if(pset%grtn.ne.zero) then
                    FRVT(ia,it)%XPP = FRVT(ia,it)%XPP + denv(ib,it)%gf*CVTT(ka,kb)%XPP + denv(ib,it)%fg*CVTS(ka,kb)%XPP
                    FRVT(ia,it)%XPM = FRVT(ia,it)%XPM + denv(ib,it)%gg*CVTT(ka,kb)%XPM + denv(ib,it)%ff*CVTS(ka,kb)%XPM
                    FRVT(ia,it)%XMP = FRVT(ia,it)%XMP + denv(ib,it)%ff*CVTT(ka,kb)%XMP + denv(ib,it)%gg*CVTS(ka,kb)%XMP
                    FRVT(ia,it)%XMM = FRVT(ia,it)%XMM + denv(ib,it)%fg*CVTT(ka,kb)%XMM + denv(ib,it)%gf*CVTS(ka,kb)%XMM

                    FRTV(ia,it)%XPP = FRTV(ia,it)%XPP + denv(ib,it)%fg*CTVT(ka,kb)%XPP + denv(ib,it)%gf*CTVS(ka,kb)%XPP
                    FRTV(ia,it)%XPM = FRTV(ia,it)%XPM + denv(ib,it)%ff*CTVT(ka,kb)%XPM + denv(ib,it)%gg*CTVS(ka,kb)%XPM
                    FRTV(ia,it)%XMP = FRTV(ia,it)%XMP + denv(ib,it)%gg*CTVT(ka,kb)%XMP + denv(ib,it)%ff*CTVS(ka,kb)%XMP
                    FRTV(ia,it)%XMM = FRTV(ia,it)%XMM + denv(ib,it)%gf*CTVT(ka,kb)%XMM + denv(ib,it)%fg*CTVS(ka,kb)%XMM

                    FRTN(ia,it)%XPP = FRTN(ia,it)%XPP + denv(ib,it)%ff*CRTT(ka,kb)%XPP + denv(ib,it)%gg*CRTS(ka,kb)%XPP
                    FRTN(ia,it)%XPM = FRTN(ia,it)%XPM + denv(ib,it)%fg*CRTT(ka,kb)%XPM + denv(ib,it)%gf*CRTS(ka,kb)%XPM
                    FRTN(ia,it)%XMP = FRTN(ia,it)%XMP + denv(ib,it)%gf*CRTT(ka,kb)%XMP + denv(ib,it)%fg*CRTS(ka,kb)%XMP
                    FRTN(ia,it)%XMM = FRTN(ia,it)%XMM + denv(ib,it)%gg*CRTT(ka,kb)%XMM + denv(ib,it)%ff*CRTS(ka,kb)%XMM
                end if

            end do  !--- end loop of ib

        !--- Multiply with the density-dependent coupling contants
            do i = 1, well%npt
                FSIG(ia,it)%XPP(:,i)    = FSIG(ia,it)%XPP(:,i)*cct%gsig(i)
                FSIG(ia,it)%XPM(:,i)    = FSIG(ia,it)%XPM(:,i)*cct%gsig(i)
                FSIG(ia,it)%XMP(:,i)    = FSIG(ia,it)%XMP(:,i)*cct%gsig(i)
                FSIG(ia,it)%XMM(:,i)    = FSIG(ia,it)%XMM(:,i)*cct%gsig(i)
        
                FOME(ia,it)%XPP(:,i)    = FOME(ia,it)%XPP(:,i)*cct%gome(i)
                FOME(ia,it)%XPM(:,i)    = FOME(ia,it)%XPM(:,i)*cct%gome(i)
                FOME(ia,it)%XMP(:,i)    = FOME(ia,it)%XMP(:,i)*cct%gome(i)
                FOME(ia,it)%XMM(:,i)    = FOME(ia,it)%XMM(:,i)*cct%gome(i)

                FRHO(ia,it)%XPP(:,i)    = FRHO(ia,it)%XPP(:,i)*cct%grho(i)
                FRHO(ia,it)%XPM(:,i)    = FRHO(ia,it)%XPM(:,i)*cct%grho(i)
                FRHO(ia,it)%XMP(:,i)    = FRHO(ia,it)%XMP(:,i)*cct%grho(i)
                FRHO(ia,it)%XMM(:,i)    = FRHO(ia,it)%XMM(:,i)*cct%grho(i)

                FPIO(ia,it)%XPP(:,i)    = FPIO(ia,it)%XPP(:,i)*cct%fpio(i)
                FPIO(ia,it)%XPM(:,i)    = FPIO(ia,it)%XPM(:,i)*cct%fpio(i)
                FPIO(ia,it)%XMP(:,i)    = FPIO(ia,it)%XMP(:,i)*cct%fpio(i)
                FPIO(ia,it)%XMM(:,i)    = FPIO(ia,it)%XMM(:,i)*cct%fpio(i)

                FRTN(ia,it)%XPP(:,i)    = FRTN(ia,it)%XPP(:,i)*cct%grtn(i)
                FRTN(ia,it)%XPM(:,i)    = FRTN(ia,it)%XPM(:,i)*cct%grtn(i)
                FRTN(ia,it)%XMP(:,i)    = FRTN(ia,it)%XMP(:,i)*cct%grtn(i)
                FRTN(ia,it)%XMM(:,i)    = FRTN(ia,it)%XMM(:,i)*cct%grtn(i)

                FRVT(ia,it)%XPP(:,i)    = FRVT(ia,it)%XPP(:,i)*cct%grtn(i)
                FRVT(ia,it)%XPM(:,i)    = FRVT(ia,it)%XPM(:,i)*cct%grtn(i)
                FRVT(ia,it)%XMP(:,i)    = FRVT(ia,it)%XMP(:,i)*cct%grtn(i)
                FRVT(ia,it)%XMM(:,i)    = FRVT(ia,it)%XMM(:,i)*cct%grtn(i)

                FRTV(ia,it)%XPP(:,i)    = FRTV(ia,it)%XPP(:,i)*cct%grho(i)
                FRTV(ia,it)%XPM(:,i)    = FRTV(ia,it)%XPM(:,i)*cct%grho(i)
                FRTV(ia,it)%XMP(:,i)    = FRTV(ia,it)%XMP(:,i)*cct%grho(i)
                FRTV(ia,it)%XMM(:,i)    = FRTV(ia,it)%XMM(:,i)*cct%grho(i)
            end do

        end do  !--- end loop of ia
    !$OMP END PARALLEL Do

    end do  !--- end loop of it
    return

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine Fock                                                                               !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine Rearrange                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- This subroutine is used to calculate the rearrangement term which are
!--- original from the density-dependence of meson-nucleon coupling.
!
    integer :: it, ib, n, i
    double precision :: vv
    double precision, dimension(MSD) :: r2, TmpR, dsig, dome, drho, drtn, dpio 

!--- The contributions from the direct parts of meson-nucleon couplings
!$OMP PARALLEL DO PRIVATE(i) SCHEDULE(STATIC)
    do i = 1, well%npt
        SigR(i) = den%rs(i) *cct%dsig(i)*dself%sig(i) + den%rv(i) *cct%dome(i)*dself%ome(i) + &
                & den%rv3(i)*cct%drho(i)*dself%rho(i) + den%rt3(i)*cct%drtn(i)*dself%rtn(i) + &
                & den%rv3(i)*cct%drho(i)*dself%rvt(i) + den%rt3(i)*cct%drtn(i)*dself%rtv(i) 
    end do
!$OMP END PARALLEL DO

    if(pset%IE.eq.1) go to 33

!--- The contributions from the exchange parts
    r2      = well%xr**2*4.0*pi
    dsig    = zero;         dome    = zero;         drho    = zero
    dpio    = zero;         drtn    = zero 
    do it = 1, 2
        do ib = 1, chunk(it)%nb

        !$OMP PARALLEL DO PRIVATE(n, i) SCHEDULE(STATIC)
            do n = 1, well%npt
                do i = 1, well%npt
                    dsig(n) = dsig(n) + denf(ib,it)%gg(n,i)*FSIG(ib,it)%XPP(n,i) &
                            &         + denf(ib,it)%gf(n,i)*FSIG(ib,it)%XPM(n,i) &
                            &         + denf(ib,it)%fg(n,i)*FSIG(ib,it)%XMP(n,i) &
                            &         + denf(ib,it)%ff(n,i)*FSIG(ib,it)%XMM(n,i)
                    dome(n) = dome(n) + denf(ib,it)%gg(n,i)*FOME(ib,it)%XPP(n,i) &
                            &         + denf(ib,it)%gf(n,i)*FOME(ib,it)%XPM(n,i) &
                            &         + denf(ib,it)%fg(n,i)*FOME(ib,it)%XMP(n,i) &
                            &         + denf(ib,it)%ff(n,i)*FOME(ib,it)%XMM(n,i)
                    drho(n) = drho(n) + denf(ib,it)%gg(n,i)*FRHO(ib,it)%XPP(n,i) &
                            &         + denf(ib,it)%gf(n,i)*FRHO(ib,it)%XPM(n,i) &
                            &         + denf(ib,it)%fg(n,i)*FRHO(ib,it)%XMP(n,i) &
                            &         + denf(ib,it)%ff(n,i)*FRHO(ib,it)%XMM(n,i) &
                            &         + denf(ib,it)%gg(n,i)*FRVT(ib,it)%XPP(n,i) &
                            &         + denf(ib,it)%gf(n,i)*FRVT(ib,it)%XPM(n,i) &
                            &         + denf(ib,it)%fg(n,i)*FRVT(ib,it)%XMP(n,i) &
                            &         + denf(ib,it)%ff(n,i)*FRVT(ib,it)%XMM(n,i)
                    dpio(n) = dpio(n) + denf(ib,it)%gg(n,i)*FPIO(ib,it)%XPP(n,i) &
                            &         + denf(ib,it)%gf(n,i)*FPIO(ib,it)%XPM(n,i) &
                            &         + denf(ib,it)%fg(n,i)*FPIO(ib,it)%XMP(n,i) &
                            &         + denf(ib,it)%ff(n,i)*FPIO(ib,it)%XMM(n,i)
                    drtn(n) = drtn(n) + denf(ib,it)%gg(n,i)*FRTN(ib,it)%XPP(n,i) &
                            &         + denf(ib,it)%gf(n,i)*FRTN(ib,it)%XPM(n,i) &
                            &         + denf(ib,it)%fg(n,i)*FRTN(ib,it)%XMP(n,i) &
                            &         + denf(ib,it)%ff(n,i)*FRTN(ib,it)%XMM(n,i) &
                            &         + denf(ib,it)%gg(n,i)*FRTV(ib,it)%XPP(n,i) &
                            &         + denf(ib,it)%gf(n,i)*FRTV(ib,it)%XPM(n,i) &
                            &         + denf(ib,it)%fg(n,i)*FRTV(ib,it)%XMP(n,i) &
                            &         + denf(ib,it)%ff(n,i)*FRTV(ib,it)%XMM(n,i) 
                end do
            end do
        !$OMP END PARALLEL DO
        end do
    end do
    
!$OMP PARALLEL DO private(i) SCHEDULE(STATIC)
    do i = 1, well%npt
        TmpR(i) = dsig(i)*cct%dsig(i) + dome(i)*cct%dome(i) + drho(i)*cct%drho(i) + &
                & dpio(i)*cct%dpio(i) + drtn(i)*cct%drtn(i)

        SigR(i) = SigR(i) + TmpR(i)/r2(i)
    end do
!$OMP END PARALLEL DO

 33 SigR    = SigR/pset%rvs 

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine Rearrange                                                                          !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!*************************************************************************************************!
!
End Module Meanfield
!
!*************************************************************************************************!
