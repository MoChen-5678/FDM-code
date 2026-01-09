!*************************************************************************************************!
!                                                                                                 !
Module density                                                                                    !
!                                                                                                 !
!*************************************************************************************************!
Use Define
Use rhflib
Use omp_lib
implicit none

!--- DEFINE THE DENSITIES
    type density1
        double precision, dimension(MSD) :: rs, rv, rt, rc  !总的密度
        double precision, dimension(MSD) :: rs3, rv3, rt3   !rho_n-rho_p
    end type density1
    type (density1) :: den  

!--- NEUTRON AND PROTON DENSITIES
    type density2
        double precision, dimension(MSD) :: rs, rv, rt
    end type density2
    
    type (density2), dimension(2) :: dens 

!--- Non-local density
    type Nonlocal
        double precision, dimension(MSD, MSD) :: fg, ff, gg, gf
    end type Nonlocal
    type (Nonlocal), dimension(NBX, 2) :: denf, denv


Contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine densit                                                                                 !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- To calcualte the densities

    double precision :: h, rel, xvv
    double precision :: ftemp, gtemp, zeta, gftemp
    integer :: it, itx, i, npt, i0, ib, n, nd

    npt = well%npt
    h   = well%h

!--- Calculate rhon and rhop
    do it = 1, 2

    !$OMP PARALLEL DO private(i, i0, xvv) SCHEDULE(STATIC)
        do i = 2, npt
            dens(it)%rs(i)  = zero
            dens(it)%rv(i)  = zero
            dens(it)%rt(i)  = zero
            do i0 = 1, lev(it)%nt
                xvv = lev(it)%vv(i0)*lev(it)%mu(i0)
                dens(it)%rs(i)  = dens(it)%rs(i) + (wav(i0,it)%xg(i)**2 - wav(i0,it)%xf(i)**2)*xvv
                dens(it)%rv(i)  = dens(it)%rv(i) + (wav(i0,it)%xg(i)**2 + wav(i0,it)%xf(i)**2)*xvv
                dens(it)%rt(i)  = dens(it)%rt(i) +  wav(i0,it)%xg(i)*     wav(i0,it)%xf(i)*two*xvv
            end do
            dens(it)%rs(i)  = dens(it)%rs(i)/(4.*pi*well%xr(i)**2)
            dens(it)%rv(i)  = dens(it)%rv(i)/(4.*pi*well%xr(i)**2)
            dens(it)%rt(i)  = dens(it)%rt(i)/(4.*pi*well%xr(i)**2)
        end do
    !$OMP END PARALLEL DO

        dens(it)%rs(1)  = 3.*(dens(it)%rs(2) - dens(it)%rs(3)) + dens(it)%rs(4)  !r=0的点，
        dens(it)%rv(1)  = 3.*(dens(it)%rv(2) - dens(it)%rv(3)) + dens(it)%rv(4)
        dens(it)%rt(1)  = 3.*(dens(it)%rt(2) - dens(it)%rt(3)) + dens(it)%rt(4)
    end do

!--- To calculate rho and rho3
!$OMP PARALLEL DO SCHEDULE(STATIC)
    do i = 1, npt
        den%rs(i)   = dens(1)%rs(i) + dens(2)%rs(i);        den%rs3(i)  = dens(1)%rs(i) - dens(2)%rs(i)
        den%rv(i)   = dens(1)%rv(i) + dens(2)%rv(i);        den%rv3(i)  = dens(1)%rv(i) - dens(2)%rv(i)
        den%rt(i)   = dens(1)%rt(i) + dens(2)%rt(i);        den%rt3(i)  = dens(1)%rt(i) - dens(2)%rt(i)
    end do
!$OMP END PARALLEL DO

!--- Calculate the non-local density for Fock terms
    if(pset%IE.eq.2) then
        do it = 1, 2
            do ib = 1, chunk(it)%nb
                denf(ib,it)%gg  = zero;             denf(ib,it)%ff  = zero
                denf(ib,it)%gf  = zero;             denf(ib,it)%fg  = zero
                
                do n = 1, chunk(it)%id(ib)
                    i0  = chunk(it)%ia(ib) + n
                    xvv = lev(it)%vv(i0)*lev(it)%mu(i0)
                !$OMP PARALLEL DO private(i) SCHEDULE(STATIC)
                    do i = 1, npt
                        denf(ib,it)%gg(:,i) = denf(ib,it)%gg(:,i) + wav(i0,it)%xg(:)*wav(i0,it)%xg(i)*xvv
                        denf(ib,it)%ff(:,i) = denf(ib,it)%ff(:,i) + wav(i0,it)%xf(:)*wav(i0,it)%xf(i)*xvv

                        denf(ib,it)%gf(:,i) = denf(ib,it)%gf(:,i) + wav(i0,it)%xg(:)*wav(i0,it)%xf(i)*xvv
                        denf(ib,it)%fg(:,i) = denf(ib,it)%fg(:,i) + wav(i0,it)%xf(:)*wav(i0,it)%xg(i)*xvv
                    end do 
                !$OMP END PARALLEL DO
                end do
            end do
        end do
    
    !--- non-local density in the isovector channels
        do it = 1, IBX
            itx = 3 - it
        !--- It should be noted that neutron and proton share the same configuration
            do ib = 1, chunk(it)%nb
                denv(ib,it)%gg  = denf(ib,it)%gg + denf(ib,itx)%gg*two
                denv(ib,it)%gf  = denf(ib,it)%gf + denf(ib,itx)%gf*two
                denv(ib,it)%fg  = denf(ib,it)%fg + denf(ib,itx)%fg*two
                denv(ib,it)%ff  = denf(ib,it)%ff + denf(ib,itx)%ff*two
            end do
        end do
    end if

 !--- CALCULATE THE COUPLING CONSTANTS WITH RESPECT TO THE BARYONIC DENSITY
 !$OMP PARALLEL DO private(zeta, gtemp, ftemp, gftemp) SCHEDULE(STATIC)
    do n = 1, npt
        zeta        =  den%rv(n)/pset%rvs
    
        gtemp       =  (one + pset%bsig*(zeta + pset%dsig)**2)
        ftemp       =  (one + pset%csig*(zeta + pset%dsig)**2)
        gftemp      =  (pset%bsig - pset%csig)*(zeta + pset%dsig)

        cct%gsig(n) =  pset%gsig*pset%asig* gtemp/ftemp
        cct%dsig(n) =  pset%gsig*pset%asig*gftemp/ftemp**2*two

        gtemp       =  (one + pset%bome*(zeta + pset%dome)**2)
        ftemp       =  (one + pset%come*(zeta + pset%dome)**2)
        gftemp      =  (pset%bome - pset%come)*(zeta + pset%dome)

        cct%gome(n) =  pset%gome*pset%aome* gtemp/ftemp
        cct%dome(n) =  pset%gome*pset%aome*gftemp/ftemp**2*two

        gtemp       =  dexp(pset%arho*zeta)
        cct%grho(n) =  pset%grho/gtemp
        cct%drho(n) = -pset%arho*cct%grho(n)

        gtemp       =  dexp(pset%artn*zeta)
        cct%grtn(n) =  pset%grtn/gtemp
        cct%drtn(n) = -pset%artn*cct%grtn(n)

        gtemp       =  dexp(pset%apio*zeta)
        cct%fpio(n) =  pset%fpio/gtemp
        cct%dpio(n) = -pset%apio*cct%fpio(n)
    end do
!$OMP END PARALLEL DO
    return

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine densit                                                                             ! 
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine distrib(cm2)                                                                           ! 
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Calculate the charge distributions combined with the center-of-mass correction
!--- and proton charge distribution
!
!--- Expectation of the square of center-of-mass momentum
    double precision, intent(in) :: cm2

!--- Coefficients for the Guassian transformation
    double precision :: xlam, xalp, xbet, r, rp, xlam2
    double precision, dimension(MSD) :: fun
    integer :: i, j
                                                                    
!----ref: appendix of doctoral thesis, with the center-of-mass correction
    xalp    = 2./3.*(0.862**2 - 0.336**2*nuc%npr(1)/nuc%npr(2))
    xbet    = 1.5/cm2

    if(xalp.le.xbet) stop ' The results is bad in Distrib!'
    xlam2   = one/(xalp - xbet)
    xlam    = dsqrt(xlam2)

    do i = 1, well%npt
        r   = well%xr(i)
        do j = 1, well%npt
            rp      = well%xr(j)
            fun(j)  = dexp(-xlam2*(r-rp)**2) - dexp(-xlam2*(r+rp)**2)
            fun(j)  = fun(j)*dens(2)%rv(j)*rp
        end do
        Call simps(fun, well%npt, well%h, den%rc(i))
        den%rc(i)   = den%rc(i)*xlam/(dsqrt(pi)*r)
    end do
    return

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine distrib                                                                            !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!*************************************************************************************************!
!                                                                                                 !
End Module density                                                                                !
!                                                                                                 !
!*************************************************************************************************!