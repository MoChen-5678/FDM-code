!*************************************************************************************************!
!                                                                                                 !
Module DiracR                                                                                     !
!                                                                                                 !
!*************************************************************************************************!
!
Use Define
Use rhflib
Use PotelHF, only: epotl
implicit none

!--- The informations of states
    type states
        double precision :: eig
        integer :: i0, it, kappa, node, matc
        double precision, dimension(2,2) :: gf
        character(len = 8) :: lab
    end type states
    type (states), private :: sta

!--- The wave functions
    double precision, dimension(MSD), private :: sg, sf

!--- The X and Y components
    double precision, dimension(MSD), private :: XG, XF, YG, YF, XGh, YGh, XFh, YFh

Contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine anaori(eigfm, ka, it, i0, h, gmesh, fmesh)                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- This subroutine is used to the original behavior of wave functions
!
    double precision, intent(in) :: eigfm, h
    integer, intent(in) :: ka, it, i0
    double precision, intent(out) :: gmesh, fmesh
    double precision :: alph, beta, original

    original    = -1.d-6

    if(ka.gt.0) then
        alph    = ((two*ka + one)/h + dpotl(it)%vtt(2) + XG(2))*original 
        beta    = eigfm - dpotl(it)%vms(2) - XF(2)
        gmesh   = original
        fmesh   = alph/beta
    else
        alph    = (eigfm - dpotl(it)%vps(2) - YG(2))*original 
        beta    = (two*ka - 1)/h  + dpotl(it)%vtt(2) + YF(2)
        gmesh   = original
        fmesh   = alph/beta
    end if
    return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End subroutine anaori                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Function detgf(eig)                                                                               !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
    double precision, intent(in) :: eig
    double precision :: detgf, eigfm, h, h2, gmesh, fmesh

    double precision, dimension(4) :: ag, af
    double precision, dimension(2, 2) :: aa
    double precision :: r, r1, r2, r4, sg1, sf1, sg2, sf2, u1g, u1f, u2g, u2f, u1x, u2y
    integer :: i0, it, ka, i, mat, jtem, j, jk, item, npt

!--- Initialization of state's quantities
    i0      = sta%i0;       it      = sta%it;       ka      = sta%kappa;        mat     = sta%matc
    h       = well%h;       h2      = well%h*half;  npt     = well%npt;         eigfm   = eig/hbc

!--- Determine the wave functions on the original and infinity position
    sg(1)   = zero;                     sf(1)   = zero
    Call anaori(eigfm, ka, it, i0, h, gmesh, fmesh)
    sg(2)   = gmesh;                    sf(2)   = fmesh
    sg(npt) = zero;                     sf(npt) = one

!--- Integrate inward from (npt - 1) *mesh
    jtem      = npt - mat
    do j = 1, jtem
        i   = npt - j + 1   ! when j=1, i= npt-1+1=npt;  when j=jtem, i = npt-(npt-mat)+1 = mat+1
        r   = well%xr(i)
        r1  = ka/ r
        r2  = ka/(r - h2)  ! back half step
        r4  = ka/(r - h)       ! back one step

        sg1 = sg(i);            sf1 = sf(i)

        do  jk = 1, 4
            go to (36,37,37,38) , jk     ! go to(36) is for jk=1, go to(37) for jk=2 and 3, go to (38) for jk=4
 36         sg2     = sg1;                                      sf2     = sf1
            u1g     = r1  + dpotl(it)%vtt(i) + XG(i);           u1f     = eigfm - dpotl(it)%vms(i) - XF(i)
            u2f     = r1  + dpotl(it)%vtt(i) + YF(i);           u2g     = eigfm - dpotl(it)%vps(i) - YG(i)
            go to 35

 37         sg2     = sg1 + ag(jk-1);                           sf2     = sf1 + af(jk-1)
            u1g     = r2  + dpotl(it)%vtth(i-1) + XGh(i-1);     u1f     = eigfm - dpotl(it)%vmsh(i-1) - XFh(i-1)
            u2f     = r2  + dpotl(it)%vtth(i-1) + YFh(i-1);     u2g     = eigfm - dpotl(it)%vpsh(i-1) - YGh(i-1)
            go to 35

 38         sg2     = sg1 + two*ag(3);                          sf2     = sf1 + two*af(3)
            u1g     = r4  + dpotl(it)%vtt(i-1) + XG(i-1);       u1f     = eigfm - dpotl(it)%vms(i-1) - XF(i-1)
            u2f     = r4  + dpotl(it)%vtt(i-1) + YF(i-1);       u2g     = eigfm - dpotl(it)%vps(i-1) - YG(i-1)
         
 35         ag(jk)  = - h2*(-u1g*sg2 + u1f*sf2);                af(jk)  = - h2*( u2f*sf2 - u2g*sg2)
        end do
        sg2     = (ag(1) + two*(ag(2)+ag(3)) + ag(4))*third;    sf2     = (af(1) + two*(af(2)+af(3)) + af(4))*third
        sg(i-1) = sg1 + sg2;                                    sf(i-1) = sf1 + sf2
    end do
    sta%gf(1,1) = sg(mat);                                  sta%gf(2,1) = sf(mat)

!--- Integrate out to match point
    item    = mat - 1
    do i  = 2, item
        r   = well%xr(i)
        r1  = ka/ r
        r2  = ka/(r + h2)
        r4  = ka/(r + h)
      
        sg1 = sg(i)
        sf1 = sf(i)
     
        do  jk = 1, 4
            go to (46,47,47,48) , jk
 46         sg2     = sg1;                                      sf2     = sf1
            u1g     = r1  + dpotl(it)%vtt(i) + XG(i);           u1f     = eigfm - dpotl(it)%vms(i) - XF(i)
            u2f     = r1  + dpotl(it)%vtt(i) + YF(i);           u2g     = eigfm - dpotl(it)%vps(i) - YG(i)
            go to 45

 47         sg2     = sg1 + ag(jk-1);                           sf2     = sf1 + af(jk-1)
            u1g     = r2  + dpotl(it)%vtth(i) + XGh(i);         u1f     = eigfm - dpotl(it)%vmsh(i) - XFh(i)
            u2f     = r2  + dpotl(it)%vtth(i) + YFh(i);         u2g     = eigfm - dpotl(it)%vpsh(i) - YGh(i)
            go to 45

 48         sg2     = sg1 + two*ag(3);                          sf2     = sf1 + two*af(3)
            u1g     = r4  + dpotl(it)%vtt(i+1) + XG(i+1);       u1f     = eigfm - dpotl(it)%vms(i+1) - XF(i+1)
            u2f     = r4  + dpotl(it)%vtt(i+1) + YF(i+1);       u2g     = eigfm - dpotl(it)%vps(i+1) - YG(i+1)

 45         ag(jk)  = h2*(-u1g*sg2 + u1f*sf2);                  af(jk)  = h2*( u2f*sf2 - u2g*sg2)
        end do
        sg2     = (ag(1) + two*(ag(2)+ag(3)) + ag(4))*third;    sf2     = (af(1) + two*(af(2)+af(3)) + af(4))*third
        sg(i+1) = sg1 + sg2;                                    sf(i+1) = sf1 + sf2
    end do ! end iteration outward to match point.
    sta%gf(1,2) = sg(mat);                                  sta%gf(2,2) = sf(mat)

!--- Return the value of matching determinant
    aa      = sta%gf
    detgf   = DET(2, 2, aa)  !function DET(ma, n, aa) in module RHFlib.f90, which is used to calculate the determinant of matrix aa(mn, n)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End function detgf                                                                                !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine normal(wave, nc, lprx)                                                                 !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Normalize the match matrix and determine the match factors to build continous wave functions
!
    type (waves), intent(out) :: wave
    integer, intent(out) :: nc
    logical, intent(out) :: lprx

    integer :: i1, i2, inode, npt, mat, i, j
    double precision :: s, s1, xstep

    double precision, dimension(2, 2) :: aa, bb
    double precision, dimension(2) :: eem, ae, zz, g
    double precision, dimension(well%npt) :: fun

    xstep   = well%h
    npt     = well%npt
    mat     = sta%matc

!--- Normalize the matching matrix and rescale the wave function
    do i = 1, 2
        s   = zero
        do j = 1, 2
            s           = s + sta%gf(j,i)**2
        end do
        s   = one/dsqrt(s)

        sta%gf(1:2, i)  = s*sta%gf(1:2,i)

        if(i.eq.1) then
            i1          = mat + 1
            i2          = npt
        else
            i1          = 2
            i2          = mat
        end if
        sg(i1:i2)   = s*sg(i1:i2)
        sf(i1:i2)   = s*sf(i1:i2)
    end do

!--- Calculate MT*M
    do i = 1, 2
        do j = 1, 2
            s   = zero
            do i1 = 1, 2
                s   = s + sta%gf(i1,i)*sta%gf(i1,j)
            end do
            aa(i,j) = s
        end do
    end do

!--- The eigenvalues and eigenvectors for aa
    Call sdiag(2,2, aa, eem, aa, zz, 1)

!--- Rescale the wave function with the eigenvectors of aa
    i1              = mat + 1;                      i2              = npt
    wave%xg(i1:i2)  = -sg(i1:i2)*aa(1, 1);          wave%xf(i1:i2)  = -sf(i1:i2)*aa(1, 1)
    
    i1              = 1;                            i2              = mat
    wave%xg(i1:i2)  =  sg(i1:i2)*aa(2, 1);          wave%xf(i1:i2)  =  sf(i1:i2)*aa(2, 1)

!--- Normalization
    fun(1:npt)      = wave%xg(1:npt)**2 + wave%xf(1:npt)**2

    Call simps(fun, npt, xstep, s);                 s               = one/dsqrt(s)

    wave%xg(1:npt)  = wave%xg(1:npt)*s;             wave%xf(1:npt)  = wave%xf(1:npt)*s

!--- Node check
    inode   = 0
    nc      = 0
    lprx    = .false.
    do i = 3, npt

        s = wave%xg(i)*wave%xg(i-1)

        if(s.lt.zero) then
            inode   = inode + 1
            if(abs(i-mat).le.1) then
                lprx        = .true.
            end if
        end if
    end do
    nc      = inode
    return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine normal                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine Inner(ep0, eig, wave, n, ip)                                                           !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Determine the wave functions and single particle energy for given states
!

    integer, intent(in) :: n
    integer, intent(out) :: ip
    double precision, intent(inout) :: eig
    double precision :: epoint, ep0
    type (waves), intent(out) :: wave

    double precision :: es, esi, f1, f2, eig1, eig2, estep, ff
    double precision, parameter :: Etop = 100.0d0 
    double precision, parameter :: tol = 1.d-13
    integer :: nc, i
    logical :: lprx

    data es/5.0d0/, esi/0.001d0/

    estep   = es
    ip      = 0
    epoint  = ep0

!--- Find the wavefunctions and eigenvalues for the I'th state
21  eig     = epoint
    eig1    = eig  + esi 
    eig2    = eig1 + estep 
    f1      = detgf(eig1)
    f2      = detgf(eig2)

    do while( f1*f2.gt.zero .and. eig2.lt.Etop )
        f1      = f2
        eig1    = eig2
        eig2    = eig1 + estep
        f2      = detgf(eig2)
    end do

!--- If the trial value is larger than the top value
    if( eig2 .gt. Etop ) then
        estep   = estep*half
        if( dabs(estep) .lt. dabs(esi) ) then
            ip  = 1
            write(*,'(a, a)') sta%lab, ' In Detgff: exceed the max value of Eig!'
            return
        end if
        go to 21
    end if

!--- Find the eigenvalue
    eig     = rtbrent(detgf, eig1, eig2, f1, f2, tol)
    ff      = detgf(eig)

!--- Normalization and Node check
    Call normal(wave, nc, lprx)

    if(nc .ne. n-1) then
        estep       = estep*half
        if(nc .lt. n-1) epoint   = eig + esi
        if( dabs(estep) .lt. dabs(esi) ) then
            ip  = 1
            write(*,'(a, f12.6, i3, a)') sta%lab, eig, nc, ' In Detgff: search alternation is too small!'
            return
            do i = 1, 100
                write(*, '(f9.3, 2f12.6)') well%xr(i), wave%xg(i), wave%xf(i)
            end do
            return
        end if
        go to 21
    end if

!--- If the match point is close to the node of wave function
    if(lprx) then
        eig         = epoint
        sta%matc    = sta%matc + 15
!        write(*,*)  'The match point is changed from', sta%matc -5, 'to', sta%matc
        go to 21
    end if

    return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine Inner                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine Detgff(ip, lpr)                                                                        !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Solve the Dirac equation with local and non-local potentials in coordinate space
!--- Integro-differential Dirac equation is transformed into equivalent differential one
!--- Shoot method is applied in sovling the equations
!
    logical, intent(in) :: lpr
    integer, intent(out) :: ip
    integer :: ib, i0, n0, nmax, n, it, i, nps, Lxy, ide
    double precision :: sf, eig, epoint, rel, sxy, esi, eps, dex
    double precision, dimension(MSD) :: fun, X0, Y0
    type (waves) :: wave
    
    double precision :: Ebot = -80.d0 
    data esi/0.001d0/

    eps = epsi*0.1d0
    
!--- Loop over neutron and proton
    do it = 1, 2
        sta%it  = it

    !--- Loop over the single particle levels
        do ib = 1, chunk(it)%nb
    
        !--- to set up the quantum numbers for the I0'th state
            n0          = chunk(it)%ia(ib)
            nmax        = chunk(it)%id(ib)

        !--- For positive energy states
        !--- Starting value of the eigenvalues
            eig         = Ebot + 0.5d0

        !--- Epoint: current starting value of the energy
            epoint      = eig
            do n = 1, nmax
                i0              = n0 + n
                sta%i0          = i0
                sta%kappa       = lev(it)%nk(i0)
                sta%node        = lev(it)%nr(i0) - 1
                sta%matc        = lev(it)%ma(i0)
                sta%lab         = lev(it)%tb(i0)

                dex = one
                ide = 1
                do while(ide.le.50 .and. dex.gt.eps)

                !--- Calculate XG, XF, YG, YF
                    Call DetXY(i0, it)

                !--- Finding the eigenvalue
                    Call Inner(epoint, eig, wave, n, ip)
                    if(ip.eq.1) return

                    dex             = abs( eig - lev(it)%ee(i0) )
                    wav(i0,it)%xg   = wave%xg
                    wav(i0,it)%xf   = wave%xf
                    lev(it)%ee(i0)  = eig
                    lev(it)%ma(i0)  = sta%matc
                    
                  !  write(6, '(2x,a, 3i5, f12.6)') lev(it)%tb(i0), i0, it, ide, eig
                    ide = ide + 1

                end do

            !--- New starting point for the next iteration
                epoint          = eig + esi
               ! write(6, '(2x,a, 3i5, f12.6)') lev(it)%tb(i0), i0, it, ide, eig

            end do
            
        end do

!--- Calculate the inner part of wavefunctions: the values of localized possibility density can be used to determine the resonances.
        nps     = well%npt
        do i0 = 1, lev(it)%nt
            
            fun = ( wav(i0,it)%xg**2 + wav(i0,it)%xf**2 )*well%xr**2
            Call simps(fun, nps, well%h, rel)
            wav(i0,it)%r2  = dsqrt(rel)
            
            wav(i0,it)%rho  = (wav(i0,it)%xg**2 + wav(i0,it)%xf**2)/(4.0*pi*well%xr**2)
            
            fun = wav(i0,it)%xg**2;         Call simps(fun, nps, well%h, wav(i0,it)%gg)
            fun = wav(i0,it)%xf**2;         Call simps(fun, nps, well%h, wav(i0,it)%ff)

        end do


    end do
    ip  = 0
    return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine Detgff                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine DetXY(i0, it)                                                                          !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Calculate the X and Y Components, and determine equivalent localized XG, XF, YG, and YF terms
!

    integer, intent(in) :: i0, it

    integer :: i, j, n, ib, ka
    double precision, dimension(MSD) :: S, XX, YY

    n   = well%npt
    ka  = lev(it)%nk(i0)
    ib  = chunk(it)%ib(ka)

!$OMP PARALLEL DO private(i, j) SCHEDULE(STATIC)
    do j = 2, n
        XX(j)  = zero
        YY(j)  = zero
        
        do i = 1, n
            YY(j)   = YY(j) + epotl(ib,it)%YG(j,i)*wav(i0,it)%xg(i) + epotl(ib,it)%YF(j,i)*wav(i0,it)%xf(i)
            XX(j)   = XX(j) + epotl(ib,it)%XG(j,i)*wav(i0,it)%xg(i) + epotl(ib,it)%XF(j,i)*wav(i0,it)%xf(i)
        end do
    end do
        
!$OMP END PARALLEL DO
    YY(1)   = YY(4) + 3.0*( YY(2) - YY(3) )
    XX(1)   = XX(4) + 3.0*( XX(2) - XX(3) )
    
    S(2:n)  = one/( wav(i0,it)%xg(2:n)**2 + wav(i0,it)%xf(2:n)**2 )
    
    XG(2:n) = ( XX(2:n) )*wav(i0,it)%xg(2:n)*S(2:n)        
    XF(2:n) = ( XX(2:n) )*wav(i0,it)%xf(2:n)*S(2:n)
    YG(2:n) = ( YY(2:n) )*wav(i0,it)%xg(2:n)*S(2:n)        
    YF(2:n) = ( YY(2:n) )*wav(i0,it)%xf(2:n)*S(2:n)

    XG(1)   = 3.0*(XG(2) - XG(3)) + XG(4);        XF(1)   = 3.0*(XF(2) - XF(3)) + XF(4)
    YG(1)   = 3.0*(YG(2) - YG(3)) + YG(4);        YF(1)   = 3.0*(YF(2) - YF(3)) + YF(4)
    
    Call intpol6(XG, XGh, n);               Call intpol6(YG, YGh, n)
    Call intpol6(XF, XFh, n);               Call intpol6(YF, YFh, n)

    return
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine DetXY                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!*************************************************************************************************!
!                                                                                                 !
End Module DiracR                                                                                 !
!                                                                                                 !
!*************************************************************************************************!