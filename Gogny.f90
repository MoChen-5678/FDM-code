!*************************************************************************************************!
!                                                                                                 !
Module Gogny                                                                                      !
!                                                                                                 !
!*************************************************************************************************!
!
Use Define
Use rhflib
Use Combination, only: C3J, C3L, C3S, Folding
Use density
Use DBESIK
Implicit None

!--- Combination of the factors and radial part
    Type GognyC
        double precision, dimension(MSD, MSD) :: gg, gf, ff
    end type GognyC
    
    Type (GognyC), dimension(NBX, NBX), private :: GoV

!--- Multipole expansion of the Gogny force
    Type GognyM
        double precision, dimension(MSD, MSD) :: vlam
    End type GognyM
    
    Type (GognyM), dimension(0:NBX+1, 2), private :: VGL

    Type GognyT
        double precision, dimension(2*MSD) :: vlam
    End type GognyT
    
    Type (GognyT), dimension(0:NBX+1, 2), private :: GoFC

    Type GognyI
        double precision, dimension(NTX, NTX) :: vg
    End type GognyI
    
    Type (GognyI), dimension(2) :: GME

!--- Definition for Gogny force
    type GognyF
        character(len=10) :: txtgog
        double precision, dimension(2) :: gr, gw, gb, gh, gm
        double precision :: w0, t3, vg
    end type GognyF

    type (GognyF) :: D1, DSA1, GPF, D1N, D1M

    data D1%txtgog/'Gogny D1'/
    data D1%gr/   0.7d0,  1.20d0/, D1%gw/-402.40d0, -21.30d0/, D1%gb/-100.d0,-11.77d0/ 
    data D1%gh/-496.2d0, 37.27d0/, D1%gm/ -23.56d0, -68.81d0/
    data D1%w0/-115.0d0/, D1%t3/1350.d0/, D1%vg/1.0d0/

    data DSA1%txtgog/'Gogny DSA1'/
    data DSA1%gr/0.7d0,1.2d0/, DSA1%gw/-1720.30d0, 103.639d0/, DSA1%gb/ 1300.00d0,-163.483d0/
    data DSA1%gh/-1813.53d0, 162.812d0/, DSA1%gm/ 1397.60d0,-223.934d0/
    data DSA1%w0/-130.d0/, DSA1%t3/1350.d0/, DSA1%vg/1.0d0/

    data D1N%txtgog/'Gogny D1N'/
    data D1N%gr/0.8d0,1.2d0/, D1N%gw/-2047.61d0, 293.02d0/, D1N%gb/ 1700.00d0,-300.78d0/
    data D1N%gh/-2414.93d0, 414.59d0/, D1N%gm/ 1519.35d0,-316.84d0/
    data D1N%w0/-115.d0/, D1N%t3/1609.46d0/, D1N%vg/1.0d0/

    data D1M%txtgog/'Gogny D1M'/
    data D1M%gr/0.5d0,1.0d0/, D1M%gw/-12797.57d0, 490.95d0/, D1M%gb/ 14048.85d0,-752.27d0/
    data D1M%gh/-15144.43d0, 675.12d0/, D1M%gm/ 11963.89d0,-693.57d0/
    data D1M%w0/-115.36d0/, D1M%t3/1562.22d0/, D1M%vg/1.0d0/

!--- Parameters for pairing with density-dependent delta force
    type DeltaF
        double precision :: V0, Gamma
    end type DeltaF
    type (DeltaF), private :: YKT, GSNL, DBCS

!--- International Journal of Modern E:, Toki 2004
    data YKT%V0/ 320.0d0/, YKT%Gamma/0.0d0/

!--- Phys. Rev. C 64, 064321(2001).
    data GSNL%V0/ 500.00d0/, GSNL%Gamma/1.0d0/

    Contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine MediaG                                                                                 !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
    Integer :: i, nd1, nd2, L, NZ, KODE, KM, nd, k1, k2, lu1, lu2, ld1, ld2, LI, LM
    double precision :: r1, r2, xrr(2), ALPH, Y(NBX+2), h6, gg, gf, ff, ex, er
    double precision, dimension(MSD) :: pgt
    double precision, dimension(2*MSD) :: rh2

!--- Coefficient of the Gogny force with/without spin exchange operator: GD/GA
    double precision, dimension(2) :: GA, GD
    
!--- Calculate the factors in the pairing interacting matrix elements
    do i = 1, 2
        GA(i)   = GPF%vg*( GPF%gw(i) - GPF%gh(i) )
        GD(i)   = GPF%vg*( GPF%gb(i) - GPF%gm(i) )
    end do

!--- Calculate the radial part of the Gogny force
    nd  = well%npt;                 h6      = well%h/6.0d0      !--- mesh information

    do nd1 = 2, 2*nd - 1            !--- Prepare for folding the integral weights
        rh2(nd1)    = (nd1 - 1)*half*well%h
    end do
    rh2(1)  = well%xr(1)
    
!--- Variables for Calling Bessel functions
    KM  = NBX + 2;                  ALPH    = half;                 KODE    = 2
    
!--- Calculate the multipole expansion: V(lambda, i, r, r')
    do nd1 = 1, nd          !--- Loop over r
        r1  = well%xr(nd1)

        do nd2  = 1, 2*nd - 1       !--- Loop over r' while divided by 2
            r2  = rh2(nd2)
            xrr = two*r1*r2/GPF%gr**2

            
            do i = 1, 2
                er  = (r1 - r2)**2/GPF%gr(i)**2
                ex  = dexp(-er)*dsqrt(two*pi/xrr(i))
                Call DBESI(xrr(i), ALPH, KODE, KM, Y, NZ)
                if(NZ.gt.0) Stop ' DBESI in Gogny: NZ > 0!'
                    
                do L = 0, NBX+1
                    GoFC(L,i)%vlam(nd2) = Y(L+1)*ex
                end do
            end do
        end do

    !--- Folding with integral weights for r'
        do i = 1, 2
            do L = 0, NBX+1
                Call Folding( GoFC(L,i)%vlam, nd, pgt )
                VGL(L,i)%vlam(nd1, 1:nd)    = pgt(1:nd)*h6
            end do
        end do
    end do

!--- Calculate the combination of the factors and radial parts of Gogny force
    do k1 = 1, NBX

        lu1 = chunk(1)%l(k1);           ld1 = chunk(1)%lp(k1)
        do k2 = 1, NBX
            lu2 = chunk(1)%l(k2);           ld2 = chunk(1)%lp(k2)
        
        !--- Initialize
            GoV(k1, k2)%gg  = zero
            GoV(k1, k2)%gf  = zero
            GoV(k1, k2)%ff  = zero

            LI  = min( abs(lu1-lu2), abs(ld1-ld2) )
            LM  = max(     lu1+lu2,      ld1+ld2  )
            do L = LI, LM, 2
                do i = 1, 2
                    gf  = (two*L + one)*C3J(L,k1,k2)*( GA(i) + GD(i) )*half
                    gg  = gf - (two*L + one)*C3L(L,k1, k2)*GD(i)*half
                    ff  = gf - (two*L + one)*C3S(L,k1, k2)*GD(i)*half
                    
                    GoV(k1, k2)%gg  = VGL(L,i)%vlam*gg + GoV(k1, k2)%gg
                    GoV(k1, k2)%gf  = VGL(L,i)%vlam*gf + GoV(k1, k2)%gf
                    GoV(k1, k2)%ff  = VGL(L,i)%vlam*ff + GoV(k1, k2)%ff
                end do
            end do
        end do
    end do 
    
    return
                    
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine MediaG                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine PairI(it)                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!---- Calculate the interacting matrix element in pairing channel
!
    Integer, intent(in) :: it
    
    Integer :: i1, i2, i, nd, k1, k2, kp1, kp2
    double precision, dimension(MSD, MSD) :: gg, gf, ff
    double precision, dimension(MSD) :: fun
    double precision :: h, rel, Gamma, V0
    
    nd  = well%npt
    h   = well%h
    
!--- Pairing effects with zero-range pairing force   ! ref paper of Toki, 
    !--- no Density-dependence
    if(pair%is(it).eq.2) then
        Gamma   = YKT%Gamma
        V0      = YKT%V0/(8.0d0*pi)
    else if(pair%is(it).eq.3) then !--- with density-dependence
        Gamma   = GSNL%Gamma
        V0      = GSNL%V0/(8.0d0*pi)
    end if

    GME(it)%vg  = zero
    !$OMP PARALLEL DO PRIVATE(k1, kp1, i2, k2, kp2, i, fun, gg, gf, ff) SCHEDULE(STATIC)
    do i1 = 1, lev(it)%nt
        
        kp1 = lev(it)%nk(i1)
        if( kp1.lt.0 ) k1 = iabs(kp1);      if( kp1.gt.0 ) k1 = KMX(1) + kp1
        
        do i2 = i1, lev(it)%nt

            kp2 = lev(it)%nk(i2)
            if( kp2.lt.0 ) k2 = iabs(kp1);      if( kp2.gt.0 ) k2 = KMX(1) + kp2
            
            fun = zero
            if(pair%is(it).le.3) then
                fun =  (wav(i1,it)%xg*wav(i2,it)%xg + wav(i1,it)%xf*wav(i2,it)%xf)**2 
                fun = -V0*fun/well%xr**2
                
                if(pair%is(it).eq.3) fun = fun*(one - (dens(it)%rv/pset%rvs)**Gamma)

            else if(pair%is(it).eq.4) then
                do i = 1, nd
                    gg(:,i) = wav(i1,it)%xg(:)*wav(i1,it)%xg(i)*wav(i2,it)%xg(i)*wav(i2,it)%xg(:)
                    gf(:,i) = wav(i1,it)%xg(:)*wav(i1,it)%xf(i)*wav(i2,it)%xf(i)*wav(i2,it)%xg(:) + &
                            & wav(i1,it)%xf(:)*wav(i1,it)%xg(i)*wav(i2,it)%xg(i)*wav(i2,it)%xf(:)
                    ff(:,i) = wav(i1,it)%xf(:)*wav(i1,it)%xf(i)*wav(i2,it)%xf(i)*wav(i2,it)%xf(:)
                            
                    fun(:)  = fun(:) + GoV(k1,k2)%gg(:,i)*gg(:,i) + GoV(k1,k2)%ff(:,i)*ff(:,i) &
                            &        + GoV(k1,k2)%gf(:,i)*gf(:,i)
                end do
            end if
            Call simps(fun(2), nd-1, h, GME(it)%vg(i1,i2))
            GME(it)%vg(i2,i1)   =  GME(it)%vg(i1,i2)

        end do
    end do
    !$OMP END PARALLEL DO
    
    return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine PairI                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!*************************************************************************************************!
!                                                                                                 !
End Module Gogny                                                                                  !
!                                                                                                 !
!*************************************************************************************************!
