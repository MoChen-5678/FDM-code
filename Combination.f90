!*************************************************************************************************!
!
Module Combination
!
!*************************************************************************************************!
!
!--- In the module, we calculate the expansion of the propagators, the CG coefficient, 
!    as well as the combinations of these two for various channels
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
Use Define
use rhflib
use jsymbols
Use DBESIK
Implicit None
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!
!--- Number of meson and photon fiels
    Integer, parameter :: MPX = 5

!--- Multipole expansions for the propagators: Bessel functions 
    TYPE Bessel
        double precision, dimension(2*MSD-1) :: xjl, xhl
    END TYPE Bessel
    TYPE (Bessel), dimension(0:LYX+1, MPX) :: bes


!--- CG coefficient appearing in the Fock terms
    DOUBLE PRECISION, DIMENSION(0:LDX, NBX, NBX) :: C3J, C3L, C3S
    DOUBLE PRECISION, DIMENSION(0:LDX, 0:LDX, 2) :: TFJ 
    TYPE BabJL_T
        double precision, dimension(0:LDX, -1:1) :: XB
    END TYPE BabJL_T
    TYPE (BabJL_T), DIMENSION( -KPX:KPX, -KPX:KPX ) :: TBJ
        

!--- Radial part of the propagators in the sigma, omega, rho, and coulomb fields
    TYPE BINTS
        double precision, dimension(MSD, MSD) :: wsig, wome, wrho, wcou 
    END TYPE BINTS
    TYPE (BINTS), DIMENSION(0:LYX) :: XBIS

!--- Radial part of the space gradient of the propagators
    TYPE BINTT
        double precision, dimension(MSD, MSD) :: wrtn, wpio, wsig, wome 
    END TYPE BINTT
    TYPE (BINTT), DIMENSION(0:LYX+1, 4) :: XBIT

!--- Radial part of the propagator for the rho vector-tensor couplings
    TYPE BINTV
        double precision, dimension(MSD, MSD) :: wrtv, wrvt
    END TYPE BINTV
    TYPE (BINTV), DIMENSION(0:LYX+1, 2) :: XBIV

!--- Combination of the propagator and CG coefficients
!--- The sigma scalar coupling and the time component of the vector couplings
    TYPE NLMF_S0
        DOUBLE PRECISION, DIMENSION(MSD,MSD) :: XPP 
    END TYPE NLMF_S0
    TYPE (NLMF_S0), DIMENSION(NBX, NBX) :: CSIG
    TYPE (NLMF_S0), DIMENSION(NBX, NBX) :: COVT, CRVT , CAVT

!--- The space component of the vector couplings
    TYPE NLMF_V
        DOUBLE PRECISION, DIMENSION(MSD,MSD) :: XPP, XPM, XMM
    END TYPE NLMF_V
    TYPE (NLMF_V), DIMENSION(NBX, NBX) :: COVS, CRVS, CAVS 

!--- The pion-PV, rho-T, rho-VT, rho-TV couplings
    TYPE NLMF_T
        DOUBLE PRECISION, DIMENSION(MSD,MSD) :: XPP, XPM, XMP, XMM
    END TYPE NLMF_T
    TYPE (NLMF_T), DIMENSION(NBX, NBX) :: CPIO
    TYPE (NLMF_T), DIMENSION(NBX, NBX) :: CRTT, CRTS
    TYPE (NLMF_T), DIMENSION(NBX, NBX) :: CVTT, CVTS
    TYPE (NLMF_T), DIMENSION(NBX, NBX) :: CTVT, CTVS
    
!--- Contact terms of the rho-T coupling, including the compenstation terms
    DOUBLE PRECISION, DIMENSION(NBX, NBX) :: DRT1, DRT2, DRT3, DRT4
    
!--- Denote the index of the gradient term 
    INTEGER, DIMENSION(4) :: idi1, idi2
    data idi1/-1,-1, 1, 1/, idi2/-1, 1, -1, 1/
    
    INTEGER, DIMENSION(2) :: idi3
    data idi3/-1, 1/

!--- for interpolation
    double precision, dimension(5,6), private :: AI
    data AI(1, 1:6) / 0.24609375d0,   1.23046875d0,  -0.8203125d0,   0.4921875d0,  -0.17578125d0,   0.02734375d0/   !--- p=-1.50
    data AI(2, 1:6) /-0.02734375d0,   0.41015625d0,   0.8203125d0,  -0.2734375d0,   0.08203125d0,  -0.01171875d0/   !--- p=-0.50
    data AI(3, 1:6) / 0.01171875d0,  -0.09765625d0,   0.5859375d0,   0.5859375d0,  -0.09765625d0,   0.01171875d0/   !--- p=+0.50
    data AI(4, 1:6) /-0.01171875d0,   0.08203125d0,  -0.2734375d0,   0.8203125d0,   0.41015625d0,  -0.02734375d0/   !--- p=+1.50
    data AI(5, 1:6) / 0.02734375d0,  -0.17578125d0,   0.4921875d0,  -0.8203125d0,   1.23046875d0,   0.24609375d0/   !--- p=+2.50

Contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine Media                                                                                  !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Calculate the multiple expansions of the propagators with radial variables
!
!
    integer :: im, id, L, KODE, NZ 
    DOUBLE PRECISION :: X, Z, ALPHA, h2

    double precision, dimension(0:LYX+1) :: YI, YK
    double precision, dimension(MPX-1) :: amuf

!--- Change the unit of the meson masses
    amuf(1) = pset%amsig/hbc;       amuf(2) = pset%amome/hbc
    amuf(3) = pset%amrho/hbc;       amuf(4) = pset%ampio/hbc 

!--- Variables for calling bessel function
    ALPHA   = half;                 h2      = well%h*half;          KODE    = 1

!--- Rank of the bessel functions from 0 to KM-1.

!--- LOOP For meson fields
    do im = 1, MPX - 1
        do id = 2, 2*MSD-1
            Z   = (id - 1)*h2;             X   = amuf(im)*Z

        !--- Calculate I_{L+1/2}, K_{L+1/2}
        !--- SUBROUTINE DBESI (X, ALPHA, KODE, N, Y, NZ) --> I/SUB(ALPHA+K-1)/(X), K=1,...,N
            CALL DBESI(X, ALPHA, KODE, LYX+2, YI(0), NZ);       if(NZ.gt.0) Stop ' DBESI: NZ > 0'

        !--- SUBROUTINE DBESK (X, FNU, KODE, N, Y, NZ) --> K/SUB(FNU+I-1)/(X), I=1,...,N
            CALL DBESK(X, ALPHA, KODE, LYX+2, YK(0), NZ);       if(NZ.gt.0) Stop ' DBESK: NZ > 0'

            bes(:,im)%xjl(id)   = YI(:)/dsqrt(Z)
            bes(:,im)%xhl(id)   = YK(:)/dsqrt(Z)
        end do

    !--- For original point 
        Z   = well%xr(1);                                   X   = amuf(im)*Z
        CALL DBESI(X, ALPHA, KODE, LYX+2, YI, NZ);          if(NZ.gt.0) Stop ' DBESI: NZ > 0'
        CALL DBESK(X, ALPHA, KODE, LYX+2, YK, NZ);          if(NZ.gt.0) Stop ' DBESK: NZ > 0'

        bes(:,im)%xjl(1)    = YI(:)/dsqrt(Z)
        bes(:,im)%xhl(1)    = YK(:)/dsqrt(Z)
    end do

!--- Photon field (Coulomb field)
    im  = MPX
    do L = 0, LYX+1
        do id = 2, 2*MSD
            X                   = (id - 1)*h2
            bes(L,im)%xjl(id)   = X**L;                     bes(L,im)%xhl(id)   = one/X**(L+1)
        end do

    !--- For original point
        X   = well%xr(1)
        bes(L,im)%xjl(1)    = X**L;                     bes(L,im)%xhl(1)    = one/X**(L+1)
    end do
    
    return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine Media                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine Folding(fun, npt, pgt)                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Interpolation and Integration concerning the propagators
!
!    1. For given r, the integration over r' shall be separated into two parts: [0, r] and [r, r_c]
!       where the basic integration formula is three-point simpson rule. Regarding the spacing of 
!       the interval, one needs to provide mid value for each interval [i, i+1], namely i+0.5. Here
!       six-point interpolation is applied.
!
!    2. For the analytic propagator terms, the mid values of each interval are calculated precisely. 
!       While for the wave functions, one only has the values at the spacing points. By applying the
!       the interpolation, the mid values of the intervals are obtained for the wave functions. Namely 
!       the interpolation is introduced for wave functions, and the relavant weights are integrated 
!       with the propagators.
!
!-------------------------------------------------------------------------------------------------!
!
    double precision, dimension(2*MSD-1), intent(in) :: fun
    double precision, dimension(MSD), intent(out) :: pgt

    integer :: j, npt, k, id, ik

!--- INTERPOLATION POINTS
    pgt = zero
    do j = 1, 6
        pgt(j)  = zero
        do k = 1, 3
            pgt(j)  = fun(2*k)*AI(k, j) + pgt(j)
        end do
    end do
    
    do j = 2, 6
        do k = 1, j-1
            id      = 2*k + 6
            pgt(j)  = fun(id)*AI(3, j-k) + pgt(j)
        end do
    end do

    do j = 7, npt - 6
        pgt(j)  = zero
        do k = 1, 6
            id      = 2*(j - 4 + k)
            pgt(j)  = fun(id)*AI(3,k) + pgt(j)
        end do
    end do

    do j = npt -5, npt
        pgt(j)  = zero
        do k = 3, 5
            id      = 2*( npt - 6 + k )
            ik      = j - ( npt - 6 )
            pgt(j)  = fun(id)*AI(k,ik) + pgt(j)
        end do
    end do
    
    do j = npt-5, npt-1
        do k = 1, npt - j
            id      = 2*( npt - 3 - k )
            ik      = j - ( npt - 6 ) + k
            pgt(j)  = fun(id)*AI(3,ik) + pgt(j)
        end do
    end do

!--- EXTENDED SIMPSON'S RULE
    pgt(1)      = pgt(1)*4.0d0 + fun(1)
    do j = 2, npt-1
        pgt(j)  = pgt(j)*4.0d0 + fun(2*j-1)*2.0d0
    end do
    pgt(npt)    = pgt(npt)*4.0d0 + fun(2*npt-1)

    return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine Folding                                                                         !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
Subroutine PreMedia
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
    integer :: i, j, npt, L, i1, i2, L1, L2, idw
    double precision, dimension(2*MSD) :: fun
    double precision, dimension(MSD) :: r2, pgt

    npt     = well%npt
    Call Media

!--- Prepare propagator matrix for sigma, ome, rho, coulomb fields
!$OMP PARALLEL DO PRIVATE(L, i, pgt, fun) SCHEDULE(STATIC)
    do L = 0, LYX

        do i = 1, npt
        !--- BESSEL FUNCTIONS OF SIGMA FIELD
            fun(  1:2*i-1  )        = bes(L,1)%xjl(  1:2*i-1  )*bes(L,1)%xhl(2*i-1)
            fun(2*i:2*npt-1)        = bes(L,1)%xhl(2*i:2*npt-1)*bes(L,1)%xjl(2*i-1)

            Call Folding(fun, npt, pgt)

        !--- CONTRIBUTIONS FROM THE L'th COMPONENT
            XBIS(L)%wsig(i, 1:npt)  = pgt(1:npt)

        !--- BESSEL FUNCTIONS OF OMEGA FIELD
            fun(  1:2*i-1  )        = bes(L,2)%xjl(  1:2*i-1  )*bes(L,2)%xhl(2*i-1)
            fun(2*i:2*npt-1)        = bes(L,2)%xhl(2*i:2*npt-1)*bes(L,2)%xjl(2*i-1)

            Call Folding(fun, npt, pgt)

        !--- CONTRIBUTIONS FROM THE L'th COMPONENT
            XBIS(L)%wome(i, 1:npt)  = pgt(1:npt)

        !--- BESSEL FUNCTIONS OF RHO FIELD
            fun(  1:2*i-1  )        = bes(L,3)%xjl(  1:2*i-1  )*bes(L,3)%xhl(2*i-1)
            fun(2*i:2*npt-1)        = bes(L,3)%xhl(2*i:2*npt-1)*bes(L,3)%xjl(2*i-1)

            Call Folding(fun, npt, pgt)

        !--- CONTRIBUTIONS FROM THE L'th COMPONENT
            XBIS(L)%wrho(i, 1:npt)  = pgt(1:npt)

        !--- BESSEL FUNCTIONS OF COULOMB FIELD
            fun(  1:2*i-1  )        = bes(L,MPX)%xjl(  1:2*i-1  )*bes(L,MPX)%xhl(2*i-1)
            fun(2*i:2*npt-1)        = bes(L,MPX)%xhl(2*i:2*npt-1)*bes(L,MPX)%xjl(2*i-1)

            Call Folding(fun, npt, pgt)

        !--- CONTRIBUTIONS FROM THE L'th COMPONENT
            XBIS(L)%wcou(i, 1:npt)  = pgt(1:npt)
        end do
    end do
!$OMP END PARALLEL DO

!--- Pion pseudo-vector, Rho tensor and vector-tensor coupling
!$OMP PARALLEL DO PRIVATE(L, idw, i1, i2, L1, L2, i, pgt, fun) SCHEDULE(STATIC)
    do L = 0, LYX+1
        
        do idw = 1, 2
            i1  = idi3(idw);        L1  = L + i1;       if(L1.lt. 0 .or. L1.gt. LYX) cycle

        !--- Radial part of the propagator in the vector-tensor couplings
            do i = 1, npt
            !--- RVT: vector tensor couplings
                fun(  1:2*i-1  )            = -bes(L1,3)%xjl(  1:2*i-1  )*bes(L,3)%xhl(2*i-1)       !--- r > r'
                fun(2*i:2*npt-1)            =  bes(L1,3)%xhl(2*i:2*npt-1)*bes(L,3)%xjl(2*i-1)       !--- r < r'

                fun(2*i-1)                  = -half*(bes(L1,3)%xjl(2*i-1)*bes(L,3)%xhl(2*i-1) - &   !--- r = r'
                                            &        bes(L1,3)%xhl(2*i-1)*bes(L,3)%xjl(2*i-1))

                Call Folding(fun, npt, pgt)

                XBIV(L, idw)%wrvt(i, 1:npt) =  pgt(1:npt)           !--- S_{\lambda\lambda_1}

            !--- RTV: tensor-vector couplings
                fun(  1:2*i-1  )            = -bes(L,3)%xjl(  1:2*i-1  )*bes(L1,3)%xhl(2*i-1)       !--- r > r'
                fun(2*i:2*npt-1)            =  bes(L,3)%xhl(2*i:2*npt-1)*bes(L1,3)%xjl(2*i-1)       !--- r < r'

                fun(2*i-1)                  = -half*(bes(L,3)%xjl(2*i-1)*bes(L1,3)%xhl(2*i-1) - &   !--- r = r'
                                            &        bes(L,3)%xhl(2*i-1)*bes(L1,3)%xjl(2*i-1))

                Call Folding(fun, npt, pgt)

                XBIV(L, idw)%wrtv(i, 1:npt) =  pgt(1:npt)           !--- S_{\lambda_1\lambda}
            end do
        end do
        
    !--- Radial part of the propagator in pseudo-vector and tensor couplings
        do idw = 1, 4
            i1  = idi1(idw);        L1  = L + i1;       if(L1.lt. 0 .or. L1.gt. LYX) cycle
            i2  = idi2(idw);        L2  = L + i2;       if(L2.lt. 0 .or. L2.gt. LYX) cycle
            do i = 1, npt

            !--- Rho-Tensor coupling
                fun(  1:2*i-1  )            =  bes(L2,3)%xjl(  1:2*i-1  )*bes(L1,3)%xhl(2*i-1)      !--- r > r'
                fun(2*i:2*npt-1)            =  bes(L2,3)%xhl(2*i:2*npt-1)*bes(L1,3)%xjl(2*i-1)      !--- r < r'

                fun(2*i-1)                  =  half*(bes(L2,3)%xjl(2*i-1)*bes(L1,3)%xhl(2*i-1) + &  !--- r = r'
                                            &        bes(L2,3)%xhl(2*i-1)*bes(L1,3)%xjl(2*i-1))

                Call Folding(fun, npt, pgt)

                XBIT(L, idw)%wrtn(i, 1:npt) = pgt(1:npt)        !--- R_{\lambda_1\lambda_2}

            !--- Pion-pseudo-vector coupling
                fun(  1:2*i-1  )            =  bes(L2,4)%xjl(  1:2*i-1  )*bes(L1,4)%xhl(2*i-1)
                fun(2*i:2*npt-1)            =  bes(L2,4)%xhl(2*i:2*npt-1)*bes(L1,4)%xjl(2*i-1)

                fun(2*i-1)                  =  half*(bes(L2,4)%xjl(2*i-1)*bes(L1,4)%xhl(2*i-1) + &
                                            &        bes(L2,4)%xhl(2*i-1)*bes(L1,4)%xjl(2*i-1))

                Call Folding(fun, npt, pgt)

                XBIT(L, idw)%wpio(i, 1:npt) = pgt(1:npt)
                    
            !--- Tensor force component in sigma coupling
                fun(  1:2*i-1  )            =  bes(L2,1)%xjl(  1:2*i-1  )*bes(L1,1)%xhl(2*i-1)
                fun(2*i:2*npt-1)            =  bes(L2,1)%xhl(2*i:2*npt-1)*bes(L1,1)%xjl(2*i-1)

                fun(2*i-1)                  =  half*(bes(L2,1)%xjl(2*i-1)*bes(L1,1)%xhl(2*i-1) + &
                                            &        bes(L2,1)%xhl(2*i-1)*bes(L1,1)%xjl(2*i-1))

                Call Folding(fun, npt, pgt)

                XBIT(L, idw)%wsig(i, 1:npt) = pgt(1:npt)
                    
            !--- Tensor force component in omega coupling
                fun(  1:2*i-1  )            =  bes(L2,2)%xjl(  1:2*i-1  )*bes(L1,2)%xhl(2*i-1)
                fun(2*i:2*npt-1)            =  bes(L2,2)%xhl(2*i:2*npt-1)*bes(L1,2)%xjl(2*i-1)

                fun(2*i-1)                  =  half*(bes(L2,2)%xjl(2*i-1)*bes(L1,2)%xhl(2*i-1) + &
                                            &        bes(L2,2)%xhl(2*i-1)*bes(L1,2)%xjl(2*i-1))

                Call Folding(fun, npt, pgt)

                XBIT(L, idw)%wome(i, 1:npt) = pgt(1:npt)
            end do
        end do

!--- End loop over L
    end do
!$OMP END PARALLEL DO

    return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine PreMedia                                                                           !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Function BJLab(kpa, kpb, J, L)                                                                      !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Calcualate the symbols B appearing in the space component of the rho-T and rho-VT couplings
!
!-------------------------------------------------------------------------------------------------!
!
    integer, intent(in) :: kpa, kpb, J, L
    double precision :: BJLab

    integer :: k1, k2, la, lb
    double precision :: hJ, tmp, tmr, xks, tmpj

!--- Calculate sqrt(2J + 1)*BJLab
    k1  = abs(kpa);         k2      = abs(kpb)  
    
    BJLab   = zero
    
    if(J.eq.0 .and. L.eq.0) return

!--- Determine the orbital angular momentum associated with kpb
    la  = k1;               if(kpa.lt.0) la = la - 1
    
    Select CASE ( J.eq.L )
    Case (.true.) 
        hJ  = two*J + one;          hJ  = dsqrt(hJ)
        tmp = J*(J + one);          tmp = dsqrt(tmp)
        
        BJLab   = hJ*( k1 + k2*(-1.0)**(k1+k2+J-1) )/tmp
    
    Case ( .false. )
        tmp = (L-J)*(L+J+1)*half
        tmr = dabs(tmp);            tmr = dsqrt(tmr)
        
        BJLab   = (-1.0)**(k1 + la)*( kpa + kpb + tmp )/tmr
    End Select

    return

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
End function BJLab
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
Function FLJLp(L, J, L1)
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Calculate the symbol F_{LJ}^{L1}
!
!-------------------------------------------------------------------------------------------------!
!
    Integer, intent(in) :: L, J, L1
    Double precision :: FLJLp, hL, hL1
    
    FLJLp   = zero
    If ( mod(L+L1+1,2).ne.0 ) return
    
    hL  = two*L + one;          hL  = dsqrt(hL)
    hL1 = two*L1 + one;         hL1 = dsqrt(hL1)

    FLJLp   = hL*hL1*wiglll(L, 1, L1)*racah0(L1,L,1, 1,1,J)
    
    return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
End Function FLJLp
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine MediaX                                                                                 !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Calculate the CG coefficient appearing in the Fock terms
!
!-------------------------------------------------------------------------------------------------!
!
    Integer :: ka, kb, k1, lu1, ld1, k2, lu2, ld2, L, J, kp1, kp2, L1, i1
    Double precision :: X, Z

!--- Determine the kappa quantities
    do ka = 1, NBX

        select case ( ka.le.KMX(1) )
    !--- For the case j = l+1/2
        case (.true.)
            KLP(ka)%kp  = -ka
            KLP(ka)%lu  = ka - 1;       KLP(ka)%ld  = ka

    !--- For the case j = l-1/2
        case (.false.)
            KLP(ka)%kp  = ka - KMX(1)
            KLP(ka)%lu  = KLP(ka)%kp;   KLP(ka)%ld  = KLP(ka)%kp - 1
        end select

    !--- absolute value of kappa and parity (0 for even, 1 for odd)
        KLP(ka)%ka  = iabs(KLP(ka)%kp)
        KLP(ka)%ip  = ( 1 - (-1)**KLP(ka)%lu )/2
    end do

!--- Calculate the squares of the C-G coefficients
!--- Prepreparation
    Call gfv
    
!--- loop over ka and kb
    C3J = zero;     C3L = zero;     C3S = zero
    do ka = 1, NBX
        k1  = KLP(ka)%ka;       lu1     = KLP(ka)%lu;       ld1     = KLP(ka)%ld
        
        do kb = 1, NBX
            k2  = KLP(kb)%ka;       lu2     = KLP(kb)%lu;       ld2     = KLP(kb)%ld

            do L = 0, k1+k2
            !
            !--- In calculationg the CG coefficients, FUNCTION wiglll(l1,l2,l3) and 
            !    FUNCTION racslj(k,l1,l2,j2,j1) are called.
            !--- For the obtained results, the coefficients (2L+1) are not included.
            !
                Select Case ( mod(lu1+lu2+L,2).eq.0 )
                    
                Case (.true.)
                !--- For the CG coefficients: |ja 1/2, jb -1/2; L0>
                    X   = wiglll(lu1, lu2, L)*racslj(L,lu1, lu2, k2, k1)
                    C3J(L, ka, kb)  = X**2*(two*lu1 + 1)*(two*lu2 + 1)

                !--- For the CG coefficients: |lu1 0, lu2 0; L0>
                    C3L(L, ka, kb)  = wiglll(lu1,lu2,L)**2

                !--- For the CG coefficients: |ld1 0, ld2 0; L0>
                    C3S(L, ka, kb)  = wiglll(ld1,ld2,L)**2

                Case (.false.) 
                !--- For the CG coefficients: |ja 1/2, jb -1/2; L0>
                    Z   = wiglll(ld1,lu2, L)*racslj(L,ld1,lu2, k2, k1)
                    C3J(L, ka, kb)  = Z**2*(two*ld1+ 1)*(two*lu2 + 1)

                !--- For the CG coefficients: |ld1 0, lu2 0; L0>
                    C3L(L, ka, kb)  = wiglll(ld1,lu2,L)**2

                !--- For the CG coefficients: |lu1 0, ld2 0; L0>
                    C3S(L, ka, kb)  = wiglll(lu1,ld2,L)**2
                End Select

            end do
        end do
    end do

!--- Calculate BabJL
    do kp1 = -KPX, KPX
        
        k1  = iabs(kp1)
        do kp2 = -KPX, KPX
            TBJ(kp1,kp2)%XB = zero

            if( kp1*kp2.eq.0 ) cycle
            
            k2  = iabs(kp2)
            do J = iabs(k1-k2), k1 + k2 - 1
                do L = iabs(J-1), J+1
                    TBJ(kp1,kp2)%XB(J,L-J)  = BJLab(kp1, kp2, J, L)
                end do
            end do

        end do
    end do
    
!--- Calculate the FLJL1
    TFJ = zero
    do J = 0, LDX
        do L = 0, LYX
            if ( iabs(J-L).gt.1 ) cycle
            
            do i1 = 1, 2
                L1  = L + idi3(i1)
                if( L1.lt.0 ) cycle
                if( iabs(J-L1).gt.1 ) cycle

                TFJ(L,J,i1) = FLJLp( L, J, L1)
            end do
        end do
    end do

    return

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine MediaX                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine PrePotels                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Calcuate the two body interactions: XMIS, XMIV, XMIT, which will be stored during the iteration
!
!

    double precision :: h, h6, xka, xkb 
    double precision :: xmrho, xmpio, hLS, tmp, APP, APM, AMP, AMM

    double precision, dimension(-1:1) :: alp
    
    integer :: ka, kb, k1, k2, kp1, kp2, lu1, lu2, ld1, ld2, idw, id1, id2
    integer :: LI, LM, LIX, LMX, JMI, JMX, L, L1, L2, i1, i2, J
    integer :: i, npt
    double precision, dimension(MSD) :: r2, xdel

!--- Calculate the CG coefficients
    Call MediaX
    
!--- Preparation
    npt     = well%npt;             h       = well%h;               h6      = well%h/6.0d0
    xmrho   = pset%amrho/hbc;       xmpio   = pset%ampio/hbc;       r2      = well%xr**2

!--- Determine the quantum numbers
!$OMP PARALLEL PRIVATE(ka,xka,kp1,k1,lu1,ld1, kb,xkb,kp2,k2,lu2,ld2, APP,APM,AMP,AMM, &
!$OMP    & J,L,hLS,LI,LM,LIX,LMX,JMI,JMX, i,i1,i2,L1,L2, alp,tmp,xdel, idw,id1,id2 )
!$OMP DO SCHEDULE(STATIC)
    do ka = 1, NBX
        k1  = KLP(ka)%ka;       kp1 = KLP(ka)%kp;       xka = kp1
        lu1 = KLP(ka)%lu;       ld1 = KLP(ka)%ld

        do kb = 1, NBX
            k2  = KLP(kb)%ka;       kp2 = KLP(kb)%kp;       xkb = kp2
            lu2 = KLP(kb)%lu;       ld2 = KLP(kb)%ld

            JMI = iabs(k1-k2);      JMX = k1 + k2 - 1

        !--- sigma-scalar, time-component vector couplings: lu1 + lu2 + L to be even
            CSIG(ka,kb)%XPP = zero;     COVT(ka,kb)%XPP = zero
            CRVT(ka,kb)%XPP = zero;     CAVT(ka,kb)%XPP = zero

            LI  = min( iabs(lu1-lu2), iabs(ld1-ld2) )
            LM  = max( iabs(lu1+lu2), iabs(ld1+ld2) )
            do L = LI, LM, 2
                hLS = ( two*L + one )*h6/(4.0d0*pi)

                CSIG(ka,kb)%XPP = CSIG(ka,kb)%XPP + XBIS(L)%wsig*C3J(L,ka,kb)*hLS
                COVT(ka,kb)%XPP = COVT(ka,kb)%XPP - XBIS(L)%wome*C3J(L,ka,kb)*hLS
                CRVT(ka,kb)%XPP = CRVT(ka,kb)%XPP - XBIS(L)%wrho*C3J(L,ka,kb)*hLS
                                
                hLS = h6/alphi
                CAVT(ka,kb)%XPP = CAVT(ka,kb)%XPP - XBIS(L)%wcou*C3J(L,ka,kb)*hLS
            end do

        !--- Space part for the vector fields: omega (OMEx), rho (RHOx) and Coulomb (COUx) fields (x=1, 2,3)
        !--- L + lu1 + ld2 to be even
            COVS(ka,kb)%XPP = zero;     CRVS(ka,kb)%XPP = zero;     CAVS(ka,kb)%XPP = zero
            COVS(ka,kb)%XPM = zero;     CRVS(ka,kb)%XPM = zero;     CAVS(ka,kb)%XPM = zero
            COVS(ka,kb)%XMM = zero;     CRVS(ka,kb)%XMM = zero;     CAVS(ka,kb)%XMM = zero
            
            LI  = min( iabs(lu1-ld2), iabs(ld1-lu2) )
            LM  = max( iabs(lu1+ld2), iabs(ld1+lu2) )
            do L = LI, LM, 2
                
                hLS = h6*(two*L + one)/(4.0d0*pi)
                
                APP = C3S(L,ka,kb)*two - C3J(L,ka,kb)
                APM = C3J(L,ka,kb)
                AMM = C3L(L,ka,kb)*two - C3J(L,ka,kb)
                
                COVS(ka,kb)%XPP = COVS(ka,kb)%XPP + XBIS(L)%wome*APP*hLS
                COVS(ka,kb)%XPM = COVS(ka,kb)%XPM + XBIS(L)%wome*APM*hLS
                COVS(ka,kb)%XMM = COVS(ka,kb)%XMM + XBIS(L)%wome*AMM*hLS

                CRVS(ka,kb)%XPP = CRVS(ka,kb)%XPP + XBIS(L)%wrho*APP*hLS
                CRVS(ka,kb)%XPM = CRVS(ka,kb)%XPM + XBIS(L)%wrho*APM*hLS
                CRVS(ka,kb)%XMM = CRVS(ka,kb)%XMM + XBIS(L)%wrho*AMM*hLS
                
                hLS = h6/alphi
                CAVS(ka,kb)%XPP = CAVS(ka,kb)%XPP + XBIS(L)%wcou*APP*hLS
                CAVS(ka,kb)%XPM = CAVS(ka,kb)%XPM + XBIS(L)%wcou*APM*hLS
                CAVS(ka,kb)%XMM = CAVS(ka,kb)%XMM + XBIS(L)%wcou*AMM*hLS
            end do

        !--- Pion (pseudo-vector couplings)
            if(pset%fpio.eq.zero) go to 20
            CPIO(ka,kb)%XPP = zero;             CPIO(ka,kb)%XPM = zero
            CPIO(ka,kb)%XMP = zero;             CPIO(ka,kb)%XMM = zero

            LI  = min( iabs(lu1-lu2), iabs(ld1-ld2) )
            LM  = max( iabs(lu1+lu2), iabs(ld1+ld2) )
            Do L = LI-1, LM+1, 2
                alp(-1) = -L;           alp(1)  = L + one
                
                hLS = h6/(two*L + one)/(4.0d0*pi)*C3J(L,ka,kb)
                do idw = 1, 4
                    i1  = idi1(idw);        L1  = L + i1;       if( L1.lt.LI .or. L1.gt.LM ) cycle
                    i2  = idi2(idw);        L2  = L + i2;       if( L2.lt.LI .or. L2.gt.LM ) cycle
                    
                    APP = -( xka+xkb + alp(i1) )*( xka+xkb + alp(i2) )*i1*i2*hLS
                    APM = +( xka+xkb + alp(i1) )*( xka+xkb - alp(i2) )*i1*i2*hLS
                    AMP = +( xka+xkb - alp(i1) )*( xka+xkb + alp(i2) )*i1*i2*hLS
                    AMM = -( xka+xkb - alp(i1) )*( xka+xkb - alp(i2) )*i1*i2*hLS
                    
                    CPIO(ka,kb)%XPP = CPIO(ka,kb)%XPP + XBIT(L,idw)%wpio*APP
                    CPIO(ka,kb)%XPM = CPIO(ka,kb)%XPM + XBIT(L,idw)%wpio*APM
                    CPIO(ka,kb)%XMP = CPIO(ka,kb)%XMP + XBIT(L,idw)%wpio*AMP
                    CPIO(ka,kb)%XMM = CPIO(ka,kb)%XMM + XBIT(L,idw)%wpio*AMM

                end do
            end do

        !--- Remove the contact term in pion pseudo-vector
            xdel    = third*two/(xmpio**2*r2*4.0*pi)
            do i = 1, npt
                CPIO(ka,kb)%XPM(i,i)    = CPIO(ka,kb)%XPM(i,i) + xdel(i)
                CPIO(ka,kb)%XMP(i,i)    = CPIO(ka,kb)%XMP(i,i) + xdel(i)
            end do

!--- Calculate the 0 component of rho-Tensor coupling (T1) and (VT1)
        20  if(pset%grtn.eq.zero) cycle
            if(pset%txtfor.eq.'PKA1') then
                LI  = abs(lu1-lu2);                   LM  = lu1 + lu2
                LMX = min(ld1+lu2, lu1+ld2);          LIX = max(abs(ld1-lu2), abs(lu1-ld2))
            else
                LI  = min( iabs(lu1-lu2), iabs(ld1-ld2) );  LM  = max( lu1+lu2, ld1+ld2 )
                LIX = min( iabs(ld1-lu2), iabs(lu1-ld2) );  LMX = max( ld1+lu2, lu1+ld2 )
            end if

            CTVT(ka,kb)%XPP = zero;     CVTT(ka,kb)%XPP = zero;     CRTT(ka,kb)%XPP = zero
            CTVT(ka,kb)%XPM = zero;     CVTT(ka,kb)%XPM = zero;     CRTT(ka,kb)%XPM = zero
            CTVT(ka,kb)%XMP = zero;     CVTT(ka,kb)%XMP = zero;     CRTT(ka,kb)%XMP = zero
            CTVT(ka,kb)%XMM = zero;     CVTT(ka,kb)%XMM = zero;     CRTT(ka,kb)%XMM = zero
            do L = LI, LM, 2
                alp(-1) = -L;               alp(1)  = L + one

                hLS = h6*C3J(L,ka,kb)/(4.0d0*pi)*xmrho                
                do idw = 1, 2
                    i1  = idi3(idw);        L1  = L + i1;       if( L1.lt.LIX .or. L1.gt.LMX ) cycle

                    APP = (xka - xkb + alp(i1))*i1*hLS
                    AMM = (xka - xkb - alp(i1))*i1*hLS
                    CVTT(ka,kb)%XPP = CVTT(ka,kb)%XPP - XBIV(L,idw)%wrvt*APP      !--- G gf G
                    CVTT(ka,kb)%XPM = CVTT(ka,kb)%XPM + XBIV(L,idw)%wrvt*AMM      !--- G gg F
                    CVTT(ka,kb)%XMP = CVTT(ka,kb)%XMP - XBIV(L,idw)%wrvt*APP      !--- F ff G
                    CVTT(ka,kb)%XMM = CVTT(ka,kb)%XMM + XBIV(L,idw)%wrvt*AMM      !--- F fg F

                    CTVT(ka,kb)%XPP = CTVT(ka,kb)%XPP + XBIV(L,idw)%wrtv*APP      !--- G fg G
                    CTVT(ka,kb)%XPM = CTVT(ka,kb)%XPM + XBIV(L,idw)%wrtv*APP      !--- G ff F
                    CTVT(ka,kb)%XMP = CTVT(ka,kb)%XMP - XBIV(L,idw)%wrtv*AMM      !--- F gg G
                    CTVT(ka,kb)%XMM = CTVT(ka,kb)%XMM - XBIV(L,idw)%wrtv*AMM      !--- F gf F
                end do
            end do

            do L = LI, LM, 2
                alp(-1) = -L;               alp(1)  = L + one
                hLS = h6*C3J(L,ka,kb)/(two*L + one)/(4.0d0*pi)*xmrho**2
                do idw = 1, 4
                    i1  = idi1(idw);    L1  = L + i1;       if(L1.lt. LIX .or. L1.gt.LMX) cycle
                    i2  = idi2(idw);    L2  = L + i2;       if(L2.lt. LIX .or. L2.gt.LMX) cycle
                    
                    APP = (xka - xkb + alp(i1))*(xka - xkb + alp(i2))*i1*i2*hLS
                    APM = (xka - xkb + alp(i1))*(xka - xkb - alp(i2))*i1*i2*hLS
                    AMP = (xka - xkb - alp(i1))*(xka - xkb + alp(i2))*i1*i2*hLS
                    AMM = (xka - xkb - alp(i1))*(xka - xkb - alp(i2))*i1*i2*hLS
                    
                    CRTT(ka,kb)%XPP = CRTT(ka,kb)%XPP + XBIT(L,idw)%wrtn*APP    !--- G ff G
                    CRTT(ka,kb)%XPM = CRTT(ka,kb)%XPM - XBIT(L,idw)%wrtn*APM    !--- G fg F
                    CRTT(ka,kb)%XMP = CRTT(ka,kb)%XMP - XBIT(L,idw)%wrtn*AMP    !--- F gf G
                    CRTT(ka,kb)%XMM = CRTT(ka,kb)%XMM + XBIT(L,idw)%wrtn*AMM    !--- F gg F
                end do
            end do

     !--- Rho tensor couplings (T2) and (VT2)
            LI  = min( iabs(lu1-ld2), iabs(ld1-lu2) )
            LM  = max( iabs(lu1+ld2), iabs(ld1+lu2) )
            LIX = min( iabs(lu1-lu2), iabs(ld1-ld2) )
            LMX = max( iabs(lu1+lu2), iabs(ld1+ld2) )

            CTVS(ka,kb)%XPP = zero;     CVTS(ka,kb)%XPP = zero;     CRTS(ka,kb)%XPP = zero;     DRT1(ka,kb) = zero
            CTVS(ka,kb)%XPM = zero;     CVTS(ka,kb)%XPM = zero;     CRTS(ka,kb)%XPM = zero;     DRT2(ka,kb) = zero
            CTVS(ka,kb)%XMP = zero;     CVTS(ka,kb)%XMP = zero;     CRTS(ka,kb)%XMP = zero;     DRT3(ka,kb) = zero
            CTVS(ka,kb)%XMM = zero;     CVTS(ka,kb)%XMM = zero;     CRTS(ka,kb)%XMM = zero;     DRT4(ka,kb) = zero
            do L = LI, LM, 2

                hLS = h6/(4.0d0*pi)*dsqrt(6.0d0)*xmrho
               
                do idw = 1, 2
                    i1  = idi3(idw);    L1  = L + i1;   if(L1.lt. LIX .or. L1 .gt. LMX) cycle
                    
                    APP = zero;     APM = zero;     AMP = zero;     AMM = zero
                    do J = abs(L-1), L+1
                        if( J.lt.JMI .or. J.gt.JMX .or. J.lt. abs(L1-1) .or. J.gt. L1+1 ) cycle

                        tmp = TFJ(L,J,idw)*C3J(J,ka,kb)*hLS*(-1)**J
                        APP = APP + TBJ( kp1,-kp2)%XB(J,L-J)*TBJ( kp1, kp2)%XB(J,L1-J)*tmp ! +-++
                        APM = APM + TBJ( kp1,-kp2)%XB(J,L-J)*TBJ(-kp1,-kp2)%XB(J,L1-J)*tmp ! +---
                        AMP = AMP + TBJ(-kp1, kp2)%XB(J,L-J)*TBJ( kp1, kp2)%XB(J,L1-J)*tmp ! -+++
                        AMM = AMM + TBJ(-kp1, kp2)%XB(J,L-J)*TBJ(-kp1,-kp2)%XB(J,L1-J)*tmp ! -+--
                    end do

                    CVTS(ka,kb)%XPP = CVTS(ka,kb)%XPP - XBIV(L,idw)%wrvt*APP      !--- G fg G
                    CVTS(ka,kb)%XPM = CVTS(ka,kb)%XPM - XBIV(L,idw)%wrvt*APM      !--- G ff F
                    CVTS(ka,kb)%XMP = CVTS(ka,kb)%XMP - XBIV(L,idw)%wrvt*AMP      !--- F gg G
                    CVTS(ka,kb)%XMM = CVTS(ka,kb)%XMM - XBIV(L,idw)%wrvt*AMM      !--- F gf F
    
                    CTVS(ka,kb)%XPP = CTVS(ka,kb)%XPP + XBIV(L,idw)%wrtv*APP      !--- G gf G
                    CTVS(ka,kb)%XPM = CTVS(ka,kb)%XPM + XBIV(L,idw)%wrtv*AMP      !--- G gg F
                    CTVS(ka,kb)%XMP = CTVS(ka,kb)%XMP + XBIV(L,idw)%wrtv*APM      !--- F ff G
                    CTVS(ka,kb)%XMM = CTVS(ka,kb)%XMM + XBIV(L,idw)%wrtv*AMM      !--- F fg F
                end do
                
                hLS = h6/(4.0d0*pi)*6.0d0*xmrho**2
                do idw = 1, 4
                    i1  = idi1(idw);    L1  = L + i1;   if(L1.lt. LIX .or. L1 .gt. LMX) cycle
                    i2  = idi2(idw);    L2  = L + i2;   if(L2.lt. LIX .or. L2 .gt. LMX) cycle

                    id1 = 1 + (i1+1)/2;                 id2 = 1 + (i2+1)/2
                    APP = zero;     APM = zero;     AMP = zero;     AMM = zero
                    do J = abs(L-1), L+1
                        if( J.lt.JMI .or. J.gt.JMX ) cycle
                        if(J.lt. abs(L1-1) .or. J.gt. L1+1 .or. J.lt. abs(L2-1) .or. J.gt. L2+1) cycle
                        
                        tmp = TFJ(L,J,id1)*TFJ(L,J,id2)*C3J(J,ka,kb)*hLS
                        APP = APP + TBJ( kp1, kp2)%XB(J,L1-J)*TBJ( kp1, kp2)%XB(J,L2-J)*tmp
                        AMP = AMP + TBJ( kp1, kp2)%XB(J,L1-J)*TBJ(-kp1,-kp2)%XB(J,L2-J)*tmp
                        APM = APM + TBJ(-kp1,-kp2)%XB(J,L1-J)*TBJ( kp1, kp2)%XB(J,L2-J)*tmp
                        AMM = AMM + TBJ(-kp1,-kp2)%XB(J,L1-J)*TBJ(-kp1,-kp2)%XB(J,L2-J)*tmp
                    end do
                    
                    CRTS(ka,kb)%XPP = CRTS(ka,kb)%XPP - XBIT(L,idw)%wrtn*APP      !--- G gg G
                    CRTS(ka,kb)%XPM = CRTS(ka,kb)%XPM - XBIT(L,idw)%wrtn*AMP      !--- G gf F
                    CRTS(ka,kb)%XMP = CRTS(ka,kb)%XMP - XBIT(L,idw)%wrtn*APM      !--- F fg G
                    CRTS(ka,kb)%XMM = CRTS(ka,kb)%XMM - XBIT(L,idw)%wrtn*AMM      !--- F ff F
                    
                    DRT1(ka,kb) = DRT1(ka,kb) + APP/h6
                    DRT2(ka,kb) = DRT2(ka,kb) + AMP/h6
                    DRT3(ka,kb) = DRT3(ka,kb) + APM/h6
                    DRT4(ka,kb) = DRT4(ka,kb) + AMM/h6
                end do
            end do

            xdel    = one/(r2*4.0d0*pi)
            do i = 1, npt
                CRTS(ka,kb)%XPP(i, i)   = CRTS(ka,kb)%XPP(i, i) - xdel(i)
                CRTS(ka,kb)%XPM(i, i)   = CRTS(ka,kb)%XPM(i, i) - xdel(i)
                CRTS(ka,kb)%XMP(i, i)   = CRTS(ka,kb)%XMP(i, i) - xdel(i)
                CRTS(ka,kb)%XMM(i, i)   = CRTS(ka,kb)%XMM(i, i) - xdel(i)
            end do

            xdel    = one/(r2*xmrho**2)
            do i = 1, npt
                CRTS(ka,kb)%XPP(i, i)   = CRTS(ka,kb)%XPP(i, i) + xdel(i)*DRT1(ka,kb)
                CRTS(ka,kb)%XPM(i, i)   = CRTS(ka,kb)%XPM(i, i) + xdel(i)*DRT2(ka,kb)
                CRTS(ka,kb)%XMP(i, i)   = CRTS(ka,kb)%XMP(i, i) + xdel(i)*DRT3(ka,kb)
                CRTS(ka,kb)%XMM(i, i)   = CRTS(ka,kb)%XMM(i, i) + xdel(i)*DRT4(ka,kb)
            end do

        end do  ! kb
    end do   ! ka
!$OMP END DO
!$OMP END PARALLEL
    return

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine PrePotels                                                                          !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!*************************************************************************************************!
!
End Module Combination
!
!*************************************************************************************************!
