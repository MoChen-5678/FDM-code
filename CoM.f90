!*************************************************************************************************!
!                                                                                                 !
Module Center                                                                                     !
!                                                                                                 !
!*************************************************************************************************!
!
!--- Calulate the center-of-mass correction on binding energy and radius
!
use Define
use rhflib
Use Combination, only: C3J
implicit none

    double precision, dimension(0:IBX) :: PCM, RCM, ECM

Contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine CoM(is, lpr)                                                                           !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!--- Center of mass correction
!--- Projection-after-variation in first-order approximation
!--- see Eur. Phys. J. A 7,467-478(2000) while corrected

    logical, intent(in) :: lpr
    integer, intent(in) :: is

    integer :: i1, i2, kap, ka, lua, lda, kbp, kb, lub, ldb, it
    double precision :: hom, tem, cmr, vv, uv, APA, AAP, ALP(-1:1), a1, a2, hja, hjb
    double precision, dimension(MSD) :: fun, fup, df2, dg2, dga, dgb, dfa, dfb
     
    if(lpr) write(l6,*) '********************** Begin Center of Mass correction ******************'
    if(is.eq.0) then
        hom     = 41.0d0*nuc%amas**(-third)
        ECM(0)  = -0.75d0*hom
        return
    else if(is.eq.1) then ! Use Microscopic methods
    !--- Coefficient and unit
        hom =  two*(pset%amu(1)*nuc%npr(1) + pset%amu(2)*nuc%npr(2))/hbc
        hom = -one/hom*hbc
        
    !--- Initialize
        PCM = zero;             RCM = zero;             ECM = zero
    
    !--- Loop over isospin
        do it = 1, IBX
            
        !--- Loop over single-particle states
            do i1 = 1, lev(it)%nt
                
            !--- kappa quantity and angular momentum for given states
                kap = lev(it)%nk(i1);                       hja = two*abs(kap)
                ka  = iabs(kap);                            if(kap.gt.0) ka = KMX(1) + kap
                lua = lev(it)%lu(i1);                       lda = lev(it)%ld(i1)

            !--- Probabilities scaled with the maximum occupation number
                vv  = lev(it)%vv(i1)*lev(it)%mu(i1)

            !--- Calculate the second derivative term: 
                Call deriv2(wav(i1,it)%xf(2), df2(2), well%npt-1, well%h)
                Call deriv2(wav(i1,it)%xg(2), dg2(2), well%npt-1, well%h)

                fun = ( wav(i1,it)%xf*df2 - lda*(lda + 1.)*wav(i1,it)%xf**2/well%xr**2 + &
                    &   wav(i1,it)%xg*dg2 - lua*(lua + 1.)*wav(i1,it)%xg**2/well%xr**2 )*vv 

                call simps(fun(2), well%npt-1, well%h, tem)
                ECM(it) = ECM(it) - tem*hom
                PCM(it) = PCM(it) - tem

            !--- Calculate the second part: squre of the first derivate terms
                tem = zero;                                 cmr = zero

            !--- Radial derivative of the radial wave function
                Call deriv( wav(i1,it)%xg(2), dga(2), well%npt-1, well%h )
                Call deriv( wav(i1,it)%xf(2), dfa(2), well%npt-1, well%h )

            !--- Loop over alpha'
                do i2 = 1, lev(it)%nt

                !--- kappa quantity and angular momentum for state alpha'
                    kbp = lev(it)%nk(i2);                       hjb = two*abs(kbp)
                    kb  = iabs(kbp);                            if(kbp.gt.0) kb = KMX(1) + kbp
                    lub = lev(it)%lu(i2);                       ldb = lev(it)%ld(i2)

                !--- angular momentum should differ with 1
                    if( abs(abs(kap) - abs(kbp)).gt.2 ) cycle
                    if( abs(lua - lub).ne.1 .and. abs(lda - ldb).ne.1 ) cycle

                !--- Coefficients for different terms in P square
                    vv  = lev(it)%vv(i1)*lev(it)%vv(i2)
                    uv  = vv*( one - lev(it)%vv(i1) )*( one - lev(it)%vv(i2) )
                    uv  = dsqrt(uv)

                !--- Radial derivative of the radial wave function
                    Call deriv( wav(i2,it)%xg(2), dgb(2), well%npt-1, well%h )
                    Call deriv( wav(i2,it)%xf(2), dfb(2), well%npt-1, well%h )
                            
                !--- Calculate the AAP & APA
                    fun = zero;                     fup = zero
                    if( abs(lua-lub).eq.1 ) then
                    !--- for AAP
                        ALP(-1) = lub;                      ALP(1)  = -lub - one                        
                        fun = wav(i1, it)%xg*( dgb + ALP(lua-lub)*wav(i2, it)%xg/well%xr )

                    !--- for APA
                        ALP(-1) = lua;                      ALP(1)  = -lua - one
                        fup = wav(i2, it)%xg*( dga + ALP(lub-lua)*wav(i1, it)%xg/well%xr )
                    end if
                    
                    if( abs(lda-ldb).eq.1 ) then
                    !--- for AAP
                        ALP(-1) = ldb;                      ALP(1)  = -ldb - one
                        fun = fun + wav(i1, it)%xf*( dfb + ALP(lda-ldb)*wav(i2, it)%xf/well%xr )
                        
                        ALP(-1) = lda;                      ALP(1)  = -lda - one
                        fup = fup + wav(i2, it)%xf*( dfa + ALP(ldb-lda)*wav(i1, it)%xf/well%xr )
                    end if

                    call simps( fun(2), well%npt-1, well%h, AAP )
                    call simps( fup(2), well%npt-1, well%h, APA )

                    tem = tem + hja*hjb*C3J(1,ka,kb)*AAP*( APA*vv- AAP*uv*two )
                    
                !--- Correction over radius
                    a1  = one;          if( abs(lua-lub).ne.1 ) a1 = zero
                    a2  = one;          if( abs(lda-ldb).ne.1 ) a2 = zero
                    fun = ( wav(i1,it)%xg*wav(i2,it)%xg*a1 + wav(i1,it)%xf*wav(i2,it)%xf*a2 )*well%xr

                    Call simps( fun(2), well%npt-1, well%h, AAP )

                    cmr = cmr - hja*hjb*AAP**2*C3J(1,ka,kb)*( vv  - two*uv ) 

                end do
                ECM(it) = ECM(it) + tem*hom
                PCM(it) = PCM(it) + tem
                RCM(it) = RCM(it) + cmr
            end do
            
            RCM(it) = RCM(it)/nuc%npr(it)
                
        !--- The total results: radius, density and binding energy corrections 
            RCM(0)  = RCM(it)*nuc%npr(it)/nuc%amas + RCM(0)
            PCM(0)  = PCM(it) + PCM(0)
            ECM(0)  = ECM(it) + ECM(0)
        end do

    end if
    
    if(lpr) write(l6,*) '********************* End Center of Mass Correction ************************'
    return
    
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine CoM                                                                                !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!*************************************************************************************************!
!                                                                                                 !
End module Center                                                                                 !
!                                                                                                 !
!*************************************************************************************************!
