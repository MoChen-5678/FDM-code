!*************************************************************************************************!
!                                                                                                 !
Module PotelHF                                                                                    !
!                                                                                                 !
!*************************************************************************************************!
!
Use Define
Use Meanfield
Use rhflib, only: intpol6
Use BASE, only: BASIS, NBD
Use omp_lib
Implicit None

!--- Non-local (Fock terms) potentials
    type epotel
        double precision, dimension(MSD, MSD) :: XG, XF, YG, YF
    end type epotel
    type (epotel), dimension(NBX, IBX) :: epotl

    type StiffM
        Integer :: NAA
        Integer,          allocatable, dimension(:) :: i1, i2 
        double precision, allocatable, dimension(:) :: AH
    end type StiffM
    type (StiffM), dimension(NBX, IBX) :: STM
    
    logical, private :: lstm = .true.

Contains    

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine Potel(lpr)                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- Determine the potential of nuclear system
!
!     This subroutine generates VPS and VMS
!
    logical, intent(in) :: lpr
    double precision :: emcc(2)
    integer :: i, it, ib 
    double precision :: sv
    double precision, dimension(MSD) :: SS, VS, VT, TT, svp, svm, svt
    double precision, dimension(MSD, MSD) :: sxg, sxf, syg, syf
    type (epotel) :: epot0
    type (dpotel) :: dpot0

    if(lpr)  write(l6,*) ' ****** Begin POTEL **********************************'

!--- Calculate the self-energies from all channels.
    Call Hartree
    if(pset%IE.eq.2) Call Fock
    Call Rearrange

!--- The direct part of meson field
    emcc    = pset%amu/hbc*two
    si      = zero

!$OMP PARALLEL DO private(i) SCHEDULE(STATIC)
    do i = 1, well%npt
        SS(i)   = cct%gsig(i)* dself%sig(i)
        VS(i)   = cct%gome(i)* dself%ome(i) + SigR(i)
        VT(i)   = cct%grho(i)*(dself%rho(i) + dself%rvt(i))
        TT(i)   = cct%grtn(i)*(dself%rtn(i) + dself%rtv(i))
    end do
!$OMP END PARALLEL DO

    do it = 1, 2
        dpot0   = dpotl(it)
    !$OMP PARALLEL DO private(i) SCHEDULE(STATIC)
        do i = 1, well%npt
            dpotl(it)%vps(i)    = VS(i) + SS(i) + tauz(it)*VT(i) + tauc(it)*dself%cou(i)
            dpotl(it)%vms(i)    = VS(i) - SS(i) + tauz(it)*VT(i) + tauc(it)*dself%cou(i) - emcc(it)
            dpotl(it)%vtt(i)    = TT(i)*tauz(it)

        !--- Find the deviation
            svp(i)              = dpot0%vps(i) - dpotl(it)%vps(i)
            svm(i)              = dpot0%vms(i) - dpotl(it)%vms(i)
            svt(i)              = dpot0%vtt(i) - dpotl(it)%vtt(i)
    
        !--- Folding
            if(ii.eq.1) cycle
            dpotl(it)%vps(i)    = dpot0%vps(i) - xmix*svp(i)
            dpotl(it)%vms(i)    = dpot0%vms(i) - xmix*svm(i)
            dpotl(it)%vtt(i)    = dpot0%vtt(i) - xmix*svt(i)
        end do 
    !$OMP END PARALLEL DO
        sv  = maxval( dabs(svp) );          si  = max(si, dabs(sv))
        sv  = maxval( dabs(svm) );          si  = max(si, dabs(sv))
        sv  = maxval( dabs(svt) );          si  = max(si, dabs(sv))
    end do

!---- interpolation of potentials between mesh points
    do it = 1, 2
        Call intpol6(dpotl(it)%vms, dpotl(it)%vmsh, well%npt)
        Call intpol6(dpotl(it)%vps, dpotl(it)%vpsh, well%npt)
        Call intpol6(dpotl(it)%vtt, dpotl(it)%vtth, well%npt)
    end do

    if(pset%IE.eq.1) go to 201

!--- The exchange part
    do it = 1, 2
        do ib = 1, chunk(it)%nb

            epot0   = epotl(ib,it)

        !$OMP PARALLEL DO PRIVATE(i) SCHEDULE(STATIC)
            do i = 1, well%npt
            !--- Calculate the non-local potentials
                epotl(ib,it)%XG(:,i)    = cct%gsig(:)*FSIG(ib,it)%XMP(:,i) + cct%gome(:)*FOME(ib,it)%XMP(:,i) + &
                                        & cct%grho(:)*FRHO(ib,it)%XMP(:,i) + cct%grtn(:)*FRTN(ib,it)%XMP(:,i) + &
                                        & cct%grho(:)*FRVT(ib,it)%XMP(:,i) + cct%grtn(:)*FRTV(ib,it)%XMP(:,i) + &
                                        &             FCOU(ib,it)%XMP(:,i) + cct%fpio(:)*FPIO(ib,it)%XMP(:,i) 

                epotl(ib,it)%XF(:,i)    = cct%gsig(:)*FSIG(ib,it)%XMM(:,i) + cct%gome(:)*FOME(ib,it)%XMM(:,i) + &
                                        & cct%grho(:)*FRHO(ib,it)%XMM(:,i) + cct%grtn(:)*FRTN(ib,it)%XMM(:,i) + &
                                        & cct%grho(:)*FRVT(ib,it)%XMM(:,i) + cct%grtn(:)*FRTV(ib,it)%XMM(:,i) + &
                                        &             FCOU(ib,it)%XMM(:,i) + cct%fpio(:)*FPIO(ib,it)%XMM(:,i) 

                epotl(ib,it)%YG(:,i)    = cct%gsig(:)*FSIG(ib,it)%XPP(:,i) + cct%gome(:)*FOME(ib,it)%XPP(:,i) + &
                                        & cct%grho(:)*FRHO(ib,it)%XPP(:,i) + cct%grtn(:)*FRTN(ib,it)%XPP(:,i) + &
                                        & cct%grho(:)*FRVT(ib,it)%XPP(:,i) + cct%grtn(:)*FRTV(ib,it)%XPP(:,i) + &
                                        &             FCOU(ib,it)%XPP(:,i) + cct%fpio(:)*FPIO(ib,it)%XPP(:,i) 
                    
                epotl(ib,it)%YF(:,i)    = cct%gsig(:)*FSIG(ib,it)%XPM(:,i) + cct%gome(:)*FOME(ib,it)%XPM(:,i) + &
                                        & cct%grho(:)*FRHO(ib,it)%XPM(:,i) + cct%grtn(:)*FRTN(ib,it)%XPM(:,i) + &
                                        & cct%grho(:)*FRVT(ib,it)%XPM(:,i) + cct%grtn(:)*FRTV(ib,it)%XPM(:,i) + &
                                        &             FCOU(ib,it)%XPM(:,i) + cct%fpio(:)*FPIO(ib,it)%XPM(:,i) 

            !--- Check for convergence
                sxg(:,i)    = epot0%XG(:,i) - epotl(ib,it)%XG(:,i)      
                sxf(:,i)    = epot0%XF(:,i) - epotl(ib,it)%XF(:,i)      
                syg(:,i)    = epot0%YG(:,i) - epotl(ib,it)%YG(:,i)      
                syf(:,i)    = epot0%YF(:,i) - epotl(ib,it)%YF(:,i)   
                    
            !--- folding
                if(ii.eq.1) cycle
                epotl(ib,it)%XG(:,i)    = epot0%XG(:,i) - xmix*sxg(:,i)
                epotl(ib,it)%XF(:,i)    = epot0%XF(:,i) - xmix*sxf(:,i)
                epotl(ib,it)%YG(:,i)    = epot0%YG(:,i) - xmix*syg(:,i)
                epotl(ib,it)%YF(:,i)    = epot0%YF(:,i) - xmix*syf(:,i) 
            end do
        !$OMP END PARALLEL DO
            sv  = maxval( dabs(sxf) );              si  = max(dabs(sv),  si)
            sv  = maxval( dabs(syg) );              si  = max(dabs(sv),  si)
            sv  = maxval( dabs(syf) );              si  = max(dabs(sv),  si)
            sv  = maxval( dabs(sxg) );              si  = max(dabs(sv),  si)
        end do
    end do

201 continue
    if (lpr) then
        write(l6,101)
    101 format(/,6x,'x (fm)', 6x,'V+S n',7x,'V-S n',7x,' VT n',7x,'V+S p',7x,'V-S p',7x,' VT p')
        do i = 1, well%npt
            write(l6,100) well%xr(i), dpotl(1)%vps(i), dpotl(1)%vms(i)+emcc(1), dpotl(1)%vtt(i), &
            &                         dpotl(2)%vps(i), dpotl(2)%vms(i)+emcc(2), dpotl(2)%vtt(i)
        enddo
    100 format(f12.4, 6e12.4)    
        write(l6,*) ' ****** END POTEL **********************************'
    endif

    return

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine Potel                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine Stiff                                                                                  !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!

!
!--- This subroutine is used to calculate the stiff matrix
!
!-------------------------------------------------------------------------------------------------!
!
    integer :: ib, it
    double precision, dimension(MSD) :: fun
    double precision :: rel
    
    integer :: i1, i2, i12, i, ka, NAA, NA, ND

!--- Bulid the link between the stiff matrix and base
    If( lstm ) then
        do it = 1, IBX
            do ib = 1, chunk(it)%nb
                ka  = chunk(it)%ka(ib)
                NA  = BASIS(ka,it)%NDF
                NAA = (NA + 1)*NA/2;            STM(ib,it)%NAA  = NAA
                if( .not.allocated( STM(ib,it)%AH ) ) allocate( STM(ib,it)%AH(NAA) )
                if( .not.allocated( STM(ib,it)%i1 ) ) allocate( STM(ib,it)%i1(NAA) )
                if( .not.allocated( STM(ib,it)%i2 ) ) allocate( STM(ib,it)%i2(NAA) )
            
                i12 = 0;            ND  = BASIS(ka,it)%ND
                do i1 = 1, NA
                    do i2 = i1, NA
                        i12 = i12 + 1
                        STM(ib,it)%i1(i12)  = i1;       if( i1.gt.ND ) STM(ib,it)%i1(i12) = NBD + i1 - ND
                        STM(ib,it)%i2(i12)  = i2;       if( i2.gt.ND ) STM(ib,it)%i2(i12) = NBD + i2 - ND
                    end do
                end do
            end do
        end do

        lstm    = .false.
    end if

    do it = 1, IBX
        do ib = 1, chunk(it)%nb
            
            ka  = chunk(it)%ka(ib)

            NAA = STM(ib,it)%NAA

        !--- Calculate the direct part
        !$OMP PARALLEL DO PRIVATE(fun, rel, i12, i1, i2, i) 
            do i12 = 1, NAA
                i1  = STM(ib,it)%i1(i12);       i2  = STM(ib,it)%i2(i12)
                
                fun(:)  = BASIS(ka,it)%XG(:,i1)*BASIS(ka,it)%XG(:,i2)*dpotl(it)%vps(:) + &
                        & BASIS(ka,it)%XF(:,i1)*BASIS(ka,it)%XF(:,i2)*dpotl(it)%vms(:) + &
                        & BASIS(ka,it)%XG(:,i1)*BASIS(ka,it)%XF(:,i2)*dpotl(it)%vtt(:) + &
                        & BASIS(ka,it)%XF(:,i1)*BASIS(ka,it)%XG(:,i2)*dpotl(it)%vtt(:)

                Call simps(fun, well%npt, well%h, rel)

                STM(ib,it)%AH(i12)  = BASIS(ka,it)%XKT(i1,i2) + rel
    
                if( pset%IE.eq.1 ) cycle

                fun = zero
                do i = 1, well%npt
                    fun(:)  = BASIS(ka,it)%XG(:,i1)*epotl(ib,it)%YG(:,i)*BASIS(ka,it)%XG(i,i2) + &
                            & BASIS(ka,it)%XG(:,i1)*epotl(ib,it)%YF(:,i)*BASIS(ka,it)%XF(i,i2) + &
                            & BASIS(ka,it)%XF(:,i1)*epotl(ib,it)%XG(:,i)*BASIS(ka,it)%XG(i,i2) + &
                            & BASIS(ka,it)%XF(:,i1)*epotl(ib,it)%XF(:,i)*BASIS(ka,it)%XF(i,i2) + fun(:)
                end do

                Call simps(fun, well%npt, well%h, rel)

                STM(ib,it)%AH(i12)  = STM(ib,it)%AH(i12) + rel

            end do
        !$OMP END PARALLEL DO

        end do
    end do

    return


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine Stiff                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!*************************************************************************************************!
!                                                                                                 !
End Module PotelHF                                                                                !
!                                                                                                 !
!*************************************************************************************************!