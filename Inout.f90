!*************************************************************************************************!
!                                                                                                 !
Module INOUTS                                                                                     !
!                                                                                                 !
!*************************************************************************************************!
Use Define
Use rhflib
Use Gogny
Use Expectation
implicit none

contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine Reader                                                                                 !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- This subroutine reads the input parameters from the file rhf.dat
!
!----------------------------------------------------------------------c

!--- Input file name
    integer :: i, Nmin

!--- Open the input file
    if (lin.ne.0) open(lin,file='rhf.dat',status='old')
    read(lin,100) l6

!--- Parameters for the radial Mesh:
    read(lin,101) well%R, r0   ! Box-Radius and r0 
    read(lin,101) well%h       ! Step-size
    well%npt   = (well%R + 0.0001)/well%h + 1     ! No. of steps
    do i = 1, well%npt
        well%xr(i)   = (i-1)*well%h  !value of radial coordinate of each mesh
    end do
    well%xr(1)   = well%h*1.d-2  ! 1.d-2 equals 1.e-2 !上一语句使xmesh(1)=0,该语句对其重新赋值为0.001.

    if (well%npt.gt.MSD) Stop ' in READER: npt too large'
    if (well%npt.lt.MSD) write(*,*) ' Notice that npt is less than MSD in Reader!'

!--- Parameter for the iteration:
    read(lin,100) maxi          ! maximal number of iterations
    read(lin,101) xmix          ! mixing parameter
    xmix0  = xmix

!--- Initialization of wavefunctions:
    read(lin,100) inin          ! initialization of wavefunctions

!--- Nucleus information: element name, mass number, nucleon numbers
    read(lin,'(2a1,3i4)') nuc%nucnam, nuc%npr 
    nuc%amas    = nuc%npr(0)
    Call NUCLEI(1, nuc%npr(2), nuc%nucnam)
    
!--- Cutoff energy for pairing correlations: 
    read(lin,101) pair%ECT      
!---- Chemical Potentials
    read(lin,101) fermi%ala

!---- Block information: need to check
    read(lin,102) pair%ikb, pair%inb
    
!--- Input the pairing choice
    read(lin,102) pair%is
    
!--- Input the factor for Gogny force
    if(pair%is(1).eq.4 .or. pair%is(2).eq.4) then
        GPF = DSA1
        read(lin,101) GPF%vg
    else
        read(lin,*)
    end if

    close(lin)
100 format(10x,i8)
101 format(10x,2f10.3)
102 format(10x,6i4)


!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine reader                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine output                                                                                 !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- For output

!--- Determine OUTPUT filename
    if(nuc%npr(0).lt.10) then
        if(nuc%nucnam(1).eq.' ') then
            write(name,'(a1,i1)') nuc%nucnam(2), nuc%npr(0)
            IFN = 2
        else
            write(name,'(2a1,i1)') nuc%nucnam, nuc%npr(0)
            IFN = 3
        end if
    else if(nuc%npr(0).lt.100) then
        if(nuc%nucnam(1).eq.' ') then
            write(name,'(a1,i2)') nuc%nucnam(2), nuc%npr(0)
            IFN = 3
        else
            write(name,'(2a1,i2)') nuc%nucnam, nuc%npr(0)
            IFN = 4
        end if
    else 
        if(nuc%nucnam(1).eq.' ') then
            write(name,'(a1,i3)') nuc%nucnam(2), nuc%npr(0)
            IFN = 4
        else
            write(name,'(2a1,i3)') nuc%nucnam, nuc%npr(0)
            IFN = 5
        end if
    end if

!--- Open the output file: 'element name//mass number//.//parameter set name'
    Open(l6,file=name(1:IFN)//'.'//pset%txtfor,status='unknown')
    write(l6,*) ' *************** BEGIN READER ***************************'
    write(l6,201) ' Box-size  Rmax and r0 : ', well%R, r0
    write(l6,201) ' Step-size: ', well%h
    write(l6,200) ' Number of Mesh-points: ', well%npt
    write(l6,200) ' Maximal numb. of iterations: ',maxi
    write(l6,201) ' Mixing parameter     : ',xmix
    write(l6,200) ' inin: initial wavefunctions:  ', inin
    write(l6,'(a,20x,2a1,3i4)') ' Nucleus:  ',nuc%nucnam, nuc%npr
    
    write(l6,*)
    write(l6,'(/,a,a)') 'The effective interaction Used: ', pset%txtfor
    write(l6,'(a,2f11.6)')      '   M =', pset%amu
    write(l6,'(4(a,f11.6,2x))') 'msig =', pset%amsig,' mome =', pset%amome,' mrho =', pset%amrho,' mpi  =', pset%ampio
    write(l6,'(5(a,f11.6,2x))') 'gsig =', pset%gsig, ' asig =', pset%asig, ' bsig =', pset%bsig, ' csig =', pset%csig, &
                            &   'dsig =', pset%dsig
    write(l6,'(5(a,f11.6,2x))') 'gome =', pset%gome, ' aome =', pset%aome, ' bome =', pset%bome, ' come =', pset%come, &
                            &   'dome =', pset%dome
    write(l6,'(6(a,f11.6,2x))') 'grho =', pset%grho, ' arho =', pset%arho, ' frho =', pset%grtn, ' artn =', pset%artn, &
                            &   ' fpi =', pset%fpio, ' api  =', pset%apio    
    write(l6,'(a, f11.6)' )     'rho0 =', pset%rvs

    write(6,'(/,a,a)') 'The effective interaction Used: ', pset%txtfor
    write(6,'(1x,a,2f10.4)')   '   M = ', pset%amu
    write(6,'(4(1x,a,f10.4))') 'msig = ', pset%amsig,' mome =', pset%amome,' mrho =', pset%amrho,' mpi  =', pset%ampio
    write(6,'(5(1x,a,f10.4))') 'gsig = ', pset%gsig, ' asig =', pset%asig, ' bsig =', pset%bsig, ' csig =', pset%csig, &
                            &  'dsig =', pset%dsig
    write(6,'(5(1x,a,f10.4))') 'gome = ', pset%gome, ' aome =', pset%aome, ' bome =', pset%bome, ' come =', pset%come, &
                            &  'dome =', pset%dome
    write(6,'(6(1x,a,f10.4))') 'grho = ', pset%grho, ' arho =', pset%arho, ' frho =', pset%grtn, ' artn =', pset%artn, &
                            &  ' fpi =', pset%fpio, ' api  =', pset%apio
    write(6,'(1x,a,f10.4)' )   'rho0 = ', pset%rvs

!---- Frozen Gap-Parameters
    write(l6,201) ' Frozen Gap Parameters:      ',pair%dec

!---- Pairing-Constants
    write(l6,201) ' Pairing-Constant:           ',pair%ga

!---- Initial Gap-Parameters
    write(l6,201) ' Initial Gap Parameters:     ',pair%del

!--- Cutoff energy for pairing correlations
    write(l6,201) ' Cutoff energy for d-BCS:    ',pair%ECT

!---- Chemical Potentials
    write(l6,201) ' Chemical Potentials:        ', fermi%ala

!---- Chemical Potentials
    write(l6,202) ' Block levels and numbers:   ', pair%ikb, pair%inb

200 format(a,i8)
201 format(a,2f10.3)
202 format(a,6i4)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine output                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine inout(is)                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
    integer :: is, it, i, i0, npt
    double precision :: rtmp, fun(MSD)

!--- Input the initial potential  !!(new: Input the initial wave function)
    npt = well%npt
    if(is.eq.1) then   
        if(inin.eq.1) return

!--- INPUT the wave functions and occupation probabilitiess
        do it = 1, 2
            open(20, file = 'WAV/'//name(1:IFN)//'.G-'//tex%tit(it)//'.'//pset%txtfor, status = 'old')
            open(21, file = 'WAV/'//name(1:IFN)//'.F-'//tex%tit(it)//'.'//pset%txtfor, status = 'old')
            read(20, 113) (lev(it)%vv(i0), i0 = 1, lev(it)%nt)
            read(21,*)
            read(20,*)
            read(21,*)
            do i = 1, npt
                read(20, 112) rtmp, (wav(i0,it)%xg(i), i0 = 1, lev(it)%nt)
                read(21, 112) rtmp, (wav(i0,it)%xf(i), i0 = 1, lev(it)%nt)
            end do
            do i0 = 1, lev(it)%nt
                fun = wav(i0,it)%xg**2 + wav(i0,it)%xf**2
                call simps(fun, npt, well%h, rtmp)
                rtmp    = one/dsqrt(rtmp)
                wav(i0,it)%xg   = wav(i0,it)%xg*rtmp
                wav(i0,it)%xf   = wav(i0,it)%xf*rtmp
            end do
            
            close(20)
            close(21)
        end do
              
    else

!--- OUTPUT the wave functions and occupation probabilities
        Call system("mkdir WAV")
        do it = 1, 2
            open(20, file = 'WAV/'//name(1:IFN)//'.G-'//tex%tit(it)//'.'//pset%txtfor, status = 'unknown')
            open(21, file = 'WAV/'//name(1:IFN)//'.F-'//tex%tit(it)//'.'//pset%txtfor, status = 'unknown')
            write(20, 113) (lev(it)%vv(i0), i0 = 1, lev(it)%nt)
            write(21, 113) (lev(it)%vv(i0), i0 = 1, lev(it)%nt)
            write(20, 114) '    r   ', (lev(it)%tb(i0), i0 = 1, lev(it)%nt)
            write(21, 114) '    r   ', (lev(it)%tb(i0), i0 = 1, lev(it)%nt)
            
            do i = 1, npt
                write(20, 112) well%xr(i), ( wav(i0,it)%xg(i), i0 = 1, lev(it)%nt )
                write(21, 112) well%xr(i), ( wav(i0,it)%xf(i), i0 = 1, lev(it)%nt )
            end do
            close(20);          close(21)
        end do
112     format(f8.4, 100e20.12)
113     format(8x, 100e20.12)
114     format(a, 100(8x,a8, 4x) )
    end if

    if(is.eq.1) return

!--- Print out the single particle levels to other files
    Call system("mkdir LEV")
    do it = 1, 2
    !--- OUTPUT the positive energy states
        open(30, file = 'LEV/'//name(1:IFN)//'.psp-'//tex%tit(it)//'.'//pset%txtfor, status = 'unknown')
        write(30, 77) 'i', 'State', 'mu', 'vv', 'eig', 'del', 'spk', 'rr', 'gg', 'ff', 'fp'
    77  format(2x, a, 4x, a, 3x, a, 6x, a, 4x, 3(5x, a, 4x), 4(6x, a, 4x) )         
        do i  = 1, lev(it)%nt
            write(30,107) i, lev(it)%tb(i), lev(it)%mu(i), lev(it)%vv(i), lev(it)%ee(i), lev(it)%del(i), &
                & lev(it)%spk(i), wav(i,it)%r2, wav(i,it)%gg, wav(i,it)%ff, wav(i,it)%ff/(wav(i,it)%gg + wav(i,it)%ff)
        end do
        close(30)
    107 format(i3, 2x, a, i4, 8f12.6)

    end do

!--- output the bulk quantities to special file
    if(nuc%nucnam(1).eq.' ') then
        open(3,file=name(1:1)//'.'//pset%txtfor, status='unknown')
    else
        open(3,file=name(1:2)//'.'//pset%txtfor, status='unknown')
    end if
12     read(3,'(a)',end = 13)  
    go to 12
13     write(3, 108)nuc%nucnam, nuc%amas, ene%etot(0), ene%rms,ene%rch, ene%rcd, ene%ecom(0), pair%del
    close(3)
108 format(2a1,f8.1, f14.6, 9f12.6)

!--- Output the densities
    Call system("mkdir DEN")
    open(33, file = 'DEN/'//name(1:IFN)//'.den'//'.'//pset%txtfor, status = 'unknown')
    write(33,'(4x,a,10x,a, 6(10x,a))') 'r', 'rhos(n)', 'rhos(p)', 'rhob(n)', 'rhob(p)', &
        &                                   'rhot(n)', 'rhot(p)', 'rhoc'
    do i = 1, well%npt
        write(33, 109) well%xr(i), dens(1)%rs(i), dens(2)%rs(i), dens(1)%rv(i), dens(2)%rv(i), &
            &                      dens(1)%rt(i), dens(2)%rt(i), den%rc(i)
    end do
109 format(f8.3, 7E17.8)
    close(33)

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine inout                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!*************************************************************************************************!
!                                                                                                 !
End Module INOUTS                                                                                 !
!                                                                                                 !
!*************************************************************************************************!
