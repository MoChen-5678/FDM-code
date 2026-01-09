!*************************************************************************************************!
!                                                                                                 !
Module Configuration                                                                              !
!                                                                                                 !
!*************************************************************************************************!
!
Use Define
Use Density, only: dens 
Use rhflib
Use Powell_hybrid
Use Gogny
implicit none


!--- Quantities for Delta
    integer, dimension(NTX), private :: kdx
    integer, private :: ist, isp, ncl
    double precision, dimension(NTX, NTX), private :: RI
    
    double precision, dimension(MSD, MSD), private :: gg, gf, fg, ff

Contains
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine Config(lpr)                                                                            !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!---  determines the levels in spherical Dirac program
!
!----------------------------------------------------------------------c
    logical, intent(in) :: lpr
    character(len=8) :: tb1
    character(len=1) :: tp, tl,tis

    integer, dimension(KMX(1)) :: nmax
    integer, dimension(2) :: mt, ks, ls

    integer :: i, j, it, ip, ib, l, i0, n, itz, nn, n0, kappa, nlk
    double precision :: rmatch, ee1

    data ks/-1, 1/, ls/-1, 0/
    
    if (lpr) write(l6,*) ' ********************** Initializing Configurations *********************************'

!--- loop over neutron, proton  (it=1: neutron, it=2: proton)
    do it = 1,2                !loops over isospin
        ib              = 0    !ib is number of blocks
        mt(it)          = 0
        lev(it)%nt      = 0
        chunk(it)%ia(1) = 0

    !--- loop over spin direction (ip=1: j=l+1/2, ip=2: j=l-1/2)
        do ip = 1,2 !loops of spin direction

        !--- Detemrine the quantum numbers for each orbit, prepare the index for single particle levels
            do j = 1, KMX(ip)

                nlk = NMX - (j+ip-1)/2 + mod(j+ip,2)
                if( nlk.le.0 ) cycle

                ib                  = ib + 1                            !----Number of kappa-chunk
                lev(it)%nt          = lev(it)%nt + nlk                  !--- Number of single particle levels
                chunk(it)%id(ib)    = nlk                               !--- Number of orbits in a kappa chunk

            !--- Number of single particle levels before current chunk
                if (ib.gt.1) chunk(it)%ia(ib)   = chunk(it)%ia(ib-1) + chunk(it)%id(ib-1)

                if (ib.gt.NBX)          stop ' in Config: NBX too small'
                if (lev(it)%nt.gt.NTX)  stop ' in Config: NTX too small'

            !--- Determine kappa values: kappa = -(j+1/2) if j=l+1/2, kappa = j+1/2 if j=l-1/2
                chunk(it)%kpa(ib)   = j*ks(ip)  ! ks/-1, 1/; ip=1 for kappa<0, ip=2 for kappa>0
                    
            !--- build the link from kappa to ib: there is an one-to-one correspondence between kappa and ib 
                kappa               = chunk(it)%kpa(ib)
                chunk(it)%ib(kappa) = ib
                
            !--- Determine the index of the block ib in the series of KMX
                chunk(it)%ka(ib)    = iabs(kappa)
                if( kappa.gt.0 ) chunk(it)%ka(ib) = KMX(1) + kappa

            !--- orbital angular momentum l and l'
                chunk(it)%l(ib)     = j + ls(ip)    
                chunk(it)%lp(ib)    = 2*j - 1 - chunk(it)%l(ib) 
                l                   = chunk(it)%l(ib)

            !--- Determine the radial number, degeneration, initializing the energy, matching point for the single particle states
                i0                  = chunk(it)%ia(ib)
                do n = 1,chunk(it)%id(ib)
                    lev(it)%nk(i0+n)    = chunk(it)%kpa(ib)
                    lev(it)%ib(i0+n)    = ib
                    lev(it)%ic(i0+n)    = 2*j   ! 
                    lev(it)%mu(i0+n)    = 2*j   ! \hat{j}^2=2j+1
                    lev(it)%nr(i0+n)    = n - 1
                    lev(it)%lu(i0+n)    = l
                    lev(it)%ld(i0+n)    = chunk(it)%lp(ib)
                    lev(it)%ma(i0+n)    = 1.2*nuc%amas**third/well%h      ! nuc%amas is mass number

                !--- Determine the notations of the single particle states
                    if(n.lt.10 .and. 2*j-1.lt.10) then
                        write(lev(it)%tb(i0+n), 200) tex%tit(it), '.', n, tex%tl(l), '.', 2*j-1   ! tb is character,
                    else if(n.lt.10 .and. 2*j-1.ge.10) then
                        write(lev(it)%tb(i0+n), 201) tex%tit(it), '.', n, tex%tl(l), 2*j-1
                    else if(n.ge.10 .and. 2*j-1.lt.10) then
                        write(lev(it)%tb(i0+n), 202) tex%tit(it), n, tex%tl(l), '.', 2*j-1
                    else
                        write(lev(it)%tb(i0+n), 203) tex%tit(it), n, tex%tl(l), 2*j-1
                    end if
                        
                !--- In the case of block
                    lev(it)%lpb(i0+n)   = .false.
                    if( (pair%ikb(it).eq.kappa) .and. (n.eq.pair%inb(it)) ) then
                        lev(it)%lpb(i0+n)   = .true.
                        pair%i0b(it)        = i0 + n
                    end if
                        
                end do
                    
 200            format(2a1, i1, 2a1, i1, '/2')
 201            format(2a1, i1,  a1, i2, '/2')
 202            format( a1, i2, 2a1, i1, '/2')
 203            format( a1, i2,  a1, i2, '/2')
                
            end do !--- end loops of kappa
        end do    ! end loops of spin direction

    !--- Total number of blocks
        chunk(it)%nb    = ib

    enddo   ! end loop of isospin
        
!---- OUTPUT the initial configurations
10    if (lpr) then
        write(l6,*) ' Nr    IOC    StateEig.      kappa ',&
        &        ' nodes     NMatch rad.'

        do it = 1,2
            itz = 0
            do i = 1,lev(it)%nt
                kappa = lev(it)%nk(i)
                if (kappa.gt.0) then
                    l   = kappa
                else
                    l   = -kappa-1
                endif
                nn      = 2*(lev(it)%nr(i)-1) + l
                write(l6,100) i,int(lev(it)%ic(i)),lev(it)%tb(i),lev(it)%ee(i), lev(it)%nk(i),lev(it)%nr(i),nn
100             format(i4,i7,4x,a8,8x,f6.1,3i7)
                itz     = itz + lev(it)%ic(i)
            enddo
            write(l6,101) ' Number of levels     nt  = ', lev(it)%nt, NTX 
            write(l6,101) ' Number of particles  tz  = ', itz 
!      write(l6,101) ' Blocked levels:    nrbl= ',nrbl
101 format(a,2i4)
        enddo
    endif
  
    if (lpr) write(l6,*) ' ****** END BASE ***********************************'
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
end subroutine Config                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine Occup(lpr)                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- To determine the occupation of each level
!
!     IS=1  occupation in fixed j-blocks
!        2  occupation by lambda-iteration (from the bottom)
!        3  BCS method with a delta force     
!        4  BCS method with a density-dependent delta force
!
!     nt:    number of levels for each particle type
!     ibl:   number of blocked level
!     tz:    particle number
!     al:    chemical potential 
!     gg:    pairing strength
!     dec:   fixed external pairing field
!     del:   gap parameter 
!     pwi:   pairing wIndow
!
!     ee:    single particle energies
!     vv:    occupation probabilities (0 < vv < 2*mu)
!     mu:    multiplicities of j-shells (j+1/2)
!----------------------------------------------------------------------c
!
    logical, intent(in) :: lpr

    integer, dimension(NTX) :: jx
    double precision, dimension(NTX) ::  eetemp,mutemp
    

    integer, parameter :: MAXD=500
    double precision, parameter :: epsd=1.d-6

    integer :: i, j, k, nstran, m, kt, numnuc, it, is, j1, j2, i1, i2, npt, nd, kt0, ks, k1, k2, mu
    integer, dimension(2) :: nz
    double precision :: etemp, mtemp, de, e1, e2, uv, xmu
    double precision :: al, g, xnz, spk, epair, es, sn, bx, ax, efermi
    double precision :: Gamma, V0, dix, sdi, top
    double precision, dimension(NTX+1) ::  del, del0, fvec
    integer, dimension(NTX) :: kdx0
    double precision, dimension(MSD) :: funs

!--- Variables for DNSQE
    integer :: IOPT, NPRINT, INFO
    double precision :: TOL

    TOL     = 1.d-10

!--- If inin = 0, the configurations is imported from the files
    if(inin.eq.0 .and. ii.eq.1) return

    if (lpr) write(l6,*) ' ***************************** BEGIN OCCUP ********************************'

    ks    = 0
!---- loop over neutron-proton
    do it = 1,2
        is  = pair%is(it)
        
!---- fixed occupation within the differten j-orbits  ! ��is>1ʱ�����е��ܼ���Ҫ����
        if (is.ge.1) then
            nz(it)  = nuc%npr(it)            
            eetemp  = lev(it)%ee 
            mutemp  = lev(it)%mu 
            
        !--- In the case of blocking
            if(pair%ikb(it).ne.0) then
                nz(it)      = nz(it) - 1
                j           = pair%i0b(it)
                mutemp(j)   = mutemp(j) - 2
            end if

       !--- Ordering the orbits by the values of single particle energy  �����ݵ�����������������
            do j = 1,lev(it)%nt - 2
                do k = j+1, lev(it)%nt - 1
                    if( eetemp(k) .LT. eetemp(j) ) then
                        etemp       = eetemp(j)
                        eetemp(j)   = eetemp(k)
                        eetemp(k)   = etemp

                        mtemp       = mutemp(j)
                        mutemp(j)   = mutemp(k)
                        mutemp(k)   = mtemp
                    end if
                end do
            end do

        !--- Determine the fermi energy: in the case of blocking, the blocked one could be above the Fermi level
            numnuc  = 0
            j       = 1
            do while( (numnuc.le.nz(it)) .and. (j.lt.lev(it)%nt) )
                efermi  = eetemp(j)
                numnuc  = numnuc + mutemp(j)
                j       = j + 1
            end do

        !--- If there exist two or more levels which have same energies
            nstran  = 0
            do j = 1, lev(it)%nt - 1
                if (eetemp(j) .EQ. efermi) nstran = nstran + 1
            end do
            if (nstran.GT.1) STOP '2 levels with same energies in occup'

        !--- The orbits below Fermi level have full occupation probability (vv=one)
            numnuc  = 0
            do j = 1, lev(it)%nt

                if ( lev(it)%ee(j).lt.efermi ) then
                    lev(it)%vv(j)   = one
                    mu              = lev(it)%mu(j)
                    
                !--- In the case of blocking and blocking level under the Fermi one
                    if( lev(it)%lpb(j) ) then
                        mu              = mu - 1
                        lev(it)%vv(j)   = dble(mu)/lev(it)%mu(j)
                    end if

                    numnuc          = numnuc + mu
                end if

            end do

        !--- Fermi level and levels above
            do j = 1, lev(it)%nt

                if (lev(it)%ee(j) .EQ. efermi) then

                    mu              = nz(it) - numnuc
                    mu              = min( mu, lev(it)%mu(j) )

                    numnuc          = numnuc + mu
                    lev(it)%vv(j)   = dble(mu)/dble(lev(it)%mu(j))
                    
                end if

                if (lev(it)%ee(j).GT.efermi) then
                    lev(it)%vv(j) = zero 
                    
                !--- if the blocking level is above the Fermi one
                    if( lev(it)%lpb(j) ) then
                        lev(it)%vv(j)   = one/lev(it)%mu(j)
                        numnuc          = numnuc + 1
                    end if

                end if
            end do

            if (numnuc.LT.nuc%npr(it)) STOP 'in occup, less levels'
        end if

    !--- Determine the Hartree-Fock Fermi level
        efermi  = lev(it)%ee(1)
        do j = 1, lev(it)%nt
            if( (lev(it)%vv(j).gt.zero) .and. (efermi.le.lev(it)%ee(j)) )  efermi = lev(it)%ee(j)
        end do

    !--- Fermi energy
        if(is.eq.1) fermi%ala(it)   = efermi


    !---- BCS pairing
        if (is.gt.1) then
            al      = efermi
            de      = pair%del(it)
            g       = pair%gg(it)

        !--- Determine states involved in pairing-window, initializing the pairing gaps
            
            nz(it)  = nuc%npr(it)
            if( pair%ikb(it).ne.0 ) nz(it) = nz(it) - 1       !--- For blocking

        !--- With Gogny pairing force, all the states are taken into account in pairing 
            del     = 2.5d0
            nd      = lev(it)%nt + 1
            del(nd) = al
            ncl     = nz(it)
            ist     = it
            
        !--- Calculate the pairing matrix  elements
            if(is.ge.2)  Call PairI(it)

!--- Solve the gap equations and determine the Fermi level
!--- Solve the gap equations and determine the Fermi level
!--- SUBROUTINE DNSQE (FCN, JAC, IOPT, N, X, FVEC, TOL, NPRINT, INFO) !, WA, LWA)
            IOPT    =  1
            NPRINT  = -1

            Call DNSQE(Delta, DeltaX, IOPT, nd, del, fvec, TOL, NPRINT, INFO) 
            
            If (INFO.eq.0) write(*,*) 'Improper input parameters'
!            If (INFO.eq.1) write(*,*) 'Algorithm estimates that the relative error between X and the solution is at most TOL'
            If (INFO.eq.2) write(*,*) 'Number of calls to FCN has reached or exceeded 100*(N+1) for IOPT=1 or 200*(N+1) for IOPT=2'
            If (INFO.eq.3) write(*,*) 'Number of calls to FCN has reached or exceeded 100*(N+1) for IOPT=1 or 200*(N+1) for IOPT=2'
            If (INFO.eq.4) write(*,*) 'Iteration is not making good progress'

        !--- Effective pairing gap: del(nd) stores the fermi energy, del(1:nd-1) store the pairing gaps
            al              = del(nd)
            de              = zero
            dix             = zero

        !--- Determine the occupation probablities and pairing energy
            do j = 1, lev(it)%nt
                e1              = lev(it)%ee(j) - al
                e2              = one/dsqrt(e1**2 + del(j)**2)
                uv              = half*del(j)*e2

                lev(it)%del(j)  = del(j)
                lev(it)%vv(j)   = half*dabs(1.0d0 - e1*e2)                
                lev(it)%spk(j)  = uv*half
                
                xmu             = lev(it)%mu(j)
                if( lev(it)%lpb(j) ) then
                    xmu             = xmu - 2
                    lev(it)%vv(j)   = ( xmu*lev(it)%vv(j) + 1 )/lev(it)%mu(j)
                end if
                de              = de + lev(it)%del(j)*uv*xmu
                dix             = dix + uv*xmu
            end do
            de              = de/dix

        !--- Effective pairing gaps and Fermi energy
            pair%del(it)    = de
            fermi%ala(it)   = al
        end if


!---- end of neutron-proton loop
    end do

    if (lpr) write(l6,*) ' ****************************** END OCCUP **********************************'
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
end subroutine occup                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
                                                                                                  !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine Delta(nd, xd, fvec, IFlag)                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- 
    integer, intent(in) :: nd
    integer :: IFlag
    double precision, dimension(nd), intent(in) :: xd
    double precision, dimension(nd), intent(out) :: fvec
!
!--- Use kdx and ist to transfer the states included in the pairings.
!
    double precision :: al, e1, e2, dn, en, xmu
    integer :: i1, i2, N, it

!--- The input pairing gap and fermi level

!--- Neutron or Proton
    it          = ist
    N           = nd - 1
    al          = xd(nd)

    fvec(nd)    = -ncl
    do i1 = 1, N

!--- Calculate the matrix elements and new gap parameters
        fvec(i1)        = xd(i1)
        do i2 = 1, N
            e1          = lev(it)%ee(i2) - al
            e2          = one/dsqrt(e1**2 + xd(i2)**2)

            xmu         = lev(it)%mu(i2);               if( lev(it)%lpb(i2) ) xmu = xmu - 2
            dn          = xmu*xd(i2)*e2
            en          = xmu*e1*e2**3

!--- The derivative with respect to delta and lambda for gap equations
            fvec(i1)    = GME(it)%vg(i1,i2)*dn*half + fvec(i1)
        end do
        e1          = lev(it)%ee(i1) - al
        e2          = one/dsqrt(e1**2 + xd(i1)**2)

        xmu         = lev(it)%mu(i1);               if( lev(it)%lpb(i1) ) xmu = xmu - 2    

        fvec(nd)    = fvec(nd) + xmu*(one-e1*e2)*half
    end do
    return
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine Delta                                                                              !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
Subroutine DeltaX(nd, xd, fvec, fjac, Jnd, IFlag)                                                 !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!
!--- 
    integer, intent(in) :: nd, Jnd
    integer :: IFlag
    double precision, dimension(nd), intent(in) :: xd
    double precision, dimension(nd), intent(in) :: fvec
    double precision, dimension(Jnd, nd), intent(out) :: fjac
!
!--- Use kdx and ist to transfer the states included in the pairings.
!
    double precision :: al, e1, e2, dn, en, xmu
    integer :: it, i1, i2, N


!--- The input pairing gap and fermi level

!--- Neutron or Proton
    it          = ist
    N           = nd - 1
    al          = xd(nd)

    fjac(nd,nd) = zero
    do i1 = 1, N

!--- Calculate the matrix elements and new gap parameters
        fjac(i1, nd)    = zero
        do i2 = 1, N
            e1          = lev(it)%ee(i2) - al
            e2          = one/dsqrt(e1**2 + xd(i2)**2)

            xmu         = lev(it)%mu(i2);               if( lev(it)%lpb(i2) ) xmu = xmu - 2

            dn          = xmu*xd(i2)*e2
            en          = xmu*e1*e2**3

!--- The derivative with respect to delta and lambda for gap equations
            fjac(i1,i2) = GME(it)%vg(i1,i2)*en*half*e1
            fjac(i1,nd) = GME(it)%vg(i1,i2)*en*half*xd(i2) + fjac(i1,nd)
        end do
        fjac(i1,i1) = fjac(i1,i1) + one

        e1          = lev(it)%ee(i1) - al
        e2          = one/dsqrt(e1**2 + xd(i1)**2)

        xmu         = lev(it)%mu(i1);               if( lev(it)%lpb(i1) ) xmu = xmu - 2    

        fjac(nd,i1) = xmu*e1*e2**3*xd(i1)*half
        fjac(nd,nd) = fjac(nd,nd) + xmu*e2**3*xd(i1)**2*half
    end do
    return

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!                                                                                                 !
End Subroutine DeltaX                                                                             !
!                                                                                                 !
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++!
!*************************************************************************************************!
!                                                                                                 !
End module Configuration                                                                          !
!                                                                                                 !
!*************************************************************************************************!
