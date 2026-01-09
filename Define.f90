!*************************************************************************************************!
!                                                                                                 !
Module Define                                                                                     !
!                                                                                                 !
!*************************************************************************************************!
!
!--- Define the grobal variables: parameter sets, radial dimension, etc.
!
implicit none

!--- Mesh point number
    integer, parameter :: MSD = 201
    Integer, parameter :: IBX = 2

!---- number of blocks: KMX(1) for j=l+1/2 orbits and KMX(2) for j=l-1/2 orbits
    integer, parameter, dimension(2) :: KMX = (/4, 4/)
    integer, parameter :: NBX = KMX(1) + KMX(2), KPX = max( KMX(1),KMX(2) )
    Integer, parameter :: LYX = 2*KPX, LDX = LYX

!--- Index for kappa quantities
    TYPE KAPPA_T
        INTEGER :: kp, lu, ld, ip, ka
    END TYPE KAPPA_T
    TYPE (KAPPA_T), DIMENSION(NBX) :: KLP

!---- Number of total levels
    INTEGER, PARAMETER :: NMX = 6
    INTEGER, PARAMETER, private :: NM1 = NMX, NM2 = NMX - 1, KM2 = KMX(2)-1, KM1 = KMX(1)
    INTEGER, PARAMETER :: NTX = ( 2*NM1 + 1 - KM1/2 )*KM1/2 + (NM1 - KM1/2)*MOD(KM1, 2) + &
                              & ( 2*NM2 + 1 - KM2/2 )*KM2/2 + (NM2 - KM2/2)*MOD(KM2, 2) + NMX    
    
!--- Definition of the effective interactions
    type Lagrangian
        character(len=8) :: txtfor                        !--- Interction name
        Integer :: IE                                     !--- IE = 1 for RMF, and IE = 2 for RHF
        double precision :: rvs                           !--- Saturation density determined by the interaction
        double precision, dimension(IBX) :: amu           !--- Nucleon masses
        double precision :: amsig, amome, amrho, ampio    !--- Meson masses
        double precision :: gsig, gome, grho, grtn, fpio  !--- Coupling constants in all channels
        double precision :: arho, artn, apio              !--- Density dependent parameters in isovector channels
        double precision :: asig, bsig, csig, dsig        !--- Density dependent parameters for sigma-scalar
        double precision :: aome, bome, come, dome        !--- Density dependent parameters for omega-vector
    end type Lagrangian

    type (Lagrangian) :: pset                              !--- Variable used through the whole routine
    type (Lagrangian) :: PKO1, PKO2, PKO3, PKA1            !--- Effective interactions in DDRHF
    type (Lagrangian) :: PKDD, DDME1, DDME2, TW99, DDLZ1   !--- Effective interactions in DDRMF

!--- Density dependent parameters and their derivatives with respect to the density
    type couplings
        double precision, dimension(MSD) :: gsig, gome, grho, fpio, grtn
        double precision, dimension(MSD) :: dsig, dome, drho, dpio, drtn
    end type couplings
    type (couplings) :: cct

!--- PKA1: Physcical Reviews C 76 (2007) 034314
    data PKA1%txtfor/'PKA1'/, PKA1%rvs/0.159996/, PKA1%amu/938.9, 938.9/, PKA1%IE/2/
    data PKA1%amsig/488.227904/, PKA1%amome/783.000000/, PKA1%amrho/769.000000/, PKA1%ampio/138.000000/
    data PKA1%gsig /  8.372672/, PKA1%gome / 11.270457/, PKA1%grho /  3.649857/, PKA1%fpio /  1.030722/
    data PKA1%asig /  1.103589/, PKA1%bsig / 16.490109/, PKA1%csig / 18.278714/, PKA1%dsig /  0.135041/
    data PKA1%aome /  1.126166/, PKA1%bome /  0.108010/, PKA1%come /  0.141251/, PKA1%dome /  1.536183/
    data PKA1%arho /  0.544017/, PKA1%apio /  1.200000/, PKA1%grtn /  3.199491/, PKA1%artn /  0.820583/

!--- PKO1: Physics Letters B 640 (2006) 150¨C154
    data PKO1%txtfor/'PKO1'/, PKO1%rvs/0.151989/, PKO1%amu/938.9, 938.9/, PKO1%IE/2/
    data PKO1%amsig/525.769084/, PKO1%amome/783.000000/, PKO1%amrho/769.000000/, PKO1%ampio/138.000000/
    data PKO1%gsig /  8.833239/, PKO1%gome / 10.729933/, PKO1%grho /  2.629000/, PKO1%fpio /  1.000000/
    data PKO1%asig /  1.384494/, PKO1%bsig /  1.513190/, PKO1%csig /  2.296615/, PKO1%dsig /  0.380974/
    data PKO1%aome /  1.403347/, PKO1%bome /  2.008719/, PKO1%come /  3.046686/, PKO1%dome /  0.330770/ 
    data PKO1%arho /  0.076760/, PKO1%apio /  1.231976/, PKO1%grtn /  0.000000/, PKO1%artn /  0.000000/
    
!--- PKO2: EPL, 82 (2008) 12001
    data PKO2%txtfor/'PKO2'/, PKO2%rvs/0.151021/, PKO2%amu/938.9, 938.9/, PKO2%IE/2/
    data PKO2%amsig/534.461766/, PKO2%amome/783.000000/, PKO2%amrho/769.000000/, PKO2%ampio/138.000000/
    data PKO2%gsig /  8.920597/, PKO2%gome / 10.550553/, PKO2%grho /  4.068299/, PKO2%fpio /  0.000000/
    data PKO2%asig /  1.375772/, PKO2%bsig /  2.064391/, PKO2%csig /  3.052417/, PKO2%dsig /  0.330459/
    data PKO2%aome /  1.451420/, PKO2%bome /  3.574373/, PKO2%come /  5.478373/, PKO2%dome /  0.246668/
    data PKO2%arho /  0.631605/, PKO2%apio /  0.000000/, PKO2%grtn /  0.000000/, PKO2%artn /  0.000000/ 

!--- PKO3: EPL, 82 (2008) 12001
    data PKO3%txtfor/'PKO3'/, PKO3%rvs/0.153006/, PKO3%amu/938.9, 938.9/, PKO3%IE/2/
    data PKO3%amsig/525.667686/, PKO3%amome/783.000000/, PKO3%amrho/769.000000/, PKO3%ampio/138.000000/
    data PKO3%gsig /  8.895635/, PKO3%gome / 10.802690/, PKO3%grho /  3.832480/, PKO3%fpio /  1.000000/
    data PKO3%asig /  1.244635/, PKO3%bsig /  1.566659/, PKO3%csig /  2.074581/, PKO3%dsig /  0.400843/
    data PKO3%aome /  1.245714/, PKO3%bome /  1.645754/, PKO3%come /  2.177077/, PKO3%dome /  0.391293/ 
    data PKO3%arho /  0.635336/, PKO3%apio /  0.934122/, PKO3%grtn /  0.000000/, PKO3%artn /  0.000000/ 

!--- TW99: Nuclear Physics A 656 (1999) 331
    data TW99%txtfor/'TW99'/, TW99%amu/939.0, 939.0/, TW99%rvs/0.153004/, TW99%IE/1/
    data TW99%amsig/550.0     /, TW99%amome/783.0     /, TW99%amrho/763.0     /, TW99%ampio/138.0     /
    data TW99%gsig / 10.72854 /, TW99%gome / 13.29015 /, TW99%grho/   6.127157/, TW99%fpio /  0.0     /
    data TW99%asig /  1.365469/, TW99%bsig /  0.226061/, TW99%csig/   0.409704/, TW99%dsig /  0.901995/
    data TW99%aome /  1.402488/, TW99%bome /  0.172577/, TW99%come/   0.344293/, TW99%dome /  0.983955/
    data TW99%arho /  0.515   /, TW99%apio /  0.000000/, TW99%grtn /  0.000000/, TW99%artn /  0.000000/ 

!--- DDME1: Phys. Rev. C66, 024306(2002)
    data DDME1%txtfor/'DD-ME1'/, DDME1%amu/938.5, 938.5/, DDME1%rvs/0.151962/, DDME1%IE/1/
    data DDME1%amsig/549.5255/, DDME1%amome/783.0   /, DDME1%amrho/763.0   /, DDME1%ampio/138.0   /
    data DDME1%gsig / 10.4434/, DDME1%gome / 12.8939/, DDME1%grho /  6.2789/, DDME1%fpio /  0.0   /
    data DDME1%asig /  1.3854/, DDME1%bsig /  0.9781/, DDME1%csig /  1.5342/, DDME1%dsig /  0.4661/
    data DDME1%aome /  1.3879/, DDME1%bome /  0.8525/, DDME1%come /  1.3566/, DDME1%dome /  0.4957/
    data DDME1%arho /  0.5008/, DDME1%apio /  0.0   /, DDME1%grtn /  0.0   /, DDME1%artn /  0.0   / 

!--- DDME2: Phys. Rev. C71, 024312(2005)
    data DDME2%txtfor/'DD-ME2'/, DDME2%amu/938.9, 938.9/, DDME2%rvs/0.152/, DDME2%IE/1/
    data DDME2%amsig/550.1238/, DDME2%amome/783.0   /, DDME2%amrho/763.0   /, DDME2%ampio/138.0   /
    data DDME2%gsig / 10.5396/, DDME2%gome / 13.0189/, DDME2%grho /  6.4792/, DDME2%fpio /  0.0   /
    data DDME2%asig /  1.3881/, DDME2%bsig /  1.0943/, DDME2%csig /  1.7057/, DDME2%dsig /  0.4421/
    data DDME2%aome /  1.3892/, DDME2%bome /  0.9240/, DDME2%come /  1.4620/, DDME2%dome /  0.4775/
    data DDME2%arho /  0.5647/, DDME2%apio /  0.0   /, DDME2%grtn /  0.0   /, DDME2%artn /  0.0   / 

!--- PKDD: Phys. Rev. C69, 0034319(2004)
    data PKDD%txtfor/'PKDD'/, PKDD%amu/939.5731, 938.2796/, PKDD%rvs/0.149552/, PKDD%IE/1/
    data PKDD%amsig/555.511236  /, PKDD%amome/783.0       /, PKDD%amrho/763.0       /, PKDD%ampio/138.0       /
    data PKDD%gsig / 10.738508  /, PKDD%gome / 13.147623  /, PKDD%grho /  5.164857  /, PKDD%fpio/  0.0       /
    data PKDD%asig /  1.32742274/, PKDD%bsig /  0.43512557/, PKDD%csig /  0.69166629/, PKDD%dsig/  0.69421032/
    data PKDD%aome /  1.34217027/, PKDD%bome /  0.37116653/, PKDD%come /  0.61139691/, PKDD%dome/  0.73837631/
    data PKDD%arho /  0.18330476/, PKDD%apio /  0.0       /, PKDD%grtn /  0.0       /, PKDD%artn/  0.0       / 
    
!--- DDLZ1
    data DDLZ1%txtfor/'DD-LZ1'/, DDLZ1%amu/938.9, 938.9/, DDLZ1%rvs/0.158100/, DDLZ1%IE/1/
    data DDLZ1%amsig/538.619216/, DDLZ1%amome/783.000000/, DDLZ1%amrho/769.000000/, DDLZ1%ampio/138.0     /
    data DDLZ1%gsig / 12.001429/, DDLZ1%gome / 14.292525/, DDLZ1%grho /  7.575467/, DDLZ1%fpio/  0.0     /
    data DDLZ1%asig /  1.062748/, DDLZ1%bsig /  1.763627/, DDLZ1%csig /  2.308928/, DDLZ1%dsig/  0.379957/
    data DDLZ1%aome /  1.059181/, DDLZ1%bome /  0.418273/, DDLZ1%come /  0.538663/, DDLZ1%dome/  0.786649/
    data DDLZ1%arho /  0.776095/, DDLZ1%apio /  0.000000/, DDLZ1%grtn /  0.000000/, DDLZ1%artn/  0.0     / 
    
!--- Classical values of the radius factor
    double precision :: r0
    data r0/1.2/   
    
!--- Index of isospin and charge
    double precision, dimension(IBX), parameter :: tauz = (/1.0d0, -1.0d0/), tauc = (/0.0d0, 1.0d0/)

!--- Variables of the self-consistent iteration
    double precision :: si, siold, epsi, xmix, xmix0, xmax, sie
    integer :: maxi, ii, inxt, iaut
    data si/1.0d0/, epsi/0.1d-5/, xmix/.50d0/, xmix0/.50d0/, xmax/1.0d0/
    data maxi/500/, inxt/1/, iaut/0/

!--- Variables of the input and output
    integer :: l6, lin 
    data l6/10/, lin/30/ 

!--- Initial choices: if inin = 0, get the wavefunctions from the files, otherwise, determined by self-consistent iteration
    integer :: inin
    data inin/1/

!--- File name for output
    character*8 :: name
    integer :: IFN

!--- Variables for nucleus: element name, mass number, nucleon numbers, etc.
    type nuclide 
        double precision :: amas, etot, rch
        integer, dimension(0:IBX) :: npr
        character, dimension(2) :: nucnam
    end type nuclide
    type (nuclide) :: nuc

!--- Notations of the physical quantities: nucleon kinds, orbits, etc.
    type text
        character, dimension(:) :: tp(2), tis(2), tit(2)
        character, dimension(:) :: tl(0:20)
    end type text

    type (text) :: tex
    data tex%tp/'+', '-'/, tex%tis/'n', 'p'/, tex%tit/'N', 'P'/
    data tex%tl/'s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o', 'P', 'q', 'r', &
             &  'S', 't', 'u', 'v', 'w'/

!--- Radial distance, which is equally spaced by h and cut off at R
    type mesh
        integer :: npt  !No. of steps
        double precision :: h, R
        double precision, dimension(MSD) :: xr
    end type
    type (mesh) :: well
    data well%R/20.0d0/, well%h/0.10/

!--- Variables of the kappa chunk
    type bloblo
        integer :: nb
        integer, dimension(NBX) :: kpa, ka, id, ia
        Integer, dimension(-KMX(1):KMX(2)) :: ib
        integer, dimension(NBX) :: l, lp
    end type bloblo
    type (bloblo), dimension(IBX) :: chunk

!--- Variables of the single-particle states
    type level
        integer :: nt
        logical, dimension(NTX) :: lpb
        integer, dimension(NTX) :: nk, nr, ma, mu, ic, lu, ld, ib
        double precision, dimension(NTX) :: ee, vv, del, spk
        character(len=8), dimension(NTX) :: tb
    end type level
    type (level), dimension(IBX) :: lev 

!--- Upper and lower components of the Wave functions
    type waves
        double precision :: r2, gg, ff
        double precision, dimension(MSD) :: xg, xf 
        double precision, dimension(MSD) :: rho
    end type waves

    type (waves), dimension(NTX, IBX) :: wav 

!--- Local (Hartree and rearrangement terms) potentials
    type dpotel
        double precision, dimension(MSD) :: vps, vms, vtt, vpsh, vmsh, vtth
    end type dpotel

    type (dpotel), dimension(IBX) :: dpotl 

!!--- Non-local (Fock terms) potentials
!    type epotel
!        double precision, dimension(MSD, MSD) :: XG, XF, YG, YF
!    end type epotel
!    type (epotel), dimension(NBX, IBX) :: epotl
!    
!!--- Direct part of the mean fields 
!    type dselfs
!        double precision, dimension(MSD) :: sig, ome, rho, cou, rtn, rvt, rtv 
!    end type dselfs
!
!    type (dselfs) :: dself
!
!!--- Exchange part of the mean fields
!    type eselfs
!        double precision, dimension(MSD, MSD) :: sigxg, omexg, rhoxg, couxg, pioxg, rtnxg, rvtxg, rtvxg
!        double precision, dimension(MSD, MSD) :: sigyg, omeyg, rhoyg, couyg, pioyg, rtnyg, rvtyg, rtvyg
!        double precision, dimension(MSD, MSD) :: sigxf, omexf, rhoxf, couxf, pioxf, rtnxf, rvtxf, rtvxf
!        double precision, dimension(MSD, MSD) :: sigyf, omeyf, rhoyf, couyf, pioyf, rtnyf, rvtyf, rtvyf
!    end type eselfs
!
!    type (eselfs), dimension(NBX, IBX) :: eself
!
!!--- Rearrangement term
!    double precision, dimension(MSD) :: SigR 

!--- Information for Fermi levels
    type fermis
        double precision, dimension(IBX) :: ala 
    end type fermis
    type (fermis) :: fermi
    data fermi%ala/-7.0,-7.0/

!--- Pairing Information
    type pairs
        integer, dimension(IBX) :: ikb, inb, is, i0b
        double precision, dimension(IBX) :: pwi, apwi, ECT
        double precision, dimension(IBX) :: ga, gg, del, spk, dec
    end type pairs
    type (pairs) :: pair
    data pair%ga/0.0,0.0/,pair%del/0.0,0.0/,pair%dec/0.0,0.0/,pair%pwi/3.0d0,  3.0d0/
    data pair%spk/0.0,0.0/, pair%apwi/0.5d0, 1.0d0/, pair%ECT/10.0, 10.0/, pair%is/1,1/

!--- Integral parameters
    double precision, dimension(6) :: oint, eint
    data oint/1.921875d0, 8.218750d0, 4.000000d0, 7.218750d0, 5.578125d0, 6.062500d0/
    data eint/6.062500d0, 5.578125d0, 7.218750d0, 4.000000d0, 8.218750d0, 1.921875d0/

!--- Define math and physics constants
    double precision :: one, two, half, third, zero, hbc, pi, alphi

    data one/1.0d0/, two/2.0d0/, half/0.5d0/, zero/0.0d0/
    data third/0.3333333333333333333333333333333333333333333333333333333333d0/
    data pi/3.141592653589793d0/, hbc/197.328284d0/, alphi/137.03602d0/


!*************************************************************************************************!
!                                                                                                 !
End Module Define                                                                                 !
!                                                                                                 !
!*************************************************************************************************!
