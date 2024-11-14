module common_arrays

  implicit none

!******************************************************************************
  ! unit of time      (to)    = 1.000000000 x 10**(-12) seconds   (pico-seconds)
  ! unit of length    (lo)    = 1.000000000 x 10**(-10) metres    (angstroms)
  ! unit of mass      (mo)    = 1.660540200 x 10**(-27) kilograms (daltons)
  ! unit of charge    (qo)    = 1.602177330 x 10**(-19) coulombs  (electrons)
  ! unit of energy    (eo)    = 1.660540200 x 10**(-23) joules    (10 j mol^-1)
  ! unit of pressure  (po)    = 1.660540200 x 10**(  7) pascals  (163.882576
  ! atm).
  !                           = 1.638825760 x 10**(  2) atmospheres
!*****************************************************************************

!************************************************************************
!*
!*          index of parameters and common variables
!*          ========================================
!*
  integer, parameter :: nom=10     ! max number of molecules (32,108,256,500,864)
  integer, parameter :: nsite=30   ! max number of sites 
  integer,parameter  :: nclus=3500
  integer, parameter :: npimax=50  ! number of pair interactions
  integer, parameter :: ncsmax=3   ! max coreshell pairs 
  integer, parameter :: npolc=3  ! number of polarization sites per molecule
  integer, parameter :: nhist = 1000 !dimension of rdf histogram 
  integer, parameter :: ksqmax = 200 ! max squred vlues of k vectors indices in FOurier sum 
  real(8), parameter :: pi=dacos(-1.d0)  ! pi
  real(8), parameter :: boltz=8.31451115d-1  ! Boltzman Const eo/(k*mol) (dlpoly)
  real(8), parameter :: avsno=6.0222e+23    ! Avogadro number
  real(8), parameter :: r4pie0=138935.4835d0 ! 1/4pi*epsilon in internal units
  real(8), parameter :: ewp = 1.0d-6   ! precision for Ewald sum
  real(8), parameter :: prsunt=0.163882576d0  ! internal units to katm
  real(8), parameter :: bohr2a = 0.529177249d0
  real(8), parameter :: h2kcal = 627.510d0   ! Hartree to kcal/mol conversion
  real(8), parameter :: twosqpi= 2.0d0/sqrt(pi)  

  integer  :: nm                   ! number of molecule in the system
  integer  :: ns                   ! sites per molecules
  integer  :: natom                ! number  of massive sites per molecules  
  integer :: npi                   ! number of pair interactions in field file
  character(len=15) :: system_name ! name of the system in the field file
  character(len=3) :: an(nsite)    ! name of atoms in a molecule
  character(len=3) :: ann(nsite,nom)    ! name of atoms in a molecule

  character(len=15) :: initial_config

  real(8) :: box,hbox,diameter ! molecular diameter
  character(len=3) :: an1(npimax),an2(npimax) ! site pair names in vdw interactions  
  real(8) :: massites(nsite) ! mass of sites 
  real(8) :: coords(3,nsite) ! coordinates of given molecular geometry 
  real(8) :: pcoord(3,nsite) ! principal corrdiates of each sites in molecule
  real(8) :: pmi(3)          ! principal moment of inertia
  real(8) :: charge(nsite)   ! charge on sites
  real(8) :: volm   ! volume of simulation box
  real(8) :: totm   ! mass of single molecule

 !  arrays for molecular centers  
  real(8) :: x(nom),y(nom),z(nom) ! position of com of molecule
  real(8) :: q0(nom),q1(nom),q2(nom),q3(nom) ! quaternions   
  real(8) :: jx(nom),jy(nom),jz(nom)  ! angular momentum of com of molecule

  real(8) :: xs(nsite,nom), ys(nsite,nom), zs(nsite,nom) ! sites corrds
  real(8) :: xxst(nsite,nom),yyst(nsite,nom),zzst(nsite,nom) ! for intermediate saving 
  real(8) :: xxt(nclus,nsite,nom),yyt(nclus,nsite,nom),zzt(nclus,nsite,nom) ! for intermediate saving
  real(8) :: chgs(nsite,nom)                              ! charges on sites

  integer :: integrator ! integer to select integrator
  character(len=15) :: pot_name   ! water potential
  character(len=10) :: ensemble   ! ensemble for simulation

  real(8) :: vdwpe,vir,tke,rke ! PE, virial, trans. KE, rot. KE
  real(8) :: rcpe,kcpe,scpe ! real, kspace, self coulumb pe
  real(8) :: vir_vdw,vir_rc,vir_kc,vir_self
  real(8) :: indpe,point_indpe,phys_indpe          ! induction energy (point dipole polarization model)
  real(8) :: fs2pe,fs3pe          ! 3b energy  
  real(8) :: vdwlrc,virlrc  ! long corrections for VDW and virial
  real(8) :: sqf,sqt        
  real(8) :: stpvdwpe ! step vdw pe
  real(8) :: stpcpe   ! sum of real, kspace, self coulumb pe
  real(8) :: stpindpe,stpindpe_point,stpindpe_phys ! step induction energy
  real(8) :: stp3bpe ! step 3b energy
  real(8) :: stppe          ! step total pe
  real(8) :: stptke ! step inst trans. ke
  real(8) :: stpttp ! step inst trans. temp
  real(8) :: stprke ! inst rot. ke
  real(8) :: stprtp ! inst rot. temp
  real(8) :: stpvir ! step vir
  real(8) :: stpte  ! step total energy
  real(8) :: stpprs ! step presurre 
  real(8) :: stptotke ! step tot kinetic energy
  real(8) :: stptmp ! step total temperature

  real(8) :: sumttp,sumrtp,sumprs
  real(8) :: sumpe,sumvdwpe,sumcpe,sumvir,sumtke,sumrke,sumte,sumtmp
  real(8) :: sumindpe, sum3bpe, sumtotke

  real(8) :: ssqvdwpe,ssqcpe,ssqpe,ssqvir,ssqtke,ssqrke,ssqte,ssqprs
  real(8) :: ssqttp,ssqrtp,ssqtmp,ssqindpe,ssq3bpe,ssqtotke

  real(8) :: avvdwpe,avcpe,avindpe,av3bpe,avtotke
  real(8) :: avpe,avvir,avtke,avrke,avte,avprs,avttp,avrtp,avtmp  
  real(8) :: flcvdwpe,flccpe,flcindpe,flc3bpe,flctotke 
  real(8) :: flcpe,flcvir,flctke,flcrke,flcte,flcprs,flcttp,flcrtp,flctmp

  real(8) :: t, tr ! step time

  ! sapt5s,  ccpol5 or ccpol8s parameters variables
  real(8) :: beta_ccpol(nsite,nsite),alp_ccpol(nsite,nsite),a1_ccpol(nsite,nsite),a2_ccpol(nsite,nsite)
  real(8) :: a3_ccpol(nsite,nsite),a12_ccpol(nsite,nsite),a0_ccpol(nsite,nsite),ae_ccpol(nsite,nsite)
  real(8) :: c0_ccpol(nsite,nsite),c1_ccpol(nsite,nsite),c2_ccpol(nsite,nsite),c3_ccpol(nsite,nsite)
  real(8) :: c6_ccpol(nsite,nsite),c8_ccpol(nsite,nsite),c10_ccpol(nsite,nsite),d1_ccpol(nsite,nsite)
  real(8) :: d6_ccpol(nsite,nsite),d8_ccpol(nsite,nsite),d10_ccpol(nsite,nsite),dp_ccpol(nsite,nsite)
  real(8) :: domo_ccpol(nsite,nsite),expalp_ccpol(nsite,nsite)

  integer :: styp(nsite),stypa(npimax),stypb(npimax),index,nint0,nminter

  real(8) :: rBuf(50),sptype(1000,1000),interdist0(nsite,nsite),aMinter(2,nsite,nsite,nsite,nsite)
  integer :: styp1,styp2,styp3,styp4

  ! induction model parameter

  logical :: induction      ! true or false (wants to include induction model) 
  character(len=20) :: induction_type         ! (damped or undamped)
  character(len=20) :: ind_model              ! physical dipole or point dipole  
  character(len=20) :: method       ! method for long range correctin to avoid total energy drift 

  real(8) :: E0x(nsite,nom),E0y(nsite,nom),E0z(nsite,nom)    ! efield of static charges
  real(8) :: Eidmx(nsite,nom),Eidmy(nsite,nom),Eidmz(nsite,nom) ! efld of ind dipoles
  real(8) :: Etx(nsite,nom),Ety(nsite,nom),Etz(nsite,nom)

  real(8) :: idmx(nsite,nom),idmy(nsite,nom),idmz(nsite,nom) ! induced dipoles

  real(8) :: polarizability(nsite)
  real(8) :: apol(nsite,nom)
  integer :: npol
 
end module common_arrays
