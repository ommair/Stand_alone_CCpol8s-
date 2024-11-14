module initialization

  use, intrinsic :: iso_fortran_env
  use common_arrays
  use utility
  use molecular_sites
  use vdw_ccpol8s_omo
  use coulomb
  use nb_induction_model 

  implicit none  

  private
  public :: control_parameters, molecule, shift_to_com
  public :: principle_moment_inertia,fcc_positions
  public :: ini_quaternions,R_Euler_to_xyz
contains

  subroutine control_parameters

    implicit none
  
    open(unit=10, file='CONTROL', form='formatted')
 
    read(10,*) pot_name   ! name of potential required 
    read(10,*) initial_config ! selects initial config (fcc or random) 
    read(10,*) induction_type ! damped or undamped induction model
    read(10,*) ind_model      ! physical dipole or mathematical dipole

  close(10)

  end subroutine control_parameters

  subroutine molecule

    implicit none

    integer :: i,j,is,js,k,ics 
    integer :: is3,nlin0,ii3

    if (induction_type .eq. 'undamped_ind') then

       open(unit=20, file='FIELD_ccpol8s_omo',form='formatted')

    else if (induction_type .eq. 'damped_ind') then

       open(unit=20, file='FIELD_ccdpol8s_omo',form='formatted')

    end if 

    read(20,*) system_name
    read(20,*) nm  ! number of molecules
    read(20,*) ns ! number of sites in a molecule
    read(20,*) npol ! number of pol centers in a molecule
    read(20,*)

    do is=1,ns
       read(20,*)an(is),(coords(j,is),j=1,3),massites(is),charge(is),polarizability(is),styp(is)
    end do

    read(20,*)
    read(20,*) npi
    read(20,*)

    if (pot_name.eq.'ccpol8s_omo') then

       do k=1,npi
         read(20,*) stypa(k),stypb(k),an1(k),an2(k),rBuf(1:7)

         index = 1
         beta_ccpol(stypa(k),stypb(k)) = rBuf(index)
         beta_ccpol(stypb(k),stypa(k)) = beta_ccpol(stypa(k),stypb(k))

         index = index + 1
         c0_ccpol(stypa(k),stypb(k)) = rBuf(index)*418.4d0
         c0_ccpol(stypb(k),stypa(k)) = c0_ccpol(stypa(k),stypb(k))

         index = index + 1
         c1_ccpol(stypa(k),stypb(k)) = rBuf(index)*418.4d0
         c1_ccpol(stypb(k),stypa(k)) = c1_ccpol(stypa(k),stypb(k))

         index = index + 1
         c2_ccpol(stypa(k),stypb(k)) = rBuf(index)*418.4d0
         c2_ccpol(stypb(k),stypa(k)) = c2_ccpol(stypa(k),stypb(k))

         index = index + 1
         c6_ccpol(stypa(k),stypb(k)) = rBuf(index)*418.4d0
         c6_ccpol(stypb(k),stypa(k)) = c6_ccpol(stypa(k),stypb(k))

         index = index + 1
         c8_ccpol(stypa(k),stypb(k)) = rBuf(index)*418.4d0
         c8_ccpol(stypb(k),stypa(k)) = c8_ccpol(stypa(k),stypb(k))

         index = index + 1
         domo_ccpol(stypa(k),stypb(k)) = rBuf(index)
         domo_ccpol(stypb(k),stypa(k)) = domo_ccpol(stypa(k),stypb(k))

      end do 
 
    else
       stop
       write(100,*) 'no potential is selected'
    end if

    close(20)


    totm = 0.d0     ! initaiting total mass of the molecule
    natom = 0      ! initiating number of massive centers

    do i=1,ns
       if (massites(i).ne.0.d0) then
         natom = natom + 1            ! total number of massive sites
         totm  = totm + massites(i)    ! total mass of a molecule
       end if
    end do

!    write(*,*) totm

  end subroutine molecule

  subroutine shift_to_com

    implicit none

    integer :: i,j,is,js
    real(8) :: rx, ry, rz  ! com of coords

    rx = 0.d0
    ry = 0.d0
    rz = 0.d0

    do j=1, ns  ! loop on sites 
          rx = rx + massites(j)*coords(1,j)/totm
          ry = ry + massites(j)*coords(2,j)/totm
          rz = rz + massites(j)*coords(3,j)/totm
    end do

    ! sites coordinates relative to com coords

    do j=1, ns
       coords(1,j) = coords(1,j) - rx
       coords(2,j) = coords(2,j) - ry
       coords(3,j) = coords(3,j) - rz
    end do

  end subroutine shift_to_com

  subroutine principle_moment_inertia

    implicit none

    integer :: i,j,is,js
    real(8) :: rmi(3,3), rot(3,3), size
    real(8) :: rxx,ryy,rzz
    real(8) :: rotxyz,aa1

    do i=1,3
       do j=1,3
          rmi(i,j) = 0.d0   ! zero all components
       end do
    end do


   do i=1,3
       do j=1,3
          do is=1, ns
             rmi(i,j) = rmi(i,j) - massites(is)*coords(i,is)*coords(j,is)
             if (i.eq.j) rmi(i,j) = rmi(i,j) + massites(is)* & 
                         (coords(1,is)**2+coords(2,is)**2+coords(3,is)**2 )
          end do
       end do
    end do

   call jacobi(rmi,rot,3)  

   rmi(1,2)=rmi(2,1)
   rmi(1,3)=rmi(3,1)
   rmi(2,3)=rmi(3,2) 

   size=sqrt(rmi(1,1)**2+rmi(2,2)**2+rmi(3,3)**2)/3.0
      if(rmi(1,1)/size.lt.1.0e-8.or.   &
        rmi(2,2)/size.lt.1.0e-8.or.    &
        rmi(3,3)/size.lt.1.0e-8) then 
        stop
        write(100,*)'molecule off-digonal inertia error '
      end if
!  assigning principle moment of inertia

   pmi(1) = rmi(1,1)
   pmi(2) = rmi(2,2)
   pmi(3) = rmi(3,3) 
  
   do is=1,ns
      rxx = coords(1,is)
      ryy = coords(2,is)
      rzz = coords(3,is)

!  multiplication of rotation matrix with site corrdinates to transform into
!  principle corrdinates

     coords(1,is) = rxx*rot(1,1)+ryy*rot(2,1)+rzz*rot(3,1) 
     coords(2,is) = rxx*rot(1,2)+ryy*rot(2,2)+rzz*rot(3,2)
     coords(3,is) = rxx*rot(1,3)+ryy*rot(2,3)+rzz*rot(3,3)

   end do 

!   writee(*,*) pmi(1),pmi(2),pmi(3) 

   do is=1,ns

     pcoord(1,is) = coords(1,is)
     pcoord(2,is) = coords(2,is)
     pcoord(3,is) = coords(3,is) 

   end do  

  
  end subroutine principle_moment_inertia
  
  subroutine fcc_positions

    implicit none

    integer :: i,j,is,js
    integer :: ix,iy,iz,aa
    integer :: nc ! ! number of fcc unit cells in each coordinate direction; n=4*nc**3
    real(8) :: ndens ! number density =  na* mass_density / molarmass[pera^3]
    real(8) :: cell ! ! unit cell
    real(8) :: r(3,nm)   ! components and no of particles
    real(8), parameter :: mdens = 1.d0 ! mass density       ! g/cm^3
    real(8), parameter :: na = 6.02214076e+23 ! per mol

    ! ! sets up the fcc lattice: four molecules per unit cell
    real, dimension(3,4), parameter :: r_fcc = reshape ( [ &
         & 0.25, 0.25, 0.25, &
         & 0.25, 0.75, 0.75, &
         & 0.75, 0.75, 0.25, &
         & 0.75, 0.25, 0.75 ], [3,4] ) ! positions in unit cell
    ndens = (na * mdens)/totm
    ndens = ndens* 1.e-24
    box  = ( real(nm) / ndens ) ** ( 1.d0/3.d0 )
    nc = nint(real(nm/4)**(1.0/3.0) )
    !write(*,*) nm,4*nc**3    
    if (nm .ne. 4*nc**3) then
        write(*,*)' invalid number of molecules for fcc lattice '
        stop
    end if
    cell = box / real(nc) ! unit cell
    hbox = box / 2.0      ! half box length
    i = 0
    do iz = 0, nc-1  ! begin triple loop over unit cell indices
         do iy = 0, nc-1
            do ix = 0, nc-1
               do aa=1,4  ! begin loop over atoms in unit cell
                  i=i+1
                  r(1:3,i) = r_fcc(:,aa) +real([ix,iy,iz])  ! in range 0..real(nc)
                  r(1:3,i) = r(1:3,i) * cell            ! in range 0..real(nc)
               end do   ! end loop over atoms in unit cell
            end do
         end do
    end do    ! end triple loop over unit cell indices
   
    do i=1,nm
         x(i) = r(1,i)
         y(i) = r(2,i)
         z(i) = r(3,i)
    end do

  end subroutine fcc_positions
      
  subroutine ini_quaternions
    
    implicit none
   
    integer :: i,j,is,js
!    real(8) :: theta,phi,psi
    real(8) :: conorm 
    integer :: seed
    real(8) :: qrr(3) 
   
    seed = 96873683

  ! set quaternions to give random orientiation
    call random_seed(seed)

    do j=1,nm   ! loop over molecules

       call random_number(q0(j))
       call random_number(q1(j))
       call random_number(q2(j))
       call random_number(q3(j))

       conorm=sqrt(q0(j)**2+q1(j)**2+q2(j)**2+q3(j)**2)
       q0(j)=q0(j)/conorm
       q1(j)=q1(j)/conorm
       q2(j)=q2(j)/conorm
       q3(j)=q3(j)/conorm

    end do  

   end subroutine ini_quaternions


   subroutine R_Euler_to_xyz

      implicit none
 
      integer, parameter :: nmax=1000
      integer, parameter :: maxlin=2000,maxpar=1000,maxdat=50000
! COM coordinates (bohr) and z-y-z Euler angles (in LAB system, degree)
      real(8) :: R(3,nmax), E(3,nmax)
      real(8) :: ElA(3),ElB(3)
      real(8) :: RA(3),RB(3)
      integer, parameter :: nsite=26,nmol=2
      integer, parameter :: nsites=2*nsite*maxdat
      integer, parameter :: Nmax8=Nmax*nsite
      real(8) :: sites(3,nsite)
      real(8) :: sitesAB(3,nsites)
      real(8) :: T(3,3),vec(3)
      real(8) :: coorN(3,Nmax8)

      integer :: i,ndat,nd,indB,indA
      real(8) :: pol,sig,plen,dmpfct,zero,pi180,dist

      data pol /9.6278d0/         ! bohr^3
      data sig /0.367911875040999981d0/ ! position of polarizable center
      data plen /1.1216873242d0/
      data dmpfct /1.d0/
      data zero /0.d0/
!      data bohr2a /0.529177249d0/
!      data h2kcal /627.510d0/

!     cartesian coordinates of the sites (in bohr) CCpol8s ref geometry

      data sites / &
      0.d0, 0.d0, 0.1255334885d0, &  
     -1.45365196228170d0, 0.d0,-0.9961538357d0, &  
      1.45365196228170d0, 0.d0,-0.9961538357d0, &  
      0.d0,-0.2067213d0,-0.2462589474d0, &  
      0.d0, 0.2067213d0,-0.2462589474d0, &  
      0.90d0,  0.053d0,  0.64d0,  &   
     -0.90d0,  0.053d0,  0.64d0,  &   
      0.90d0, -0.053d0,  0.64d0,  &   
     -0.90d0, -0.053d0,  0.64d0,  &   
      1.29d0,  0.14d0,  -0.91d0,  &   
     -1.29d0,  0.14d0,  -0.91d0,  &   
      1.29d0, -0.14d0,  -0.91d0,  &   
     -1.29d0, -0.14d0,  -0.91d0,  &   
      0.56d0,  0.24d0,  -0.48d0,  &   
     -0.56d0,  0.24d0,  -0.48d0,  &   
      0.56d0, -0.24d0,  -0.48d0,  &   
     -0.56d0, -0.24d0,  -0.48d0,  &   
      0.91d0,  0.21d0,  -0.68d0,  &   
     -0.91d0,  0.21d0,  -0.68d0,  &   
      0.91d0, -0.21d0,  -0.68d0,  &   
     -0.91d0, -0.21d0,  -0.68d0,  &   
      1.48d0,  0.26d0,  -0.62d0,  &   
     -1.48d0,  0.26d0,  -0.62d0,  &   
      1.48d0, -0.26d0,  -0.62d0,  &   
     -1.48d0, -0.26d0,  -0.62d0,  &   
      0.d0,0.d0,-0.24237838654099997d0/ ! pol site

      pi180=dacos(-1.d0)/180.d0

      open(unit=130, file='input_geo_xyz.dat')
!      open(unit=140, file='optim.dat.613.txt')
      read (*,*) ndat
      write(130,*) ndat
      write(130,*) nmol
      write(130,*) nsite
      do nd=1,ndat
!         read (*,*) nd
         read (*,*) dist, ElA(2),ElA(3),ElB(1),ElB(2),ElB(3)   ! bohr,degrees
!         write(*,*) dist, ElA(2),ElA(3),ElB(1),ElB(2),ElB(3)
         ! switch to radians
         ElA(1)=0.d0
         do i=1,3
         ElA(i)=ElA(i)*pi180
         ElB(i)=ElB(i)*pi180
         end do

         RA(1)=0.d0
         RA(2)=0.d0
         RA(3)=0.d0
         RB(1)=0.d0
         RB(2)=0.d0
         RB(3)=dist/bohr2a     

! calculate site positions (in bohr) after rotations and translations

      call rotmat_molecule (ElA(1),ElA(2),ElA(3), T)
      indA = (nd-1)*nsite*2
      indB = (nd-1)*nsite*2 + nsite
      do ns=1,nsite
         call Av (T,sites(1,ns), vec)
         sitesAB(1,indA+ns)=vec(1)+RA(1)
         sitesAB(2,indA+ns)=vec(2)+RA(2)
         sitesAB(3,indA+ns)=vec(3)+RA(3)
      end do
      call rotmat_molecule (ElB(1),ElB(2),ElB(3), T)
      do ns=1,nsite
         call Av(T,sites(1,ns), vec)
         sitesAB(1,indB+ns)=vec(1)+RB(1)
         sitesAB(2,indB+ns)=vec(2)+RB(2)
         sitesAB(3,indB+ns)=vec(3)+RB(3)
      end do

!!      write(*,*) nd,' dimer'
!      write(*,*) 'molA'
      do ns=1,nsite
      write(130,*)ns,sitesAB(1,indA+ns)*bohr2a,sitesAB(2,indA+ns)*bohr2a, &
                  sitesAB(3,indA+ns)*bohr2a
      end do

!!      write(*,*) 'molB'
      do ns=1,nsite
      write(130,*)ns,sitesAB(1,indB+ns)*bohr2a,sitesAB(2,indB+ns)*bohr2a, &
                  sitesAB(3,indB+ns)*bohr2a
      end do

      end do   ! loop ends on cluster

!      close(140)
      close(130)

!------------------------------------------------------------------------
   end subroutine R_Euler_to_xyz

      subroutine Av (A,v,b)
! calculates A*v=b
      implicit none
      integer :: i,j
      real(8) :: A(3,3),v(3),b(3)
      real(8) :: zero,sum
      data zero /0.d0/

      do 10 i=1,3
      sum=zero
      do 20 j=1,3
      sum=sum +a(i,j)*v(j)
 20   continue
      b(i)=sum
 10   continue

      end

        subroutine rotmat_molecule(ad,bd,gd,u)
!
! Construct the rotation matrix corresponding to
! the Euler angles ad,bd,gd (in radians).
!
        implicit real*8 (a-h,o-z)
        dimension u(3,3)
        sad = dsin(ad)
        sbd = dsin(bd)
        sgd = dsin(gd)
        cad = dcos(ad)
        cbd = dcos(bd)
        cgd = dcos(gd)
!----- construct the transformation matrix
        u(1,1) = cad*cbd*cgd - sad*sgd
        u(1,2) = -cad*cbd*sgd - sad*cgd
        u(1,3) = cad*sbd
        u(2,1) = sad*cbd*cgd + cad*sgd
        u(2,2) = -sad*cbd*sgd + cad*cgd
        u(2,3) = sad*sbd
        u(3,1) = -sbd*cgd
        u(3,2) = sbd*sgd
        u(3,3) = cbd
        return
        end

!----------------------------------------------------------
      subroutine xycopy (n,x,y)
      real*8 x(*),y(*)
      integer :: i,n

      do 10 i=1,n
 10   y(i)=x(i)
      end

end module initialization
