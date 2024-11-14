module  nb_induction_model

  use, intrinsic :: iso_fortran_env
  use common_arrays
  
  implicit none
  private


  public :: nb_induction
  public :: iteridm_mathematical_dipole
  public :: efield_stat,efield_idm


contains

  subroutine nb_induction(qindpe)

    implicit none

    integer :: is,js,i,j,k,isteps
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi
    real(8) :: qindpe,point_qindpe,phys_qindpe    

    if (ind_model .eq. 'mathematical_dipole') then
!   for mathematical dipole model
    call efield_stat ! efield of permanent charges

    call iteridm_mathematical_dipole(qindpe)  ! iterate to get self consistent indued dipoles

    else if (ind_model .eq. 'physical_dipole') then
!   for physical dipole model
    call iteridm_physical_dipole(qindpe)

    end if

  end subroutine nb_induction

! *********** MATHEMATICAL DIPOLE MODEL *************************  
  subroutine iteridm_mathematical_dipole(qindpe)

    implicit none

    integer :: is,js,i,j,k,isteps
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: efldxi,efldyi,efldzi,efldxj,efldyj,efldzj
    real(8) :: xpolc,ypolc,zpolc !separation between dipole centers
    real(8) :: dist2,rccsq,dist,doti,dotj
    real(8) :: qis,qjs  ! charges on sites
    real(8) :: boxi,xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i
    real(8) :: polr,d3ii,d5ii,expar2,selffac
    real(8) :: thr_iter,change
    real(8) :: energy,qindpe,engtst
    real(8), parameter :: maxit = 500
    real(8) :: polE1x,polE1y,polE1z

    do i=1,nm
       do is=1,ns

       idmx(is,i) = 0.d0 ! initiating induced dipoles
       idmy(is,i) = 0.d0
       idmz(is,i) = 0.d0

       Etx(is,i) = 0.d0 ! initiating induced dipoles
       Ety(is,i) = 0.d0
       Etz(is,i) = 0.d0

       end do
    end do

!**************************************************************************
! Now iterate to calculate converged polarization energy        starts here

    thr_iter = 1.d-20
    change = 10.d0
    isteps = 0

    do while (change.gt.thr_iter.and.isteps.lt.maxit)  ! while loop starts

    energy = 0.0d0
    change = 0.d0 

    call efield_idm 

    do i=1,nm    ! loop over ith molecule starts

       do is=ns-npol+1,ns 

          Etx(is,i) =  Eidmx(is,i) + E0x(is,i) ! total efield
          Ety(is,i) =  Eidmy(is,i) + E0y(is,i) 
          Etz(is,i) =  Eidmz(is,i) + E0z(is,i) 

          polr = apol(is,i) !polarizability(is)

          polE1x = polr* Etx(is,i) 
          polE1y = polr* Ety(is,i) 
          polE1z = polr* Etz(is,i)   
 
          change = (idmx(is,i)-polE1x)**2 + &
                   (idmy(is,i)-polE1y)**2 + &
                   (idmz(is,i)-polE1z)**2 + change

          idmx(is,i) = polE1x
          idmy(is,i) = polE1y
          idmz(is,i) = polE1z

          energy = -0.5d0*polr*(Etx(is,i)*E0x(is,i) + &
                                Ety(is,i)*E0y(is,i) + &
                                Etz(is,i)*E0z(is,i)) +  energy

         end do   ! loop over is 

      end do  ! loop over ith molecule ends  

      isteps = isteps + 1

    end do       ! while loop ends

    if (isteps.ge.maxit) then
        write (*,*) 'No convergence in indN_iter'
        write (*,'(a,g12.3)') 'energy change=',change
        write (*,'(a,g12.3)') 'thr_iter=',thr_iter
        stop
    end if

    qindpe = r4pie0*energy/418.4d0

!    write(*,*) qindpe

  end subroutine iteridm_mathematical_dipole

  subroutine efield_stat

    implicit none

    integer :: is,js,i,j,k,isteps
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: efldx,efldy,efldz
    real(8) :: dist2,rccsq,dist,doti,dotj,sigcs
    real(8) :: qis,qjs  ! charges on sites
    real(8) :: boxi,xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i
    real(8) :: polr,f2,f3  

    real(8) :: qindpe
    real(8) :: expdr,detr

    do i=1,nm
       do is=1,ns

       E0x(is,i)= 0.d0 
       E0y(is,i)= 0.d0
       E0z(is,i)= 0.d0

!       xxs(is,i) = xs(is,i)
!       yys(is,i) = ys(is,i)
!       zzs(is,i) = zs(is,i)

       end do
    end do

! ************* Calculating Electric from permanant charges ****************
! starts here

    do i=1,nm       ! loop over molecule i   starts

        do j=1,nm    ! loop over molecule j   starts

           if (i.ne.j) then 
 
             do is=ns-npol+1,ns

                do js=1,ns

                   qis = charge(is)
                   qjs = charge(js)

                   ddx = xxst(is,i) - xxst(js,j)
                   ddy = yyst(is,i) - yyst(js,j)
                   ddz = zzst(is,i) - zzst(js,j) 

                   dist2 = ddx**2 + ddy**2 + ddz**2

                   dist = sqrt(dist2)
                   d1i = 1.d0/dist
                   d2i = d1i**2
                   d3i = d2i*d1i

                   f2 = 1.d0

                   if (induction_type .eq. 'damped_ind') then

                   sigcs = 0.59d0**2.d0 

                   f2 = (erf(dist/(sqrt(2.d0*sigcs)))-twosqpi*dist/(sqrt(2.d0*sigcs))*exp(-dist2/(2.d0*sigcs)))  
 
                   else 
 
                   f2 = 1.d0
               
                   end if
 
                   ! components of efield in internal units (multiply
                   ! r4pie0)! Ei= qi* vecr_rij/rij^3

                   E0x(is,i)= E0x(is,i) + f2*qjs*ddx*d3i
                   E0y(is,i)= E0y(is,i) + f2*qjs*ddy*d3i
                   E0z(is,i)= E0z(is,i) + f2*qjs*ddz*d3i     

                end do

             end do

           end if     
 
        end do        ! loop over molecule j ends

    end do           ! loop over molecule i ends

! ************* Calculating Electric from permanant charges ends ***********
 
  end subroutine efield_stat


! *************** Short Range Damped Induction terms ***************************

  subroutine efield_idm

    implicit none

    integer :: is,js,i,j,k,isteps
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: efldx,efldy,efldz
    real(8) :: dist2,rccsq,dist,doti,dotj
    real(8) :: qis,qjs  ! charges on sites
    real(8) :: boxi,xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i
    real(8) :: polr,f2,f3,sigcs  

    real(8) :: expdr,detr

    real(8) :: qindpe

    do i=1,nm
       do is=1,ns

       Eidmx(is,i)= 0.d0
       Eidmy(is,i)= 0.d0
       Eidmz(is,i)= 0.d0

!       xxs(is,i) = xs(is,i)
!       yys(is,i) = ys(is,i)
!       zzs(is,i) = zs(is,i)

       end do
    end do

!**************************************************************************
! Now Calculating Electric from point induced dipoles

    do i=1,nm   ! loop over ith molecule starts

       do j=1,nm   ! loop over jth molecule starts

          if (i.ne.j) then
 
             do is=ns-npol+1,ns

                do js=ns-npol+1,ns    ! loop over sites of molecule j  Starts

                   ddx = xxst(is,i) - xxst(js,j)
                   ddy = yyst(is,i) - yyst(js,j)
                   ddz = zzst(is,i) - zzst(js,j)

                   dist2 = ddx**2 + ddy**2 + ddz**2

                   dist = sqrt(dist2)
                   d1i = 1.d0/dist
                   d2i = d1i**2
                   d3i = d2i*d1i
                   d5i = d3i*d2i
                  
                   f2 = 1.d0
                   f3 = 1.d0

                   if (induction_type .eq. 'damped_ind') then

                   sigcs = 0.59d0**2.d0
 
                   f2 =erf(dist/(sqrt(4.d0*sigcs)))-(twosqpi*dist/(sqrt(4.d0*sigcs)) + &
                        (2.d0/3.d0)*twosqpi*dist**3.d0/(sqrt(4.d0*sigcs))**3.d0)*exp(-dist2/(4.d0*sigcs))  
                             
                   f3 =(erf(dist/(sqrt(4.d0*sigcs)))-twosqpi*dist/(sqrt(4.d0*sigcs))*exp(-dist2/(4.d0*sigcs)))
               
                   else

                   f2 = 1.d0
                   f3 = 1.d0

                   end if

                      ! dot product of induced dipole and separation vector
                   doti = idmx(is,i)*ddx+idmy(is,i)*ddy+idmz(is,i)*ddz
                   dotj = idmx(js,j)*ddx+idmy(js,j)*ddy+idmz(js,j)*ddz  
     
                   ! calculate the  efield of inducd dipole of molecule j at
                   ! induced dipole of molecule i

                   efldx = f2*3.d0*d5i*dotj*ddx - f3*idmx(js,j)*d3i
                   efldy = f2*3.d0*d5i*dotj*ddy - f3*idmy(js,j)*d3i
                   efldz = f2*3.d0*d5i*dotj*ddz - f3*idmz(js,j)*d3i

                   ! calculate total efield at induced dipole at molecule i
                   ! because of parmenant charges and induced dipoles at
                   ! molecule j

                   Eidmx(is,i)= Eidmx(is,i) + efldx
                   Eidmy(is,i)= Eidmy(is,i) + efldy
                   Eidmz(is,i)= Eidmz(is,i) + efldz 

                   end do     ! loop over sites js  ends

                end do     ! loop over sites is  ends

          end if

       end do      ! loop over jth molecule ends

    end do  ! loop over ith molecule ends


  end subroutine efield_idm


! *********** PHYSICSL DIPOLE MODEL *********************
  subroutine iteridm_physical_dipole(qindpe) ! or qindpe

    implicit none

    integer :: is,js,i,j,k,isteps
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    !real(8) :: xxst(nsite,nom),yyst(nsite,nom),zzst(nsite,nom) ! for intermediate saving
    !real(8) :: Etx(nsite,nom),Ety(nsite,nom),Etz(nsite,nom)  ! total efield
    real(8) :: Ex,Ey,Ez      ! tot efield
    real(8) :: efldx,efldy,efldz
    real(8) :: dist2,rccsq,dist,f2,rhois,rhojs
    real(8) :: qis,qjs  ! charges on sites
    real(8) :: boxi,xcc,ycc,zcc,dx,dy,dz,ddx,ddy,ddz,d1i,d2i,d3i,d5i,d7i
    real(8) :: thr_iter,change
    real(8) :: energy,qindpe,engtst,engtst2,qqindpe
    real(8) :: elst1,elst2,sigcs
    real(8) :: fdp(nsite,nsite)
    real(8) :: dipolesize(nom)
    real(8) :: separation(nom)
    real(8) :: qc(nsite,nom),qs(nsite,nom)

    real(8) :: xcor(nsite,nom), ycor(nsite,nom), zcor(nsite,nom) ! core corrds
    real(8) :: xshl(nsite,nom), yshl(nsite,nom), zshl(nsite,nom) ! shell corrds
    real(8) :: xshlt(nsite,nom), yshlt(nsite,nom), zshlt(nsite,nom)
 
    integer, parameter :: maxit = 100

    if (induction_type .eq. 'damped_ind') then

       sigcs = 0.59d0**2.d0  ! A^-1

    else

       sigcs=0.d0

    end if

    do i=1,nm
       do is=1,ns

       E0x(is,i)= 0.d0
       E0y(is,i)= 0.d0
       E0z(is,i)= 0.d0

       qc(is,i)= 0.d0
       qs(is,i)= 0.d0

       xcor(is,i) = 0.d0
       ycor(is,i) = 0.d0
       zcor(is,i) = 0.d0

       xshl(is,i) = 0.d0
       yshl(is,i) = 0.d0
       zshl(is,i) = 0.d0

       end do
    end do


    do i=1,nm
       qc(26,i)=-8.d0
       qs(26,i)=8.d0

       xcor(26,i) = xxst(26,i)
       ycor(26,i) = yyst(26,i)
       zcor(26,i) = zzst(26,i)

       xshl(26,i) = xxst(26,i)
       yshl(26,i) = yyst(26,i)
       zshl(26,i) = zzst(26,i)

    end do

    ! electric field due t permanent charges of the system

    do i=1,nm ! loop over ith molecule

       do j=1,nm   ! loop over jth molecule

          if (i.ne.j) then

             do is=ns-npol+1,ns   ! core-shell sites of molecules i

                do js=1,ns-npol ! loop over sites of molecule j

                   ddx = xxst(is,i) - xxst(js,j)
                   ddy = yyst(is,i) - yyst(js,j)
                   ddz = zzst(is,i) - zzst(js,j)

                   dist2 = ddx**2 + ddy**2 + ddz**2

                   dist = sqrt(dist2)
                   d1i = 1.d0/dist
                   d2i = d1i**2
                   d3i = d2i*d1i
                   d5i = d3i*d2i

                   if (induction_type .eq. 'damped_ind') then
                       sigcs = (0.59d0)**2.d0  ! A^-1
                    f2 = (erf(dist/(sqrt(2.d0*sigcs)))-twosqpi*dist/(sqrt(2.d0*sigcs))*exp(-dist2/(2.d0*sigcs)))
                   else
                    f2=1.d0
                   end if

                   qjs = chgs(js,j) !*erf(dist/(sqrt(2.d0*sigcs))) 

                   efldx = f2*r4pie0*qjs*ddx*d3i
                   efldy = f2*r4pie0*qjs*ddy*d3i
                   efldz = f2*r4pie0*qjs*ddz*d3i

                   E0x(is,i) = E0x(is,i) + efldx
                   E0y(is,i) = E0y(is,i) + efldy
                   E0z(is,i) = E0z(is,i) + efldz

                end do    
             end do

          end if 
       end do

    end do


    ! iterate electric field on the core positions to adjust shells

    thr_iter = 1.d-10
    change = 10.d0
    isteps = 0

    do while (isteps.lt.maxit)   ! while loop starts here

    change = 0.d0

    do i=1,nm       ! loop over molecule i   starts

        do is=ns-npol+1,ns    ! loop over molecule j   starts

           Etx(is,i)= 0.d0
           Ety(is,i)= 0.d0
           Etz(is,i)= 0.d0

           do j=1,nm

           if (i.ne.j) then

                ! Efield due to permanent charges
                do js=1,ns-npol

                   ddx = xxst(is,i)-xxst(js,j)
                   ddy = yyst(is,i)-yyst(js,j)
                   ddz = zzst(is,i)-zzst(js,j)

                   dist2 = ddx**2 + ddy**2 + ddz**2

                   dist = sqrt(dist2)
                   d1i = 1.d0/dist
                   d2i = d1i**2
                   d3i = d2i*d1i 

                   if (induction_type .eq. 'damped_ind') then
                       sigcs = (0.59d0)**2.d0  ! A^-1
                   f2 = (erf(dist/(sqrt(2.d0*sigcs)))-twosqpi*dist/(sqrt(2.d0*sigcs))*exp(-dist2/(2.d0*sigcs)))
                   else
                    f2=1.d0
                   end if

                   qjs = charge(js) !*erf(dist/(sqrt(2.d0*sigcs)))

                   Etx(is,i)= Etx(is,i) + f2*r4pie0*qjs*ddx*d3i 
                   Ety(is,i)= Ety(is,i) + f2*r4pie0*qjs*ddy*d3i 
                   Etz(is,i)= Etz(is,i) + f2*r4pie0*qjs*ddz*d3i

                end do

                ! Efield due to core charges
                do js=ns-npol+1,ns

                   ddx = xxst(is,i)-xcor(js,j)
                   ddy = yyst(is,i)-ycor(js,j)
                   ddz = zzst(is,i)-zcor(js,j)

                   dist2 = ddx**2 + ddy**2 + ddz**2

                   dist = sqrt(dist2)
                   d1i = 1.d0/dist
                   d2i = d1i**2
                   d3i = d2i*d1i

                   if (induction_type .eq. 'damped_ind') then
                       sigcs = (0.59d0)**2.d0  ! A^-1
                   f2 = (erf(dist/(sqrt(4.d0*sigcs)))-twosqpi*dist/(sqrt(4.d0*sigcs))*exp(-dist2/(4.d0*sigcs)))
                   else
                    f2=1.d0
                   end if

                   rhojs = qc(js,j) !*erf(dist/(sqrt(4.d0*sigcs)))

                   ! components of efield in internal units (multiply
                   ! r4pie0)! Ei= qi* vecr_rij/rij^3

                   Etx(is,i)= Etx(is,i) + f2*r4pie0*rhojs*ddx*d3i
                   Ety(is,i)= Ety(is,i) + f2*r4pie0*rhojs*ddy*d3i
                   Etz(is,i)= Etz(is,i) + f2*r4pie0*rhojs*ddz*d3i

                end do

                ! Efield due to shell charges
                do js=ns-npol+1,ns

                   ddx = xxst(is,i)-xshl(js,j)
                   ddy = yyst(is,i)-yshl(js,j)
                   ddz = zzst(is,i)-zshl(js,j)

                   dist2 = ddx**2 + ddy**2 + ddz**2

                   dist = sqrt(dist2)
                   d1i = 1.d0/dist
                   d2i = d1i**2
                   d3i = d2i*d1i

                   if (induction_type .eq. 'damped_ind') then
                       sigcs = (0.59d0)**2.d0  ! A^-1
                   f2 = (erf(dist/(sqrt(4.d0*sigcs)))-twosqpi*dist/(sqrt(4.d0*sigcs))*exp(-dist2/(4.d0*sigcs)))
                   else
                    f2=1.d0
                   end if

                   rhojs = qs(js,j) !*erf(dist/(sqrt(4.d0*sigcs)))

                   ! components of efield in internal units (multiply
                   ! r4pie0)! Ei= qi* vecr_rij/rij^3

                   Etx(is,i)= Etx(is,i) + f2*r4pie0*rhojs*ddx*d3i
                   Ety(is,i)= Ety(is,i) + f2*r4pie0*rhojs*ddy*d3i
                   Etz(is,i)= Etz(is,i) + f2*r4pie0*rhojs*ddz*d3i

                end do

           end if

           end do     ! loop over molecule j ends

        end do        ! loop over is ends

    end do           ! loop over molecule i ends


    do i=1,nm

       do is=ns-npol+1,ns

       xshl(is,i) = xcor(is,i) + apol(is,i)*Etx(is,i)/(qs(is,i)*r4pie0)
       yshl(is,i) = ycor(is,i) + apol(is,i)*Ety(is,i)/(qs(is,i)*r4pie0)
       zshl(is,i) = zcor(is,i) + apol(is,i)*Etz(is,i)/(qs(is,i)*r4pie0)

       end do
    end do    
  
    isteps = isteps + 1

!    write(*,*) isteps

    do i=1,nm
       do is=ns-npol+1,ns

          separation(i) = sqrt((xshl(is,i)-xcor(is,i))**2+ &
                         (yshl(is,i)-ycor(is,i))**2+ &
                         (zshl(is,i)-zcor(is,i))**2)

!          write(*,*) i,separation(i),separation(i)*qs(is,i)

       end do
    end do

    end do       ! while loop ends

! ************* Calculating Electric from permanant charges ends ***********

! ******************************************************************************
! e_ind = -0.5d0*mu*E0 = -0.5d0*(polarizability * E_total)

    engtst = 0.d0

    do i=1,nm
       do is=ns-npol+1,ns

          engtst = engtst - 0.5d0*apol(is,i)/r4pie0*(Etx(is,i)*E0x(is,i) + &
                                             Ety(is,i)*E0y(is,i) + &
                                             Etz(is,i)*E0z(is,i) )

       end do
    end do

    qqindpe = engtst/418.4d0

!   write(*,*) qqindpe

! ****************************************************************************
! ****************************************************************************

  ! calulating coulomb energy

      energy = 0.d0

      do i=1,nm-1       ! loop over molecule i   starts

          do j=i+1,nm

             do is=1,ns

                do js=1,ns

                   ddx = xxst(is,i) - xxst(js,j)
                   ddy = yyst(is,i) - yyst(js,j)
                   ddz = zzst(is,i) - zzst(js,j)

                   dist2 = ddx**2 + ddy**2 + ddz**2

                   rhois = qc(is,i)*(2.d0*pi*sigcs)**(-3/2)* &
                              exp(-dist2/(2.d0*sigcs))

                   rhojs = qc(js,j)*(2.d0*pi*sigcs)**(-3/2)* &
                              exp(-dist2/(2.d0*sigcs))

                   qis = chgs(is,i)
                   qjs = chgs(js,j)

                   dist = sqrt(dist2)
                   d1i = 1.d0/dist
                   d2i = d1i**2
                   d3i = d2i*d1i

                   !energy = energy + (qis+rhois)*(qjs+rhojs)*d1i*r4pie0
                   energy = energy + (qis*qjs+qis*qc(js,j)*erf(dist/(sqrt(2.d0*sigcs))) + &
                                      qc(is,i)*qjs*erf(dist/(sqrt(2.d0*sigcs))) + &
                                      qc(is,i)*qc(js,j)*erf(dist/(sqrt(4.d0*sigcs))))*d1i*r4pie0

                end do

             end do

             do is=1,ns

                do js=ns-npol+1,ns

                   ddx = xxst(is,i) - xshl(js,j)
                   ddy = yyst(is,i) - yshl(js,j)
                   ddz = zzst(is,i) - zshl(js,j)

                   dist2 = ddx**2 + ddy**2 + ddz**2

                   rhois = qc(is,i)*(2.d0*pi*sigcs)**(-3/2)* &
                              exp(-dist2/(2.d0*sigcs))

                   rhojs = qs(js,j)*(2.d0*pi*sigcs)**(-3/2)* &
                              exp(-dist2/(2.d0*sigcs))
                   qis = chgs(is,i)

                   dist = sqrt(dist2)
                   d1i = 1.d0/dist
                   d2i = d1i**2
                   d3i = d2i*d1i

!                   energy = energy + (qis+rhois)*rhojs*d1i*r4pie0
                   energy = energy + (qis*qs(js,j)*erf(dist/(sqrt(2.d0*sigcs))) + &
                                      qc(is,i)*qs(js,j)*erf(dist/(sqrt(4.d0*sigcs))))*d1i*r4pie0
                end do

             end do

             do is=ns-npol+1,ns

                do js=1,ns

                   ddx = xshl(is,i) - xxst(js,j)
                   ddy = yshl(is,i) - yyst(js,j)
                   ddz = zshl(is,i) - zzst(js,j)

                   dist2 = ddx**2 + ddy**2 + ddz**2

                   rhois = qs(is,i)*(2.d0*pi*sigcs)**(-3/2)* &
                              exp(-dist2/(2.d0*sigcs))

                   rhojs = qc(js,j)*(2.d0*pi*sigcs)**(-3/2)* &
                              exp(-dist2/(2.d0*sigcs))

                   qjs = chgs(js,j)

                   dist = sqrt(dist2)
                   d1i = 1.d0/dist
                   d2i = d1i**2
                   d3i = d2i*d1i

!                   energy = energy + rhois*(qjs+rhojs)*d1i*r4pie0
                   energy = energy + (qjs*qs(is,i)*erf(dist/(sqrt(2.d0*sigcs)))+ &
                                     qc(js,j)*qs(is,i)*erf(dist/(sqrt(4.d0*sigcs))))*d1i*r4pie0

                end do

             end do

             do is=ns-npol+1,ns

                do js=ns-npol+1,ns

                   ddx = xshl(is,i) - xshl(js,j)
                   ddy = yshl(is,i) - yshl(js,j)
                   ddz = zshl(is,i) - zshl(js,j)

                   dist2 = ddx**2 + ddy**2 + ddz**2

                   rhois = qs(is,i)*(2.d0*pi*sigcs)**(-3/2)* &
                              exp(-dist2/(2.d0*sigcs))

                   rhojs = qs(js,j)*(2.d0*pi*sigcs)**(-3/2)* &
                              exp(-dist2/(2.d0*sigcs))

                   dist = sqrt(dist2)
                   d1i = 1.d0/dist
                   d2i = d1i**2
                   d3i = d2i*d1i

!                   energy = energy + rhois*rhojs*d1i*r4pie0
                   energy = energy + (qs(js,j)*qs(is,i)*erf(dist/(sqrt(4.d0*sigcs))))*d1i*r4pie0

                end do

             end do

          end do

      end do

      ! dipole creation energy

      do i=1,nm
         do is=ns-npol+1,ns

         ddx = xxst(is,i) - xshl(is,i)
         ddy = yyst(is,i) - yshl(is,i)
         ddz = zzst(is,i) - zshl(is,i)

         dist2 = ddx**2 + ddy**2 + ddz**2

         rhois = qs(is,i) !*(2.d0*pi*sigcs)**(-3/2)* &
                          !   exp(-dist2/(2.d0*sigcs))

         dist = sqrt(dist2)

         energy = energy + 0.5d0*r4pie0*rhois**2.d0*dist2/apol(is,i)

         end do
      end do

    qindpe = energy/418.4d0 - rcpe

!    write(*,*) qindpe,qqindpe

  end subroutine iteridm_physical_dipole

end module nb_induction_model
