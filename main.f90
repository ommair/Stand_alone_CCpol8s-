program main

 use common_arrays
 use utility
 use initialization
 use molecular_sites
 use vdw_ccpol8s_omo
 use coulomb 
 use nb_induction_model

  implicit none

  integer :: i,j,is,js,k,m,ics,is3,nlin0,ii3,ii,iis,nd,ndimers

!  call R_Euler_to_xyz

  open(100, file='OUTPUT')

    call control_parameters     ! simulation parameters
    call molecule               ! reads molecular geometery(site'smasses,charges,coords)

    write(*,*) 'runninng water cluster with ', system_name, ' calculation using ',induction_type,' ',ind_model
    write(*,*)
    write(*,*) "        dimer","     VDW(kcal/mol) ","           COUL(kcal/mol) ", &
                 "           IND(kcal/mol)", "      CCpol8s'+NB(ind)(kcal/mol)"


    write(100,*) 'runninng water cluster with ', system_name, ' calculation using ',induction_type,' ',ind_model
    write(100,*)
    write(100,*) "        dimer","     VDW(kcal/mol) ","           COUL(kcal/mol) ", & 
                 "           IND(kcal/mol)", "      CCpol8s'+NB(ind)(kcal/mol)"

    if (initial_config .eq. 'fcc_config') then

       call shift_to_com           ! shifiting molecule to its com
       call principle_moment_inertia  ! calculaing principle moment of inertia

       write(*,*)
       write(*,*) '*****Reading fcc configuration XYZ in Angstroms*****'

       call fcc_positions
       call ini_quaternions ! assigns initial quaternions to com
       call sites           ! creates sites for each com using quaternions


       do i=1,nm
          do is=1,ns

             xxst(is,i) = x(i)+xs(is,i) ! site position in the box
             yyst(is,i) = y(i)+ys(is,i)
             zzst(is,i) = z(i)+zs(is,i)

          end do
       end do

    else if (initial_config .eq. 'user_config') then

         call R_Euler_to_xyz  ! converting R/Euler to xyz        

!         write(*,*)
!         write(*,*) '*****Reading User input configuration XYZ in Angstroms*****'
!         open(unit=30, file='config.inp')
         open(unit=30, file='input_geo_xyz.dat')
             read(30,*) ndimers ! no of dimers in the data
!             write(*,*) ndimers,' number of dimers'
!             write(*,*) '*****check OUTPUT files for energies*****'
             read(30,*) nm      ! no of molecules in a single dimer
             read(30,*) ns      ! no of sites in the molecule
             do nd=1,ndimers
                do i=1,nm
                   do is=1,ns

                      read(30,*)iis,xxst(is,i),yyst(is,i),zzst(is,i)

                   end do
                end do

                do i=1,nm
                   do is=1,ns

                      chgs(is,i)=charge(is)
                      apol(is,i)=polarizability(is)*bohr2a**(3.d0)

                   end do
               end do


               if (pot_name.eq.'ccpol8s_omo') then
                   call ccpol8s_omo(vdwpe)
               end if

               !  Ewald summation for charges

               call electrostatics(rcpe)

               !  calling N-Body Iterative Induction model

               call nb_induction(indpe)


               stpvdwpe = vdwpe
               stpcpe = scpe+rcpe+kcpe
               stpindpe = indpe

               write(100,*)nd, stpvdwpe,stpcpe,stpindpe,stpvdwpe+stpcpe+stpindpe

               write(*,*)nd, stpvdwpe,stpcpe,stpindpe,stpvdwpe+stpcpe+stpindpe


   end do

   close(30)

   end if

   close(100)
end program main
