module coulomb

  use, intrinsic :: iso_fortran_env
  use common_arrays

  implicit none
  private

  public :: electrostatics


contains

  subroutine electrostatics(qrcpe)

    implicit none

    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi

    integer :: i,j,is,js,k
    real(8) :: qvir_rc,qvirial
    real(8) :: wij,alp
    real(8) :: xcc,ycc,zcc
    real(8) :: dx,dy,dz
    real(8) :: rsx,rsy,rsz
    real(8) :: rccsq,rssq,drr
    real(8) :: qis,qjs,chgprd
    real(8) :: qrcpe,fc,cfij
    real(8) :: fxsij,fysij,fzsij

    qvirial = 0.d0
    qrcpe =0.d0

     do i=1,nm-1
        do j=i+1, nm

              do is=1,ns

                  do js=1,ns
                 
                  qis = chgs(is,i)
                  qjs = chgs(js,j)  
                  
                  rsx = xxst(is,i) - xxst(js,j) 
                  rsy = yyst(is,i) - yyst(js,j) 
                  rsz = zzst(is,i) - zzst(js,j) 
                 
                  rssq = rsx**2 + rsy**2 + rsz**2
                          
                  drr = sqrt(rssq)

                  qrcpe = qrcpe + r4pie0*qis*qjs/drr

 
                 end do 

              end do
           
        end do
     end do

     qrcpe = qrcpe/418.4d0       ! real space energy for coulomb 
 
  end subroutine electrostatics

end module coulomb
