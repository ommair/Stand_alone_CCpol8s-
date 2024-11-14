module vdw_ccpol8s_omo

  use, intrinsic :: iso_fortran_env
  use common_arrays
  use utility 
  
  implicit none
  private

  public :: ccpol8s_omo

contains

  subroutine ccpol8s_omo(qqvdwe)

    implicit none
 
    real(8) :: xx(nom),yy(nom),zz(nom)
    real(8) :: xxs(nsite,nom),yys(nsite,nom),zzs(nsite,nom)
    real(8) :: boxi   

    integer :: i,j,is,js,k
    real(8) :: qqvdwe,qvdwe, qvirial
    real(8) :: wij,qis,qjs
    real(8) :: xcc,ycc,zcc
    real(8) :: dx,dy,dz
    real(8) :: rccsq
    real(8) :: rssq, rsx,rsy,rsz
    real(8) :: vfij
    real(8) :: sqr
    real(8) :: eks,u_exp,du_exp
    real(8) :: u_ele, f1,df1, du_ele
    real(8) :: u_ind,f6,f8,f10,df6,df8,df10
    real(8) :: r6,r8,r10,du_ind
!    real(8) :: tt

    qvdwe=0.d0

     do i=1,nm-1
        do j=i+1, nm

             do is=1,ns

                do js=1,ns
                 
                  rsx = xxst(is,i) - xxst(js,j) 
                  rsy = yyst(is,i) - yyst(js,j) 
                  rsz = zzst(is,i) - zzst(js,j)  
                  
                  qis = chgs(is,i)
                  qjs = chgs(js,j)

                  rssq = rsx**2 + rsy**2 + rsz**2
 
                  sqr = sqrt(rssq)
          
                  ! exponential (excahnge - repulsion) 
                  eks = exp(-beta_ccpol(styp(is),styp(js))*sqr) 
                  u_exp = eks*(c0_ccpol(styp(is),styp(js)) + c1_ccpol(styp(is),styp(js))*sqr &
                              + c2_ccpol(styp(is),styp(js))*sqr**2)

                  ! induction-dispersion
                  f6 = tt(6, domo_ccpol(styp(is),styp(js)), sqr)
                  f8 = tt(8, domo_ccpol(styp(is),styp(js)), sqr)
                  r6 = rssq*rssq*rssq
                  r8 = rssq*r6
                  u_ind = f6*c6_ccpol(styp(is),styp(js))/r6+f8*c8_ccpol(styp(is),styp(js))/r8  

                  qvdwe = qvdwe + u_exp + u_ind
 
                end do 

             end do
           
        end do
     end do

     qqvdwe = qvdwe/418.4d0

!     write(*,*) qqvdwe
 
  end subroutine ccpol8s_omo


end module vdw_ccpol8s_omo 
