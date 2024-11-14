module utility

  use, intrinsic :: iso_fortran_env
  use common_arrays

  implicit none

  private

  public :: jacobi, gauss
  public :: tt,dtt
  public :: matvec1

contains

! scalar product
  function scalp(a,b)  
        implicit none
        
        real(8) ::  a(3),b(3)
        real(8) :: scalp

        scalp = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)
  end


! Multiply vetor v by the matrix a. Store in u.

  subroutine matvec1(a,v,u)
    implicit none
    integer :: i,j,nnn
    real(8) ::  a(3,3),v(3),u(3)

    nnn = 3
    do i=1,nnn
       u(i) = 0.d0
       do j=1,nnn
          u(i) = u(i) + a(i,j)*v(j)
       end do
    end do
  
  end


! compute the factorial of n
!
  function fact(nmax)
    implicit none
    integer :: nmax
    integer :: i
    real(8) :: fact
    fact = 1.d0
    do i=1,nmax
    fact = fact*i
    end do
    
!    write(200,*)'fact', fact    
  end function fact

  ! tang-toennies damoing function
  function tt(n, delta, r)
    implicit none
    integer :: i,ncn
    integer, intent(in) :: n
    real(8), intent(in) :: delta, r
    real(8) :: tt
    real(8) :: term,ssum,deltar

    deltar = delta * r
    term = 1.0d0
    ssum = 1.0d0
    ncn = n
    do i=1, ncn
       term = term*deltar/i
       ssum = ssum + term
    end do
    tt = 1.0 - exp(-deltar)*ssum

  end function tt

  function dtt(n, delta, r)
    implicit none
    integer :: i,ncn
    integer, intent(in) :: n
    real(8), intent(in) :: delta, r
    real(8) :: dtt
    real(8) :: delr, ffact

    ffact = 1.0d0
    ncn = n
    do i=1, ncn
       ffact = ffact * i
    end do
    delr = delta * r
    dtt = delta*exp(-delr)*delr**n/ffact

  end function dtt

   subroutine jacobi(a,v,n)

!c***********************************************************************
!c
!c     diagonalisation of real symmetric matices by jacobi method
!c
!c     input parameters:
!c
!c     a(n,n) is the matrix to be diagonalised
!c     v(n,n) is the eigenvector matrix
!c     n   is the dimension of the matrices
!c
!c     jacobi processes lower triangle only (upper triangle unchanged)
!c
!c     variable rho sets absolute tolerance on convergence
!c     variable tes is a moving tolerance that diminishes
!c     on each pass until at true convergence tes<rho
!c
!c     author w.smith 1993
!c
!c***********************************************************************

      implicit none

      logical :: pass
      integer :: n,i,j,k
      real(8) :: a,v,rho,tes,scl,v1,v2,v3,u,omg,s,c,tem

      dimension a(n,n),v(n,n)

      rho=1.0d-16
      tes=0.0d0
      scl=0.0d0

!c     initialize eigenvectors

      do i=1,n
        do j=1,n
          v(i,j)=0.0d0
        end do
        v(i,i)=1.0d0
      end do

!c     rescale matrix for optimal accuracy

      do i=1,n
        if(abs(a(i,i)).gt.scl) scl=abs(a(i,i))
      end do
      do i=1,n
        do j=1,i
          a(i,j)=a(i,j)/scl
        end do
      end do

!c     set initial value of moving tolerance

      do i=2,n
        do j=1,i-1
          tes=tes+2.0d0*a(i,j)*a(i,j)
        enddo
      enddo
      tes=sqrt(tes)

!c     recycle until absolute tolerance satisfied

      do while(tes.gt.rho)

        tes=tes/dble(n)
        if(tes.lt.rho)tes=rho

!c     jacobi diagonalisation

        pass=.true.

!c     recycle until moving tolerance satisfied

        do while(pass)

          pass=.false.

          do i=2,n

            do j=1,i-1

              if(abs(a(i,j)).ge.tes)then
                pass=.true.
                v1=a(j,j)
                v2=a(i,j)
                v3=a(i,i)
                u=0.5d0*(v1-v3)
                if(abs(u).lt.rho)then
                  omg=-1.0d0
                else
                  omg=-v2/sqrt(v2*v2+u*u)
                  if(u.lt.0.0d0)omg=-omg
                endif
                s=omg/sqrt(2.0d0*(1.0d0+sqrt(1.0d0-omg*omg)))
                c=sqrt(1.0d0-s*s)
                do k=1,n
                  if(k.ge.i)then
                    tem=a(k,j)*c-a(k,i)*s
                    a(k,i)=a(k,j)*s+a(k,i)*c
                    a(k,j)=tem
                  else if(k.lt.j)then
                    tem=a(j,k)*c-a(i,k)*s
                    a(i,k)=a(j,k)*s+a(i,k)*c
                    a(j,k)=tem
                  else
                    tem=a(k,j)*c-a(i,k)*s
                    a(i,k)=a(k,j)*s+a(i,k)*c
                    a(k,j)=tem
                  endif
                  tem=v(k,j)*c-v(k,i)*s
                  v(k,i)=v(k,j)*s+v(k,i)*c
                  v(k,j)=tem
                enddo
                a(j,j)=v1*c*c+v3*s*s-2.0d0*v2*s*c
                a(i,i)=v1*s*s+v3*c*c+2.0d0*v2*s*c
                a(i,j)=(v1-v3)*s*c+v2*(c*c-s*s)
              endif

            enddo

          enddo

        enddo

      enddo

!c     rescale matrix

      do i=1,n
        do j=1,i
          a(i,j)=scl*a(i,j)
        enddo
      enddo

      return
      end subroutine jacobi

      function duni()

!c*********************************************************************
!c     
!c     dl_poly random number generator based on the universal
!c     random number generator of marsaglia, zaman and tsang
!c     (stats and prob. lett. 8 (1990) 35-39.) it must be
!c     called once to initialise parameters u,c,cd,cm
!c     
!c     copyright daresbury laboratory 1992
!c     author -  w.smith         july 1992
!c     
!c*********************************************************************

      implicit none

      logical new
      integer ir,jr,i,j,k,l,m,ii,jj
      real(4) s,t,u,c,cd,cm,uni
      real(8) duni
      dimension u(97)
      save u,c,cd,cm,uni,ir,jr,new
      data new/.true./

      if(new)then

!c     initial values of i,j,k must be in range 1 to 178 (not all 1)
!c     initial value of l must be in range 0 to 168.

        i=12
        j=34
        k=56
        l=78
!c     
        ir=97
        jr=33
        new=.false.

        do 200 ii=1,97
          s=0.0
          t=0.5
          do 100 jj=1,24
            m=mod(mod(i*j,179)*k,179)
            i=j
            j=k
            k=m
            l=mod(53*l+1,169)
            if(mod(l*m,64).ge.32)s=s+t
            t=0.5*t
  100     continue
          u(ii)=s
  200   continue
        c =  362436.0/16777216.0
        cd= 7654321.0/16777216.0
        cm=16777213.0/16777216.0
      else

!c     calculate random number
        uni=u(ir)-u(jr)
        if(uni.lt.0.0)uni=uni+1.0
        u(ir)=uni
        ir=ir-1
        if(ir.eq.0)ir=97
        jr=jr-1
        if(jr.eq.0)jr=97
        c=c-cd
        if(c.lt.0.0)c=c+cm
        uni=uni-c
        if(uni.lt.0.0)uni=uni+1.0
        duni=dble(uni)
      endif
      
      return
      end function duni


     subroutine gauss(nnn,vxx,vyy,vzz)

!c*********************************************************************
!c     
!c     dl_poly subroutine for constructing velocity arrays
!c     with a gaussian distribution of unit variance.
!c     
!c     based on the box-muller method
!c     
!c     note - this version uses a universal random number 
!c     generator, which generates pseudo-random numbers between
!c     0 and 1. it is based on the algorithm of marsaglia, zaman
!c     and tsang in: stats and prob. lett. 8 (1990) 35-39.
!c     
!c     copyright daresbury laboratory 2007
!c     author - w. smith         nov  2007
!c     
!c*********************************************************************
      
!      use setup_module
      
      implicit none

      integer nnn,ii
      real(8) vxx,vyy,vzz,rrr,rr1,rr2
      
      dimension vxx(nnn),vyy(nnn),vzz(nnn)
      
!c     initialise random number generator
      
      rrr=duni()
      
!c     calculate gaussian random numbers
      
      do ii=1,2*(nnn/2),2
        
        rr1=sqrt(-2.d0*log(duni()))
        rr2=2.d0*pi*duni()
        vxx(ii)=rr1*cos(rr2)
        vxx(ii+1)=rr1*sin(rr2)

        rr1=sqrt(-2.d0*log(duni()))
        rr2=2.d0*pi*duni()
        vyy(ii)=rr1*cos(rr2)
        vyy(ii+1)=rr1*sin(rr2)

        rr1=sqrt(-2.d0*log(duni()))
        rr2=2.d0*pi*duni()
        vzz(ii)=rr1*cos(rr2)
        vzz(ii+1)=rr1*sin(rr2)
        
      enddo
      if(mod(nnn,2).ne.0)then
        
        rr1=sqrt(-2.d0*log(duni()))
        rr2=2.d0*pi*duni()
        vxx(nnn)=rr1*cos(rr2)
        vyy(nnn)=rr1*sin(rr2)
        rr1=sqrt(-2.d0*log(duni()))
        rr2=2.d0*pi*duni()
        vzz(nnn)=rr1*cos(rr2)
        
      endif
      
      return
      end subroutine gauss

end module utility  
