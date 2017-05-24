
! ********************************************************************
SUBROUTINE FLUX(t,tout)
      
  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE
      
!---------------------------------------------------------------------
! Dummy Variables
  REAL*8, INTENT(IN) :: t,tout
! Local Variables
  INTEGER ::  i, k
  REAL*8  :: dtt,eps,mymax,mymin
  REAL*8  :: fx,fz,is1,is2,ap1,ap2,phix1,phix2,phiz1,phiz2
!---------------------------------------------------------------------
  
! compute time step
  dtt = dt
  
  eps = 1.0d-8
  do k=drange(3),drange(4)
     do i=drange(1),drange(2)+1
        !           X-DERIVATIVE CALCULATION
        fx = 0.d0
        if(phi(i-1,k).le.phi(i+1,k).and.&
             phi(i-1,k).lt.phi(i,k)) then
           is1 = (phi(i-2,k)-2.*phi(i-1,k)+phi(i,k)  )**2
           is2 = (phi(i-1,k)-2.*phi(i,k)  +phi(i+1,k))**2
           ap1 = 0.5/(eps+is1)
           ap2 = 1.0/(eps+is2)
           phix1 = (phi(i-2,k)-4.*phi(i-1,k)+3.*phi(i,k))&
                /(2.*dx)
           phix2 = (phi(i+1,k)-phi(i-1,k))/(2.*dx)
           fx = (ap1*phix1+ap2*phix2)/(ap1+ap2)
        endif
        
        if(phi(i+1,k).le.phi(i-1,k).and.&
             phi(i+1,k).lt.phi(i,k)) then
           is1 = (phi(i-1,k)-2.*phi(i,k)  +phi(i+1,k))**2
           is2 = (phi(i,k)  -2.*phi(i+1,k)+phi(i+2,k))**2
           ap1 = 1.0/(eps+is1)
           ap2 = 0.5/(eps+is2)
           phix1 = (phi(i+1,k)-phi(i-1,k))/(2.*dx)
           phix2 = (-3.*phi(i,k)+4.*phi(i+1,k)-phi(i+2,k))/(2.*dx)
           fx = (ap1*phix1+ap2*phix2)/(ap1+ap2)
        endif

        phit(i,k) = rb(i,k)*sqrt(one+fx*fx)
     enddo
  enddo
  
  mymax = maxval(phit(drange(1):drange(2),drange(3)))
  mymin = minval(phit(drange(1):drange(2),drange(3)))

  CALL MPI_ALLREDUCE(mymax,surfvel_max,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
  CALL MPI_ALLREDUCE(mymin,surfvel_min,1,MPI_DOUBLE_PRECISION,MPI_MIN,comm3d,ierr)

 
  !     Update phi
  do k=drange(3),drange(4)
     do i=drange(1),drange(2)

           phi(i,k) = phi(i,k) - dtt*phit(i,k)

     enddo
  enddo
  
!---------------------------------------------------------------------------
  RETURN 
END SUBROUTINE FLUX
!***************************************************************************
