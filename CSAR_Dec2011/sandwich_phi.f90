! ********************************************************************
SUBROUTINE SANDWICH_PHI(t,tout)

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

!---------------------------------------------------------------------
  REAL*8, INTENT(IN) :: t
  REAL*8, INTENT(OUT) :: tout
! Local variables
  INTEGER ::  i, k,imax,kmax,xs,xe,zs,ze, itime,maxtime
  REAL*8  :: h1,h2,h3,h4
  REAL*8  :: h1a,h3a,maxvec,allmaxvec,philoc
!---------------------------------------------------------------------
!
! Compute new surface shape
!
  tout = t + dt
  CALL PYROLYSIS

  maxtime = merge(50, 10, ncyc > 5000)
  do itime = 1,maxtime

     oldphi=phi
     CALL FILL_PHI_GHOST
     CALL SANDWICH_FLUX(t,tout)
 
  enddo

  CALL FILL_PHI_GHOST
!
! Compute derivatives
!
  h1 = 1.0d0/(2.0d0*dx)
  h1a = 1.0d0/dx
  h2 = 1.0d0/(dx*dx)
  h3 = 1.0d0/(2.0d0*dz)
  h3a = 1.0d0/dz
  h4 = 1.0d0/(dz*dz)

  xs = drange(1);xe = drange(2);zs=drange(3);ze = drange(4)
  maxvec=0.0
  do k = drange(3)-1, drange(4)+1
     do i = drange(1)-1, drange(2)+1
        dphidx(i,k) = (phi(i+1,k)-phi(i-1,k))*h1
        dphi2dx2(i,k) = (phi(i+1,k)&
             -2.0d0*phi(i,k)+phi(i-1,k))*h2
        dphidz(i,k) = (phi(i,k+1)-phi(i,k-1))*h3
        dphi2dz2(i,k) = (phi(i,k+1)&
             -2.0d0*phi(i,k)+phi(i,k-1))*h4
        dphidxa(i,k) = (phi(i+1,k)-phi(i,k))*h1a
        dphidza(i,k) = (phi(i,k+1)-phi(i,k))*h3a
        vec(i,k)=sqrt(1.0d0+dphidx(i,k)**2+dphidz(i,k)**2)
     end do
  end do

  if(num_mg_meshes > 1)&
       CALL UPDATE_PHI_MG

!---------------------------------------------------------------------
  RETURN
END SUBROUTINE SANDWICH_PHI
!*********************************************************************

! ********************************************************************
SUBROUTINE SANDWICH_FLUX(t,tout)
      
  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE
      
!---------------------------------------------------------------------
! Dummy Variables
  REAL*8, INTENT(IN) :: t,tout
! Local Variables
  INTEGER ::  i, k
  REAL*8  :: eps,mymax,mymin
  REAL*8  :: fx,fz,is1,is2,ap1,ap2,phix1,phix2,phiz1,phiz2
!---------------------------------------------------------------------
  
  ! compute time step
  
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

!!!  if(mod(ncyc,2*writemod) == 0) write(*,*)ncyc,t,tout,'MAX MIN SURF VELOCITIES',surfvel_max,surfvel_min
  
  
  !     Update phi
  do k=drange(3),drange(4)
     do i=drange(1),drange(2)
        phi(i,k) = phi(i,k) - dt*phit(i,k)*3.0
     enddo
  enddo
  
!---------------------------------------------------------------------------
  RETURN 
END SUBROUTINE SANDWICH_FLUX
!***************************************************************************
