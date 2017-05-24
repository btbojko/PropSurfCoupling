! ********************************************************************
SUBROUTINE UPDATE_PHI(t,tout)

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

!---------------------------------------------------------------------
  REAL*8, INTENT(IN) :: t,tout
! Local variables
  INTEGER ::  i, k,imax,kmax,xs,xe,zs,ze,itemp
  INTEGER ::  done_period
  REAL*8  ::  h1,h2,h3,h4
  REAL*8  ::  h1a,h3a,maxvec,allmaxvec
!---------------------------------------------------------------------

  
  if (tout == zero) GOTO 100

  oldphi = phi

  CALL PYROLYSIS

  CALL FLUX(t,tout)

  CALL FILL_PHI_GHOST
  itemp = 6   !Gross
  CALL DIRICHLET_BC(itemp,finest_mesh)

  if(ncyc < ncyc_constant) then
    phi = zero  !this is for a flat surface
  endif

100 CONTINUE

!--- COMPUTE shape derivatives
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

!
! NOTE, each processor will have a different phi_per, t_per,
! but we will consider only the one for processor (0)
!

  done_period = 0
  if(phi(0,0) < phi_per + period .AND. &         ! [phi < 0]
       oldphi(0,0) >= phi_per + period) then

     phi_per = oldphi(0,0)
     t_per = t

     done_period = 1
     IF(myid == 0) then
        write(*,*) 'NEW PERIOD', phi_per,phi(0,0),oldphi(0,0),t_per,period
       
     ENDIF

  endif  
  
  CALL MPI_BCAST(done_period,1,MPI_INTEGER,0,comm3d,ierr)

  !!if(done_period == 1) &
  !!     CALL DUMP(t)

!---------------------------------------------------------------------
  RETURN
END SUBROUTINE UPDATE_PHI
!*********************************************************************


!*********************************************************************
SUBROUTINE UPDATE_PHI_MG

  USE GLOBAL_DATA
  IMPLICIT NONE

  include "parallel.h"

!--------------------------------------------------------------------
! local variables
  INTEGER ::  i,k,j, nxc,nzc, im,km,jm,f2c,col(6),msh,imsh
  REAL *8 :: h1a, h1
  TYPE(mesh_pointers), POINTER :: current_mesh
!---------------------------------------------------------------------

! just inject the value of phi in wrk vector !
! note, im,km are the indexes  on the fine mesh that correspond to i,k
! on mesh number <mesh>
! 
  DO j=2,iex
     jm = j-1
     msh = lx_mg(j)
     nxc=nx_mg(msh)
     nzc=nz_mg(msh)
     f2c = 2**(jm)
!
! Find the mesh corresponding to mesh number <msh>. this way of doing 
! is a little crude but should not be time consuming
!
     current_mesh => finest_mesh
     do imsh = 2,msh
        current_mesh => current_mesh%coarse_mesh
     enddo

     if(current_mesh%blank) EXIT

     do i=0,nxc
        im=f2c*(i+ib_mg(msh))-ib_mg(1)
     do k=0,nzc
        km=f2c*(k+kb_mg(msh))-kb_mg(1)
        wrk(i,k,jm) = 0 !phi(im,km)
     enddo
     enddo


!
! update the ghost cells for each level
! now that we pass also the mesh pointer maybe us the case to avoid 
! passing the wrk vecto as a separate pointer but we can pass 
! current_mesh%f which would be included in current_mesh
!
     col=0; col(5:6)=jm; col(2) = nxc-1; col(4)=nzc-1
     call parallel_swap(wrk,col,current_mesh)

!
! divide by dx  (h1a==h3a) to simplify the implementation
!
     h1a = 1.0d0/(dx*dble(f2c))
     h1 = h1a/two

     do i= -1,nxc-1
     do k= -1,nzc-1
        dphidxa_mg(i,k,jm) = (wrk(i+1,k,jm)-wrk(i,k,jm))*h1a 
        dphidza_mg(i,k,jm) = (wrk(i,k+1,jm)-wrk(i,k,jm))*h1a
     enddo
     enddo

     do i= 0,nxc-1
     do k= 0,nzc-1
        dphidx_mg(i,k,jm) = (wrk(i+1,k,jm)-wrk(i-1,k,jm))*h1
        dphidz_mg(i,k,jm) = (wrk(i,k+1,jm)-wrk(i,k-1,jm))*h1
     enddo
     enddo
    
  ENDDO  !do j=2,iex

  if(is2D) THEN
    dphidza_mg = zero
    dphidz_mg = zero
  ENDIF

!------------------------------------------------------------------------
  RETURN
END SUBROUTINE UPDATE_PHI_MG
!*************************************************************************

! ********************************************************************
SUBROUTINE PROPAGATE_SURF(t,tout,maxtime,done)

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

!---------------------------------------------------------------------
  REAL*8, INTENT(INOUT) :: t,tout
  INTEGER, INTENT(IN) :: maxtime
  LOGICAL, INTENT(OUT) :: done
! Local variables
  INTEGER :: itime
!---------------------------------------------------------------------

  do itime = 1,maxtime
     
     tout = t + dt
     CALL UPDATE_PHI(t,tout)
     CALL UPDATE_PSI(t,tout)

     done = abs(t-tstop) <= 1.0d-10

     if(done) EXIT

  enddo


!---------------------------------------------------------------------
  RETURN
END SUBROUTINE PROPAGATE_SURF
!*********************************************************************
