!!***************************************************************
SUBROUTINE BCQ

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE
  include "parallel.h"

!----------------------------------------------------------------
! Local variables
  INTEGER :: i, k, j, eqn,xs,zs,xe,ze,col(7)
  REAL *8 :: poo,coe,usurf,wsurf,vsurf
  REAL *8 :: pressure, rratio, my_flux, dya, rhogas, fact1
!----------------------------------------------------------------

  CALL UPDATE_CONDITIONS_BC

  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
  poo = PRESSURE(tcyc)
  coe = two*dim_fact*poo
  col(1)=1; col(2)=2; col(3)=1; col(4) = ndim;
  col(5)=drange(5); col(6)=col(5)      !j=0 ONLY!
  col(7) = 0

!
!  IMPOSE velocity boundary conditions
!
  do i=xs-1,nxv(3)+1
     do k=zs-1,ze+1
!
!SOUTH
        rhogas=coe/(f(i,k,0,1)+f(i,k,1,1)) !rhogas on the surf
        rratio = rhos(i,k,0)/rhogas-1.0d0
        vsurf=rb(i,k)*rratio/vec(i,k)    
        if(vsurf<0)then 
           write(*,*)i,k,'ERROR vsur <0 at wall',vsurf,myid,f(i,k,0,1),f(i,k,1,1),ncyc
           write(*,*)'pressure and density',rhogas,rhos(i,k,0),rratio,coe,poo,ncyc
        endif
        q(i,k,0,3) = vsurf 
        vec_dirichlet(3,3,i) = q(i,0,0,3)
        vvel(i,k,0) = rhos(i,k,0)*rb(i,k)*vec(i,k)    !rhogas*vvel
!NORTH
        vvel(i,k,ny) = vvel(i,k,ny-1)
!........................................................................
! These formulae preerve the derivative at the boundary j=ny
!........................................................................
!      q(i,k,ny+1,1) = q(i,k,ny-2,1)-3.0*q(i,k,ny-1,1)+3.0*q(i,k,ny,1)
!      q(i,k,ny+1,2) = q(i,k,ny-2,2)-3.0*q(i,k,ny-1,2)+3.0*q(i,k,ny,2)
!
     enddo
  enddo

  do i=xs-1,nxv(1)
     do k=zs-1,ze
        usurf = -half*(q(i+1,k,0,3)+q(i,k,0,3))*dphidxa(i,k)  !-vsurf*fx
        wsurf = -half*(q(i,k,0,3)+q(i,k+1,0,3))*dphidza(i,k)  !-vsurf*fz
        q(i,k,0,1)=usurf
        q(i,k,0,2)=wsurf
     enddo
     vec_dirichlet(1,3,i) = q(i,0,0,1)
     vec_dirichlet(2,3,i) = q(i,0,0,2)
  enddo

!--- EXTRA COMMUNICATION
  finest_mesh%w => q
  var_updated = 1
  call parallel_swap4(col,finest_mesh)  

!--- DIRICHLET BOUNDARY CONDITIONS
  do i = 1,ndim
     CALL  DIRICHLET_BC(i,finest_mesh)
  enddo
  

  do i=xs-1,nxv(1)+1
     do k=zs-1,ze+1
        uvel(i,k,0) = q(i,k,0,1)
        wvel(i,k,0) = q(i,k,0,2)
        qsurf(i,k,1:2) =  q(i,k,0,1:2)  !this is on the surface
        
!
!   2nd order extrapolation at the wall, file BCUW.maple
!
        q(i,k,0,1)=(8.0d0*uvel(i,k,0)-6.0d0*q(i,k,1,1)+q(i,k,2,1))*third
        q(i,k,0,2)=(8.0d0*wvel(i,k,0)-6.0d0*q(i,k,1,2)+q(i,k,2,2))*third
  
     enddo
     vec_dirichlet(1,3,i) = q(i,0,0,1)
     vec_dirichlet(2,3,i) = q(i,0,0,2)
  enddo

!--- UPDATE NORTH BOUNDARY CONDITIONS -- these will have to be changed according to the problem
  q(:,:,ny,3) = q(:,:,ny-2,3)   !zero gradient for v
  q(:,:,ny,2) = zero            !2d conditions for  w
  if(ALL(type_bc(1,1,:) == 1)) then
     q(:,:,ny,1) = zero                  !conservation of mass for u
  elseif(ALL(type_bc(1,1,:) == 0)) then  !conservation of mass for u
     do k=zs-1,ze+1
        do i=xs-1,nxv(1)+1
           fact1 = (f(i,k,ny,1)+f(i+1,k,ny,1))/&
                &  (f(ibc(1,1,1),k,ny,1)+f(ibc(1,1,1)+1,k,ny,1))
           q(i,k,ny,1) =  vec_dirichlet(1,1,ny)*fact1
        enddo
     enddo
  endif

!-------------------------------------------------------------------------------
  RETURN
END SUBROUTINE BCQ
!*******************************************************************************
