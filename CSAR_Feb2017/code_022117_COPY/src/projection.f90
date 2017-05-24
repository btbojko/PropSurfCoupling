!*********************************************************************
SUBROUTINE PROJECTION(flagValidate)

  USE GLOBAL_DATA
  USE mg_solver
  USE MYMPI

  IMPLICIT NONE

  include "parallel.h"

!-----------------------------------------------------------------------
!  Local variables
  INTEGER :: i, k, j, eqn, col(7), error
  INTEGER :: xs,zs,xe,ze,ys,ye,flagValidate
  REAL*8  :: gx,gy,gz,dxa,dya,dza,scoef(3),sourceterm,myerr_divt,err_divt
!-----------------------------------------------------------------------
!
  if(SkipPressure.or. .not. doprojection) then
     p = 0d0;pold = 0d0;divt = 0d0;rate = 0d0;
     return
  endif
!
  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
  ys = drange(5)+1
  ye = drange(6)-1
!
!  the array col contains information about the color of the pints to be updated
!  col(1:2), the dimension to be updated, col(3:4), and how many point in the
!  j direction need to be updated col(5:6)
!
  col(1)=1; col(2)=2; col(3)=1; col(4) = ndim; 
  col(5)=drange(5); col(6)=drange(6); col(7) = 0

  dxa = 1.0d0/dx
  dza = 1.0d0/dz

  CALL VELOCITY( merge(5,4,flagValidate > 0) )

!
!  evaluate the rhs of the pressure Poisson equation (PPE)
!
  qpick = (/1, 15, 21 , 1/)
  sourceterm = zero;
  do j=ys,ye
     dya = detady(j)/dy
     do k=zs,ze
        do i=xs,xe
           quartuple =(/myid,i,k,j/)
           MQterm(i,k,0) = (one +  MQchipr(j,1)*MQphi(i,k))
           if( .not. skipContSource) then
              scoef(1) = (two*dtCont(1)+dtCont(2))/(dtCont(1)+dtCont(2))/dtCont(1)
              scoef(2) = -(dtCont(1)+dtCont(2))/dtCont(1)/dtCont(2)
              scoef(3) = dtCont(1)/dtCont(2)/(dtCont(1)+dtCont(2))
              sourceterm = sum(scoef*rhogasCont(i,k,j,1:3)) 
           endif
           gx = (uvel(i,k,j)-uvel(i-1,k,j))*dxa
           gz = (wvel(i,k,j)-wvel(i,k-1,j))*dza
           gy = (vvel(i,k,j)-vvel(i,k,j-1))*dya
           if(flagValidate >= 1) then
              Postdivt(i,k,j)= (gx+gy+gz) + sourceterm*MQterm(i,k,0) 
              if(j<JRFLO(1) .and. abs(Postdivt(i,k,j)) > 1d2 .and. flagValidate == 1) then
                 write(*,'(4i4,1p2e12.4,a,15e12.4,a,15e12.4)')myid,i,k,j,Postdivt(i,k,j),divt(i,k,j),'MAX DIV excdeed'
              endif
           else
              divt(i,k,j)= (gx+gy+gz) + sourceterm*MQterm(i,k,0)
           endif
        enddo
     enddo
  enddo

!---set-up convenience boundary conditions
  divt(:,:,ny) = divt(:,:,ny-1)  !This value should never be used
  divt(:,:,0) = zero
  !call FourierExtVelUpdate(3,1)

  if(flagValidate == 1) then
     myerr_divt = maxval(abs(Postdivt(:,:,1:min(JRFLO(1)-1,ny-1)) ))
     CALL MPI_ALLREDUCE(myerr_divt,err_divt,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
     if(myid ==  0)write(*,'(i6,a,1p4e12.4)') ncyc,'DIV_ERROR',err_divt,Uimposed
  endif
  if(flagValidate >0) RETURN 

!---discretize
  CALL get_coeff_ppe

!
!  SHIFT the vectors to be used in the MGsolver (easy to coarsen)
!
  divt=cshift(divt,1,3)
  p=cshift(p,1,3)

  finest_mesh%x => p
  finest_mesh%f => divt

  CALL fmg_solve   !CALL THE MULTIGRID SOLVER

  divt=cshift(divt,-1,3)
  p=cshift(p,-1,3)


  CALL UPDATE_Q

  CALL parallel_swap4(col,finest_mesh)

!-----------------------------------------------------------------------
  RETURN
END SUBROUTINE PROJECTION
!*************************************************************************


!*************************************************************************
SUBROUTINE UPDATE_Q

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

  include "parallel.h"
!-------------------------------------------------------------------------   
! Local variables
  INTEGER :: i, k, j, xs,zs,xe,ze,ys,ye
  REAL*8  :: gx,gy,gz,dx1a,dya,dz1a,dya4
  REAL*8  :: pressure,poo,coe,rhogas
!-------------------------------------------------------------------------   

  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
  ys = drange(5)+1
  ye = drange(6)-1

  dx1a = 1.0d0/dx
  dz1a = 1.0d0/dz    

  poo = PRESSURE(tcyc(2))
  coe = two*dim_fact*poo

!...........................................................................
! These BCs for the pressure are pplied so update_q is valid at j=1 & j=ny-1
!...........................................................................
  !!p(:,:,0) = two*p(:,:,1)-p(:,:,2)  !!OLD
  p(:,:,0) = p(:,:,1)
  p(:,:,ny) = zero

  do i=xs-1,xe+1
     do k=zs-1,ze+1   
        do j =0,ny
           pold(i,k,j) = pold(i,k,j) + p(i,k,j)
        enddo
     Enddo
  Enddo


!
!  UPDATE VELOCITY AND PRESSURE GRADIENT
!
  DO j=ys,ye
     dya4 = quarter*detady(j)/dy


     do k = zs-1, ze+1
        do i = xs-1, xe+1
           MQterm(i,k,0) = (one +  MQchipr(j,1)*(MQphi(i,k)))
           MQterm(i,k,1) = (one +  MQchipr(j,1)*(MQphi(i+1,k)+MQphi(i,k))*half)
           MQterm(i,k,2) = (one +  MQchipr(j,1)*(MQphi(i,k+1)+MQphi(i,k))*half)
           MQterm(i,k,3) = (one +  MQchipr(j,2)*(MQphi(i,k)))
        enddo
     enddo

     do i=xs-1,xe+1
        do k=zs-1,ze+1
           if(skipVarDensity) then
              rhogas = rhoRef
           else
              rhogas = coe/(f(i,k,j,1)+f(i+1,k,j,1))
           endif
           gy = (MQchi(j+1,1)*(p(i,k,j+1)+p(i+1,k,j+1)) - MQchi(j-1,1)*(p(i,k,j-1)+p(i+1,k,j-1)))*dya4
           gx = (MQterm(i+1,k,0)*p(i+1,k,j)-MQterm(i,k,0)*p(i,k,j))*dx1a
           q(i,k,j,1) = q(i,k,j,1) + (dphidxa(i,k)*gy - gx)/rhogas/MQterm(i,k,1)


           if(skipVarDensity) then
              rhogas = rhoRef
           else
              rhogas = coe/(f(i,k,j,1)+f(i,k+1,j,1))
           endif
           gy = (MQchi(j+1,1)*(p(i,k,j+1)+p(i,k+1,j+1))-MQchi(j-1,1)*(p(i,k,j-1)+p(i,k+1,j-1)))*dya4
           gz = (MQterm(i,k+1,0)*p(i,k+1,j)-MQterm(i,k,0)*p(i,k,j))*dz1a
           q(i,k,j,2) = q(i,k,j,2) + (dphidza(i,k)*gy - gz)/rhogas/MQterm(i,k,2)
        enddo
     enddo

     dya = detadya(j)/dy
     do i=xs-1,xe+1
        do k=zs-1,ze+1

           if(skipVarDensity) then
              rhogas = rhoRef
           else
              rhogas = coe/(f(i,k,j+1,1)+f(i,k,j,1))
           endif
           gy = (p(i,k,j+1)-p(i,k,j))*dya
           q(i,k,j,3) = q(i,k,j,3) - gy/rhogas/MQterm(i,k,3)

           gy = (MQchi(j+1,1)*p(i,k,j+1)- MQchi(j,1)*p(i,k,j))*dya
           gx = (MQterm(i+1,k,3)*(p(i+1,k,j+1)+p(i+1,k,j)) &
                -MQterm(i-1,k,3)*(p(i-1,k,j+1)+p(i-1,k,j)))*dx1a/4d0
           gz = (MQterm(i,k+1,3)*(p(i,k+1,j+1)+p(i,k+1,j)) &
                -MQterm(i,k-1,3)*(p(i,k-1,j+1)+p(i,k-1,j)))*dz1a/4d0
           rate(i,k,j,1) = rate(i,k,j,1) + ( dphidx(i,k)*gy - gx)/rhogas/MQterm(i,k,3)
           rate(i,k,j,2) = rate(i,k,j,2) + ( dphidz(i,k)*gy - gz)/rhogas/MQterm(i,k,3)

        enddo
     enddo
  ENDDO
  !!$  CALL MGTEST

!-------------------------------------------------------------------------
  RETURN
END SUBROUTINE UPDATE_Q
!*************************************************************************


!*************************************************************************
SUBROUTINE get_coeff_ppe

  USE data_types
  USE global_data

  IMPLICIT NONE

!-------------------------------------------------------------------------
! Local variables:
  TYPE(mesh_pointers), POINTER :: current_mesh
  TYPE(cppevec), DIMENSION(:,:), POINTER :: CPPE
  REAL(KIND=double), DIMENSION(:), POINTER :: ey,eya
  REAL(KIND=double), DIMENSION(:,:), POINTER :: chi_mg,chipr_mg

  INTEGER :: xe,ze,ye,i,k,j,mnx,mnz,m

  REAL(KIND=double) :: dxf,dzf,dyf,dxc,dzc,dyc,dxa,dza,dya
  REAL(KIND=double) :: dxa2,dza2,dy2a,dxdy4,dzdy4
  REAL(KIND=double) :: term,c379
  REAL(KIND=double) :: dya4,dyxz4
!--------------------------------------------------------------------------

  current_mesh => finest_mesh
  dxf = finest_mesh%dx ; dzf = finest_mesh%dz; dyf = finest_mesh%dy

  DO

     xe = current_mesh%xe ; ze = current_mesh%ze; ye = current_mesh%ye
     dxc = current_mesh%dx ; dzc = current_mesh%dz; dyc = current_mesh%dy
     mnx = anint(dxc/dxf); mnz = anint(dzc/dzf);
     m= anint( log10 (dble(mnx)) / log10 (2.0) )

     dxa = 1.0/dxc ; dza = 1.0/dzc; dya = 1.0/dyc
     dxa2 = dxa*dxa; dza2 = dza*dza; dy2a = dya*dya
     dxdy4 = dxa*dya/4.0d0;dzdy4 = dza*dya/4.0d0
     dya4 = dya/4.0d0;


     mcomp   = 1
     mgradjM = 0
     mgradj  = 2

     ey       => current_mesh%ey 
     eya      => current_mesh%eya
     chi_mg   => current_mesh%chi
     chipr_mg => current_mesh%chipr
     CPPE     =>  current_mesh%cppe


!
!   b1=c1
!
     current_mesh%c1 = -dxa2


     DO j = 0,ye

        do i=-1,xe+1
           do k=-1,ze+1
              MQterm(i,k,mcomp)   = (one +  chipr_mg(j,1)  *  phi_mg(i,k,m) )
              MQterm(i,k,mgradjM) = (one +  chipr_mg(j-1,2)*  phi_mg(i,k,m) )
              MQterm(i,k,mgradj)  = (one +  chipr_mg(j,2)  *  phi_mg(i,k,m) )
           enddo
        enddo

        dyterm = ey(j)*  (/ eya(j-1), eya(j) /)  *dya**2
        dyxz4 = ey(j)*dya4*dxa

        quartuple =(/myid,i,k,j+1/)
        if(j == 0) then

           do i=0,xe
              do k=0,ze
                 CPPE(i,k)%ccmg(2,j) = current_mesh%c1*MQterm(i-1,k,mcomp) - &
                      &chi_mg(j,mgradj)*MQterm(i-1,k,mgradj)/MQterm(i,k,mgradj)*dphidx_mg(i,k,m)*dyxz4 &
                      &+ dphidxa_mg(i-1,k,m)*chi_mg(j-1,mcomp)*dyxz4

                 CPPE(i,k)%ccmg(3,j) = current_mesh%c1*MQterm(i+1,k,mcomp) + &
                      &chi_mg(j,mgradj)*MQterm(i+1,k,mgradj)/MQterm(i,k,mgradj)*dphidx_mg(i,k,m)*dyxz4 &
                      &- dphidxa_mg(i,k,m)*chi_mg(j-1,mcomp)*dyxz4

                 CPPE(i,k)%ccmg(4,j) = current_mesh%c1*MQterm(i,k-1,mcomp) - &
                      &chi_mg(j,mgradj)*MQterm(i,k-1,mgradj)/MQterm(i,k,mgradj)*dphidz_mg(i,k,m)*dyxz4 &
                      &+ dphidza_mg(i,k-1,m)*chi_mg(j-1,mcomp)*dyxz4

                 CPPE(i,k)%ccmg(5,j) = current_mesh%c1*MQterm(i,k+1,mcomp) + &
                      &chi_mg(j,mgradj)*MQterm(i,k+1,mgradj)/MQterm(i,k,mgradj)*dphidz_mg(i,k,m)*dyxz4 &
                      &- dphidza_mg(i,k,m)*chi_mg(j-1,mcomp)*dyxz4



                 CPPE(i,k)%ccmg(6,j) = zero
                 CPPE(i,k)%ccmg(7,j) = ( (-dphidx_mg(i,k,m)**2-dphidz_mg(i,k,m)**2)/MQterm(i,k&
                      &,mgradj)*chi_mg(j,mgradj)*chi_mg(j+1,mcomp) - one/MQterm(i&
                      &,k,mgradj)  )*dyterm(2) + (dphidxa_mg(i,k,m)-dphidxa_mg(i-1,k,m) +&
                      & dphidza_mg(i,k,m)-dphidza_mg(i,k-1,m) ) *chi_mg(j+1,mcomp)*dyxz4


                 CPPE(i,k)%ccmg(8,j) = zero
                 CPPE(i,k)%ccmg(9,j) = zero
                 CPPE(i,k)%ccmg(10,j) =( dphidxa_mg(i,k,m)*chi_mg(j+1,mcomp)+chi_mg(j,mgradj)/&
                      &MQterm(i,k,mgradj)*MQterm(i+1,k,mgradj)*dphidx_mg(i,k,m))*dyxz4
                 CPPE(i,k)%ccmg(11,j) =(-dphidxa_mg(i-1,k,m)*chi_mg(j+1,mcomp)-chi_mg(j,mgradj)&
                      &/MQterm(i,k,mgradj)*MQterm(i-1,k,mgradj)*dphidx_mg(i,k,m))*dyxz4

                 CPPE(i,k)%ccmg(12,j) = zero
                 CPPE(i,k)%ccmg(13,j) = zero
                 CPPE(i,k)%ccmg(14,j) =(chi_mg(j,mgradj)/MQterm(i,k,mgradj)*MQterm(i,k+1,mgradj)&
                      &*dphidz_mg(i,k,m)+dphidza_mg(i,k,m)*chi_mg(j+1,mcomp))*dyxz4
                 CPPE(i,k)%ccmg(15,j) =(-chi_mg(j,mgradj)/MQterm(i,k,mgradj)*MQterm(i,k-1,mgradj)&
                      &*dphidz_mg(i,k,m)-dphidza_mg(i,k-1,m)*chi_mg(j+1,mcomp))*dyxz4

                 CPPE(i,k)%ccmg(1,j) = ( (-dphidx_mg(i,k,m)**2-dphidz_mg(i,k,m)**2)*chi_mg(j,mgradj)*&
                      &chi_mg(j,mcomp) - one  )*dyterm(2)/MQterm(i,k,mgradj) &
                      &- 4d0* MQterm(i, k, mcomp)*dxa**2 &
                      &-(-dphidxa_mg(i,k,m)+dphidxa_mg(i-1,k,m)-dphidza_mg(i,k,m)+dphidza_mg(i,k-1,m))*&
                      &chi_mg(j-1,mcomp)*dyxz4  
              enddo
           enddo
        else
           do i=0,xe
              do k=0,ze
                 CPPE(i,k)%ccmg(2,j) = current_mesh%c1*MQterm(i-1,k,mcomp) - &
                      &( chi_mg(j,mgradj)  *MQterm(i-1,k,mgradj) /MQterm(i,k,mgradj) &
                      &- chi_mg(j-1,mgradj)*MQterm(i-1,k,mgradjM)/MQterm(i,k,mgradjM))&
                      &*dphidx_mg(i,k,m)*dyxz4
                 CPPE(i,k)%ccmg(3,j) = current_mesh%c1*MQterm(i+1,k,mcomp) + &
                      &( chi_mg(j,mgradj)  *MQterm(i+1,k,mgradj) /MQterm(i,k,mgradj) &
                      & -chi_mg(j-1,mgradj)*MQterm(i+1,k,mgradjM)/MQterm(i,k,mgradjM))&
                      &*dphidx_mg(i,k,m)*dyxz4
                 CPPE(i,k)%ccmg(4,j) = current_mesh%c1*MQterm(i,k-1,mcomp) - &
                      &( chi_mg(j,mgradj)  *MQterm(i,k-1,mgradj) /MQterm(i,k,mgradj) &
                      & -chi_mg(j-1,mgradj)*MQterm(i,k-1,mgradjM)/MQterm(i,k,mgradjM))&
                      &*dphidz_mg(i,k,m)*dyxz4

                 CPPE(i,k)%ccmg(5,j) = current_mesh%c1*MQterm(i,k+1,mcomp) + &
                      &( chi_mg(j,mgradj)  *MQterm(i,k+1,mgradj) /MQterm(i,k,mgradj) -&
                      &  chi_mg(j-1,mgradj)*MQterm(i,k+1,mgradjM)/MQterm(i,k,mgradjM))&
                      &*dphidz_mg(i,k,m)*dyxz4
                 CPPE(i,k)%ccmg(6,j) = ( (-dphidx_mg(i,k,m)**2-dphidz_mg(i,k,m)**2)/MQterm(i,k&
                      &,mgradjM)*chi_mg(j-1,mgradj)*chi_mg(j-1,mcomp) - one/MQterm(i&
                      &,k,mgradjM)  )*dyterm(1) + (-dphidxa_mg(i,k,m)+dphidxa_mg(i-1,k,m) -&
                      & dphidza_mg(i,k,m)+dphidza_mg(i,k-1,m))*chi_mg(j-1,mcomp)*dyxz4
                 CPPE(i,k)%ccmg(7,j) = ( (-dphidx_mg(i,k,m)**2-dphidz_mg(i,k,m)**2)/MQterm(i,k&
                      &,mgradj)*chi_mg(j,mgradj)*chi_mg(j+1,mcomp) - one/MQterm(i&
                      &,k,mgradj)  )*dyterm(2) + (dphidxa_mg(i,k,m)-dphidxa_mg(i-1,k,m) +&
                      & dphidza_mg(i,k,m)-dphidza_mg(i,k-1,m) ) *chi_mg(j+1,mcomp)*dyxz4
                 CPPE(i,k)%ccmg(8,j) = (dphidxa_mg(i-1,k,m)*chi_mg(j-1,mcomp) + chi_mg(j-1,mgradj)&
                      &/MQterm(i,k,mgradjM)*MQterm(i-1,k,mgradjM)*dphidx_mg(i,k,m))*dyxz4
                 CPPE(i,k)%ccmg(9,j) = (-dphidxa_mg(i,k,m)*chi_mg(j-1,mcomp)-chi_mg(j-1,mgradj)&
                      &/MQterm(i,k,mgradjM)*MQterm(i+1,k,mgradjM)*dphidx_mg(i,k,m))*dyxz4
                 CPPE(i,k)%ccmg(10,j) =( dphidxa_mg(i,k,m)*chi_mg(j+1,mcomp)+chi_mg(j,mgradj)/&
                      &MQterm(i,k,mgradj)*MQterm(i+1,k,mgradj)*dphidx_mg(i,k,m))*dyxz4
                 CPPE(i,k)%ccmg(11,j) =(-dphidxa_mg(i-1,k,m)*chi_mg(j+1,mcomp)-chi_mg(j,mgradj)&
                      &/MQterm(i,k,mgradj)*MQterm(i-1,k,mgradj)*dphidx_mg(i,k,m))*dyxz4
                 CPPE(i,k)%ccmg(12,j) =(chi_mg(j-1,mgradj)/MQterm(i,k,mgradjM)*&
                      &MQterm(i,k-1,mgradjM)*dphidz_mg(i,k,m)+dphidza_mg(i,k-1,m)*chi_mg(j-1,mcomp))*dyxz4
                 CPPE(i,k)%ccmg(13,j) =(-chi_mg(j-1,mgradj)/MQterm(i,k,mgradjM)*MQterm(i,k+1,mgradjM)&
                      &*dphidz_mg(i,k,m)-dphidza_mg(i,k,m)*chi_mg(j-1,mcomp))*dyxz4
                 CPPE(i,k)%ccmg(14,j) =(chi_mg(j,mgradj)/MQterm(i,k,mgradj)*MQterm(i,k+1,mgradj)&
                      &*dphidz_mg(i,k,m)+dphidza_mg(i,k,m)*chi_mg(j+1,mcomp))*dyxz4
                 CPPE(i,k)%ccmg(15,j) =(-chi_mg(j,mgradj)/MQterm(i,k,mgradj)*MQterm(i,k-1,mgradj)&
                      &*dphidz_mg(i,k,m)-dphidza_mg(i,k-1,m)*chi_mg(j+1,mcomp))*dyxz4
                 CPPE(i,k)%ccmg(1,j) = ( (-dphidx_mg(i,k,m)**2-dphidz_mg(i,k,m)**2)&
                      &*chi_mg(j-1,mgradj)*chi_mg(j,mcomp)- one  )*dyterm(1)/MQterm(i,k,mgradjM)&
                      &+ ( (-dphidx_mg(i,k,m)**2-dphidz_mg(i,k,m)**2)*chi_mg(j,mgradj)*&
                      &chi_mg(j,mcomp) - one  )*dyterm(2)/MQterm(i,k,mgradj) &
                      &- 4d0* MQterm(i, k, mcomp)*dxa**2
              enddo
           enddo
        endif
     ENDDO
     if(.false.) then
        do j = 0,0
           do i=0,xe
              do k=0,ze
                 CPPE(i,k)%ccmg(1,j) = sum(CPPE(i,k)%ccmg(2:15,j))
              enddo
           enddo
        enddo
     endif

     IF ( .NOT. ASSOCIATED(current_mesh%coarse_mesh) ) EXIT
     current_mesh => current_mesh%coarse_mesh
  END DO

!---------------------------------------------------------------------
  RETURN
END SUBROUTINE get_coeff_ppe
!**********************************************************************


!**********************************************************************
SUBROUTINE RHSPTEST

  USE GLOBAL_DATA
  USE data_types
  USE MYMPI

  IMPLICIT NONE

  include "parallel.h"

!-----------------------------------------------------------------------
!  Local variables
  INTEGER :: i, k, j, eqn, col(7)
  INTEGER :: xs,zs,xe,ze,ys,ye
  REAL*8  :: gx,gy,gz,dxa,dya,dza
  REAL*8  :: pressure,poo,coe,maxq(2),allmaxq(2),rhogas
!-----------------------------------------------------------------------

  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
  ys = drange(5)+1
  ye = drange(6)-1
  col(1)=1; col(2)=2; col(3)=1; col(4) = ndim; 
  col(5)=drange(5); col(6)=drange(6); col(7) = 0

  dxa = 1.0d0/dx
  dza = 1.0d0/dz

  poo = PRESSURE(tcyc(2))
  coe = two*dim_fact*poo

  do j=ys,ye
     dya = detadya(j)/dy
     do i=xs,xe
        do k=zs,ze
           rhogas = coe/(f(i,k,j+1,1)+f(i,k,j,1))
           gy = (p(i,k,j+1)-p(i,k,j))*dya
           gx = 0.25*(p(i+1,k,j+1)-p(i-1,k,j+1)+p(i+1,k,j)-p(i-1,k,j))*dxa
           gz = 0.25*(p(i,k+1,j+1)-p(i,k-1,j+1)+p(i,k+1,j)-p(i,k-1,j))*dza
           uvel(i,k,j) = uvel(i,k,j) + (dphidx(i,k)*gy - gx)/rhogas
           wvel(i,k,j) = wvel(i,k,j) + (dphidz(i,k)*gy - gz)/rhogas
        enddo
     enddo
  enddo


  CALL parallel_swap4(col,finest_mesh)

  CALL VELOCITY(6)

!
!  TEST the projection
!
  wrk = zero
  do j=ys,ye
     dya = detady(j)/dy
     do k=zs,ze
        do i=xs,xe
           gx = (uvel(i,k,j)-uvel(i-1,k,j))*dxa
           gz = (wvel(i,k,j)-wvel(i,k-1,j))*dza
           gy = (vvel(i,k,j)-vvel(i,k,j-1))*dya
           wrk(i,k,j)= gx+gy+gz 
        enddo
     enddo
  enddo
!
  slambda = zero
  do j=0,ny
     do k=zs,ze
        do i=xs,xe
           slambda(j) = slambda(j) + vvel(i,k,j)
        enddo
     enddo
  enddo

  maxq(1) = maxval(abs(vvel(0:xe,0:ze,0:ye))); maxq(2) = maxval(abs(wrk))
  CALL MPI_ALLREDUCE(maxq,allmaxq,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
  if(myid == 0)&
       write(*,*) ncyc,'RHPTEST,maxvvel,maxresid', allmaxq
  CALL MPI_ALLREDUCE(slambda(0),rlambda(0),50,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
  if(myid == 0) then
     do j=0,ny
        write(*,*) j,'VELSUM',rlambda(j)
     enddo
  endif

  !!write(*,*) 'RHPTEST', maxval(abs(wrk)),maxloc(abs(wrk)),ye,finest_mesh%ye,ubound(wrk)
!
!-----------------------------------------------------------------------
  RETURN
END SUBROUTINE RHSPTEST
!**********************************************************************



!**********************************************************************
SUBROUTINE MGTEST

  USE data_types
  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!-----------------------------------------------------------------------
!  Local variables
  INTEGER :: i, k, j,jmax,ierrWG,jj
  INTEGER :: xs,zs,xe,ze,ys,ye
  REAL*8  :: gx,gy,gz,dxa,dya,dza, myerr_mg(2),err_mg(2)
  REAL(KIND=double) :: dxa2,dza2,dy2a,dxdy4,dzdy4
  REAL(KIND=double) :: cd(4),cj(9),c3,c7,ey,c1,c9,term,c379,vv(15),cc(15)
  LOGICAL :: outcome
!-----------------------------------------------------------------------

  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)
  ys = drange(5)+1
  ye = drange(6)-1

  if(JRFLO(1) >=ny) then
     jmax = ye+1
  else
     jmax = JRFLO(1)
  endif


!
  wrk= zero

  do i=0,xe
     do k=0,ze

        do j=1,jmax-2
           jj = j+1

           quartuple =(/myid,i,k,jj/)
           qpick = (/6,27,13, 4/)

           vv = (/-p(i,k,jj),p(i-1,k,jj),p(i+1,k,jj),p(i,k-1,jj),p(i,k+1,jj),p(i,k,jj-1),p(i,k,jj+1),&
                &p(i-1,k,jj-1),p(i+1,k,jj-1),p(i+1,k,jj+1),p(i-1,k,jj+1),&
                &p(i,k-1,jj-1),p(i,k+1,jj-1),p(i,k+1,jj+1),p(i,k-1,jj+1) /)
           cc(1:15) = finest_mesh%cPPE(i,k)%ccmg(1:15,j)



           wrk(i,k,j)= divt(i,k,jj)  + sum(cc*vv)

           !!$           if(wrk(i,k,j) > 1.5d2) then
           !!$              write(*,'(4i4,1p2e12.4,a,15e12.4,a,15e12.4)')myid,i,k,j,wrk(i,k,j),divt(i,k,jj),'CC',cc,'VV',vv
           !!$              call MPI_FINALIZE(ierr)
           !!$              stop
           !!$           endif
           !!$           if(all(quartuple == qpick)) then
           !!$              print*,'MGTEST---',sum(cc*vv)
           !!$           endif

           if(ncyc_run < ncyc_debug) then
              call NAN_VAL(q(i,k,j,1:3),outcome,3)
              if(outcome) then
                 write(6,*)'============================================'
                 write(6,*)'DEBUG SECTION FOR PROJECTION SOLVER'
                 write(6,*)'N_CYCLE',ncyc
                 write(6,*)'Myid',myid
                 write(6,*)'POS REL This proc',i,k,j
                 write(6,*)'grid limits',drange(2),drange(4),drange(6)
                 write(6,*)'POS ABS',i+ib_mg(1),k+kb_mg(1),j
                 write(6,'(a,1p50e12.4)')'NEW SOLN',q(i,k,j,1:3)
                 write(6,'(a,1p50e12.4)')'OLD SOLN',oldsoln(i,k,j,maxsize+1:maxsize+3)
                 write(6,'(a,1p50e12.4)')'SOURCE TERM',divt(i,k,j),rhogasCont(i,k,j,1:3),dtCont
                 write(6,'(a,1p50e12.4)')'PRESSURE i varing',p(i-1:i+1,k,j)
                 write(6,'(a,1p50e12.4)')'NEIGH SOL,i-',q(i-2:i,k,j,1:3)
                 write(6,'(a,1p50e12.4)')'NEIGH SOL,i+',q(i:i+2,k,j,1:3)
                 write(6,'(a,1p50e12.4)')'NEIGH SOL,k-',q(i,k-2:k,j,1:3)
                 write(6,'(a,1p50e12.4)')'NEIGH SOL,k+',q(i,k:k+2,j,1:3)
                 write(6,*)'############################################'
                 call MPI_Abort(MPI_COMM_WORLD,ierrWG)
                 call MPI_finalize(MPI_COMM_WORLD)
                 stop
              endif
           endif



        enddo
     enddo
  enddo

  myerr_mg= (/ maxval(abs(wrk)), maxval(abs(divt)) /)

  CALL MPI_ALLREDUCE(myerr_mg,err_mg,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
  if(myid ==  0) then
     write(*,'(i6,1p2e12.4,a,12e12.4)')ncyc,tcyc,'MGtest',err_mg,Uimposed(1:2)
  endif
!-----------------------------------------------------------------------
  RETURN
END SUBROUTINE MGTEST
!**********************************************************************


!**********************************************************************
SUBROUTINE Test_MGcoeffs

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE
  include "parallel.h"


!-------------------------------------------------------------------------   
! Local variables
  INTEGER :: i, k, j, m,col(7),ii,kk,jj
  REAL(KIND=double) :: vv(15),cc(15),rhogas,pressure
  INTEGER ::  status(MPI_STATUS_SIZE)
!-------------------------------------------------------------------------   

  i = qpick(2);k=qpick(3);j=qpick(4);
  do m = -1,16
     if(m>=0) then
        if(m>0) p=zero;
        q = zero;rate = zero;
        jj = j
        do ii = max(Lbound(vvel,1),Lbound(vec,1)), min(Ubound(vvel,1),Ubound(vec,1))
           do kk = max(Lbound(vvel,2),Lbound(vec,2)), min(Ubound(vvel,2),Ubound(vec,2))
              if(skipVarDensity) then
                 rhogas = rhoRef
              else
                 rhogas = two*dim_fact*PRESSURE(tcyc(2)) / (f(ii,kk,jj+1,1)+f(ii,kk,jj,1))
                 if(jj == ny) rhogas = two*dim_fact*PRESSURE(tcyc(2)) /(two*f(ii,kk,jj,1))
              endif
              vvel(ii,kk,0) =  rhogas*MQvelocity(i,k,j,3)
           enddo
        enddo
     endif
     if(myid == qpick(1)) then
        if(m==1)p(i,k,j) = one
        if(m==2)p(i-1,k,j) = one
        if(m==3)p(i+1,k,j) = one
        if(m==4)p(i,k-1,j) = one
        if(m==5)p(i,k+1,j) = one
        if(m==6)p(i,k,j-1) = one
        if(m==7)p(i,k,j+1) = one
        if(m==8)p(i-1,k,j-1) = one
        if(m==9)p(i+1,k,j-1) = one
        if(m==10)p(i+1,k,j+1) = one
        if(m==11)p(i-1,k,j+1) = one
        if(m==12)p(i,k-1,j-1) = one
        if(m==13)p(i,k+1,j-1) = one
        if(m==14)p(i,k+1,j+1) = one
        if(m==15)p(i,k-1,j+1) = one
     endif

     if(k==0) then
        CALL MPI_SENDRECV(p(i,-1,j),1,MPI_DOUBLE_PRECISION,dnbr(3),0,&
             &p(i,drange(4),j),1,MPI_DOUBLE_PRECISION,dnbr(4),0,comm3d,status,ierr)
        CALL MPI_SENDRECV(p(i,-1,j+1),1,MPI_DOUBLE_PRECISION,dnbr(3),0,&
             &p(i,drange(4),j+1),1,MPI_DOUBLE_PRECISION,dnbr(4),0,comm3d,status,ierr)
        CALL MPI_SENDRECV(p(i,-1,j-1),1,MPI_DOUBLE_PRECISION,dnbr(3),0,&
             &p(i,drange(4),j-1),1,MPI_DOUBLE_PRECISION,dnbr(4),0,comm3d,status,ierr)
     endif

     if(m>=0) then
        col(1:7) = (/drange(1:4),0,ny,1/)
        call parallel_swap(p,col(1:6),finest_mesh)
        call UPDATE_Q   
        col = (/1,2,1,ndim,drange(5),drange(6),0/)
        CALL parallel_swap4(col,finest_mesh)
     endif

     if(m==4) then
        do ii = lbound(p,1),ubound(p,1)
           do kk = lbound(p,2),ubound(p,2)
              do jj = lbound(p,3),ubound(p,3)
                 if(abs(p(ii,kk,jj)) > 0d0) print*,myid,ii,kk,jj,phi_mg(ii,kk,0)
              enddo
           enddo
        enddo
     endif

     if(myid == qpick(1)) then

        call  PROJECTION(2)
        if(m>0)  then
           write(*,'(i3,1p34e14.6)')M,finest_mesh%CPPE(i,k)%ccmg(m,j-1),PostDivt(i,k,j),q(i,k-1,j,:)
        else
           vv = (/-p(i,k,j),p(i-1,k,j),p(i+1,k,j),p(i,k-1,j),p(i,k+1,j),p(i,k,j-1),p(i,k,j+1),&
                &p(i-1,k,j-1),p(i+1,k,j-1),p(i+1,k,j+1),p(i-1,k,j+1),&
                &p(i,k-1,j-1),p(i,k+1,j-1),p(i,k+1,j+1),p(i,k-1,j+1) /)
           cc(1:15) = finest_mesh%cPPE(i,k)%ccmg(1:15,j-1)
           print*,M,'SUM',sum(cc*vv),PostDivt(i,k,j)
        endif
     endif
  enddo

  call MPI_FINALIZE(i)
  stop

!-------------------------------------------------------------------------
  RETURN
END SUBROUTINE TEST_MGCOEFFS
!******************************************************************************
