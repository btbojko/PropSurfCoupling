! ********************************************************************
SUBROUTINE SOLID_PROP

USE GLOBAL_DATA

IMPLICIT NONE
  
!---------------------------------------------------------------------
! Local Variables
INTEGER :: i,j,k
REAL*8  :: tmp_x, tmp_h, tmp_t
!--------------------------------------------------------------------

!!! OLD  DENSITY

rhos_old = rhos

do k=drange(3)-1,drange(4)+1
    do i=drange(1)-1,drange(2)+1
        xsurf(i,k,2) = 0.0d0
        xsurf(i,k,3) = 1.0d0
        xsurf(i,k,4:neqgas) = 0.0d0
        if (psi(i,k,0) < zero) then
! Binder surface
            if(premixed == 0) then
                xsurf(i,k,2) = 1.0d00 - alp_W
                xsurf(i,k,3) = alp_W
            else    
                xsurf(i,k,2) = 1.0d00
                xsurf(i,k,3) = 0.0d0
            endif
            xsurf(i,k,4:neqgas) = 0.0d0
        endif   
    enddo
enddo

do k=drange(3)-1,drange(4)+1
    do i=drange(1)-1,drange(2)+1
        do j=drange(5),drange(6)
            qheats(i,k,j) = qheat_ap
            rhos(i,k,j) = rho_ap
            lambdas(i,k,j) = lambda_ap
            da_s(i,k,j) = da_ap
            theta_s(i,k,j) = theta_ap
            if ( psi(i,k,j) < zero ) then
! BINDER 
                qheats(i,k,j) = alp_W*qheat_ap+(1.0-alp_W)*qheat_binder
                rhos(i,k,j) = alp_V*rho_ap+(1.0-alp_V)*rho_binder
                tmp_x=lambda_ap/lambda_binder
                tmp_t=alp_V
                tmp_h=tmp_x+0.5d0*(1.0d0-tmp_x)**2.0d0*(1.0d0-tmp_t)**2.0d0
                if(tmp_x.lt.1.0d0) then
                    lambdas(i,k,j)=lambda_binder*(tmp_h+sqrt(tmp_h*tmp_h-tmp_x*tmp_x))
                else   
                    lambdas(i,k,j)=lambda_binder*(tmp_h-sqrt(tmp_h*tmp_h-tmp_x*tmp_x))
                endif
                da_s(i,k,j) = da_binder
                theta_s(i,k,j) = theta_binder
 
         !       da_s(i,k,j) = da_ap**alp_V * da_binder**(1.0d0-alp_V)
         !       theta_s(i,k,j) = alp_V*theta_ap+(1.0d0-alp_V)*theta_binder

            endif
        enddo
    enddo
enddo

!-------------------------------------------------------------------
RETURN
END SUBROUTINE SOLID_PROP
! ********************************************************************
! ********************************************************************
SUBROUTINE SMOOTH

USE GLOBAL_DATA

IMPLICIT NONE
  
!---------------------------------------------------------------------
! Local Variables
INTEGER :: i,j,k,m,flag_old
REAL*8  :: tmp_x, tmp_h, tmp_t,vv(neqgas,2),limit_H
REAL*8  :: Heaviside
!--------------------------------------------------------------------


limit_H = 0.3*.0
ny_l_re = 24

call reinitia(19,1,1)

 rhos_old = rhos

vv = zero

vv(2,1) = 1.0d0
vv(3,1) =  1.0d0
do k=drange(3)-1,drange(4)+1
    do i=drange(1)-1,drange(2)+1
        do m =2,neqgas
            xsurf(i,k,m)= Heaviside(vv(m,1),vv(m,2),ff(i,0),limit_H)
        enddo
    enddo
enddo

do k=drange(3)-1,drange(4)+1
    do i=drange(1)-1,drange(2)+1
        do j=drange(5),drange(6) 
            qheats(i,k,j)  = Heaviside(qheat_binder,qheat_ap, ff(i,j),limit_H)
            rhos(i,k,j)    = Heaviside(rho_binder,rho_ap,ff(i,j),limit_H)
            lambdas(i,k,j) = Heaviside(lambda_binder,lambda_ap,ff(i,j),limit_H)
            da_s(i,k,j)    = Heaviside(da_binder,da_ap,ff(i,j),limit_H)
            theta_s(i,k,j) = Heaviside(theta_binder,theta_ap,ff(i,j),limit_H)
        enddo
    enddo
enddo

!-------------------------------------------------------------------
 RETURN
 END SUBROUTINE SMOOTH
! ********************************************************************
! ********************************************************************
SUBROUTINE PRN_TMP

USE GLOBAL_DATA

IMPLICIT NONE
  
!---------------------------------------------------------------------
! Local Variables
INTEGER :: i,j,k
REAL*8  :: xx,yy,yys
!--------------------------------------------------------------------

do k=drange(3)-1,drange(4)+1
    do i=drange(1)-1,drange(2)+1
        do j=drange(5),drange(6)
            psi(i,k,j) = rhos(i,k,j)
        enddo
    enddo
enddo

open (UNIT=25,FILE='solid.dat',STATUS='UNKNOWN')
allocate(G_psi(-1:nx+1,0:ny))
ALLOCATE(G_x(0:nx),G_phi(0:nx),G_speed(0:nx))
call gather_x
call gather_phi
call gather_psi
do j = 0, ny
    do i = 0, nx
        xx = G_x(i)
        yy = y(j) + G_phi(i)
        yys = -y(j) + G_phi(i)
        write(25,900) xx, yys, G_psi(i,j)
    enddo
enddo
   
deallocate(G_psi) 
DEALLOCATE(G_x,G_phi,G_speed)

900 FORMAT(2x,19e18.6)

stop

!-------------------------------------------------------------------
RETURN
END SUBROUTINE PRN_TMP
! ********************************************************************
