! ********************************************************************
! ********************************************************************
! Updated: TLJ; 1/6/2017
! Filename: update_psinew.f90
! ********************************************************************
! ********************************************************************
!
!     ________________________________________________________________
!
!     BRIEF CODE DESCRIPTION:
!     =======================
!
!     Subroutine updates the marker function psi(i,k,j)
!       used in the solid to distinguish between pure
!       AP and binder; also updates surface and solid
!       values
!
!     Solid level-set function to compute solid values
!        psi = +1    AP, oxidizer
!        psi = -1    Binder, fuel
!
!     subroutines: UPDATE_PSI
!                  LEVELSET_3D
!                  LEVELSET_2D
!
! ********************************************************************

SUBROUTINE UPDATE_PSI(t,tout)

  USE GLOBAL_DATA

  IMPLICIT NONE

  INTEGER ::  i ,n, j, k, m, kkk
  INTEGER :: zmin,zmax,xmin,xmax
  REAL*8  :: xcut,zcut,ycut,t,tout,ycutold
!
! computes the level set for each mode separately
!
!     t > 0; ipack = 1
!

  if (ipack == 1) THEN
     if(ndim==2) then
     DO j=drange(5),drange(6)
        kkk=0
        DO m=1,ncells
           DO n=1,ncells
              kkk=kkk+1
              zmin = nzquad(m,1)
              zmax = nzquad(m,2)
              xmin = nxquad(n,1)
              xmax = nxquad(n,2)
              DO k=zmin,zmax
              DO i=xmin,xmax
                 xcut = x(i)
                 zcut = z(k)
                 ycut = phi(i,k) - y(j)
                 ycut = dmod(ycut,period)
                 ycutold=oldphi(i,k) - y(j)
                 ycutold=mod(ycutold,period)
                 if(ycut > ycutold .OR. tout <= 0.0) psi(i,k,j)=1.0
                 CALL LEVELSET_2D(xcut,ycut,i,k,j)
              END DO
              END DO
           END DO
        END DO
     END DO
     endif
     if(ndim==3) then
     DO j=drange(5),drange(6)
        kkk=0
        DO m=1,ncells
           DO n=1,ncells
              kkk=kkk+1
              zmin = nzquad(m,1)
              zmax = nzquad(m,2)
              xmin = nxquad(n,1)
              xmax = nxquad(n,2)
              DO k=zmin,zmax
              DO i=xmin,xmax
                 xcut = x(i)
                 zcut = z(k)
                 ycut = phi(i,k) - y(j)
                 ycut = dmod(ycut,period)
                 ycutold=oldphi(i,k) - y(j)
                 ycutold=mod(ycutold,period)
                 if(ycut > ycutold .OR. tout <= 0.0) psi(i,k,j)=1.0
                 CALL LEVELSET_3D(xcut,zcut,ycut,i,k,j,kkk,nquad(kkk))
              END DO
              END DO
           END DO
        END DO
     END DO
     endif
  else
     if(ipack == -1) then
        xloc = xend
     end if
     do i = drange(1)-1,drange(2)+1
        psi(i,:,:) = one          !AP
        if (abs(x(i)) <= xloc .AND. equivalence_ratio  /= 0) &
             psi(i,:,:) = -one      !BINDER
     enddo
  endif

  CALL FILL_PSI_GHOST

!------------------------------------------------------------------------
  RETURN
END SUBROUTINE UPDATE_PSI
! ********************************************************************
SUBROUTINE LEVELSET_3D(xcut,zcut,ycut,i,k,j,qu,qumax)

  USE GLOBAL_DATA

  IMPLICIT NONE

!  Local Variables
  INTEGER :: n,i,j,k,nn,qu,qumax
  REAL*8  :: psinew, r2, psi1,xcut,ycut,zcut
  REAL*8  :: eps,delx,delz,dely,delta

! computes the level set for each mode separately

  psinew = qumax
  nn=anint(abs(psi(i,k,j)))
  eps=1.d-9
  do while(psinew.ge.eps.and.nn.le.qumax)
     n = quad(qu,nn)
     delx=xcut-xcenter(n)
     delz=zcut-zcenter(n)
     r2 = rad(n)*rad(n)
     delta=r2-delx**2-delz**2
     dely=ycut-ycenter(n)
     if(dely**2.le.delta)then !inside the sphere n
           psinew=-nn
     elseif (delta.ge.0.d0.and.dely.gt.eps) then !PROJECT ON the sphere n
           psinew=nn
           nn=qumax
     endif
     nn=nn+1
  enddo
  psi(i,k,j) = -psinew

  RETURN
END SUBROUTINE LEVELSET_3D
! ********************************************************************
! ********************************************************************
SUBROUTINE LEVELSET_2D(xcut,ycut,i,k,j)

  USE GLOBAL_DATA

  IMPLICIT NONE

!  Local Variables
  INTEGER :: n,i,j,k 
  REAL*8  :: psinew, r2, psi1
  REAL*8  :: xcut, ycut, a_1,a_2, theta,a_3,alpha_test
!-----------------------------------------------------------------

! computes the level set for each mode separately

  psinew = 100.0d0
  do n=1,ncircles
!-modi0
     r2 = rad(n)*rad(n)
     psi1 = ((xcut-xcenter(n))**2 &
          + (ycut-ycenter(n))**2) - r2
     psinew = min(psi1,psinew)
  end do
  psi(i,k,j) = -psinew

!-modi0
  if (psi(i,k,j) .gt. 0.d0) psi(i,k,j) = one
  if (psi(i,k,j) .lt. 0.d0) psi(i,k,j) = -one

  RETURN
END SUBROUTINE LEVELSET_2D
!******************************************************************
