! ********************************************************************
SUBROUTINE LAMBDA(lavg)

  USE GLOBAL_DATA

  IMPLICIT NONE
!----------------------------------------------------------------------
! Local Variables
  INTEGER :: i,j,k,ys,ye,lavg
  REAL*8  :: hx, hy, hz
  REAL*8  :: coe, tavg(3), c1,c2
!----------------------------------------------------------------------
  
!
! now update gas-phase lambda
!

  hx = 1.0d0/(2.0d0*dx)
  hy = 1.0d0/(2.0d0*dy)
  hz = 1.0d0/(2.0d0*dz)
  c1 = coe_lamb1 ; c2 = coe_lamb2;

  if (lavg.eq.0) then

    do j=drange(5),drange(6)
    do k=drange(3)-1,drange(4)+1
    do i=drange(1)-1,drange(2)+1
  
       lambdag(i,k,j) = c1*f(i,k,j,1)+c2
  
    enddo
    enddo
    enddo

  ELSEif (lavg.eq.1) then

    do j=drange(5),drange(6)
    do k=drange(3)-1,drange(4)+1
    do i=drange(1)-1,drange(2)+1
  
       lambdag(i,k,j) = c1*tmdpnt(i,k,j)+c2
  
    enddo
    enddo
    enddo
  
  
    do i=drange(1)-1,drange(2)
    do k=drange(3)-1,drange(4)
    do j=drange(5)-1,drange(6)
  
       tavg(1) = quarter*(tmdpnt(i,k,j) + tmdpnt(i+1,k,j)+&   ! i+1/2,  k+1/2  ,j
                          tmdpnt(i,k+1,j)+tmdpnt(i+1,k+1,j))
       tavg(2) = quarter*(tmdpnt(i,k,j) + tmdpnt(i+1,k,j)+&   ! i+1/2,  k      ,j+1/2
                          tmdpnt(i,k,j+1)+tmdpnt(i+1,k,j+1))
       tavg(3) = quarter*(tmdpnt(i,k,j) + tmdpnt(i,k+1,j)+&   ! i,     k+1/2   ,j+1/2
                          tmdpnt(i,k,j+1)+tmdpnt(i,k+1,j+1))
       lbdavg(i,k,j,:) = c1*tavg(:)+c2
  
    enddo
    enddo
    enddo


  endif   !if (lavg.eq.0)
     
  if (lavg.eq.0) then

    ys = drange(5)
    ye = drange(6)

    if(ddrange(5).eq.0) ys = ys + 1
    if(ddrange(6).eq.ny) ye = ye - 1
  
!
! Now compute derivatives for interior grid
!
    do j = ys, ye
    do k = drange(3), drange(4)
    do i = drange(1), drange(2)
       dlgdx(i,k,j) = (lambdag(i+1,k,j)-lambdag(i-1,k,j))*hx
       dlgdy(i,k,j) = (lambdag(i,k,j+1)-lambdag(i,k,j-1))*hy
       dlgdz(i,k,j) = (lambdag(i,k+1,j)-lambdag(i,k-1,j))*hz
    end do
    end do
    end do

  endif   !if (lavg.eq.0)
  
!---------------------------------------------------------------------
  RETURN
END SUBROUTINE LAMBDA
!**********************************************************************
