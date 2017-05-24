!****************************************************************************************
SUBROUTINE BR_CORR(num_p,num_d,pressure,diameter,rate_corr)

USE GLOBAL_DATA

IMPLICIT NONE
!---------------------------------------------------------------------
INTEGER, INTENT(OUT) :: num_p,num_d
REAL*8, INTENT(OUT) :: pressure(15),diameter(15),rate_corr(15,15)
! Local variables
INTEGER ::  i,j
CHARACTER :: comt*100

!---------------------------------------------------------------------

open(unit=1,FILE='correlations.txt',STATUS='old')

read(1,*) comt
read(1,*) num_p
do i=1,num_p
	read(1,*) pressure(i)
enddo
read(1,*)
read(1,*) num_d
do i=1,num_d
	read(1,*) diameter(i)
enddo 
read(1,*)
do i=1,num_p
	read(1,*) (rate_corr(i,j),j=1,num_d)
enddo
close(unit=1)

!---------------------------------------------------------------------
return
END SUBROUTINE BR_CORR
!*********************************************************************


SUBROUTINE POLIN2(x1a,x2a,yaa,m,n,x1,x2,yout)
!	Adapted from Numerical Recipes
!	Given arrays x1a(1:m) and x2a(1:n) of independent variables, and an m by n 
!	array of function values ya(1:m,1:n), tabulated at the grid defined by x1a and 
!	x2a; and given x1 and x2 values, returns an interpolated function value of y

USE GLOBAL_DATA

IMPLICIT NONE
!---------------------------------------------------------------------
INTEGER, INTENT(IN) :: m,n
REAL*8, INTENT(IN) :: x1a(15),x2a(15),yaa(15,15),x1,x2
REAL*8, INTENT(OUT) :: yout
! Local variables
INTEGER :: j,k
REAL*8 :: ymtmp(100),yntmp(100)
!---------------------------------------------------------------------

do j=1,m
    do k=1,n
	    yntmp(k)=yaa(j,k)
	enddo
	call lin_int(x2a,yntmp,n,x2,ymtmp(j))
enddo
call lin_int(x1a,ymtmp,m,x1,yout)

!---------------------------------------------------------------------	
return
END SUBROUTINE POLIN2
!*********************************************************************

SUBROUTINE LIN_INT(xa_in,ya_in,n,xin,yout)
!	Basic linear interpolation 

USE GLOBAL_DATA

IMPLICIT NONE
!---------------------------------------------------------------------
INTEGER, INTENT(IN) :: n
REAL*8, INTENT(IN) :: xa_in(n),ya_in(n),xin
REAL*8, INTENT(OUT) :: yout
! Local variables
INTEGER :: gg,i
!---------------------------------------------------------------------

gg=0
!	determine lower bound
do i=1,n
    if(xin-xa_in(i).ge.0.d0) gg=i
enddo
!	interpolate
yout=ya_in(gg)-(xa_in(gg)-xin)*(ya_in(gg)-ya_in(gg+1))/(xa_in(gg)-xa_in(gg+1))	

!---------------------------------------------------------------------
return
END SUBROUTINE LIN_INT
!*********************************************************************


