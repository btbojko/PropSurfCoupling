! *****************************************************************
Subroutine READ_matrix
  
  USE GLOBAL_DATA
  USE MYMPI
  
  IMPLICIT NONE

!-------------------------------------------------------------------------

! Local variables
  INTEGER :: i,j,k,l,m,n
!-------------------------------------------------------------------------
  open(789,file = 'rocfire.rflo.txt')
  read(789,*) matrix%nFourier,MImposed
  call allocate_matrix
  do i = 1,matrix%nFourier
     read(789,*)matrix%M(i),matrix%F(i),matrix%A(i)
  enddo
  close (789)


!----------------------------------------------------------------------------
  RETURN
END SUBROUTINE READ_MATRIX
!****************************************************************************
