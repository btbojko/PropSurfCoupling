! ********************************************************************
SUBROUTINE GET_PROC_DECOMP
 
  USE GLOBAL_DATA
  USE MYMPI
  
  IMPLICIT NONE

!----------------------------------------------------------------------
! Local Variables
  INTEGER ::  nxprocs, pleft
!----------------------------------------------------------------------

  pdims(3) = nyprocs
!
! divide processors between x and z directions
! for now we assume tht nyproc ==1
!
  nxprocs = int(sqrt(dble(nproc/nyprocs))+0.5)
  pleft = mod(nproc,nxprocs)
  if( pleft == 0 ) then
    pdims(1) = nxprocs
    pdims(2) = nproc/nxprocs
  else
    do while ( pleft /= 0)
       nxprocs = nxprocs - 1 
       pleft = mod(nproc,nxprocs)
    enddo
    pdims(1) = nxprocs; pdims(2) = nproc/nxprocs
  endif

!
! TWO DIMENSIONAL GRID (comment out if not 2D)
!
  nxprocs = nproc/nyprocs
  pdims(1) = nxprocs; pdims(2) = 1

  if(pdims(1).lt.0.or.pdims(2).lt.0.or.pdims(3).lt.0) then
     if(mywid.eq.0) then
        write(*,*) 'Rocfire3D :: FATAL ERROR! bad nprocs/decomposition'
     endif
     STOP
  endif

!
! DETERMINE PARALLELIZATION
!
  isparallel(:)=pdims(:)>1

!
! CREATE A CARTESIAN TOPOLOGY
!
  isperiodic(1) = .TRUE.
  isperiodic(2) = .TRUE.
  isperiodic(3) = .FALSE.
 
  call MPI_CART_CREATE(MPI_COMM_WORLD, 3, pdims, isperiodic,&
       .TRUE., comm3d, ierr)
  call MPI_COMM_RANK(comm3d,myid,ierr)
 
! PRF: Get the neighboring process ranks, and return them in dnbr.
! PRF: dnbr(1) = xleft, dnbr(2) = xright
! PRF: dnbr(3) = zleft, dnbr(4) = zright
! PRF: dnbr(5) = bottom, dnbr(6) = top
! PRF: Note that left and bottom are defined to be in the negative
! PRF: index direction.
!
  call MPI_CART_SHIFT(comm3d,0,1,dnbr(1),dnbr(2),ierr)
  call MPI_CART_SHIFT(comm3d,1,1,dnbr(3),dnbr(4),ierr)
  call MPI_CART_SHIFT(comm3d,2,1,dnbr(5),dnbr(6),ierr)

  call MPI_CART_COORDS(comm3d,myid,3,coords,ierr)

!
! WRITE OUT SOME INFORMATION
!
  write(*,*)  myid, ':: coordinates: (',coords(1),&
       ',',coords(2),',',coords(3),')'
!    write(*,*)  myid, ':: Xrange :=', drange(1) , '-', &
!         drange(2)
!    write(*,*)  myid, ':: Zrange :=', drange(3) , '-', &
!         drange(4)
  !  write(*,*)  myid, ':: Yrange :=', drange(5) , '-', &
  !       drange(6)
  !  write(*,*) myid, ':: X neigbors:',dnbr(1),dnbr(2)
  !  write(*,*) myid, ':: Z neigbors:',dnbr(3),dnbr(4)
  !  write(*,*) myid, ':: Y neigbors:',dnbr(5),dnbr(6)
  !


!----------------------------------------------------------------------------
  RETURN
END SUBROUTINE GET_PROC_DECOMP
!*****************************************************************************
