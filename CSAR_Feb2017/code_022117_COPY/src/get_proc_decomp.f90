! ********************************************************************
SUBROUTINE GET_PROC_DECOMP(flag)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!----------------------------------------------------------------------
! Local Variables
  INTEGER :: i,j,k,l,m,n,pcs,flag
  INTEGER ::  nxprocs, pleft
  integer elcount,myzsize,el,ixz, myxsize
  integer elm,elsizex, lastd,lastd1,elsizez
  integer mydrange(6),mydrange1(4),mydrange11(11)
  integer myysize,mysizeC
  integer,allocatable :: tmpdrange11(:)
  INTEGER :: error
!----------------------------------------------------------------------

  if(flag == 2) goto 110

  pdims(3) = nyprocs
!
! divide processors between x and z directions
! for now we assume tht nyproc ==1
!
  nxprocs = ceiling( sqrt(dble(nproc/nyprocs)))
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
  if(is2d) then
     nxprocs = nproc/nyprocs
     pdims(1) = nxprocs; pdims(2) = 1
  end if

  if(pdims(1).lt.0.or.pdims(2).lt.0.or.pdims(3).lt.0) then
     if(myid.eq.0) then
        write(*,*) 'Rocfire3D :: FATAL ERROR! bad nprocs/decomposition'
     endif
     STOP
  endif

!
! DETERMINE PARALLELIZATION
!
  isparallel(1:3)=pdims(1:3)>1

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


  if(is2d) then
     rowcomm = comm3d
     rowid   = myid
     rowsize = nproc
  else
! PRF: Get the communicator for the procs in my X row
     remain(1) = .true.
     remain(2) = .false.
     remain(3) = .true.
     call MPI_CART_SUB(comm3d,remain,rowcomm(1),ierr)
     call MPI_COMM_RANK(rowcomm(1),rowid(1),ierr)
     call MPI_COMM_SIZE(rowcomm(1),rowsize(1),ierr)


     remain(1) = .false.
     remain(2) = .true.
     remain(3) = .true.
     call MPI_CART_SUB(comm3d,remain,rowcomm(2),ierr)
     call MPI_COMM_RANK(rowcomm(2),rowid(2),ierr)
     call MPI_COMM_SIZE(rowcomm(2),rowsize(2),ierr)
  end if



  if(myid == 0) then
     write(*,'(80("+")/,a,/,20("-"))')'output from GetProcDecomp:'
     write(6,*)'is2d',is2d
     write(6,*)'pdims',pdims
     write(6,*)'isperiodic',isperiodic
     write(*,'(20("-"),/a/,80("+"))')'END GetProcDecomp:'
  endif

  

  if (flag == 1) return

110 continue

  allocate(displacements(nproc+1),STAT=error) 
  allocate(displacements1(nproc+1),STAT=error) 
  allocate(displacements2(nproc+1),STAT=error) 

  allocate(alldrange(4*nproc),STAT=error)
  allocate(alldrange1(4*nproc),STAT=error)
  allocate(alldrange2(4*nproc),STAT=error)

  allocate(allsize(nproc+1),STAT=error)  !2D
  allocate(allsize1(nproc+1),STAT=error)  !1D
  allocate(allsize2(nproc+1),STAT=error)  

  allocate(tmpdrange11(11*nproc),STAT=error)
  allocate(alldrange11(11,nproc),STAT=error)
  alldrange = -111


  mydrange(1) = drange(1)
  mydrange(2) = drange(2)
  if(rowid(1) == rowsize(1)-1) mydrange(2) = drange(2)+1

  mydrange(3) = drange(3)
  mydrange(4) = drange(4)
  if(rowid(2) == rowsize(2)-1) mydrange(4) = drange(4)+1

  mydrange(5:6) = rowid(1:2)


  mydrange11(1:6) =  mydrange
  mydrange11(7) = ib_mg(1)
  mydrange11(8) = kb_mg(1)
  mydrange11(9:11) = xr(1:3)

  myxsize = (mydrange(2) - mydrange(1)) + 1
  myzsize = (mydrange(4) - mydrange(3)) + 1
  mysizeB = myxsize * myzsize
  call MPI_ALLGATHER(mysizeB,1,MPI_INTEGER,allsize,1,MPI_INTEGER,&
       comm3d,ierr)
  call MPI_ALLGATHER(mydrange,6,MPI_INTEGER,alldrange,6,MPI_INTEGER,&
       comm3d,ierr)
  call MPI_ALLGATHER(mydrange11,11,MPI_INTEGER,tmpdrange11,11,MPI_INTEGER,&
       comm3d,ierr)
  alldrange11 = reshape(tmpdrange11,shape = (/11,nproc/))
  deallocate(tmpdrange11)

  displacements(1) = 0
  do i = 2,nproc
     displacements(i) = allsize(i-1) + displacements(i-1)
  enddo

  lastd =  allsize(nproc) + displacements(nproc)

 
!--- ALLOCATION
!
  allocate(yloc(0:nyloc),jloc(0:nyloc),STAT=error)
  call CHECK_ERROR(error,'G_allocation')

  allocate( myphi1(0:nx+1),STAT=error)
  allocate(allphi1(0:nx+1),STAT=error)
  call CHECK_ERROR(error,'G_allocation')

  allocate(myphi2(0:mysizeB),STAT=error)
  allocate( allphi2(0:lastd+1),STAT=error)
  call CHECK_ERROR(error,'G_allocation')

  ALLOCATE(G_x(0:nx),STAT=error)
  call CHECK_ERROR(error,'G_allocation')

  ALLOCATE(G_phi(-1:nx+1,-1:nz+1),STAT=error)
  call CHECK_ERROR(error,'G_allocation')

  ALLOCATE(G_z(0:nz),STAT=error)
  call CHECK_ERROR(error,'G_allocation')

  maxQvars = 4
  maxFvars = 5

  ALLOCATE(G_f(0:nx,0:nz,0:ny),STAT=error)
  call CHECK_ERROR(error,'G_allocation')

  ALLOCATE(G_rate(-1:nx+1,-1:max(ny,nz)+1,0:max(3,maxsize)),STAT=error)  ! here I assume ny < nz sice is use in both
  call CHECK_ERROR(error,'G_allocation')

  ALLOCATE(G_psi(-1:nx+1,-1:nz+1),STAT=error)
  call CHECK_ERROR(error,'G_allocation')



  do pcs = 1, nproc-1
     elm = 6 * (pcs-1) +1
     el =  elm + 6
     elsizex = alldrange(el+1)-alldrange(el+0)
     elsizez = alldrange(el+3)-alldrange(el+2)
     if(alldrange(el+4) >  alldrange(elm+4)) then
        alldrange(el+0) = alldrange(elm+1) + 1
     else
        alldrange(el+0) =  alldrange(elm+0+ 6* (alldrange(el+4) -&
             alldrange(elm+4)) )
     endif
     alldrange(el+1) = alldrange(el+0) + elsizex

     if(alldrange(el+5) >  alldrange(elm+5)) then
        alldrange(el+2) = alldrange(elm+3) + 1
     else
        alldrange(el+2) =  alldrange(elm+2+ 6* (alldrange(el+5) -&
             alldrange(elm+5)) )
     endif

     alldrange(el+3) = alldrange(el+2) + elsizez
  enddo


!----------------------------------------------------------------------------
  RETURN
END SUBROUTINE GET_PROC_DECOMP
!*****************************************************************************
