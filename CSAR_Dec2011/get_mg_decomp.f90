! ********************************************************************
! This subroutine is the heart of the parallelization
! here the domain decomposition is established and send receive list
! are built
! ********************************************************************
SUBROUTINE GET_MG_DECOMP
 
  USE GLOBAL_DATA
  USE data_structure
  USE LUSOLVER
  USE MYMPI
  
  IMPLICIT NONE

!----------------------------------------------------------------------
! Local Variables
  INTEGER ::  k, error, status(MPI_STATUS_SIZE)
  INTEGER ::  crs,crsx,crsy
  LOGICAL ::  ibsodd,kbsodd
  REAL*8  ::  x1,x2
  TYPE(mesh_pointers),POINTER :: current_mesh
!----------------------------------------------------------------------


!
! DETERMINE THE receving map, given the send direction, 1 left,2 right, 3 up, 4 down
! finds the neighbor from who each proc receives
!
  rmap(1)=2;rmap(2)=1;rmap(3)=4;rmap(4)=3;

!................................................................................
!- ixp== #intervals on the coarsest, ixp=2 => 3 pts
!................................................................................
  ixp = nnxx;jyq = nnyy
  iex=lvlx;jey=lvly
  isoddcoarse = 2*int(ixp/2) < ixp

  if(ixp*jyq <= 0) THEN
    write(*,*)'INCORRET DOMAIN DECOMPOSITION',myid,ixp,jyq
    call MPI_Abort (MPI_COMM_WORLD, 1, ierr)
    STOP
  endif

  num_meshes = max( iex , jey )
  num_mg_meshes = num_meshes
  mg_maxcycle = 3   !w\ iguess =2 does actually mg_maxcycle-1 MG cycles
  gamma = 1; iguess = 2; icoarse=0
  nu1=3; nu2=3; nu_fine = 3
  
!
!-SUBGRID SET >> 1  @ finest ------- num_mg_meshes @ coarsest
!
  ALLOCATE(nx_mg(num_mg_meshes),nz_mg(num_mg_meshes),STAT=error)
  ALLOCATE(ib_mg(num_mg_meshes),kb_mg(num_mg_meshes),STAT=error)
  ALLOCATE(id_mg(num_mg_meshes),kd_mg(num_mg_meshes),STAT=error)
  ALLOCATE(ny_mg(num_mg_meshes),nx_tot(num_mg_meshes),STAT=error)
  ALLOCATE(nx_tot_fix(num_mg_meshes),nx_mg_fix(num_mg_meshes),STAT=error)
  ALLOCATE(nz_tot(num_mg_meshes),STAT=error)
  ALLOCATE(dx_mg(num_mg_meshes),dz_mg(num_mg_meshes),dy_mg(num_mg_meshes),STAT=error)
  ALLOCATE(iscoarse_xz(num_mg_meshes),iscoarse_y(num_mg_meshes))
  ALLOCATE(lx_mg(iex))

!
! TWO COARSNING STRATEGIES, the first is cheaper.
!
  if(icoarse == 0) then
    do k=1,num_mg_meshes
       nx_tot_fix(k) = ixp*2**(max(iex-k,0)) !number of intervals
       ny_mg(k) = jyq*2**(max(jey-k,0)) 
       lx_mg(min(iex,k))=min(iex,k)
    enddo
  else
    do k=1,num_mg_meshes
       crs=num_mg_meshes-k
       crsx=min(crs,iex-1)
       crsy=min(crs,jey-1)
       nx_tot_fix(k) = ixp*2**crsx           !number of intervals
       ny_mg(k) = jyq*2**crsy
       lx_mg(max(k+iex-num_mg_meshes,1)) = k
    enddo
  endif


  nz_tot(:) = nx_tot_fix(:)

!
! define fine grid total dimensions: nx is later defined
!
  ny = ny_mg(1)+1
  nz_tot(:) = 1
  nz =  nz_tot(1)

!
! Local processor grid dimensions
!
  nx_mg_fix(1) = int(nx_tot_fix(1)/pdims(1)) + min(max( mod(nx_tot_fix(1),pdims(1)) - coords(1) , 0 ) , 1)
  nz_mg(1) = int(nz_tot(1)/pdims(2)) + min(max( mod(nz_tot(1),pdims(2)) - coords(2) , 0 ) , 1)
  ib_mg(1) = int(nx_tot_fix(1)/pdims(1))*coords(1)+min( mod(nx_tot_fix(1),pdims(1)) , coords(1) )
  kb_mg(1) = int(nz_tot(1)/pdims(2))*coords(2)+min( mod(nz_tot(1),pdims(2)) , coords(2) )

  do k=2,num_mg_meshes
     if(nx_tot_fix(k) == nx_tot_fix(k-1)) then
       nx_mg_fix(k) = nx_mg_fix(k-1)   ! corarsen in y only, leave x,z unchanged
       ib_mg(k) = ib_mg(k-1)
     else
       ibsodd = 2*int(ib_mg(k-1)/2) < ib_mg(k-1)
       nx_mg_fix(k) = nx_mg_fix(k-1) - merge(int((nx_mg_fix(k-1)+1)/2),int(nx_mg_fix(k-1)/2),ibsodd)
       ib_mg(k) = int((ib_mg(k-1)+1)/2)
    endif
    if(nz_tot(k) == nz_tot(k-1)) then
       nz_mg(k) = nz_mg(k-1)
       kb_mg(k) = kb_mg(k-1)
     else
       kbsodd = 2*int(kb_mg(k-1)/2) < kb_mg(k-1)
       nz_mg(k) = nz_mg(k-1) - merge(int((nz_mg(k-1)+1)/2),int(nz_mg(k-1)/2),kbsodd)
       kb_mg(k) = int((kb_mg(k-1)+1)/2)
    endif
  enddo
!
! If a point shifts grid when interpolating, flag it so to correctly interpolate
! its value
!
  id_mg=0;kd_mg=0
  do k=1,num_mg_meshes-1
     id_mg(k) = 2*ib_mg(k+1)-ib_mg(k)  
     kd_mg(k) = 2*kb_mg(k+1)-kb_mg(k)
  enddo

!
! set logical flags that keep track of weather there is coarsening in xz,y or both
!
  iscoarse_xz = .FALSE. ;  iscoarse_y = .FALSE. 
  do k=1,num_mg_meshes-1
     if( nx_tot_fix(k) > nx_tot_fix(k+1) ) iscoarse_xz(k) = .TRUE.
     if( ny_mg(k) > ny_mg(k+1) ) iscoarse_y(k) = .TRUE.
  enddo

!<> set bc types

  CALL SELECT_BC_TYPE

!
! COMPUTE GRID SPACING 
!
  if(issymmetric)then
     x1 = zero
  else
     x1 = -xend
  endif

  x2 = xend
  dx = (x2-x1) / (dble(nx_tot(1)) - sum(addbnd(1:2)) )

  dy = yend/dble(ny)

  xstart_MG =  x1 + addbnd(1)*dx


  do k=1,num_mg_meshes
     dx_mg(k) = dx*dble(nx_tot_fix(1)/nx_tot_fix(k))
     dz_mg(k) = dz*dble(nz/nz_tot(k))
     dy_mg(k) = dy*dble(ny_mg(1)/ny_mg(k))
!     dy_mg(k) = (yend-dy)/dble(ny_mg(k)*pdims(3))
  enddo

!
! Global grid dimensions
!
  nx = nx_tot(1)
  nz = nx


  drange(1) = int(nx_tot_fix(1)/pdims(1))*coords(1)+min( mod(nx_tot_fix(1),pdims(1)) , coords(1) )
  drange(2) = drange(1) + nx_mg(1)-1     !note nx_mg not nx_mg_fix
  drange(3) = int(nz/pdims(2))*coords(2)+min( mod(nz,pdims(2)) , coords(2) )
  drange(4) = drange(3) + nz_mg(1)-1
  drange(5) = 0
  drange(6) = ny

  IF(is2D) THEN   !NEEDED BECAUSE I evaluate c_ijk = 4*c1+...; assuming dx=dz
     nz_tot(:) = 1
     nz =  nz_tot(1)
     drange(3) = 0
     drange(4) = 0
     dz = dx
     dz_mg = dx_mg
  ENDIF


  CALL storage_allocate

!
! allocate metric quantities used in the MG algorithm
!
  IF(num_mg_meshes == 1)THEN  !gross
    ALLOCATE(dphidx_mg(-1:0+1,-1:0+1,iex-1),STAT=error)
    ALLOCATE(dphidz_mg(-1:0+1,-1:0+1,iex-1),STAT=error)

    ALLOCATE(dphidxa_mg(-1:0,-1:0,iex-1),STAT=error)
    ALLOCATE(dphidza_mg(-1:0,-1:0,iex-1),STAT=error)
  else
    ALLOCATE(dphidx_mg(-1:nx_mg(2)+1,-1:nz_mg(2)+1,iex-1),STAT=error)
    ALLOCATE(dphidz_mg(-1:nx_mg(2)+1,-1:nz_mg(2)+1,iex-1),STAT=error)

    ALLOCATE(dphidxa_mg(-1:nx_mg(2),-1:nz_mg(2),iex-1),STAT=error)
    ALLOCATE(dphidza_mg(-1:nx_mg(2),-1:nz_mg(2),iex-1),STAT=error)
  endif
!
! define red and black point lists used in gauss_relax
!
  CALL rednblack_lists

  mg_prm(1)=ixp;mg_prm(2)=jyq;mg_prm(3)=iex;mg_prm(4)=jey


!--- ASSIGN RANGE

  xr(1) = (drange(2) - drange(1)) + 1
  ddrange(1) = drange(1)
  ddrange(2) = drange(2)
  drange(1) = 0
  drange(2) = xr(1) - 1
  xr(2) = (drange(4) - drange(3)) + 1
  ddrange(3) = drange(3)
  ddrange(4) = drange(4)
  drange(3) = 0
  drange(4) = xr(2) - 1
  xr(3) = (drange(6)-drange(5))+1
  ddrange(5) = drange(5)
  ddrange(6) = drange(6)
  drange(5) = 0
  drange(6) = xr(3) - 1

!
! MODIFICATIONS FOR NEUMANN WEST/EAST BC
!
  nxv(0:4) = drange(2)
  nyv(0:4) = ny-1 
  if(is_proc_bnd(2)) then
     nxv(0:4) = drange(2) + addptx(0:4)
     xr(1) = xr(1) + maxval(addptx(0:4))
  endif

  drange(2) = nxv(0)

!
! MAXRANGE
!
  maxrange = maxval(xr)
  call rednblack_lists_eqn  !uses nxv,nyv 


  if(ipack == -1) return
  if(.NOT. do_solve_direct) THEN
      LU_frequency = 1;
      residual_direct = 1.d-10;
      RETURN
  endif


  if(num_mg_meshes == 1) then
     current_mesh => finest_mesh
  else
     current_mesh => coarsest_mesh
  endif
!
! ALLOCATE SPARSE LU MATRIX
!
  CALL LU_ALLOC(nx_tot(num_mg_meshes), current_mesh%ye+1, 9, 150 ,&
      current_mesh%xe+1, current_mesh%ye+1 )
!----------------------------------------------------------------------------
  RETURN
END SUBROUTINE GET_MG_DECOMP
!*****************************************************************************
