! ********************************************************************
! This subroutine is the heart of the parallelization
! here the domain decomposition is established and send receive list
! are built
! ********************************************************************
SUBROUTINE GET_MG_DECOMP

  USE GLOBAL_DATA
  USE data_structure
  USE MYMPI

  IMPLICIT NONE

!----------------------------------------------------------------------
! Local Variables
  INTEGER ::  k, error, status(MPI_STATUS_SIZE)
  INTEGER ::  crs,crsx,crsy
  LOGICAL ::  ibsodd,kbsodd
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
     print '(a,1x,3i3)' ,'INCORRET DOMAIN DECOMPOSITION',myid,ixp,jyq
     call MPI_Abort (MPI_COMM_WORLD, 1, ierr)
     STOP
  endif

  num_meshes = max( iex , jey )
  num_mg_meshes = num_meshes
  mg_maxcycle = 4
  gamma = 1; iguess = 2; icoarse=0
  nu1=12; nu2=12; nu_fine = 10

!
!-SUBGRID SET >> 1  @ finest ------- num_mg_meshes @ coarsest
!
  ALLOCATE(nx_mg(num_mg_meshes),nz_mg(num_mg_meshes),STAT=error)
  ALLOCATE(ib_mg(num_mg_meshes),kb_mg(num_mg_meshes),STAT=error)
  ALLOCATE(id_mg(num_mg_meshes),kd_mg(num_mg_meshes),STAT=error)
  ALLOCATE(ny_mg(num_mg_meshes),nx_tot(num_mg_meshes),STAT=error)
  ALLOCATE(nz_tot(num_mg_meshes),STAT=error)
  ALLOCATE(dx_mg(num_mg_meshes),dz_mg(num_mg_meshes),dy_mg(num_mg_meshes),STAT=error)
  ALLOCATE(iscoarse_xz(num_mg_meshes),iscoarse_y(num_mg_meshes))
  allocate(gmres_size_mg(num_mg_meshes))
  ALLOCATE(lx_mg(iex))

!
! TWO COARSNING STRATEGIES, the first is cheaper.
!
  if(icoarse == 0) then
     do k=1,num_mg_meshes
        nx_tot(k) = ixp*2**(max(iex-k,0)) !number of intervals
        ny_mg(k) = jyq*2**(max(jey-k,0)) 
        lx_mg(min(iex,k))=min(iex,k)
        gmres_size_mg(k) = merge(2000,20,k==num_mg_meshes.and.k/=1)
     enddo
  else
     do k=1,num_mg_meshes
        crs=num_mg_meshes-k
        crsx=min(crs,iex-1)
        crsy=min(crs,jey-1)
        nx_tot(k) = ixp*2**crsx           !number of intervals
        ny_mg(k) = jyq*2**crsy
        lx_mg(max(k+iex-num_mg_meshes,1)) = k
     enddo
  endif
  nz_tot(:) = nx_tot(:)

!-- 2D only
  if(is2D) then
     nz_tot(:) = 1
  end if

!
! Global grid dimensions
!
  nx = nx_tot(1)
  nz = nz_tot(1)
  ny = ny_mg(1)+1

!
! Sub-grid dimensions
!
  nx_mg(1) = int(nx_tot(1)/pdims(1)) + min(max( mod(nx_tot(1),pdims(1)) - coords(1) , 0 ) , 1)
  nz_mg(1) = int(nz_tot(1)/pdims(2)) + min(max( mod(nz_tot(1),pdims(2)) - coords(2) , 0 ) , 1)
  ib_mg(1) = int(nx_tot(1)/pdims(1))*coords(1)+min( mod(nx_tot(1),pdims(1)) , coords(1) )
  kb_mg(1) = int(nz_tot(1)/pdims(2))*coords(2)+min( mod(nz_tot(1),pdims(2)) , coords(2) )

  do k=2,num_mg_meshes
     if(nx_tot(k) == nx_tot(k-1)) then
        nx_mg(k) = nx_mg(k-1)   ! corarsen in y only, leave x,z unchanged
        ib_mg(k) = ib_mg(k-1)
     else
        ibsodd = 2*int(ib_mg(k-1)/2) < ib_mg(k-1)
        nx_mg(k) = nx_mg(k-1) - merge(int((nx_mg(k-1)+1)/2),int(nx_mg(k-1)/2),ibsodd)
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
     if(nx_tot(k)>nx_tot(k+1)) iscoarse_xz(k) = .TRUE.
     if(ny_mg(k)>ny_mg(k+1)) iscoarse_y(k) = .TRUE.
  enddo

  drange(1) = int(nx/pdims(1))*coords(1)+min( mod(nx,pdims(1)) , coords(1) )
  drange(2) = drange(1) + nx_mg(1)-1
  drange(3) = int(nz/pdims(2))*coords(2)+min( mod(nz,pdims(2)) , coords(2) )
  drange(4) = drange(3) + nz_mg(1)-1
  drange(5) = 0
  drange(6) = ny

!
! COMPUTE GRID SPACINGS 
!

  IF(issymmetric) THEN
     dx = xend/(dble(nx))
     dz = zend/(dble(nz))
  ELSE
     dx = 2.0d0*xend/dble(nx)
     dz = 2.0d0*zend/dble(nz)
  ENDIF
  dy = yend/dble(ny)


  do k=1,num_mg_meshes
     dx_mg(k) = dx*dble(nx/nx_tot(k))
     dz_mg(k) = dz*dble(nz/nz_tot(k))
!     dy_mg(k) = dy*dfloat(ny_mg(1)/ny_mg(k))
     dy_mg(k) = (yend-dy)/dble(ny_mg(k)*pdims(3))
  enddo

  IF(is2D) THEN
     dz = dx
     dz_mg = dx_mg
  ENDIF

!    write(30+myid,*)'MG decomposition as follows'
!    do k=1,num_mg_meshes
!       write(*,*)'MESH # ',k,'NX',nx_mg(k),nx,'NZ',nz_mg(k),nz,'DX',dx_mg(k),dx,&
!                             'NY',ny_mg(k),ny,'DY',dy_mg(k),dy
!    enddo
!    write(30+myid,*)'DRANGE',drange(1),drange(2),drange(3),drange(4),drange(5),drange(6),nx_mg(:),coords
!  do k=1,num_mg_meshes
!    write(30+myid,*)myid,'ib,kb',ib_mg(k),kb_mg(k),'id,kd',id_mg(k),kd_mg(k),'MESH',k
!  enddo
!   write(30+myid,*)myid,'DNBR',dnbr(1),dnbr(2),dnbr(3),dnbr(4)

!
! ALLOCATE MEMORY ON EAH GRID USING NESTED POINTERS
!
  CALL storage_allocate

!
! allocate metric quantities used in the MG algorithm
!
  ALLOCATE(dphidx_mg(-1:nx_mg(1)+1,-1:nz_mg(1)+1,0:iex-1),STAT=error)
  ALLOCATE(dphidz_mg(-1:nx_mg(1)+1,-1:nz_mg(1)+1,0:iex-1),STAT=error)

  ALLOCATE(dphidxa_mg(-1:nx_mg(1),-1:nz_mg(1),0:iex-1),STAT=error)
  ALLOCATE(dphidza_mg(-1:nx_mg(1),-1:nz_mg(1),0:iex-1),STAT=error)
  ALLOCATE(phi_mg(-1:nx_mg(1),-1:nz_mg(1),0:iex-1),STAT=error)
  ALLOCATE(psi_mg(-1:nx_mg(1),-1:nz_mg(1),0:iex-1),STAT=error)
  ALLOCATE(temp_mg(-1:nx_mg(1),-1:nz_mg(1),0:iex-1),STAT=error)
  temp_mg = -1d0
  dphidx_mg = 0d0
  dphidz_mg = 0d0
  dphidxa_mg = 0d0
  dphidza_mg = 0d0
  phi_mg = 0d0
  psi_mg = 0d0

!
! define red and black point lists used in gauss_relax
!
  CALL rednblack_lists

  mg_prm(1)=ixp;mg_prm(2)=jyq;mg_prm(3)=iex;mg_prm(4)=jey


!-- ASSIGN RANGE
  xr(1) = (drange(2) - drange(1)) + 1
  dSerialRange(1) = drange(1)
  drange(1) = 0
  dSerialRange(2) = drange(2)
  drange(2) = xr(1) - 1
  xr(2) = (drange(4) - drange(3)) + 1
  dSerialRange(3) = drange(3)
  drange(3) = 0
  dSerialRange(4) = drange(4)
  drange(4) = xr(2) - 1
  xr(3) = (drange(6)-drange(5))+1
  dSerialRange(5) = drange(5)
  drange(5) = 0
  dSerialRange(6) = drange(6)
  drange(6) = xr(3) - 1

  dwiderange(1) = drange(1)-1
  dwiderange(2) = drange(2)+1
  dwiderange(3) = drange(3)-1
  dwiderange(4) = drange(4)+1
  dwiderange(5) = drange(5)-1
  dwiderange(6) = drange(6)+1

!
! MAXRANGE
!
  maxrange = maxval(xr)

!----------------------------------------------------------------------------
  RETURN
END SUBROUTINE GET_MG_DECOMP
!*****************************************************************************
