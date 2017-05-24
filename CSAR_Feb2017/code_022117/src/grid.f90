! ********************************************************************
SUBROUTINE GRIDSETUP

  USE data_types
  USE GLOBAL_DATA

  IMPLICIT NONE

!-------------------------------------------------------------------------
!     Local Variables

  TYPE(mesh_pointers), POINTER :: current_mesh 
  INTEGER ::  i, j, k , ii,jj,kk, m,l
  INTEGER ::  npp,nppp,n,jsliceviz
  INTEGER,allocatable :: floc(:)
  REAL*8  :: den, der1, der2,const
  REAL*8  :: eta, eta2
  REAL*8  ::  dyloc,distloc

  integer :: mnx,mnz, mny,xe,ze,ye
  REAL*8  :: xstart,xx,yy,zz,dyc,dxc,dzc,dxf,dzf,dyf,add_hdx(3),dxv(3)
!-------------------------------------------------------------------------

  xstart = merge(0d0,-xend,issymmetric)
! x-grid; uniform and with no mapping
  DO i = lbound(x,1), ubound(x,1)
     x(i) = xstart + dfloat(i+dSerialRange(1))*dx
  END DO


! z-grid; uniform and with no mapping
  DO k = lbound(z,1), ubound(z,1)
     if(is2d) then
        z(k) = zero
     else
        z(k) = -zend + dfloat(k+dSerialRange(3))*dz
     end if
  END DO
  do k=1,num_mg_meshes
     do i = -1,nx_mg(k)
        ii = i+ib_mg(k)
        x_mg(i,k) = dble(ii)*dx_mg(k) + xstart
     end do
  end do
  if(is2d) then
     z_mg = 0d0
  else 
     do k=1,num_mg_meshes
        do i = -1,nz_mg(k)
           ii = i+kb_mg(k)
           z_mg(i,k) = dble(ii)*dz_mg(k) -zend
        end do
     end do
  end if


! y-grid
  DO j = 0, ny
     do m =1,2
        eta  = dfloat(j-dSerialRange(5))*dy + (m-1)*dy/two 
        MQeta(j,m) = eta   !duplicate it
        eta2 = (eta/yend)**2
        den = 2.0d0 - eta2
        if(m == 1) then
           y(j)  = eta/den**c1y
        else
           ya(j) = eta/den**c1y
        endif
        der1 = 1.0d0/(den**c1y)+2.0d0*eta2*c1y/den**(c1y+1.0d0)
        der2 = 2.0d0*eta*c1y/den**(c1y+1.0d0)&
             + 4.0d0*eta*c1y/den**(c1y+1.0d0)&
             + 4.0d0*eta*eta2*c1y*(c1y+1.0d0)/den**(c1y+2.0d0)
        der2 = der2/(yend**2)
        if( m == 1 ) then
           detady(j) = 1.0d0/der1
           deta2dy2(j) = -der2/der1**3.0d0
        else
           detadya(j) = 1.0d0/der1
           deta2dy2a(j) = -der2/der1**3.0d0
        endif
     enddo
     if (myid==0) write(21,*) j,y(j),ya(j)
  END DO
  close(21)
   MQeta(-1,:) = MQeta(0,m)


  jstart  = ny+2  !ceiling(dble(ny)/two)
  jend    = ny+2
  ystart  = MQeta(jstart,1)
  detabar = MQeta(jend,1) - ystart;
  const=-one/detabar**3;

  do m =1,2
     DO j = 0, ny
        deta = MQeta(j,m)-ystart
        if( j < jstart .or. j > jend) then
           MQchi(j,m) = merge(one,zero,j<jend .or. j < jstart)
           MQchipr(j,m) = zero
           MQchisc(j,m)=  zero
        else
           MQchi(j,m) = one + const * deta**3.
           MQchipr(j,m) = const*three * deta**2d0 * merge(detady(j),detadya(j),m==1)
           MQchisc(j,m)=  const*6.0d0 * deta * merge(detady(j)**2,detadya(j)**2,m==1)&
                &+  const*three * deta**2d0 * merge(deta2dy2(j),deta2dy2a(j),m==1)
        endif
     ENDDO
  enddo


  current_mesh => finest_mesh
  dyf = current_mesh%dy
  do j = -1, ny-1
     jj=j+1
     current_mesh%y(j) = y(jj)
     current_mesh%ya(j) = ya(jj)
     current_mesh%eta(j) = MQeta(jj,1)
     current_mesh%ey(j) = detady(jj)
     current_mesh%eya(j) = detadya(jj)
     current_mesh%chi(j,1:2) = MQchi(jj,1:2)
     current_mesh%chipr(j,1:2) = MQchipr(jj,1:2)
  enddo
  current_mesh%chi(-1,:)=one
  current_mesh%chipr(-1,:)=zero


  DO

     IF ( .NOT. ASSOCIATED(current_mesh%coarse_mesh) ) EXIT
     current_mesh => current_mesh%coarse_mesh

     ye = current_mesh%ye
     dyc = current_mesh%dy
     mny = anint(dyc/dyf)

     do j=-1,ye
        jj=j*mny+1
        if(jj > 0) then
           current_mesh%y(j) = y(jj)
           current_mesh%ya(j) = ya(jj)
           current_mesh%eta(j) = MQeta(jj,1)
           current_mesh%ey(j) = detady(jj)
           current_mesh%eya(j) = detady(jj+mny/2)
           current_mesh%chi(j,1:2) =  (/ MQchi(jj,1), MQchi(jj+mny/2,1)/)
           current_mesh%chipr(j,1:2) =  (/ MQchipr(jj,1), MQchipr(jj+mny/2,1)/)
        else
           current_mesh%y(j) = y(0)
           current_mesh%ya(j) = y(0)
           current_mesh%eta(j) = y(0)
           current_mesh%ey(j) = detady(0)
           current_mesh%eya(j) = detady(0)
           current_mesh%chi(j,1:2) = one
           current_mesh%chipr(j,1:2) = zero
        endif
        !!current_mesh%eya(j) = dyc/(y(jj+mny)-y(jj))
        !!current_mesh%ey(j) = (two*dyc)/(y(jj+mny)-y(jj-mny))
        !!current_mesh%eyy(j) = deta2dy2(jj)
     enddo
!current_mesh%ey(0) =  detady(1)
!current_mesh%eya(0) = detady(1+mny/2)
     current_mesh%eya(0) = dyc/(y(1+mny)-y(1))
     current_mesh%ey(0) =  dyc/(y(1+mny)-y(1))
     current_mesh%eyy(0) =  deta2dy2(1)

     current_mesh%chi(-1,:)=one
     current_mesh%chipr(-1,:)=zero
  END DO


!modifications to plot values on y=const planes
  jsliceviz = ny
  allocate(floc(0:nyloc))
  npp = floor(dble(nyloc)/dble(inyloc))
  nppp = 0
  floc = 1
  do n = 0,inyloc-1
     nppp = nppp + 2**n
     floc(n*npp+1:(n+1)*npp) = 2**n
  enddo
  dyloc = y(jsliceviz)/dble(npp*nppp)
  if(inyloc*npp < nyloc) floc(inyloc*npp+1:nyloc) = 1
  floc = cshift(floc,-nyloc+inyloc*npp)

  yloc(0) = zero
  do jj = 1,nyloc
     yloc(jj) = yloc(jj-1) + floc(jj)*dyloc 
  enddo
!
! y-grid  2 loops becasue of staggering
  do jj = 0,nyloc
     distloc = 1d10
     DO j = 0, ny
        if(abs(y(j) - yloc(jj)) < distloc) then
	   distloc = abs(y(j) - yloc(jj))
           jloc(jj) = j
        endif
     ENDDO
  enddo
!
  if(myid == 0) then
     open(UNIT=18,file = 'ycoord.dat',STATUS='UNKNOWN')
     open(UNIT=19,file = 'jcoord.dat',STATUS='UNKNOWN')
     write(18,'(1p2000e12.4)')(y(jloc(l)),l=0, nyloc)  !write it out every time
     write(19,'(2000i5)')(jloc(l),l=0, nyloc)  !write it out every time
     close(18)
     close(19)
  endif
  if(myid == 0) then
     write(6,*)'lengthscale=',lengthscale
     write(20,*)lengthscale
     write(20,*)nx
     write(20,*)nz
     write(20,*)ny
     close(20)
  endif

  nxv = nx
  nzv = nz
  nyv = ny


  call mpi_barrier(comm3d,i)


  dxf = finest_mesh%dx ; dzf = finest_mesh%dz; dyf = finest_mesh%dy
  dxv = (/dx,dz,dy/)
  DO N = 0,3
     current_mesh => finest_mesh
     add_hdx = 0
     if(N>0) add_hdx(N) = dxv(N)/2d0
     DO
        dxc = current_mesh%dx ;dzc = current_mesh%dz ;
        mnx = anint(dxc/dxf); mnz = anint(dzc/dzf);
        m = current_mesh%mesh_num  !anint( log10 (dble(mnx)) / log10 (2.0) )
        xe = current_mesh%xe ; ze = current_mesh%ze; ye = current_mesh%ye
        DO j = -1,ye+1
           yy = current_mesh%eta(j) + add_hdx(3)
           jj = lbound(Mqeta,1)
           do while((yy-Mqeta(jj,1))*(yy-Mqeta(jj+1,1)) >0 .and. jj <= ny .and. yy > Mqeta(jj,1) )
              jj = jj+1
           enddo
           jj = min(jj,ubound(Mqeta,1)-1)
           do i=-1,xe+1
              xx = x_mg(i,m) + add_hdx(1)
              ii = lbound(x,1)
              do while((xx-x(ii))*(xx-x(ii+1)) >0 .and. ii <= nx_mg(1) .and. xx > x(ii) )
                 ii = ii+1
!!>                 if(myid == 0) print*,i,ii,xx,x(ii),x(ii+1)
              enddo
              ii = min(ii,ubound(x,1)-1)
!!>              if(myid == 2) print*,'------------',i,ii,xx,x(ii),x(ii+1),xe,nx_mg(m),current_mesh%mesh_num
              do k=-1,ze+1
                 zz = z_mg(k,m) + add_hdx(2)
                 kk = lbound(z,1)
                 do while((zz-z(kk))*(zz-z(kk+1)) >0 .and. kk <= nz_mg(1) .and. zz > z(kk) )
                    kk = kk+1
                 enddo
                 kk = min(kk,ubound(z,1)-1)
                 current_mesh%c2f(N,i,k,j)%indx(1:3) = (/ii,kk,jj/)
                 current_mesh%c2f(N,i,k,j)%alpha(1:3) = (/ (x(ii+1)-xx)/&
                 &max(x(ii+1)-x(ii),1d-16),(z(kk)-zz)/max(z(kk+1)-z(kk),1d-16),&
                 &(MQeta(jj+1,1)-yy)/max(MQeta(jj+1,1)-MQeta(jj,1),1d-16) /)
!!>                 write(1000*myid + 100*current_mesh%mesh_num+10*N,'(3i3,1p3e12.4,3i3,5e12.4)') i,k,j,xx,zz,yy,current_mesh%c2f(N,i,k,j)%indx(1:3),&
!!>                      &                                         current_mesh%c2f(N,i,k,j)%alpha(1:3),x(current_mesh%c2f(N,i,k,j)%indx(1)),&
!!>                      &x(current_mesh%c2f(N,i,k,j)%indx(1)+1)
              end do
           end do
        end DO

!!>        if(m>1) close(1000*myid + 100*current_mesh%mesh_num+10*N)
        IF ( .NOT. ASSOCIATED(current_mesh%coarse_mesh) ) EXIT
        current_mesh => current_mesh%coarse_mesh

     END DO
  END DO

  call mpi_barrier(comm3d,i)

  if (myid == 0) then
     write(*,'(80("%")/,a,/,20("-"))')'output from Grid:'
     write(6,'(a,3i4)')'Overall grid Dimensions',nx,nz,ny
     write(6,'(a,1p123e12.4)')'Pressure, atm ',press
     write(6,'(a,1p123e12.4)')'Smallest dy = ',y(1)-y(0)
     write(6,'(a,1p123e12.4)')'Smallest dx = ',x(1)-x(0)
     write(6,'(a,1p123e12.4)')'Smallest dz = ',z(1)-z(0)
     write(6,'(a,1p123e12.4)')'Timestep    = ',timestep
     write(*,'(20("-"),/a/,80("%"))')'END GRID'
  endif


!-------------------------------------------------------------------
  RETURN
END SUBROUTINE GRIDSETUP
!*******************************************************************
