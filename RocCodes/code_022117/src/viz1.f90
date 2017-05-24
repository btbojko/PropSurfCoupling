Module Gather_R
  USE GLOBAL_DATA
  USE MYMPI
  USE data_types
  IMPLICIT NONE

contains

!*****************************************************************************
  subroutine gatherbuff(buff,meqn)

    IMPLICIT NONE

!     Local variables
    real*8,POINTER :: buff(:,:,:)
    integer,INTENT(IN) :: meqn
    integer i, j, k, n, myxsize,pcs
    integer elcount,myysize,el
    integer mysize,elm,elsize
    integer mydrange(4)

!---------------------------------------------------------------


    if(.not. associated(myphi)) allocate(myphi(0:(nx+1)*(ny+1)+1))
    if(.not. associated(allphi)) allocate(allphi(0:(nx+1)*(ny+1)+1))
    if(.not. allocated(allsize)) allocate(allsize(nproc+1))
    if(.not. allocated(displacements)) allocate(displacements(nproc+1))
    if(.not. allocated(alldrange)) allocate(alldrange(4*nproc))
    allsize = 0
    displacements = 0

    mydrange(1) = drange(1)
    mydrange(2) = drange(2)
    if(myid == nproc-1) mydrange(2) = drange(2)+1

    mydrange(3) = 0
    mydrange(4) = ny
    if(mydrange(1).eq.0) mydrange(1) = 0
    if(mydrange(2).eq.nx) mydrange(2) = nx 
    if(mydrange(3).eq.0) mydrange(3) = 0
    if(mydrange(4).eq.ny) mydrange(4) = ny
    myxsize = (mydrange(2) - mydrange(1)) + 1
    myysize = (mydrange(4) - mydrange(3)) + 1
    mysize = myxsize * myysize

    call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allsize,1,MPI_INTEGER,comm3d,ierr)

    call MPI_ALLGATHER(mydrange,4,MPI_INTEGER,alldrange,4,MPI_INTEGER,comm3d,ierr)

    displacements(1) = 0
    do i = 2,nproc
       displacements(i) = allsize(i-1) + displacements(i-1)
    enddo
    do pcs = 1, nproc-1
       elm = 1 + 4 * (pcs-1)
       el =  1 + 4 * pcs
       elsize = alldrange(el+1)-alldrange(el+0)
       alldrange(el+0) = alldrange(elm+1) + 1
       alldrange(el+1) = alldrange(el+0) + elsize
    enddo

    k = 0
    elcount = 0
    do i = mydrange(1),mydrange(2)
       do j = mydrange(3),mydrange(4)
          myphi(elcount) = buff(i,k,j)
          elcount = elcount + 1
       enddo
    enddo

    call MPI_ALLGATHERV(myphi(0),mysize,MPI_DOUBLE_PRECISION,&
         allphi(0),allsize,displacements,MPI_DOUBLE_PRECISION,comm3d, ierr)

    elcount = 0
    do pcs = 0, nproc-1
       el = 1 + (pcs*4)
       do i = alldrange(el),alldrange(el+1)
          do j = alldrange(el+2),alldrange(el+3)
             G_rate(i,j,meqn) = allphi(elcount)
             elcount = elcount + 1 
          enddo
       enddo
    enddo

    return
  end subroutine gatherbuff
!**********************************************************

!*****************************************************************************
  subroutine gatherbuff3(buff)

	! Gathers data from all processors into a single array for output.	
	! Takes whatever is in the local buff(i,k,j) and puts it into the
	! composite G_f(i,k,j).  To use this routine, copy your local 3D
	! scalar field into buff, call the routine, and then write G_f to a
	! file.

    IMPLICIT NONE

!     Local variables
    real*8,POINTER :: buff(:,:,:)
    integer i, j, k, n, myxsize,pcs
    integer elcount,myzsize,el,lastel
    integer mysize,elm,elsize
    integer mydrange(4),mySdrange(4)

!---------------------------------------------------------------


    if(.not. associated(myphi3)) allocate(myphi3(0:(nx+1)*(nz+1)*(ny+1)+1))
    if(.not. associated(allphi3)) allocate(allphi3(0:(nx+1)*(nz+1)*(ny+1)+1))
    if(.not. allocated(allsize)) allocate(allsize(nproc+1))
    if(.not. allocated(displacements)) allocate(displacements(nproc+1))
    if(.not. allocated(alldrange)) allocate(alldrange(4*nproc))
    allsize = 0
    displacements = 0

    mydrange(1) = drange(1)
    mydrange(2) = drange(2)
    mydrange(3) = drange(3)
    mydrange(4) = drange(4)

    if(mydrange(1).eq.0) mydrange(1) = 0
    if(mydrange(2)+ib_mg(1).eq.nx-1) mydrange(2) = mydrange(2)+1
    if(mydrange(3).eq.0) mydrange(3) = 0
    if(mydrange(4)+kb_mg(1).eq.nz-1) mydrange(4) = mydrange(4)+1

    myxsize = (mydrange(2) - mydrange(1)) + 1
    myzsize = (mydrange(4) - mydrange(3)) + 1
    mysize = myxsize * myzsize*(ny-0+1)

    call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allsize,1,MPI_INTEGER,comm3d,ierr)

    mySdrange(1:2) = mydrange(1:2) + ib_mg(1)
    mySdrange(3:4) = mydrange(3:4) + kb_mg(1)
    call MPI_ALLGATHER(mySdrange,4,MPI_INTEGER,alldrange,4,MPI_INTEGER,comm3d,ierr)

    displacements(1) = 0
    do i = 2,nproc
       displacements(i) = allsize(i-1) + displacements(i-1)
    enddo
    lastel = allsize(nproc) + displacements(nproc)

    elcount = 0
    do j =0,ny
       do k = mydrange(3),mydrange(4)
          do i = mydrange(1),mydrange(2)
             myphi3(elcount) = buff(i,k,j)
             elcount = elcount + 1
          enddo
       enddo
    enddo

    call MPI_ALLGATHERV(myphi3(0),mysize,MPI_DOUBLE_PRECISION,&
         allphi3(0),allsize,displacements,MPI_DOUBLE_PRECISION,comm3d, ierr)

    elcount = 0
    do pcs = 0, nproc-1
       el = 1 + (pcs*4)
       do j =0,ny
          do k = alldrange(el+2),alldrange(el+3)
             do i = alldrange(el),alldrange(el+1)
                G_f(i,k,j) = allphi3(elcount)
                elcount = elcount + 1 
             enddo
          enddo
       enddo
    enddo

    return
  end subroutine gatherbuff3
!**********************************************************

!*****************************************************************************
subroutine gatherbuffsurf(buff,buffO)

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

!     Local variables
  real*8,POINTER :: buff(:,:)
  real*8,POINTER :: buffO(:,:)
  integer i, j,k,l,m, n,pcs,flag,nIN
  integer elcount,myzsize,el,ixz, myxsize
  integer mysize,elm,elsize, lastd
  integer mydrange(4),qplet(4),mySdrange(4)
!---------------------------------------------------------------

  if(.not. allocated(displacements)) allocate(displacements(nproc+1))
  if(.not. allocated(alldrange)) allocate(alldrange(4*nproc))


  mydrange(1) = drange(1)
  mydrange(2) = drange(2)
  mydrange(3) = drange(3)
  mydrange(4) = drange(4)

  if(mydrange(1).eq.0) mydrange(1) = 0
  if(mydrange(2)+ib_mg(1).eq.nx-1) mydrange(2) = mydrange(2)+1
  if(mydrange(3).eq.0) mydrange(3) = 0
  if(mydrange(4)+kb_mg(1).eq.nz-1) mydrange(4) = mydrange(4)+1

  mySdrange(1:2) = mydrange(1:2) + ib_mg(1)
  mySdrange(3:4) = mydrange(3:4) + kb_mg(1)
!
  myxsize = (mydrange(2) - mydrange(1)) + 1
  myzsize = (mydrange(4) - mydrange(3)) + 1
  mysize = myxsize * myzsize
  call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allsize,1,MPI_INTEGER,comm3d,ierr)
  call MPI_ALLGATHER(mySdrange,4,MPI_INTEGER,alldrange,4,MPI_INTEGER,comm3d,ierr)

!
  displacements(1) = 0
  do i = 2,nproc
     displacements(i) = allsize(i-1) + displacements(i-1)
  enddo

  lastd =  allsize(nproc) + displacements(nproc)

  if(any(alldrange < 0)) then
     print*,myid,'gatherPhiSurf NEGATIVE DRANGE',alldrange,rowid
     stop 'gatherPhiSurf '
  endif

  elcount = 0
  do i = mydrange(1),mydrange(2)
     do k = mydrange(3),mydrange(4)
        myphi2(elcount) = buff(i,k)
        elcount = elcount + 1
     enddo
  enddo
  call MPI_ALLGATHERV(myphi2(0),mysizeB,MPI_DOUBLE_PRECISION,&
       allphi2(0),allsize,displacements,MPI_DOUBLE_PRECISION,&
       comm3d, ierr)

  elcount = 0
  do pcs = 0, nproc-1
     el = 1 + (pcs*4)
     do i = alldrange(el),alldrange(el+1)
        do k = alldrange(el+2),alldrange(el+3)
           buffO(i,k) = allphi2(elcount)
           elcount = elcount + 1 
        enddo
     enddo
  enddo

!-----------------------------------------------
  return
end subroutine gatherbuffsurf
!**********************************************************

!*****************************************************************************
  subroutine gatherbuffMG(current_mesh,buff,fname,fid)

    IMPLICIT NONE

!     Local variables
    real*8,POINTER :: buff(:,:,:)
    TYPE(mesh_pointers) :: current_mesh
    integer i, j, k, n, myxsize,pcs,fid,meqn
    integer elcount,myysize,el
    integer mysize1,mysize2,mysize,elm,elsize
    integer mydrange(4),mn,mgrid
    character*(*) :: fname
    character(LEN=20) :: fn
    real*8 :: dmeanphi
!---------------------------------------------------------------

    if(current_mesh%blank) RETURN
    mn = current_mesh%mesh_num
    mgrid = anint( log10 (dble(anint(current_mesh%dx/finest_mesh%dx))) / log10 (2.0) )

    if(.not. associated(myphi1)) allocate(myphi1(0:(nx+1)))
    if(.not. associated(myphi2)) allocate(myphi2(0:(nx+1)*(nz+1)))
    if(.not. associated(myphi)) allocate(myphi(0:(nx+1)*(ny+1)))

    if(.not. associated(allphi1)) allocate(allphi1(0:(nx+1)))
    if(.not. associated(allphi2)) allocate(allphi2(0:(nx+1)*(nz+1)))
    if(.not. associated(allphi)) allocate(allphi(0:(nx+1)*(ny+1)))

    if(.not. allocated(allsize1)) allocate(allsize1(nproc+1))
    if(.not. allocated(allsize2)) allocate(allsize2(nproc+1))
    if(.not. allocated(allsize)) allocate(allsize(nproc+1))

    if(.not. allocated(displacements1)) allocate(displacements1(nproc+1))
    if(.not. allocated(displacements2)) allocate(displacements2(nproc+1))
    if(.not. allocated(displacements)) allocate(displacements(nproc+1))

    if(.not. allocated(alldrange)) allocate(alldrange(4*nproc))

    mydrange(1) = 0
    mydrange(2) = current_mesh%xe

    mydrange(3) = 0
    mydrange(4) = current_mesh%ye

    myxsize = (mydrange(2) - mydrange(1)) + 1
    myysize = (mydrange(4) - mydrange(3)) + 1
    mysize1 = myxsize                             !1 == 2
    mysize2 = myxsize                             !1 == 2
    mysize  = myxsize * myysize


    call MPI_ALLGATHER(mysize1,1,MPI_INTEGER,allsize1,1,MPI_INTEGER,current_mesh%rowcomm(1),ierr)
    call MPI_ALLGATHER(mysize2,1,MPI_INTEGER,allsize2,1,MPI_INTEGER,current_mesh%comm3d,ierr)
    call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allsize,1,MPI_INTEGER,current_mesh%comm3d,ierr)

    call MPI_ALLGATHER(mydrange,4,MPI_INTEGER,alldrange,4,MPI_INTEGER,current_mesh%comm3d,ierr)

    displacements1(1) = 0
    displacements2(1) = 0
    displacements(1) = 0
    do i = 2,nproc
       displacements1(i) = allsize1(i-1) + displacements1(i-1)
       displacements2(i) = allsize2(i-1) + displacements2(i-1)
       displacements(i) = allsize(i-1) + displacements(i-1)
    enddo
    do pcs = 1, nproc-1
       elm = 1 + 4 * (pcs-1)
       el =  1 + 4 * pcs
       elsize = alldrange(el+1)-alldrange(el+0)
       alldrange(el+0) = alldrange(elm+1) + 1
       alldrange(el+1) = alldrange(el+0) + elsize
    enddo

    do i = mydrange(1),mydrange(2)
       myphi1(i-mydrange(1)) = x_mg(i,mn)
    enddo
    k = 0
    elcount = 0
    do i = mydrange(1),mydrange(2)
       myphi2(elcount) = phi_mg(i,k,mgrid)
       elcount = elcount + 1
    enddo
    elcount = 0
    do i = mydrange(1),mydrange(2)
       do j = mydrange(3),mydrange(4)
          myphi(elcount) = buff(i,k,j)
          elcount = elcount + 1
       enddo
    enddo

    call MPI_ALLGATHERV(myphi1(0),mysize1,MPI_DOUBLE_PRECISION,allphi1(0),allsize1,displacements1,MPI_DOUBLE_PRECISION, &
         current_mesh%rowcomm(1),ierr)
    call MPI_ALLGATHERV(myphi2(0),mysize2,MPI_DOUBLE_PRECISION,allphi2(0),allsize2,displacements2,MPI_DOUBLE_PRECISION,&
         current_mesh%comm3d, ierr)
    call MPI_ALLGATHERV(myphi(0),mysize,MPI_DOUBLE_PRECISION,allphi(0),allsize,displacements,MPI_DOUBLE_PRECISION,&
         current_mesh%comm3d, ierr)
    !!                                                                                                         
    do i = 0,nx_tot(mn)
       G_x(i) = allphi1(i)
    enddo

    elcount = 0
    do pcs = 0, nproc-1
       el = 1 + (pcs*4)
       do i = alldrange(el),alldrange(el+1)
          G_phi(i,k) = allphi2(elcount)
          elcount = elcount + 1 
       enddo
    enddo
    dmeanphi = sum(G_phi(0:nx_tot(mn),0))/dble(nx_tot(mn)+1)

    meqn = 1
    elcount = 0
    do pcs = 0, nproc-1
       el = 1 + (pcs*4)
       do i = alldrange(el),alldrange(el+1)
          do j = alldrange(el+2),alldrange(el+3)
             G_rate(i,j,meqn) = allphi(elcount)
             elcount = elcount + 1 
          enddo
       enddo
    enddo

    if(myid /= 0) return

    write(fn,'(a,"_",i3.3,a)') fname,fid,'.dat'
    close(1234)
    open(1234,file = fn)

    do j = 0,current_mesh%ye
       write(1234,'(1p1230e12.4)')G_x(0:nx_tot(mn)-1)
    end do
    do j = 0,current_mesh%ye
       write(1234,'(1p1230e12.4)')G_phi(0:nx_tot(mn)-1,0) + current_mesh%y(j)-dmeanphi
    end do

    do j = 0,current_mesh%ye
       write(1234,'(1p1230e12.4)')G_rate(0:nx_tot(mn)-1,j,meqn)
    end do
    close(1234)

    if(myid == 0) print*,'Dumped', fn


    return
  end subroutine gatherbuffMG
!
!**************************************************************************************
!
  subroutine consumption_time

!  Local variables
    integer i, k,j,mystart, myend, myxsize
    integer mysize,error
    real*8 :: vel_g,stpp,stppO
    real*8,allocatable :: tcons(:),alltcons(:)
    character*30 ::fname
!
    k = 0d0
    allocate(tcons(drange(1)-1:drange(2)+1),alltcons(0:nx))
    tcons = 0d0
    do i = drange(1),drange(2)
       if(x(i) < XLOC) then
          do j = 0,ny-1
             vel_g = 2d0/(vvel(i,k,j+1)/rhog(i,k,j+1) + vvel(i,k,j)/rhog(i,k,j))
             tcons(i) = tcons(i) + vel_g*dy/detadya(j)
             stpp = f(i,k,j+1,5)*beta_al/(f(i,k,j+1,5)*beta_al+f(i,k,j+1,6))
             if(ncyc>100 .and. myid == 0.and. abs(i) <= 1) &
                  &print'(2i3,1p9e12.4)',i,j,tcons(i),vel_g,dy,detadya(j),stpp
             if(stpp < 1d-2) then
                stppO = f(i,k,j,5)*beta_al/(f(i,k,j,5)*beta_al+f(i,k,j,6))
                if(ncyc>100 .and. myid == 0.and. abs(i) <= 1) &
                     &print'(i3,1pe12.4,a,i3,6e12.4)',i,tcons(i),'Stoppage,J,stpp,vel_g',j,stpp,vel_g,y(j)
                tcons(i) = tcons(i) - merge(abs(stpp - 1d-2)/abs(stpp-stppO)*vel_g*dy/detadya(j),0d0,stpp>0d0)
                exit
             end if
          end do
       end if
    end do
    mystart = drange(1)
    myend = drange(2)
    if(rowid(1) == rowsize(1)-1) myend = drange(2)+1
    myxsize = (myend - mystart) + 1
    mysize = myxsize

    allsize1 = 0
    call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allsize1,1,MPI_INTEGER,&
         rowcomm(1),ierr)

    displacements1(1) = 0
    do i = 2,nproc
       displacements1(i) = allsize1(i-1) + displacements1(i-1)
    enddo
    do i = mystart,myend
       myphi1(i-mystart) = tcons(i)
    enddo

    call MPI_ALLGATHERV(myphi1(0),mysize,MPI_DOUBLE_PRECISION, &
         allphi1(0),allsize1,displacements1,MPI_DOUBLE_PRECISION, &
         rowcomm(1),ierr)
    !!
    do i = 0,nx
       alltcons(i) = allphi1(i)
    enddo
    if(myid /= 0) then
       return
    end if

     if(.not. NORADIATION.and. .not. NOQ4) then
    write(fname,'(a,2("_",1pe7.1),a)')trim(scase) // 'cons',KrierFACTOR,press,'.dat'
     elseif(NORADIATION) then
    write(fname,'(a,2("_",1pe7.1),a)')trim(scase) // 'consQ',KrierFACTOR,press,'.dat'
     elseif(NOQ4) then
    write(fname,'(a,2("_",1pe7.1),a)')trim(scase) // 'consR',KrierFACTOR,press,'.dat'
     end if
    open(2233,file = fname)

    write(2233,'(1p1230e12.4)')G_x(0:nx)   !cancel(uncomment)
    write(2233,'(1p1230e12.4)')alltcons(0:nx)
    close(2233)

    return
  end subroutine consumption_time

end Module Gather_R


!*****************************************************************************
subroutine VIZSAND(flag)
  USE GLOBAL_DATA
  USE MYMPI
  Use Gather_R
  IMPLICIT NONE

  LOGICAL,INTENT(IN) :: flag
  INTEGER :: i,j,k,m,idump,mcuts
  INTEGER, ALLOCATABLE,DIMENSION(:) :: icuts
  real*8,POINTER :: ptr3(:,:,:)
  real*8 :: dmeanphi
  character*30::fname
!------------------------------------------------------


  if(.not. flag) return 
  if(ipack /= 0 ) return

  if(.not. associated(G_rate)) allocate(G_rate(0:nx,0:ny,maxsize))
  allocate(ptr3(lbound(f,1):ubound(f,1),lbound(f,2):ubound(f,2),lbound(f,3):ubound(f,3)))

  call gather_xz
  call gatherPhiSurf
  call gatherbuff(radheat,1)
  do j = 0,ny
     do i = drange(1),drange(2)
        do k = drange(3),drange(4)
           ptr3(i,k,j) = f(i,k,j,1) !1d0-sum(f(i,k,j,2:neqgas))
        end do
     end do
  end do
  call gatherbuff(ptr3,2)

  ptr3(drange(1):drange(2),drange(3):drange(4),drange(5):drange(6)) = &
       &Vfrac(drange(1):drange(2),drange(3):drange(4),drange(5):drange(6),1)
  call gatherbuff(ptr3,3)
  ptr3(drange(1):drange(2),drange(3):drange(4),drange(5):drange(6)) = &
       &Vfrac(drange(1):drange(2),drange(3):drange(4),drange(5):drange(6),2)
  call gatherbuff(ptr3,4)

  do j = drange(5),drange(6)
     do k = drange(3),drange(4)
        do i = drange(1),drange(2)
           if(sum(f(i,k,j,5:6)) > 1d-6) then
              ptr3(i,k,j) = f(i,k,j,5)  !Diamp(i,k,j,1)
           else
              ptr3(i,k,j) =0d0
           end if
        enddo
     enddo
  enddo
  call gatherbuff(ptr3,5)
  do j = drange(5),drange(6)
     do k = drange(3),drange(4)
        do i = drange(1),drange(2)
           if(sum(f(i,k,j,5:6)) > 1d-6) then
              ptr3(i,k,j) = Diamp(i,k,j,1)
           else
              ptr3(i,k,j) =0d0
           end if
        enddo
     enddo
  enddo
  call gatherbuff(ptr3,6)
  
  nullify(ptr3)

  call consumption_time

  if(myid /= 0) then
     return
  end if

  dmeanphi = sum(G_phi(0:nx,0))/dble(nx+1)
  
  if(pressvary) then
     write(fname,'(a,2("_",1pe7.1),a)')trim(scase) // 'surf',KrierFACTOR,press,'.dat'
  else
     if(.not. NORADIATION.and. .not. NOQ4) then
        write(fname,'(a,"_",1pe7.1,a)')trim(scase) // 'surf',KrierFACTOR,'.dat'
     elseif(NORADIATION) then
        write(fname,'(a,"_",1pe7.1,a)')trim(scase) // 'surfQ',KrierFACTOR,'.dat'
     elseif(NOQ4) then
        write(fname,'(a,"_",1pe7.1,a)')trim(scase) // 'surfR',KrierFACTOR,'.dat'
     end if
  end if

  open(1233,file = trim(fname))

  write(1233,'(a,1p2e12.4,a,e12.4)')'% regression rate',surfvel_max, surfvel_min,'Pressure',Press
  
  write(1233,'(1p1230e12.4)')G_x(0:nx)
  write(1233,'(1p1230e12.4)')G_phi(0:nx,0)-dmeanphi
  close(1233)
  
  G_rate(nx,:,:) = G_rate(nx-1,:,:)
  if(pressvary) then
     write(fname,'(a,2("_",1pe7.1),a)')trim(scase) // 'rad',KrierFACTOR,press,'.dat'
  else
     if(.not. NORADIATION .and. .not. NOQ4) then
        write(fname,'(a,"_",1pe7.1,a)')trim(scase) // 'rad',KrierFACTOR,'.dat'
     elseif(NORADIATION) then
        write(fname,'(a,"_",1pe7.1,a)')trim(scase) // 'radQ',KrierFACTOR,'.dat'
     elseif(NOQ4) then
        write(fname,'(a,"_",1pe7.1,a)')trim(scase) // 'radR',KrierFACTOR,'.dat'
     end if
  end if


  open(1233,file = trim(fname))

  write(1233,'(a,1p2e12.4,a,e12.4)')'% regression rate',surfvel_max, surfvel_min,'Pressure',Press

  do j = 0,ny
     write(1233,'(1p1230e12.4)')G_x(0:nx)
  end do
  do j = 0,ny
     write(1233,'(1p1230e12.4)')G_phi(0:nx,0)+y(j)-dmeanphi
  end do

  do m = 1,6
     write(1233,*)'%component',m
     do j = 0,ny
        write(1233,'(1p1230e12.4)')G_rate(0:nx,j,m)
!!>        write(1233,*)G_rate(0:nx,j,m)
     end do
  enddo


  close(1233)
  
  return

!------------------------------------------------------
  RETURN
END subroutine VIZSAND
!*****************************************************************************8


!*****************************************************************************
subroutine VIZPACK3D(flag)
  USE GLOBAL_DATA
  USE MYMPI
  Use Gather_R
  Use Statistics
  IMPLICIT NONE

  LOGICAL,INTENT(IN) :: flag
  INTEGER :: i,j,k,m,idump,mcuts,error,ivar,nvar
  real*8,POINTER :: ptr3(:,:,:)
  real*8 :: dmeanphi
  character*30::fname
!------------------------------------------------------

  if(.not. flag) return 
  if(ipack /= 1 ) return

  if(myid == 0) Print*,'Beginning vizpack3d'

  if(.not. associated(G_f)) ALLOCATE(G_f(0:nx,0:nz,0:ny),STAT=error)
  allocate(ptr3(lbound(f,1):ubound(f,1),lbound(f,2):ubound(f,2),lbound(f,3):ubound(f,3)))

  call gather_xz
  call gatherbuffsurf(phi,G_phi)
  call gatherbuffsurf(BCradheat,G_psi)

  dmeanphi = sum(G_phi(0:nx,0:nz))/dble(nx+1)/dble(nz+1)
  G_psi(nx,:) = G_psi(0,:)
  G_psi(:,nz) = G_psi(:,0)


  write(fname,'(a,2("_",1pe7.1),a)')trim(scase) // 'surf',KrierFACTOR,press,'.dat'

  if(myid == 0) then
     Print*,'Opening Surf file'

     open(1233,file = trim(fname))

     write(1233,'(a,"[",1p2e12.4,"]",10(",",1x,a,e12.4))')'% regression_rate=',stat_rb(1:2),'Pressure=',Press

     write(1233,'(1p1230e12.4)')G_x(0:nx)
     write(1233,'(1p9991230e14.6)')G_phi(0:nx,0:nz)-dmeanphi
     write(1233,'(1p9991230e14.6)')G_psi(0:nx,0:nz)

     print*,'DUmping BCradheat',maxval(G_psi(0:nx,0:nz)),minval(G_psi(0:nx,0:nz))
     close(1233)
  endif

  if(myid == 0) Print*,'gatherbuff3 in vizpack3d'

  Nvar = 2
  varloop: Do ivar = 1,Nvar
     do j = 0,ny
        do i = drange(1)-1,drange(2)+1
           do k = drange(3)-1,drange(4)+1
              if(ivar == 1) then
                 ptr3(i,k,j) = f(i,k,j,1) !1d0-sum(f(i,k,j,2:neqgas))
              elseif(ivar == 2) then
                 ptr3(i,k,j) = radheat(i,k,j)
              end if
           end do
        end do
     end do
     call gatherbuff3(ptr3)
  
     if(ivar == nvar) nullify(ptr3)

     if(myid /= 0) then
        cycle varloop
     end if

     G_f(nx,:,:) = G_f(0,:,:)
     G_f(:,nz,:) = G_f(:,0,:)
     if(ivar == 1) then
        G_f(:,:,ny) = G_f(:,:,ny-1)
     elseif(ivar == 2) then
        G_f(:,:,0) = G_f(:,:,1)
        G_f(:,:,ny) = G_f(:,:,ny-1)
     end if
     if(ivar == 1) then
        write(fname,'(a,2("_",1pe7.1),a)')trim(scase) // 'rad3D',KrierFACTOR,press,'.dat'

        if(myid == 0) Print*,'Opening vizpack3d file: ',trim(fname)

        open(1233,file = trim(fname))

        write(1233,'(a,"[",1p2e12.4,"]",10(",",1x,a,"=",e12.4))')&
          '% regression_rate=',stat_rb(1:2),'Pressure',&
          Press,'nf',dble(nvar),'dmnphi',dmeanphi

        write(1233,'(1p1230e12.4)')G_x(0:nx)
        do j = 0,ny
           write(1233,'(1p9991230e14.6)')G_phi(0:nx,0:nz)+y(j)-dmeanphi
        end do
     endif

     if(myid == 0) Print*,'Writing file',&
          maxval(G_f(0:nx,0:nz,0:ny)),minval(G_f(0:nx,0:nz,0:ny))
     write(1233,'(1p999230e12.4)')G_f(0:nx,0:nz,0:ny)

  end Do varloop

  if(myid == 0) Print*,'END output  vizpack3d'

  close(1233)
  
  if(myid == 0) Print*,'END vizpack3d'
  return

!------------------------------------------------------
  RETURN
END subroutine VIZPACK3D
!*****************************************************************************8


!********************************************************************
subroutine gather_speed
  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

  integer i, mystart, myend, myxsize
  integer elcount
  integer mysize
!--------------------------------------------------------

  mystart = drange(1)
  myend = drange(2)
  if(x(myend) >= - period/two - two*dx) myend = drange(2)+1
  myxsize = (myend - mystart) + 1
  mysize = myxsize

  if(mysize.ne.(nx+1)) then

     if(mystart.eq.0) then
        mystart = 0
     endif
     if(myend.eq.nx) then 
        myend = nx
     endif
     myxsize = (myend - mystart) + 1
     mysize = myxsize         

     call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allsize1,1,MPI_INTEGER,&
          cutcomm,ierr)

     displacements1(1) = 0
     do i = 2,nproc_cut
        displacements1(i) = allsize1(i-1) + displacements1(i-1)
     enddo

     elcount = 0
     do i = mystart,myend
        myphi1(elcount)  = phit(i,kplot)
        elcount = elcount + 1
     enddo
     call MPI_ALLGATHERV(myphi1(0),mysize,MPI_DOUBLE_PRECISION,&
          allphi1(0),allsize1,displacements1,MPI_DOUBLE_PRECISION,&
          cutcomm,ierr)
     elcount = 0
     do i = 0,nx+1
        G_speed(i) = allphi1(elcount)
        elcount = elcount + 1
     enddo
  endif
  return
end subroutine gather_speed
!*****************************************************************************



!*****************************************************************************
subroutine shaneWriteFields( fileSuffix )
	! Prints field data to file.  This is a modified version of VIZPACK3D,
	! my version writes the data as x,y,z arrays rather than x,z,y.

	use GLOBAL_DATA
	use MYMPI
	use Gather_R
	implicit none
	
	integer, intent(IN) :: fileSuffix
	
	INTEGER :: i,j,k,m,error,ivar,nvar
	real*8,POINTER :: ptr3(:,:,:)
	character*100::fname


	if(.not. associated(G_f)) ALLOCATE(G_f(0:nx,0:nz,0:ny),STAT=error)
	allocate(ptr3(lbound(f,1):ubound(f,1),lbound(f,2):ubound(f,2),lbound(f,3):ubound(f,3)))

	call gather_xz							! gathers local x and z into global G_x and G_z
	call gatherbuffsurf(phi,G_phi)			! gathers local phi into G_phi

	! write out the x,y,z grids and phi
    4000 format( 1pe14.6 )
	4001 format( 1p999999e14.6 )
	
	if( myid == 0 ) then
		open(  3000, file='gridx.grid' )
		write( 3000, fmt=4000 ) G_x(0:nx)
		close( 3000 )

		open(  3000, file='gridz.grid' )
		write( 3000, fmt=4000 ) G_z(0:nz)
		close( 3000 )

		open(  3000, file='gridy.grid' )
		write( 3000, fmt=4000 ) y(0:ny)
		close( 3000 )

		write( fname, '("phii.field",BZ,I6.6)' ) fileSuffix
		open(  3000, file=fname, recl = 2147483646 )
		do i = 0, nx
				write( 3000, fmt=4001 ) G_phi(i,0:nz) ! written as column-major in x
		end do
		close( 3000 )
	end if
	
	! write out the scalar field variables
	!     f(i,k,j,1)  = T    (gas phase temperature -- using tgas)
	!     f(i,k,j,2)  = Y_O  (gas phase oxidizer X -- using Agas)
	!     f(i,k,j,3)  = Y_Z  (gas phase intermediate Z -- using Bgas)
	!     f(i,k,j,4)  = Y_F  (gas phase fuel Y -- using Cgas)
	!     f(i,k,j,5)  = T    (solid phase temperature -- using tsol)
	!     psi(i,k,j) = level set in the solid
	nvar = 6
	do ivar = 1, nvar
		
		! copy the 3D field into G_f(x,z,y)
		if( ivar < 6 ) then
			ptr3( drange(1)-1:drange(2)+1, drange(3)-1:drange(4)+1, 0:ny ) = &
			   f( drange(1)-1:drange(2)+1, drange(3)-1:drange(4)+1, 0:ny, ivar )
		elseif( ivar == 6 ) then
			ptr3( drange(1)-1:drange(2)+1, drange(3)-1:drange(4)+1, 0:ny ) = &
			 psi( drange(1)-1:drange(2)+1, drange(3)-1:drange(4)+1, 0:ny )
		end if
		call gatherbuff3(ptr3)

		! write the field to a file
		if(myid == 0) then

			! copy the periodic boundaries (I guess???)
			G_f(nx,:,:) = G_f(0,:,:)
			G_f(:,nz,:) = G_f(:,0,:)
			G_f(:,:,0) = G_f(:,:,1)
			G_f(:,:,ny) = G_f(:,:,ny-1)

			! name the file
			select case (ivar)
				case (1)
					write( fname, '("Tgas.field",BZ,I6.6)' ) fileSuffix
				case (2)
					write( fname, '("Agas.field",BZ,I6.6)' ) fileSuffix
				case (3)
					write( fname, '("Bgas.field",BZ,I6.6)' ) fileSuffix
				case (4)
					write( fname, '("Cgas.field",BZ,I6.6)' ) fileSuffix
				case (5)
					write( fname, '("Tsol.field",BZ,I6.6)' ) fileSuffix
				case (6)
					write( fname, '("psii.field",BZ,I6.6)' ) fileSuffix
			end select
			open( 3000, file = fname, recl = 2147483646 )
			
			! write the field to the file in blocks
			do i = 0, nx
				do j = 0, ny
					write( 3000, 4001 ) G_f(i,0:nz,j)
				end do
			end do
			close( 3000 )

		end if

	end do

	! clean up
	nullify(ptr3)

end subroutine shaneWriteFields
!*****************************************************************************


!*****************************************************************************
subroutine shaneWriteSlice( fileSuffix )
	! Prints the middle slice in the x-y plane to file.  This is of course far
	! less data than printing the whole field, so it should be easier to make
	! movies and such.  The way I do it here is not the most efficient, but I
	! did not feel like writing another "gatherbuff2D" routine to pull the x-y
	! plane only.

	use GLOBAL_DATA
	use MYMPI
	use Gather_R
	implicit none
	
	integer, intent(IN) :: fileSuffix
	integer, parameter :: nvar = 6
	INTEGER :: i,j,k,m,idump,mcuts,error,ivar,zmid
	integer :: varmap(nvar)
	real*8,POINTER :: ptr3(:,:,:)
	real*8 :: dmeanphi
	character*100::fname
	character*10 :: charmap(nvar)

	if(.not. associated(G_f)) ALLOCATE(G_f(0:nx,0:nz,0:ny),STAT=error)
	allocate(ptr3(lbound(f,1):ubound(f,1),lbound(f,2):ubound(f,2),lbound(f,3):ubound(f,3)))

	call gather_xz							! gathers local x and z into global G_x and G_z
	call gatherbuffsurf(phi,G_phi)			! gathers local phi into G_phi

	! find the middle of the z grid
	zmid = nz/2

	! write out the x,y,z grids and phi
    4000 format( 1pe14.6 )
	4001 format( 1p999999e14.6 )
	
	if( myid == 0 ) then
	
		! write x meshgrid
		write( fname, '("x-all",BZ,I6.6,".mesh")' ) fileSuffix
		open(  3000, file=fname, recl = 2147483646 )
		do j = 0, ny
			write( 3000, fmt=4001 ) G_x(0:nx)
		end do
		close( 3000 )

		! write y meshgrid for gas phase
		write( fname, '("y-gas",BZ,I6.6,".mesh")' ) fileSuffix
		open(  3000, file=fname, recl = 2147483646 )
		do j = 0, ny
			write( 3000, fmt=4001 ) ( y(j)+G_phi(i,zmid), i=0,nx )
		end do
		close( 3000 )

		! write y meshgrid for solid phase
		write( fname, '("y-sol",BZ,I6.6,".mesh")' ) fileSuffix
		open(  3000, file=fname, recl = 2147483646 )
		do j = ny, 0, -1
			write( 3000, fmt=4001 ) ( -y(j)+G_phi(i,zmid), i=0,nx )
		end do
		close( 3000 )

        end if
        
        ! write out the scalar field variables on the x-y plane
        !     f(i,k,j,1)  = T    (gas phase temperature -- using tgas)
        !     f(i,k,j,2)  = Y_O  (gas phase oxidizer X -- using Agas)
        !     f(i,k,j,3)  = Y_Z  (gas phase intermediate Z -- using Bgas)
        !     f(i,k,j,4)  = Y_F  (gas phase fuel Y -- using Cgas)
        !     f(i,k,j,neqmax)  = T    (solid phase temperature -- using tsol)
        !     psi(i,k,j) = level set in the solid
        varmap(1) = 1       ! Tgas
        varmap(2) = 2       ! Y_O
        varmap(3) = 3       ! Y_Z
        varmap(4) = 4       ! Y_F
        varmap(5) = neqmax  ! Tsol
        varmap(6) = -1      ! error
        charmap(1) = 'Tgas'
        charmap(2) = 'Agas'
        charmap(3) = 'Bgas'
        charmap(4) = 'Cgas'
        charmap(5) = 'Tsol'
        charmap(6) = 'psi'


        do ivar = 1, nvar

        ! copy the 3D field into G_f(x,z,y)
        if( ivar < 6 ) then
        ptr3( drange(1)-1:drange(2)+1, drange(3)-1:drange(4)+1, 0:ny ) = &
           f( drange(1)-1:drange(2)+1, drange(3)-1:drange(4)+1, 0:ny, varmap(ivar) )
        elseif( ivar == 6 ) then
        ptr3( drange(1)-1:drange(2)+1, drange(3)-1:drange(4)+1, 0:ny ) = &
         psi( drange(1)-1:drange(2)+1, drange(3)-1:drange(4)+1, 0:ny )
        end if
        call gatherbuff3(ptr3)

! write the slice to a file
if(myid == 0) then

			! copy the periodic boundaries (I guess???)
			G_f(nx,:,:) = G_f(0,:,:)
			G_f(:,nz,:) = G_f(:,0,:)
			G_f(:,:,0) = G_f(:,:,1)
			G_f(:,:,ny) = G_f(:,:,ny-1)

			! name the file
			write( fname, '(A,BZ,I6.6,".slice")' ) trim( charmap(ivar) ), fileSuffix
			open( 3000, file = fname, recl = 2147483646 )
			
			! write the slice to the file with arrangement corresponding to the
			! meshgrids above (to plot in Matlab, y is the row number, x is column)
			do j = 0, ny
				write( 3000, 4001 ) G_f(0:nx,zmid,j)
			end do
			close( 3000 )

		end if

	end do

	! clean up
	nullify(ptr3)

end subroutine shaneWriteSlice
!*****************************************************************************
