!     ****************************************************************
!     *                                                              *
!     *                     module data_structure                    *
!     *                                                              *
!     *     this module contains subroutines to manipulate the       *
!     *     data structure.                                          *
!     *                                                              * 
!     *                                                              *
!     ****************************************************************
!

MODULE data_structure 

  USE data_types
  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

  PUBLIC :: storage_allocate, rednblack_lists

CONTAINS

!***************************************************************************************************

  SUBROUTINE storage_allocate

    ! Subroutine to create the dynamic storage used by the program.

    !---------------------------------------------------------------------------------------------------

    ! Local variables:

    TYPE(mesh_pointers), POINTER :: current_mesh  ! pointer to the data for the current mesh
    TYPE(mesh_pointers), POINTER :: fine_mesh     ! pointer to the data for the fine mesh

    INTEGER :: error, imesh, r , nxx,nyy,nzz,mn,igm,i,k

    !---------------------------------------------------------------------------------------------------

    ! Allocate storage for various global arrays.

    r = zero

    ! Create the linked list of MESH_POINTERS that point to the storage required for each mesh.

    ALLOCATE( finest_mesh, STAT=error )

    IF ( error /= 0 ) THEN
       WRITE(*,*) "Allocation error for FINEST_MESH"
       STOP
    END IF

    NULLIFY( finest_mesh%fine_mesh )

    current_mesh => finest_mesh

    IF ( num_meshes == 1 ) THEN

       finest_mesh%mesh_num = 1

       finest_mesh%xe = nx_mg(1)-1
       finest_mesh%ze = nz_mg(1)-1  
       finest_mesh%ye = ny_mg(1)-1   !because we shift to 0-(ny-1) rather than (1-ny)

       finest_mesh%dx = dx_mg(1)
       finest_mesh%dz = dz_mg(1)
       finest_mesh%dy = dy_mg(1)

       NULLIFY( finest_mesh%coarse_mesh )
       NULLIFY( coarsest_mesh )

    ELSE IF ( num_meshes > 1 ) THEN

       DO imesh = 1,num_meshes


          current_mesh%mesh_num = imesh

          current_mesh%xe = nx_mg(imesh)-1   !THESE are LOOP upper LIMITS
          current_mesh%ze = nz_mg(imesh)-1
          current_mesh%ye = ny_mg(imesh)-1

          current_mesh%dx = dx_mg(imesh)
          current_mesh%dz = dz_mg(imesh)
          current_mesh%dy = dy_mg(imesh)


          IF ( imesh == 1 ) THEN  ! finest mesh

             ALLOCATE( current_mesh%coarse_mesh, STAT=error )

             IF ( error /= 0 ) THEN
                WRITE(*,*) "Allocation error for COARSE_MESH"
                STOP
             END IF

             fine_mesh    => current_mesh
             current_mesh => current_mesh%coarse_mesh

          ELSE IF ( imesh == num_meshes ) THEN  ! coarsest mesh

             current_mesh%fine_mesh => fine_mesh

             coarsest_mesh => current_mesh

             NULLIFY ( coarsest_mesh%coarse_mesh )

          ELSE  ! intermediate mesh

             ALLOCATE( current_mesh%coarse_mesh, STAT=error )

             IF ( error /= 0 ) THEN
                WRITE(*,*) "Allocation error for COARSE_MESH"
                STOP
             END IF

             current_mesh%fine_mesh => fine_mesh

             fine_mesh    => current_mesh
             current_mesh => current_mesh%coarse_mesh

          END IF

       END DO

    END IF


    ! Assign storage for the mesh vectors required in the solution algorithm on
    ! each mesh.

    current_mesh => finest_mesh
    nxx = current_mesh%xe;nzz = current_mesh%ze;nyy = current_mesh%ye+1
    ALLOCATE( wrk(-ngh:nxx+ngh,-ngh:nzz+ngh,-ngh:nyy+ngh), STAT=error )

    DO

       IF ( .NOT. ASSOCIATED(current_mesh) ) EXIT

       nxx = current_mesh%xe
       nzz = current_mesh%ze
       nyy = current_mesh%ye+1
       mn = current_mesh%mesh_num
!
!  allocate at the coefficient matrix pointers for every mesh
!
       ALLOCATE( current_mesh%cd(0:nxx,0:nzz,4), STAT=error )
       ALLOCATE( current_mesh%cj(0:nxx,0:nzz,9), STAT=error )
       ALLOCATE( current_mesh%c3(0:nxx,0:nzz,nyy), STAT=error )
       ALLOCATE( current_mesh%c7(0:nxx,0:nzz,nyy), STAT=error )
       ALLOCATE( current_mesh%c2f(0:3,-1:nxx+1,-1:nzz+1,-1:nyy+1), STAT=error )
       ALLOCATE( current_mesh%ey(-1:nyy), STAT=error )
       ALLOCATE( current_mesh%y(-1:nyy), STAT=error )
       ALLOCATE( current_mesh%ya(-1:nyy), STAT=error )
       ALLOCATE( current_mesh%eta(-1:nyy), STAT=error )
       ALLOCATE( current_mesh%eyy(-1:nyy), STAT=error )
       ALLOCATE( current_mesh%eya(-1:nyy), STAT=error )
       ALLOCATE( current_mesh%CPPE(0:nxx+1,0:nzz+1))
       do i =  lbound(current_mesh%CPPE,1), ubound(current_mesh%CPPE,1)
          do k =  lbound(current_mesh%CPPE,2), ubound(current_mesh%CPPE,2)
             ALLOCATE( current_mesh%CPPE(i,k)%ccmg(15,0:nyy+1))
          end do
       end do
       ALLOCATE( current_mesh%chi(-1:nyy,2), STAT=error )
       ALLOCATE( current_mesh%chipr(-1:nyy,2), STAT=error )
       CALL CHECK_ERROR(error,'MULTIGRID')

       ALLOCATE( current_mesh%xold(-ngh:nxx+ngh,-ngh:nzz+ngh,-ngh:nyy+ngh), STAT=error )
       ALLOCATE( current_mesh%ff(-ngh:nxx+ngh,-ngh:nzz+ngh,-ngh:nyy+ngh), STAT=error )
       ALLOCATE( current_mesh%gm(gmres_size_mg(mn)))
       do igm = lbound(current_mesh%gm,1), ubound(current_mesh%gm,1)
          ALLOCATE( current_mesh%gm(igm)%mv(-ngh:nxx+ngh,-ngh:nzz+ngh,-ngh:nyy+ngh), STAT=error )
       end do
       !NOT on FINEST*************************
       IF ( ASSOCIATED(current_mesh%fine_mesh) ) THEN  

          ALLOCATE( current_mesh%f(-ngh:nxx+ngh,-ngh:nzz+ngh,-ngh:nyy+ngh), STAT=error )
          ALLOCATE( current_mesh%x(-ngh:nxx+ngh,-ngh:nzz+ngh,-ngh:nyy+ngh), STAT=error )

          !
          !  Note on the coarse levels w has size 4 for the last dimension
          !
          !         ALLOCATE( current_mesh%w(-ngh:nxx+ngh,-ngh:nzz+ngh,-ngh:nyy+ngh,4), STAT=error )
          !         ALLOCATE( current_mesh%r(-ngh:nxx+ngh,-ngh:nzz+ngh,-ngh:nyy+ngh,4), STAT=error )
          !         ALLOCATE(current_mesh%cm(9,-ngh:nxx+ngh,-ngh:nzz+ngh,-ngh:nyy+ngh), STAT=error )
          !         ALLOCATE(current_mesh%dr(4,4, -ngh:nxx+ngh,-ngh:nzz+ngh,-ngh:nyy+ngh), STAT=error )

          current_mesh%f = zero
          current_mesh%x = zero

       END IF

       current_mesh => current_mesh%coarse_mesh

    END DO

    !---------------------------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE storage_allocate
!***************************************************************************************************

!***************************************************************************************************
  SUBROUTINE rednblack_lists

  ! Subroutine to create mapping to red n black solution scheme

!---------------------------------------------------------------------------------------------------
  ! Local variables:

    TYPE(mesh_pointers), POINTER :: current_mesh  ! pointer to the data for the current mesh
    TYPE(mesh_pointers), POINTER :: fine_mesh     ! pointer to the data for the fine mesh
    INTEGER :: xe,ze,i,k,nm,ii,kk,ired,iblck,isb(8),isr(8)
    INTEGER :: maxcolorpoints,maxnsnd,maxnrcv,error
    INTEGER :: dummy1,dummy2,color,tmpdim(3),dtmp(4),npts(4),tmpcomm
    INTEGER :: status(MPI_STATUS_SIZE),tmpcoords(3)
    LOGICAL :: iodd,kodd,iodd0,kodd0,passz3,passz4,recvz3,recvz4
    LOGICAL :: isperiodic(3)
!DELETE
    INTEGER :: icolor,indx
!---------------------------------------------------------------------------------------------------

  current_mesh => finest_mesh

  DO

    IF ( .NOT. ASSOCIATED(current_mesh) ) EXIT

    xe = current_mesh%xe;ze = current_mesh%ze;
    nm=current_mesh%mesh_num 

!
!  Create list of black and red nodes for the Gauss Seidel Scheme
!
    current_mesh%nclr(:)=0
    current_mesh%nsnd= 0;current_mesh%nrcv=0
    do i=0,xe
    do k=0,ze
       ii=i+ib_mg(nm)
       kk=k+kb_mg(nm)
       iodd = 2*int(ii/2) < ii
       kodd = 2*int(kk/2) < kk
       if(iodd.eqv.kodd) then
         current_mesh%nclr(1) = current_mesh%nclr(1)+1
       else
         current_mesh%nclr(2) = current_mesh%nclr(2)+1
       endif
    enddo
    enddo
      
! 
!  this part of the program fills up communication lists used in the red and black
!  gauss Siedel scheme. Although it is possible to implment the scheme without 
!  mapping lists, this makes life easier when the grid dimensions differ among
!  processors
! NOTATION in nsnd(*,1:2), * is:
! (1,3) is (left,up), hence send towards prcessor (left,up), and recive from
! processor (left,up), (2,4) is (right,down)
!

    do i=-1,xe+1
    do k=-1,ze+1
       ii = i+ib_mg(nm)
       kk = k+kb_mg(nm)

!
! these 4 flags are indeed always true, but in the MG solver we can avoid swapping corners
! for the finest grid saving some communication. But in the implicit solver we need to pass the 
! corners in the velocity equation. In future development, it may be the case to pass the
! corners only when actually needed
!
       passz3 = (i>=0.AND.i<=xe) .OR. nm >= 1
       passz4 = (i>=0.AND.i<=xe) .OR. nm >= 1
       recvz3 = (i>=0.AND.i<=xe) .OR. nm >= 1 
       recvz4 = (i>=0.AND.i<=xe) .OR. nm >= 1

       iodd = 2*int(abs(ii)/2) < abs(ii)
       kodd = 2*int(abs(kk)/2) < abs(kk)

       if ( ii == nx_tot(nm) ) iodd = 2*int(abs(0)/2) < abs(0)
       if ( kk == nx_tot(nm) ) kodd = 2*int(abs(0)/2) < abs(0)
       if ( ii == -1 ) iodd = 2*int(abs(nx_tot(nm)-1)/2) < abs(nx_tot(nm)-1)
       if ( kk == -1 ) kodd = 2*int(abs(nx_tot(nm)-1)/2) < abs(nx_tot(nm)-1)

       if(iodd.eqv.kodd) then
         if(i == 0.AND.(k>=0.AND.k<=ze) )  current_mesh%nsnd(1,1)=current_mesh%nsnd(1,1)+1
         if(i == xe.AND.(k>=0.AND.k<=ze))  current_mesh%nsnd(2,1)=current_mesh%nsnd(2,1)+1
         if(k == 0.AND. passz3 )           current_mesh%nsnd(3,1)=current_mesh%nsnd(3,1)+1
         if(k == ze.AND. passz4 )          current_mesh%nsnd(4,1)=current_mesh%nsnd(4,1)+1
!
         if(i ==-1  .AND.(k>=0.AND.k<=ze)) current_mesh%nrcv(1,1)=current_mesh%nrcv(1,1)+1
         if(i ==xe+1.AND.(k>=0.AND.k<=ze)) current_mesh%nrcv(2,1)=current_mesh%nrcv(2,1)+1
         if(k ==-1  .AND.recvz3 )          current_mesh%nrcv(3,1)=current_mesh%nrcv(3,1)+1
         if(k ==ze+1.AND.recvz4 )          current_mesh%nrcv(4,1)=current_mesh%nrcv(4,1)+1
       else
         if(i == 0.AND.(k>=0.AND.k<=ze) )  current_mesh%nsnd(1,2)=current_mesh%nsnd(1,2)+1
         if(i == xe.AND.(k>=0.AND.k<=ze))  current_mesh%nsnd(2,2)=current_mesh%nsnd(2,2)+1
         if(k == 0.AND. passz3 )           current_mesh%nsnd(3,2)=current_mesh%nsnd(3,2)+1
         if(k == ze.AND. passz4 )          current_mesh%nsnd(4,2)=current_mesh%nsnd(4,2)+1
!
         if(i ==-1  .AND.(k>=0.AND.k<=ze)) current_mesh%nrcv(1,2)=current_mesh%nrcv(1,2)+1
         if(i ==xe+1.AND.(k>=0.AND.k<=ze)) current_mesh%nrcv(2,2)=current_mesh%nrcv(2,2)+1
         if(k ==  -1   .AND. recvz3)       current_mesh%nrcv(3,2)=current_mesh%nrcv(3,2)+1
         if(k == ze+1  .AND. recvz4)       current_mesh%nrcv(4,2)=current_mesh%nrcv(4,2)+1
       endif
    enddo
    enddo


    maxcolorpoints=max(current_mesh%nclr(1),current_mesh%nclr(2))
    maxnsnd=maxval(current_mesh%nsnd); maxnrcv=maxval(current_mesh%nrcv);
        
    ALLOCATE( current_mesh%iclr(maxcolorpoints,2), STAT=error )
    ALLOCATE( current_mesh%kclr(maxcolorpoints,2), STAT=error )

    ALLOCATE( current_mesh%sndclr(maxnsnd,4,2), STAT=error )
    ALLOCATE( current_mesh%rcvclr(maxnrcv,4,2), STAT=error )

    current_mesh%kclr=0;current_mesh%iclr=0
    current_mesh%sndclr=0;current_mesh%rcvclr=0

    ired=0;iblck=0
    do i=0,xe
    do k=0,ze
       ii=i+ib_mg(nm)
       kk=k+kb_mg(nm)
       iodd = 2*int(ii/2) < ii
       kodd = 2*int(kk/2) < kk
       if(iodd.eqv.kodd) then
         iblck = iblck+1
         current_mesh%iclr(iblck,1) = i
         current_mesh%kclr(iblck,1) = k
       else
         ired = ired+1
         current_mesh%iclr(ired,2) = i
         current_mesh%kclr(ired,2) = k
       endif
    enddo
    enddo 

!
! Now put the points in the send-receive list
!
    isb=0;isr=0
    do i=-1,xe+1
    do k=-1,ze+1
       ii=i+ib_mg(nm)
       kk=k+kb_mg(nm)

       passz3 = (i>=0.AND.i<=xe) .OR. nm >= 1
       passz4 = (i>=0.AND.i<=xe) .OR. nm >= 1
       recvz3 = (i>=0.AND.i<=xe) .OR. nm >= 1
       recvz4 = (i>=0.AND.i<=xe) .OR. nm >= 1
!
       iodd = 2*int(abs(ii)/2) < abs(ii)
       kodd = 2*int(abs(kk)/2) < abs(kk)

       if ( ii == nx_tot(nm) ) iodd = 2*int(abs(0)/2) < abs(0)
       if ( kk == nx_tot(nm) ) kodd = 2*int(abs(0)/2) < abs(0)
       if ( ii == -1 ) iodd = 2*int(abs(nx_tot(nm)-1)/2) < abs(nx_tot(nm)-1)
       if ( kk == -1 ) kodd = 2*int(abs(nx_tot(nm)-1)/2) < abs(nx_tot(nm)-1)
!!
       IF(iodd.eqv.kodd) then
         if(i == 0.AND.(k>=0.AND.k<=ze))   then
           isb(1)=isb(1)+1; current_mesh%sndclr(isb(1),1,1)=k;
         endif 
         if(i == xe.AND.(k>=0.AND.k<=ze) ) then
           isb(2)=isb(2)+1; current_mesh%sndclr(isb(2),2,1)=k;
         endif
         if(k == 0.AND. passz3 )           then
           isb(3)=isb(3)+1; current_mesh%sndclr(isb(3),3,1)=i;
         endif
         if(k == ze.AND. passz4 )          then
           isb(4)=isb(4)+1; current_mesh%sndclr(isb(4),4,1)=i;
         endif
!
         if(i ==-1.AND.(k>=0.AND.k<=ze))   then
           isb(5)=isb(5)+1; current_mesh%rcvclr(isb(5),1,1)=k;
         endif
         if(i ==xe+1.AND.(k>=0.AND.k<=ze)) then
           isb(6)=isb(6)+1; current_mesh%rcvclr(isb(6),2,1)=k;
         endif
         if(k ==-1.AND. recvz3 )           then
           isb(7)=isb(7)+1; current_mesh%rcvclr(isb(7),3,1)=i;
         endif
         if(k ==ze+1.AND. recvz4 )         then
           isb(8)=isb(8)+1; current_mesh%rcvclr(isb(8),4,1)=i;
         endif
       ELSE
         if(i == 0.AND.(k>=0.AND.k<=ze))   then
           isr(1)=isr(1)+1; current_mesh%sndclr(isr(1),1,2)=k;
         endif 
         if(i == xe.AND.(k>=0.AND.k<=ze) ) then
           isr(2)=isr(2)+1; current_mesh%sndclr(isr(2),2,2)=k;
         endif
         if(k == 0.AND. passz3 )           then
           isr(3)=isr(3)+1; current_mesh%sndclr(isr(3),3,2)=i;
         endif
         if(k == ze.AND. passz4 )          then
           isr(4)=isr(4)+1; current_mesh%sndclr(isr(4),4,2)=i;
         endif
!
         if(i ==-1.AND.(k>=0.AND.k<=ze))   then
           isr(5)=isr(5)+1; current_mesh%rcvclr(isr(5),1,2)=k;
         endif
         if(i ==xe+1.AND.(k>=0.AND.k<=ze)) then
           isr(6)=isr(6)+1; current_mesh%rcvclr(isr(6),2,2)=k;
         endif
         if(k ==-1.AND. recvz3 )           then
           isr(7)=isr(7)+1; current_mesh%rcvclr(isr(7),3,2)=i;
         endif
         if(k ==ze+1.AND.recvz4)           then
           isr(8)=isr(8)+1; current_mesh%rcvclr(isr(8),4,2)=i;
         endif
       ENDIF
!!
    enddo
    enddo

!
!  This last part of the code has been added to handle the problems that occur
!  when due to the coarsening a processor has no points, hence it has to be
!  excluded from the communicating set. The approach is to create a new communicator
!  current_mesh%comm3d per each mesh and cache to it the cartesian topology
!
    current_mesh%dnbr(1:4) = dnbr(1:4)
    current_mesh%comm3d = comm3d
    current_mesh%myid = myid
    current_mesh%rowcomm = rowcomm
    current_mesh%rowid = rowid
    current_mesh%rowsize = rowsize
    tmpcoords(1:3)=coords(1:3)
    current_mesh%pskip(1:2) = pdims(1:2) > nx_tot(nm)
    current_mesh%send_void(1:4) = 0
    current_mesh%blank = .FALSE.

    if(any(current_mesh%pskip)) then
       current_mesh%blank = (xe+1) * (ze+1) == 0
       color = MPI_UNDEFINED
       if(xe+1 >0 .AND. ze+1  >0 ) color=1
!
!    create a separate communicator for the processors that remain, and cache cart. attribute
!
       call MPI_COMM_SPLIT(comm3d, color, myid, tmpcomm, ierr)
       if (color == 1) then
         tmpdim(1)=min(nx_tot(nm),pdims(1));tmpdim(2)=min(nx_tot(nm),pdims(2));tmpdim(3)=1
         isperiodic(1) = .true.
         isperiodic(2) = .true.
         isperiodic(3) = .false.
         CALL MPI_CART_CREATE(tmpcomm, 3, tmpdim, isperiodic,&
          .false., current_mesh%comm3d, ierr)
         CALL MPI_COMM_RANK(current_mesh%comm3d,current_mesh%myid,ierr)
         call MPI_CART_COORDS(current_mesh%comm3d,current_mesh%myid,3,tmpcoords,ierr)
         call MPI_CART_SHIFT(current_mesh%comm3d,0,1,dummy1,dummy2,ierr)
         current_mesh%dnbr(1) = dummy1
         current_mesh%dnbr(2) = dummy2
 
         call MPI_CART_SHIFT(current_mesh%comm3d,1,1,dummy1,dummy2,ierr)
         current_mesh%dnbr(3) = dummy1
         current_mesh%dnbr(4) = dummy2

         remain(1) = .true.
         remain(2) = .false.
         remain(3) = .true.
         call MPI_CART_SUB(current_mesh%comm3d,remain,current_mesh%rowcomm(1),ierr)
         call MPI_COMM_RANK(current_mesh%rowcomm(1),current_mesh%rowid(1),ierr)
         call MPI_COMM_SIZE(current_mesh%rowcomm(1),current_mesh%rowsize(1),ierr)

  
         remain(1) = .false.
         remain(2) = .true.
         remain(3) = .true.
         call MPI_CART_SUB(current_mesh%comm3d,remain,current_mesh%rowcomm(2),ierr)
         call MPI_COMM_RANK(current_mesh%rowcomm(2),current_mesh%rowid(2),ierr)
         call MPI_COMM_SIZE(current_mesh%rowcomm(2),current_mesh%rowsize(2),ierr)
      else
         current_mesh%comm3d=MPI_COMM_NULL
         current_mesh%rowcomm=MPI_COMM_NULL
       endif

!
!  IF active processors on current mesh were also not blanked on the finer mesh means that when we
!  will interpolate the cgc we need to send them values. Here create send receive lists
!  and also attempt to balance the communication
!
       if (.NOT. current_mesh%fine_mesh%blank) then
         npts=-1;
         dtmp(1:4)=current_mesh%fine_mesh%dnbr(1:4)
         tmpcomm = current_mesh%fine_mesh%comm3d
         CALL MPI_SENDRECV(nx_mg(nm),1,MPI_INTEGER,dtmp(1),0,npts(2),1,&
            MPI_INTEGER,dtmp(2),0,tmpcomm,status,ierr)
  
         CALL MPI_SENDRECV(nx_mg(nm),1,MPI_INTEGER,dtmp(2),0,npts(1),1,&
            MPI_INTEGER,dtmp(1),0,tmpcomm,status,ierr)
  
         CALL MPI_SENDRECV(nz_mg(nm),1,MPI_INTEGER,dtmp(3),0,npts(4),1,&
            MPI_INTEGER,dtmp(4),0,tmpcomm,status,ierr)
  
         CALL MPI_SENDRECV(nz_mg(nm),1,MPI_INTEGER,dtmp(4),0,npts(3),1,&
            MPI_INTEGER,dtmp(3),0,tmpcomm,status,ierr)
         do ii=1,4
           npts(ii)=merge(0,1,npts(ii) == 0)
         enddo
    
!
!   if the nodes are active on the current mesh (color ==1) fill send lists
!   the load balancing is attempted: if a processor bords both vertically
!   and horiz. with a vanishing one, tmpcoords(1) == tmpcoords(2) (** if
!   pdims(1) == pdims(2)), hence the odd even startegy should prevent
!   a processor from sending both horiz and vertically limiting the overhead.
!   If pdims(1) /= pdims(2) this is not in general true.
!
         if (color == 1) then
            indx = tmpcoords(2)
            iodd = 2*int(indx/2) < indx
            if(iodd .AND. npts(1) == 0) current_mesh%send_void(1) = 1
            if((.NOT. iodd) .AND. npts(2) == 0) current_mesh%send_void(2) = 1
            indx = tmpcoords(1)
            iodd = 2*int(indx/2) < indx
            if((.NOT. iodd) .AND. npts(3) == 0) current_mesh%send_void(3) = 1
            if(iodd .AND. npts(4) == 0) current_mesh%send_void(4) = 1
!!            if(sum(current_mesh%send_void) == sum(npts)) current_mesh%send_void(3:4) = 0
!!            if(sum(current_mesh%send_void) <= sum(npts)) current_mesh%send_void(3:4) = 0
         endif
!
!   add to the sending nodes the processors in line with those with (nx,nz) = (0,0)
!
         if(npts(1) == 0 .AND. ze+1 == 0 .AND. xe+1 > 0) current_mesh%send_void(1) = 1

!
!   now form the receiving lists
!
         CALL MPI_SENDRECV(current_mesh%send_void(1),1,MPI_INTEGER,dtmp(1),0,current_mesh%recv_void(2),1,&
            MPI_INTEGER,dtmp(2),0,tmpcomm,status,ierr)
  
         CALL MPI_SENDRECV(current_mesh%send_void(2),1,MPI_INTEGER,dtmp(2),0,current_mesh%recv_void(1),1,&
            MPI_INTEGER,dtmp(1),0,tmpcomm,status,ierr)
  
         CALL MPI_SENDRECV(current_mesh%send_void(3),1,MPI_INTEGER,dtmp(3),0,current_mesh%recv_void(4),1,&
            MPI_INTEGER,dtmp(4),0,tmpcomm,status,ierr)
  
         CALL MPI_SENDRECV(current_mesh%send_void(4),1,MPI_INTEGER,dtmp(4),0,current_mesh%recv_void(3),1,&
            MPI_INTEGER,dtmp(3),0,tmpcomm,status,ierr)
      
       endif     !not blnked on the finer mesh

    endif  !if pskip:: the proc dimensions are larger than the grid dimensions

      
    current_mesh => current_mesh%coarse_mesh

  END DO     

!---------------------------------------------------------------------------------------------------
  RETURN
  END SUBROUTINE rednblack_lists
!***************************************************************************************************

END MODULE data_structure 

!***************************************************************************************************
