! ********************************************************************
SUBROUTINE UPDATE_PHI(t,tout,iflag)

  USE GLOBAL_DATA
  USE MYMPI
  Use Restart_Dump
  IMPLICIT NONE

!---------------------------------------------------------------------
  REAL*8, INTENT(IN) :: t,tout
  integer,INTENT(IN):: iflag
! Local variables
  INTEGER :: i,k,j,imax,kmax,xs,xe,zs,ze,imax2nd,kmax2nd
  INTEGER :: imaxtime,itime,imin,kmin
  REAL*8  :: h1,h2,h3,h4,tmp1,phixAVG,phizAVG,rbxavg,rbzavg
  REAL*8  :: h1a,h3a,maxvec,allmaxvec,maxphi,maxfxxPlusxfzz
!---------------------------------------------------------------------


  if (iflag == 0) GOTO 100


  CALL PYROLYSIS

  if(ipack == 0) then
     Imaxtime = merge(60, 25, ncyc > 600)
  else
     Imaxtime = 1
  end if
  do itime = 1,Imaxtime
     oldphi = phi
     CALL FLUX(t,tout)
     CALL FILL_PHI_GHOST
  enddo

100 CONTINUE

  CALL FILL_PHI_GHOST

!--- COMPUTE shape derivatives
  h1 = 1.0d0/(2.0d0*dx)
  h1a = 1.0d0/dx
  h2 = 1.0d0/(dx*dx)
  h3 = 1.0d0/(2.0d0*dz)
  h3a = 1.0d0/dz
  h4 = 1.0d0/(dz*dz)

  xs = drange(1);xe = drange(2);zs=drange(3);ze = drange(4)

  call shape_derivatives

  do k = drange(3)-1, drange(4)+1
     do i = drange(1)-1, drange(2)+1
        vec0(i,k) = sqrt(1.0d0+dphidx(i,k)**2+dphidz(i,k)**2) 
        vec(i,k,0) = rb(i,k)*vec0(i,k)
        rbxavg = half*(rb(i,k)+rb(i+1,k)) 
        rbzavg = half*(rb(i,k)+rb(i,k+1))
        phizAVG = half*(dphidz(i,k)+dphidz(i+1,k)) 
        phixAVG = half*(dphidx(i,k)+dphidx(i,k+1))
        vec(i,k,1) = rbxavg*sqrt(1.0d0+ dphidxa(i,k)**2 + phizAVG**2) 
        vec(i,k,2) = rbzavg*sqrt(1.0d0+ phixAVG**2 + dphidza(i,k)**2)
     end do
  end do

  if (maxphi == fmax) then
     tmp1 = vec(imax,kmax,0)
  else 
     tmp1 = zero
  endif
  CALL MPI_ALLREDUCE(tmp1,velofmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)


  if(num_mg_meshes > 1)then
     CALL UPDATE_PHI_MG
  endif

  MQphi = phi - fmax
  phi_mg = phi_mg - fmax
  do k = drange(3)-1, drange(4)+1
     do i = drange(1)-1, drange(2)+1
        phi_mg(i,k,0) = MQphi(i,k) 
        dphidxa_mg(i,k,0) = dphidxa(i,k)
        dphidza_mg(i,k,0) = dphidza(i,k) 
        dphidx_mg(i,k,0) = dphidx(i,k)
        dphidz_mg(i,k,0) = dphidz(i,k)
        psi_mg(i,k,0) = psi(i,k,0)
        temp_mg(i,k,0) = f(i,k,0,1)
     enddo
  enddo


!
! NOTE, each processor will have a different phi_per, t_per,
! but we will consider only the one for processor (0)
!
  CALL GRID_SPEED

  if(iflag == 0) return

  done_period = 0
  if(phi(0,0) < phi_per + period .AND. &         ! [phi < 0]
       oldphi(0,0) >= phi_per + period) then

     phi_per = oldphi(0,0)
     t_per = t

     IF(myid == 0) then
        write(*,*) 'NEW PERIOD', phi_per,phi(0,0),oldphi(0,0),t_per,period
        OPEN (73,file='period.info',form='formatted',action="write",POSITION='APPEND')
!        OPEN (73,file='period.info',form='formatted',action="write")
        write(73,*) 'PERIOD REACHED AT ncyc', ncyc
        write(73,*) 'NEW VALUES FOR T_per and PHI_per are :',t_per, phi_per
        write(73,*) '--------------------------'
        CLOSE (73)
        done_period = 1
     ENDIF

  endif

  CALL MPI_BCAST(done_period,1,MPI_INTEGER,0,comm3d,ierr)

  if(done_period == 1) then
     CALL DUMP(t,outputfile)
     speed_timeint = 0.0d0
     dt_timeint = 0.0
  endif

!---------------------------------------------------------------------
  RETURN

contains

!--------------------------------------------------
  subroutine shape_derivatives

    integer :: i,k,it,itimesaround
    real*8  :: maxima_local(2),maxima_global(2)
    logical :: flag_local,flag_global

    itimesaround = 0
1   continue
    itimesaround = itimesaround+1

    if(itimesaround > 10) return

    maxphi = -1d99
    maxfxxPlusxfzz =  -1d99
    if (skipDeformation) then
!if a flat surface was to be considered
       dphidx = zero
       dphi2dx2 = zero
       dphidz = zero
       dphi2dz2 = zero
       dphidxa = zero
       dphidza = zero
       dphidxa_mg = zero 
       dphidza_mg = zero
       dphidx_mg = zero
       dphidz_mg = zero
       phi_mg = zero
       fmax = zero
       imax = 0
       kmax = 0
       imax2nd = 0
       kmax2nd = 0
    else
       taGsmooth = .false.
       do k = drange(3)-1, drange(4)+1
          do i = drange(1)-1, drange(2)+1
             dphidx(i,k) = (phi(i+1,k)-phi(i-1,k))*h1
             dphi2dx2(i,k) = (phi(i+1,k)&
                  -2.0d0*phi(i,k)+phi(i-1,k))*h2
             dphidz(i,k) = (phi(i,k+1)-phi(i,k-1))*h3
             dphi2dz2(i,k) = (phi(i,k+1)&
                  -2.0d0*phi(i,k)+phi(i,k-1))*h4
             dphidxa(i,k) = (phi(i+1,k)-phi(i,k))*h1a
             dphidza(i,k) = (phi(i,k+1)-phi(i,k))*h3a
             if(abs(phi(i,k) -(sum(phi(i-2:i+2,k-2:k+2))- phi(i,k))/(size(phi(i-2:i+2,k-2:k+2))-1)) > 1d-2 &
                  &.AND. (k >= zs .AND. k<= ze .AND. i>= xs .AND. i<= xe)) then
                taGsmooth(i,k) = .true.
             endif
             if(abs(dphidx(i,k) + dphidz(i,k))  > maxfxxPlusxfzz .AND. (k >= zs .AND. k<= ze .AND. i>= xs .AND. i<= xe)) then
                maxfxxPlusxfzz = dphidx(i,k) + dphidz(i,k) 
                imax2nd = i
                kmax2nd = k
             endif
             if(phi(i,k) > maxphi.AND. (k >= zs .AND. k<= ze .AND. i>= xs .AND. i<= xe)) then
                maxphi = phi(i,k)
                imax = i
                kmax = k
             endif
          end do
       end do
    endif

    maxima_local = (/ maxphi,maxfxxPlusxfzz/)
    CALL MPI_ALLREDUCE(maxima_local,maxima_global,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
    flag_local = count(tagsmooth) > 0
    CALL MPI_ALLREDUCE(flag_local,flag_global,1,MPI_LOGICAL,MPI_LOR,comm3d,ierr)

    if(is2d) then
       dphidz = 0d0
       dphi2dz2 = 0d0
       dphidza =  0d0
    end if
       

! Local surface smoothing
!!>    if(flag_global) then
!!>       if(maxima_global(2) == maxima_local(2) .and. ipack == 1) then
!!>          write(*,*)myid,'MAX(2)',imax2nd,kmax2nd,maxima_global(2),dphidx(imax2nd,kmax2nd),dphidz(imax2nd,kmax2nd)
!!>          write(*,*) ncyc,'SMOOTHING',flag_local
!!>       end if
!!>       do it = 1,5
!!>          do i = drange(1),drange(2)
!!>             do k = drange(3),drange(4)
!!>                if(tagsmooth(i,k))then
!!>                   phi(i,k) = (sum(phi(i-2:i+2,k-2:k+2))- phi(i,k))/(size(phi(i-2:i+2,k-2:k+2))-1)
!!>                endif
!!>             enddo
!!>          enddo
!!>       enddo
!!>       CALL FILL_PHI_GHOST
!!>       goto 1
!!>    end if

    return
  end subroutine shape_derivatives
!----------------------------------------------------------
END SUBROUTINE UPDATE_PHI
!*********************************************************************


!*********************************************************************
SUBROUTINE UPDATE_PHI_MG

  USE GLOBAL_DATA
  IMPLICIT NONE

  include "parallel.h"

!--------------------------------------------------------------------
! local variables
  INTEGER ::  i,k,j, nxc,nzc, im,km,jm,f2c,col(6),msh,imsh,itimes
  REAL *8 :: h1a, h1
  TYPE(mesh_pointers), POINTER :: current_mesh
!---------------------------------------------------------------------

! just inject the value of phi in wrk vector !
! note, im,km are the indexes  on the fine mesh that correspond to i,k
! on mesh number <mesh>
! 

  timesdo: Do itimes = 1,3
     meshesdo: DO j=2,iex
        jm = j-1
        msh = lx_mg(j)
        nxc=nx_mg(msh)
        nzc=nz_mg(msh)
        f2c = 2**(jm)
!
! Find the mesh corresponding to mesh number <msh>.
!
        current_mesh => finest_mesh
        do imsh = 2,msh
           current_mesh => current_mesh%coarse_mesh
        enddo
        if(current_mesh%blank) EXIT

        do i=0,nxc
           im=f2c*(i+ib_mg(msh))-ib_mg(1)
           do k=0,nzc
              km=f2c*(k+kb_mg(msh))-kb_mg(1)
              if(itimes == 1) then
                 wrk(i,k,jm) = phi(im,km)
              elseif(itimes == 2) then
                 wrk(i,k,jm) = psi(im,km,0)
              elseif(itimes == 3) then
                 wrk(i,k,jm) = f(im,km,0,1)
              end if
           enddo
        enddo

!
! update the ghost cells for each level
! now that we pass also the mesh pointer maybe us the case to avoid 
! passing the wrk vecto as a separate pointer but we can pass 
! current_mesh%f which would be included in current_mesh
!
        col=0; col(5:6)=jm; col(2) = nxc-1; col(4)=nzc-1
        call parallel_swap(wrk,col,current_mesh)


        do i= -1,nxc
           do k= -1,nzc
              if(itimes == 1) then
                 phi_mg(i,k,jm) = wrk(i,k,jm)
              elseif(itimes == 2) then
                 psi_mg(i,k,jm) = wrk(i,k,jm)
              elseif(itimes == 3) then
                 temp_mg(i,k,jm) = wrk(i,k,jm)
              end if
           enddo
        enddo

        if (itimes >= 2) cycle meshesdo

!
! divide by dx  (h1a==h3a) to simplify the implementation
!
        h1a = 1.0d0/(dx*dble(f2c))
        h1 = h1a/two

        do i= -1,nxc-1
           do k= -1,nzc-1
              dphidxa_mg(i,k,jm) = (wrk(i+1,k,jm)-wrk(i,k,jm))*h1a 
              dphidza_mg(i,k,jm) = (wrk(i,k+1,jm)-wrk(i,k,jm))*h1a
           enddo
        enddo

        do i= 0,nxc-1
           do k= 0,nzc-1
              dphidx_mg(i,k,jm) = (wrk(i+1,k,jm)-wrk(i-1,k,jm))*h1
              dphidz_mg(i,k,jm) = (wrk(i,k+1,jm)-wrk(i,k-1,jm))*h1
           enddo
        enddo

     ENDDO meshesdo  !do j=2,iex
  end Do timesdo
  if(is2D) THEN
     dphidza_mg = zero
     dphidz_mg = zero
  ENDIF

!------------------------------------------------------------------------
  RETURN
END SUBROUTINE UPDATE_PHI_MG
!*************************************************************************

! ********************************************************************
SUBROUTINE PROPAGATE_SURF(t,tout,MaxNumofTimes,done)

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

!---------------------------------------------------------------------
  REAL*8, INTENT(INOUT) :: t,tout
  INTEGER, INTENT(IN) :: MaxNumofTimes
  LOGICAL, INTENT(OUT) :: done
! Local variables
  INTEGER :: itime
!---------------------------------------------------------------------
  if(MaxNumofTimes <= 0) then  ! added to simulate non-propagating cases
     tout = t + dble(abs(MaxNumofTimes)) * dt
     done = .false.
     t = tout
     CALL UPDATE_PHI(zero,zero,1)
     CALL UPDATE_PSI(t,tout)
     CALL SOLID_PROP
     CALL PYROLYSIS
     return
  endif
!
  do itime = 1,MaxNumofTimes

     tout = t + dt
     CALL UPDATE_PHI(t,tout,1)
     CALL UPDATE_PSI(t,tout)
     CALL SOLID_PROP
     CALL PYROLYSIS   
     t = tout

     done = abs(t-tstop) <= 1.0d-10

     if(done) EXIT

  enddo


!---------------------------------------------------------------------
  RETURN
END SUBROUTINE PROPAGATE_SURF
!*********************************************************************
! ********************************************************************                                      
SUBROUTINE GRID_SPEED

  USE GLOBAL_DATA
  IMPLICIT NONE

!--------------------------------------------------------------------- 
  INTEGER :: i,k,j
!--------------------------------------------------------------------- 

  do k = drange(3)-1, drange(4)+1
     do i = drange(1)-1, drange(2)+1
        do j = 0,ny
           MQvelocity(i,k,j,0:2) = MQchi(j,1)*vec(i,k,0:2) + velofmax*(one-MQchi(j,1))
           MQvelocity(i,k,j,3)   = MQchi(j,2)*vec(i,k,0)    + velofmax*(one-MQchi(j,2))
        enddo
     enddo
  enddo


  return
END SUBROUTINE GRID_SPEED
