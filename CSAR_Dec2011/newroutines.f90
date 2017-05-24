! This subroutine is the heart of the parallelization
! here the domain decomposition is established and send receive list
! are built
! ********************************************************************
SUBROUTINE SELECT_BC_TYPE

USE GLOBAL_DATA
USE data_structure
USE MYMPI

IMPLICIT NONE

!----------------------------------------------------------------------
! Local Variables
integer :: ibnd
!----------------------------------------------------------------------

! <> ALLOCATE SOME VARIABLES
sz_bndry = max(nx_tot_fix(1)+2,ny)
ALLOCATE(type_bc(0:4,4,-2:sz_bndry))
ALLOCATE(vec_dirichlet(-maxsize+1:ndim+1,4,-2:sz_bndry))

!<> EAST-WEST PROCESSORS
is_proc_bnd(1) = coords(1) == 0
is_proc_bnd(2) = coords(1)+1 == pdims(1)
is_proc_bnd(3) = .TRUE.
is_proc_bnd(4) = .TRUE. 

addpty = 0
addbnd = zero
!.....................................

!<>  type_bc(eqn,W-E-S-N,:)   !1st_index:eqn:0temp&species,1uvel;2wvel;3vvel; 
!                              2nd_index:ibnd:1east,2west,3south;4north
!                              3rd_index:point along the boundary

!Temperature AND species
type_bc(0,1,:) =  2    !0 dirichlet   1 zero gradient  2 periodic
type_bc(0,2,:) =  2
type_bc(0,3,:) =  0  
type_bc(0,4,:) =  1
!uvelocity
type_bc(1,1,:) =  2
type_bc(1,2,:) =  2
type_bc(1,3,:) =  0
type_bc(1,4,:) =  1
!vvelocity
type_bc(3,1,:) =  2
type_bc(3,2,:) =  2
type_bc(3,3,:) =  0
type_bc(3,4,:) =  1
!Pressure !if pressure is assigned make sure to provide values
type_bc(4,1,:) =  2
type_bc(4,2,:) =  2 
type_bc(4,3,:) =  1
type_bc(4,4,:) =  0

!.........................
!----this says that if Neumann bc for the vel are assigned than the pressure
! is found from the contiuity+NS directly and is imposed at the boundary
do ibnd = 1,2
    if(ALL(type_bc(1,ibnd,:) == 1)  ) &
        type_bc(4,ibnd,:) = 0 !use pressure from y NS momentum comp
enddo 
do ibnd = 3,4
    if( ALL(type_bc(3,ibnd,:) == 1) ) &
        type_bc(4,ibnd,:) = 0 !use pressure from x NS momentum comp
enddo

!----this says that if Dirichlet bc for the vel are assigned than for pressure
! Nemann bc are assigned
do ibnd = 1,2
    if(ALL(type_bc(1,ibnd,:) == 0) ) &
        type_bc(4,ibnd,:) = 1 !use assigned grad press
enddo
do ibnd = 3,4
    if( ALL(type_bc(3,ibnd,:) == 0) ) &
        type_bc(4,ibnd,:) = 1  !use assigned grad press
enddo
!.........................

if( ALL(type_bc(1,2,:) == 0) ) then  ! dirichlet for u at east boundary
    nx_tot = nx_tot_fix + 1     
    if(is_proc_bnd(2)) then
        nx_mg = nx_mg_fix + 1
    else
        nx_mg = nx_mg_fix
    endif
else
    nx_tot = nx_tot_fix
    nx_mg = nx_mg_fix
endif

! enforce periodicity if hardwired flag are assigned
if(is_periodic_y) type_bc(:,3:4,:) = 2

if(ALL(is_periodic_xz)) then
    addptx = 0
    addbnd(1:2) = zero
    type_bc(:,1:2,:) =  2
    goto 200
endif
!
!
!
ibnd = 1  !WEST
if(  ALL(type_bc(1,ibnd,:) == 0) ) then  !Dirichlet WEST
    addbnd(ibnd) = half
    ibc(0,ibnd,1) = -1
    ibc(0,ibnd,2) = 0
    ibc(1,ibnd,1) = -1
    ibc(3,ibnd,1) = -1
    ibc(3,ibnd,2) = 0
    ibc(4,ibnd,1) = -1
elseif(  ALL(type_bc(1,ibnd,:) == 1) ) then !NEUMANN WEST
    addbnd(ibnd) = zero
    ibc(0,ibnd,1) = -1
    ibc(0,ibnd,2) = 1
    ibc(1,ibnd,1) = -1
    ibc(1,ibnd,2) = 0
    ibc(3,ibnd,1) = -1
    ibc(3,ibnd,2) = 1
    ibc(4,ibnd,1) = 0
else
    addbnd(ibnd) = 0.0d0
endif
!
ibnd = 2  !EAST
if(  ALL(type_bc(1,ibnd,:) == 0) ) then  !Dirichlet EAST
    addbnd(ibnd) = -half
    addptx(1) =  - 1
    ibc(0,ibnd,1) =  nx_mg(1)
    ibc(0,ibnd,2) = nx_mg(1)-1
    ibc(1,ibnd,1) = nx_mg(1)-1
    ibc(3,ibnd,1) =  nx_mg(1)
    ibc(4,ibnd,1) = +1         !xe+1 = nx_mg(1)
elseif(  ALL(type_bc(1,ibnd,:) == 1) ) then   !NEUMANN EAST
    addbnd(2) = zero
    addptx(0) =  + 1
    addptx(2:3) = + 1
    ibc(0,ibnd,1) =  nx_mg(1) + 1
    ibc(0,ibnd,2) = nx_mg(1) - 1
    ibc(1,ibnd,1) =  nx_mg(1)
    ibc(1,ibnd,2) = nx_mg(1) - 1
    ibc(3,ibnd,1) = nx_mg(1) + 1
    ibc(3,ibnd,2) = nx_mg(1) - 1
    ibc(4,ibnd,1) = +1         !xe+1 = nx_mg(1)
else
    addbnd(2) = 0.0d0
endif
!
200 continue
!
if(is_periodic_y)then
    goto 300
endif
ibnd = 3  !SOUTH
ibc(0,ibnd,1) = 0
ibc(0,ibnd,2) = 1
ibc(1,ibnd,1) = 0
ibc(2,ibnd,1) = 0
ibc(3,ibnd,1) = 0 

ibnd = 4  !NORTH
ibc(0,ibnd,1) = ny+1
ibc(0,ibnd,2) = ny-1
ibc(1,ibnd,1) = ny+1
ibc(1,ibnd,2) = ny-1
ibc(2,ibnd,1) = ny+1
ibc(2,ibnd,2) = ny-1
ibc(3,ibnd,1) = ny
ibc(3,ibnd,2) = ny-1
!  if(ALL(type_bc(3,4,:) == 1))&
!       addpty(0:2) = + 1
!  if(ALL(type_bc(3,3,:) == 0)) addbnd(3) = half
!  if(ALL(type_bc(3,4,:) == 0)) addbnd(4) = -half

300  continue

!use extrapolation fo higher order BOUNDARY CONDITIONS
use_extrap = .FALSE.
!  use_extrap(0,1) = .TRUE.
!  use_extrap(2:3,1) = .TRUE.
!  use_extrap(0,3) = .TRUE.
use_extrap(1:2,3) = .TRUE.

! set injection velocities. These values are used only
! if dirichlet BC are used, Dimensions are cm/sec
! first index is the velocity component and last is the  boundary
!
Q_scalar = zero
Q_scalar(1,1) = zero; !it's the inflow velocity at the west boundary

!--------------------------------------------------------
return
END SUBROUTINE SELECT_BC_TYPE
!******************************************************************************

!*******************************************************************************
SUBROUTINE DIRICHLET_BC(eqn,current_mesh)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE
  include "parallel.h"

!----------------------------------------------------------------
! Local variables
  INTEGER, INTENT(INOUT) :: eqn
  TYPE(mesh_pointers),INTENT(INOUT) :: current_mesh
  INTEGER :: i,ibnd,iref(4),n
  INTEGER :: i1,i2,i3
  INTEGER :: j1,j2,j3
  INTEGER ::  k, j, xs,zs,xe,ze,ys,ye
  INTEGER :: eqn_check,ic_in(4),ic_out(4),is(4),eqn0
  REAL*8 :: c_xtr(0:2)
  LOGICAL :: do_implicit
!----------------------------------------------------------------
!both Dirichlet and Neumann BC


  c_xtr(0) = 8.0d0/3.0d0
  c_xtr(1) = -6.0d0/3.0d0
  c_xtr(2) = one/3.0d0
  is(1) = 1
  is(2) = -1
  is(3) = 1
  is(4) = -1

  do_implicit = eqn < 0 
  eqn = mod(abs(eqn),10)
  eqn_check = abs(min(eqn,4))
  if(is_periodic_xz(eqn_check)) GOTO 200

  print*,'DOING DIRICHLET BOUND COND'

  SELECT CASE(eqn)

  case(5)

     xs = 0               ; ys = 0              ; zs = 0;    
     xe = current_mesh%xe ; ye = current_mesh%ye; ze = current_mesh%ze;
     iref(1) = 0
     iref(2) = xe 
     do ibnd = 1,2
        IF(is_proc_bnd(ibnd)) THEN
           i1 = iref(ibnd) + ibc(4,ibnd,1)
           i2 = iref(ibnd) + ibc(4,ibnd,1) - is(ibnd)
           do j = 0,ye+1
              wrk(min(i1,i2):max(i1,i2),:,j) =   zero
           enddo
        ENDIF
     enddo

  case(4)  !do west/east 

     xs = 0               ; ys = 0              ; zs = 0;    
     xe = current_mesh%xe ; ye = current_mesh%ye; ze = current_mesh%ze;
     iref(1) = 0
     iref(2) = xe


     do ibnd = 1,2
        IF(is_proc_bnd(ibnd)) THEN
           i =  iref(ibnd) + ibc(eqn,ibnd,1)
           i3 = i + is(ibnd)
           do j = 0,ye
              if(type_bc(4,ibnd,j) == 1) then
                 current_mesh%x(i,:,j) =   zero
              elseif(type_bc(4,ibnd,j) == 0) then
                 if(current_mesh%mesh_num == 1) then
                    current_mesh%x(i,:,j) = vec_dirichlet(4,ibnd,j)
                 else
                    current_mesh%x(i,:,j) =   zero
                 endif
              endif
              if(type_bc(1,ibnd,j) == 0 .AND. do_implicit)then
                 current_mesh%c9(i3,:,j) = current_mesh%c9(i3,:,j)&
                      - current_mesh%c1
              endif
           enddo
        ENDIF
     enddo



  case(0)  !temperature species surface-function  boundary conditions

     xs = drange(1)
     xe = nxv(eqn)
     zs = drange(3)
     ze = drange(4)
     ys = 0
     ye = nyv(eqn)+1
     ic_out(1) = 2
     ic_in(1) = 3
     ic_out(2) = 3
     ic_in(2) = 2

!WEST-EAST  
     do ibnd = 1,2
        do k = zs,ze
           IF(is_proc_bnd(ibnd)) THEN
              i1 = ibc(eqn,ibnd,1)
              i2 = ibc(eqn,ibnd,2)
              do j = ys,ye
                 if(type_bc(eqn,ibnd,j) == 0) then
                    do n = 1,neqmax
                       f(i1,k,j,n) = vec_dirichlet(0-n+1,ibnd,j)
                    enddo
                 elseif(type_bc(eqn,ibnd,j) == 1) then
                    f(i1,k,j,1:neqmax) = f(i2,k,j,1:neqmax)
!                     if( j == ys) phi(i1,k) = phi(i2,k)
                 endif
                 if(type_bc(eqn,ibnd,j) == 1 .AND. do_implicit) then
                    if( max(i1,i2) - min(i1,i2) == 2 ) then
                       i3 = (i1+i2)/2
                       coeff(ic_in(ibnd),i3,k,j) = coeff(ic_in(ibnd),i3,k,j) & 
                            + coeff(ic_out(ibnd),i3,k,j)
                       coeff(ic_out(ibnd),i3,k,j) = zero
                       csld(ic_in(ibnd),i3,k,j) = csld(ic_in(ibnd),i3,k,j) & 
                            + csld(ic_out(ibnd),i3,k,j)
                       csld(ic_out(ibnd),i3,k,j) = zero
                    elseif( max(i1,i2) - min(i1,i2) == 1 ) then
                       i3 = i2
                       coeff(ic_out(ibnd),i3,k,j) = zero
                       csld(ic_out(ibnd),i3,k,j) = zero
                    else
                       STOP 'ERROR BC_DIRICHLET'
                    endif
                 endif
              enddo
           ENDIF
        enddo
     enddo


  CASE(1:3)

     xs = drange(1)
     xe = nxv(eqn)
     zs = drange(3)
     ze = drange(4)
     ys = 0
     ye = nyv(eqn)+1
     ic_out(1) = 2
     ic_in(1) = 3
     ic_out(2) = 3
     ic_in(2) = 2

!WEST-EAST  
     do ibnd = 1,2
        do k = zs-1,ze+1
           IF(is_proc_bnd(ibnd)) THEN
              i1 = ibc(eqn,ibnd,1)
              i2 = ibc(eqn,ibnd,2)
              do j = ys,ye
                 if(type_bc(eqn,ibnd,j) == 0) then
                    if(.NOT. use_extrap(eqn,ibnd))then
                       q(i1,k,j,eqn) = vec_dirichlet(eqn,ibnd,j)
                    else
                       q(i1,k,j,eqn) = c_xtr(0)*vec_dirichlet(eqn,ibnd,j)&
                            +  c_xtr(1)*q(i1+is(ibnd)*1,k,j,eqn) &
                            +  c_xtr(2)*q(i1+is(ibnd)*2,k,j,eqn)                       
                    endif
                 elseif(type_bc(eqn,ibnd,j) == 1) then
                    q(i1,k,j,eqn) = q(i2,k,j,eqn)
                 endif
                 if(do_implicit) then
                    if(type_bc(eqn,ibnd,j) == 0 .AND. use_extrap(eqn,ibnd)) then
                       i3 = i1+is(ibnd)*1
                       coeff(ic_in(ibnd),i3,k,j) = coeff(ic_in(ibnd),i3,k,j)&
                            + c_xtr(2)*coeff(ic_out(ibnd),i3,k,j)
                       dqdt(i3,k,j,eqn) = dqdt(i3,k,j,eqn) &
                            +c_xtr(0)*vec_dirichlet(eqn,ibnd,j)*coeff(ic_out(ibnd),i3,k,j)
                       q(i1,k,j,eqn) = zero        !because it was included in the RHS
                       coeff(ic_out(ibnd),i3,k,j) = c_xtr(0)*coeff(ic_out(ibnd),i3,k,j)
                    elseif(type_bc(eqn,ibnd,j) == 1 ) then
                       if( max(i1,i2) - min(i1,i2) == 2 ) then
                          i3 = (i1+i2)/2
                          coeff(ic_in(ibnd),i3,k,j) = coeff(ic_in(ibnd),i3,k,j)&
                               + coeff(ic_out(ibnd),i3,k,j)
                          coeff(ic_out(ibnd),i3,k,j) = zero
                       elseif( max(i1,i2) - min(i1,i2) == 1 ) then
                          i3 = i2
                          coeff(ic_out(ibnd),i3,k,j) = zero
                       else
                          STOP 'ERROR BC_DIRICHLET Q'
                       endif
                    endif
                 endif  !do_implicit
              enddo  !j
           ENDIF  !is_proc_bnd(ibnd)
        enddo
     enddo  !ibnd 

  case(6)  ! surface-function  boundary conditions

     eqn0 = 0
     xs = drange(1)
     xe = nxv(eqn0)
     zs = drange(3)
     ze = drange(4)
     j = 0

!WEST-EAST  
     do ibnd = 1,2
        do k = zs,ze
           IF(is_proc_bnd(ibnd)) THEN
              i1 = ibc(eqn0,ibnd,1)
              i2 = ibc(eqn0,ibnd,2)
              if(type_bc(eqn0,ibnd,0) == 0) then
                 print*,'INVALID OPTION'
                 stop 'newroutines.f90'
              elseif(type_bc(eqn0,ibnd,j) == 1) then
                 phi(i1,k) = phi(i2,k)
                 phi(i1-is(ibnd),k) = phi(i2+is(ibnd),k)
              endif
           ENDIF
        enddo
     enddo

  END SELECT


200 CONTINUE

  if(is_periodic_y)goto 400


!SOUTH NORTH BOUNDARY CONDITIONS

  SELECT CASE(eqn)


  case(1:3)   !SOUTH-NORTH velocity boundary conditions

     xs = drange(1)
     xe = nxv(eqn)
     zs = drange(3)
     ze = drange(4)
     ys = 0
     ye = nyv(eqn)+1
     ic_out(3) = 4
     ic_in(3) = 5
     ic_out(4) = 5
     ic_in(4) = 4
     do ibnd = 3,3   !should be 3,4  but dont do 4 for now
        do k = zs,ze
           IF(is_proc_bnd(ibnd)) THEN
              j1 = ibc(eqn,ibnd,1)
              j2 = ibc(eqn,ibnd,2)
              do i = xs,xe
                 if(type_bc(eqn,ibnd,i) == 0) then  
                    q(i,k,j1,eqn) = vec_dirichlet(eqn,ibnd,i)
                 elseif(type_bc(eqn,ibnd,i) == 1) then
                    q(i,k,j1,eqn) = q(i,k,j2,eqn)
                 endif
                 if(do_implicit) then              
                    if(type_bc(eqn,ibnd,i) == 0 .AND. use_extrap(eqn,ibnd)) then
                       j3 = j1+is(ibnd)*1
                       coeff(ic_in(ibnd),i,k,j3) = coeff(ic_in(ibnd),i,k,j3)&
                            + c_xtr(2)*coeff(ic_out(ibnd),i,k,j3)
                       dqdt(i,k,j3,eqn) = dqdt(i,k,j3,eqn) &
                            +c_xtr(0)*qsurf(i,k,eqn)*coeff(ic_out(ibnd),i,k,j3)
!                       q(i,k,j1,eqn) = zero     !because a contrib. was added to RHS
                       coeff(ic_out(ibnd),i,k,j3) = c_xtr(0)*coeff(ic_out(ibnd),i,k,j3)
                    elseif(type_bc(eqn,ibnd,i) == 1) then
                       if( max(j1,j2) - min(j1,j2) == 2 ) then
                          j3 = (j1+j2)/2
                          coeff(ic_in(ibnd),i,k,j3) = coeff(ic_in(ibnd),i,k,j3) & 
                               + coeff(ic_out(ibnd),i,k,j3)
                          coeff(ic_out(ibnd),i,k,j3) = zero
                       elseif( max(j1,j2) - min(j1,j2) == 1 ) then
                          j3 = j2
                          coeff(ic_out(ibnd),i,k,j3) = zero
                       else
                          STOP 'ERROR BC_DIRICHLET T N-S'
                       endif
                    endif
                 endif
              enddo
           ENDIF
        enddo
     enddo


  END SELECT

400 CONTINUE
!
!HERE ENFORCE 2 DIMENSIONALITY OF THE FLOW
  if(is2D) then
     q(:,-1,:,:) = q(:,0,:,:)
     q(:,1,:,:) = q(:,0,:,:)
     f(:,-1,:,:) = f(:,0,:,:)
     f(:,1,:,:) = f(:,0,:,:)
     phi(:,-1) = phi(:,0)
     phi(:,-2) = phi(:,0)
     phi(:,1) = phi(:,0)
     phi(:,2) = phi(:,0)
  endif


!-------------------------------------------------------------------------------
  RETURN
END SUBROUTINE DIRICHLET_BC
!*******************************************************************************

!*************************************************
 SUBROUTINE SET_CONDITIONS_BC
!
!sets the boundary values
!
USE GLOBAL_DATA
IMPLICIT NONE
!--------------------------------------------------------------
! Local Variables
INTEGER :: i, k, j, eqn, ibnd
REAL*8  :: exp_u
!--------------------------------------------------------------  

exp_u = one/3.0d0
do eqn = 1,3,2
    do ibnd = 1,2
        if( ANY(type_bc(eqn,ibnd,:) == 0) ) then
            do j = 1,nyv(eqn)+1
                vec_dirichlet(eqn,ibnd,j) = Q_scalar(eqn,ibnd) *abs(ya(j)/yend)**exp_u !fixed in subroutine SELECT_BC_TYPE
            enddo
        endif
    enddo
enddo

! use extrapolation for the value at zero
vec_dirichlet(1,1,0) =(8.0d0*zero-6.0d0*vec_dirichlet(1,1,1)+vec_dirichlet(1,1,2))*third
!-----------------------------------------------------------
RETURN
END SUBROUTINE SET_CONDITIONS_BC
!*************************************************************

!*************************************************
 SUBROUTINE UPDATE_CONDITIONS_BC
!
!updates the boundary values if they vary with time
!
   USE GLOBAL_DATA
   IMPLICIT NONE
!--------------------------------------------------------------
! Local Variables
   INTEGER :: i, k, j, eqn, error, ibnd
   INTEGER :: xs,zs,xe,ze,ys,ye
   INTEGER :: ii, jj
   REAL*8  :: p_north,p_east,u_west
   REAL*8  :: coe_rho,rhogas,poo,pressure, dya, dxa
   REAL*8  :: exp_u
!--------------------------------------------------------------  

   poo = PRESSURE(tcyc)
   coe_rho = dim_fact*poo

   xs = drange(1)
   xe = drange(2)
   zs = drange(3)
   ze = drange(4)
   ys = drange(5)+1
   ye = finest_mesh%ye+1
!
! THIS is a temporary fix for periodic grids
!
   if(ALL(type_bc(3,4,:) == 1) .AND. ALL(type_bc(1,4,:) == 1) .AND. ALL(is_periodic_xz) ) then
      type_bc(1,4,:) = 0
      use_extrap(1,4) = .FALSE.
      vec_dirichlet(1,4,:)  =  zero
      type_bc(4,4,:) = 0
   endif

!
!these two are necessary if you have inflow west or east BC
!
!pressure is determined from NS + continuity
   ibnd = 1
   if( ALL(type_bc(1,ibnd,:) == 1) ) then

      ii = ibc(4,ibnd,1);
      k=0
      v_south(0) = vvel(ii,k,ibc(3,3,1))
      do j = ys,ye
         dya = detady(j)/dy
         v_south(j) = v_south(j-1) - rate(ii,k,j,4)/dya
      enddo

      p_north = zero
      vec_dirichlet(4,ibnd,ye) = p_north
      do k=zs,ze
         do j=ye,ys,-1
            dya = detadya(j)/dy
            jj = j-1
            vec_dirichlet(4,ibnd,jj) = vec_dirichlet(4,ibnd,jj+1) &
                 - (vvel(ii,k,j) - v_south(j)) / dya
         enddo
      enddo
   endif

   ibnd = 2
   if( ALL(type_bc(1,ibnd,:) == 1) ) then      
      
      ii = finest_mesh%xe + ibc(4,ibnd,1)
      k=0
      v_south(0) = vvel(ii,k,ibc(3,3,1))
      do j = ys,ye
         dya = detady(j)/dy
         v_south(j) = v_south(j-1) - rate(ii,k,j,4)/dya
      enddo
      
      p_north = zero
      vec_dirichlet(4,ibnd,ye) = p_north
      do k=zs,ze
         do j=ye,ys,-1
            dya = detadya(j)/dy
            jj = j-1
            vec_dirichlet(4,ibnd,jj) = vec_dirichlet(4,ibnd,jj+1) &
                 - (vvel(ii,k,j) - v_south(j))/dya
         enddo
      enddo
   endif

!
! now update time varing bc for u,v at the east west boundaries if present
!
   if(time_varying_inflow) then    !time_varying_inflow is set in newroutines.f90
      exp_u = one/3.0d0
      do eqn = 1,3,2
         do ibnd = 1,2
            
            if( ANY(type_bc(eqn,ibnd,:) == 0) ) then
               do j = 1,nyv(eqn)+1
                  vec_dirichlet(eqn,ibnd,j) = Q_scalar(eqn,ibnd) *abs(ya(j)/yend)**exp_u !fixed in subroutine SELECT_BC_TYPE
               enddo
            endif

         enddo
      enddo
   endif

!-----------------------------------------------------------
   RETURN
 END SUBROUTINE UPDATE_CONDITIONS_BC
!*************************************************************

!**********************************************************
SUBROUTINE NAN_VAL(fq,outcome,neq)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!-----------------------------------------------------------
!!!   Local Variables
  integer i, k, j, n,neq
  REAL*8 :: fq(neq), mx(neq),mn(neq)
  LOGICAL :: outcome
!-----------------------------------------------------------

!

  mx(1) = 4000.0d0
  mn(1) = 310.0d0
  do n=2,neq
     mn(n) = -1.0
     mx(n) = 1.d0 - mn(n)
  enddo


  outcome = .FALSE.
  do n=1,neq
     outcome = outcome .OR. (.NOT. fq(n) > mn(n)) .OR. (.NOT. fq(n) < mx(n))
  enddo

!--------------------------------------------------------
  RETURN
END SUBROUTINE NAN_VAL
!**********************************************************************
!

!**********************************************************
SUBROUTINE NAN_FIELD(outcome,doprint)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!-----------------------------------------------------------
!!!   Local Variables
  INTEGER :: i,k,j
  LOGICAL :: outcome,repeat,doprint
!-----------------------------------------------------------

!
!!check for NANs
  repeat = .FALSE.
  do i = drange(1),drange(2)
     do k = drange(3),drange(4)
        do j = drange(5)+1, drange(6)-1

           call NAN_VAL(f(i,k,j,1:neqgas),outcome,neqgas)

           repeat = repeat .or. outcome

           if(outcome.and.doprint)then
!              print*,ncyc,myid,i,j,ny,&
!                   'OUT OF RNGE',f(i,k,j,1:neqgas),'OLD',&
!                   oldsoln(i,k,j,1:neqgas),'COEFFs',coeff(:,i,k,j),&
!                   'RHS',dfdt(i,k,j,1:neqgas),'RATE',rate(i,k,j,1:neqgas),&
!                   'TRANSP',dqdt(i,k,j,1:neqgas)
              write(138,*)ncyc,myid,i,j,f(i,k,j,1:neqgas), 'CCC',coeff(1:5,i,k,j),&
                   'DDD',dfdt(i,k,j,1:neqgas)
!               
           endif
        enddo
     enddo
  enddo

  call MPI_ALLREDUCE(repeat,outcome,1,MPI_LOGICAL,MPI_LOR,comm3d,ierr)

!--------------------------------------------------------
  RETURN
END SUBROUTINE NAN_FIELD
!**********************************************************************
!

!**********************************************************************
SUBROUTINE STORAGE

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!---------------------------------------------------------------------
! Local variables
  INTEGER :: i,j,k,n,eqn 
  REAL*8 ::  heat_storage,heat_storage_old,all_storage 
  REAL*8 ::  heat_storage_der
!---------------------------------------------------------------------


  heat_storage = zero
  heat_storage_old = zero
  heat_storage_der = zero
  do j = drange(5), drange(6)
     do k = drange(3), drange(4)
        do i = drange(1), drange(2)
           heat_storage = heat_storage  &
                + rhos(i,k,j)*f(i,k,j,neqmax)*dy/detady(j)
           heat_storage_old = heat_storage_old  &
                + rhos_old(i,k,j)*oldsoln(i,k,j,neqmax)*dy/detady(j)
           heat_storage_der = heat_storage_der +  rate(i,k,j,maxsize)*dy/detady(j)
        enddo
     enddo
  enddo


  call MPI_ALLREDUCE(heat_storage,all_storage,1,&
       MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
  heat_storage = all_storage/(nx*nz)
  call MPI_ALLREDUCE(heat_storage_old,all_storage,1,&
       MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
  heat_storage_old = all_storage/(nx*nz)
  call MPI_ALLREDUCE(heat_storage_der,all_storage,1,&
       MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
  storage_term_3 = all_storage/(nx*nz)


  storage_term = -(heat_storage-heat_storage_old)/timestep


  RETURN
END SUBROUTINE STORAGE
 


  
