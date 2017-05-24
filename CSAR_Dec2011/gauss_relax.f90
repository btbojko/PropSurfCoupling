!     ****************************************************************
!     *                                                              *
!     *                      subroutine gauss_relax                  *
!     *                                                              *
!     *                             solve Ax=f                       *
!     *                                                              *
!     *     with red and black decomposition in the x &  z           *
!     *                                                              *
!     *                                                              * 
!     ****************************************************************
!



!***************************************************************************************************

SUBROUTINE gauss_relax ( current_mesh, nu, resid_flag )

  USE data_types
  USE GLOBAL_DATA 
  USE MYMPI

  IMPLICIT NONE

  include "parallel.h"

! Subroutine performs NU cycles of GS
! relaxation on mesh CURRENT_MESH.
! NOTE ASSUME NO PROCESSOR_PARTITIONING IN Y
! a 2 colors scheme works  if there is no partitioning in y direction
! there is no cross term in the x-z plane. This implementation also requires
! that the  number of pts in x and z is even
! This is necessarily true except on the coarsest grid

!---------------------------------------------------------------------------------------------------
!   Dummy variables:

  TYPE(mesh_pointers) :: current_mesh
!    TYPE(mesh_pointers), POINTER :: current_mesh

  INTEGER, INTENT(IN) :: nu                ! number of relaxations 
  LOGICAL, INTENT(IN) :: resid_flag        ! if this is true it evaluates the residual

!   Local Variables
  REAL(KIND=double), DIMENSION(:,:,:), POINTER :: px
  REAL(KIND=double), DIMENSION(:,:,:), POINTER :: fx

  REAL(KIND=double), DIMENSION(:,:,:), POINTER :: cd,c3,c7,c9,cj
  REAL(KIND=double), DIMENSION(:), POINTER :: ey

  INTEGER :: icount,i,j,k,xe,ze,ye,nuu
  INTEGER :: indx,icolor,rcolor, tmpcomm
  INTEGER :: js, jmax, imax, kmax, ib, jb

  REAL(KIND=double) :: c1,resmax,myres,omg_ppe,delta_p
!---------------------------------------------------------------------------------------------------

!
!   RETURN immediatly if the processor has no points
!
  if(current_mesh%blank) RETURN

  px => current_mesh%x ; fx => current_mesh%f
  cd => current_mesh%cd ; 


  c3 => current_mesh%c3 ; 
  c7 => current_mesh%c7 ;
  c9 => current_mesh%c9 ;

  ey => current_mesh%ey; cj => current_mesh%cj
  c1 = current_mesh%c1

  xe = current_mesh%xe ; ze = current_mesh%ze; ye = current_mesh%ye 
  tmpcomm = current_mesh%comm3d
  nuu=nu; if(current_mesh%mesh_num <= 1) nuu=nu_fine
  omg_ppe=1.00d0

  if(is_periodic_y) then
     jb = 0
     jmax = ye
  else
     jb = 1
     jmax = ye + 1
  endif
  js = 0

!
!   START ITERATIONS
!
  do  icount = 1,nuu

     do icolor=1,2

        rcolor=mod(icolor+2,2)+1     !2,1
        call color_swap(current_mesh,rcolor)

        do indx = 1 , current_mesh%nclr(icolor)

           i = current_mesh%iclr(indx,icolor)
           k = current_mesh%kclr(indx,icolor)

           if( jb == 1 ) then
              j= 0 
              delta_p = fx(i,k,0) + cj(i,k,1) * px(i+1,k,0) + cj(i,k,2) * px(i+1,k,1) &
                   + cj(i,k,4)*px(i-1,k,1)+ cj(i,k,5)*px(i-1,k,0)
              ay(j) = 0.0d0
              cy(j) = - cj(i,k,3)
              by(j) = c9(i,k,j)
              fy(j) = delta_p
           endif

           do j=jb,ye


              delta_p =  fx(i,k,j) + c1*(px(i+1,k,j) + px(i-1,k,j))+ ey(j) * ( &
                   cd(i,k,1) *(px(i+1,k,j+1)- px(i+1,k,j-1))&
                   +cd(i,k,2) *(px(i-1,k,j+1)- px(i-1,k,j-1)))
              ay(j) = - c7(i,k,j) * ey(j)
              cy(j) = - c3(i,k,j) * ey(j)
              by(j) = c9(i,k,j)
              fy(j) = delta_p
           enddo
           fy(jmax) = 0.0d0
           ay(jmax) = 0.0d0
           cy(jmax) = 0.0d0
           by(jmax) = 1.0d0

! --- THOMAS algorithm
           do j=js+1,jmax
              ay(j) = ay(j)/by(j-1)
              by(j) = by(j)-ay(j)*cy(j-1)
              fy(j) = fy(j)-ay(j)*fy(j-1)
           enddo

           fy(jmax) = fy(jmax)/by(jmax)
           px(i,k,jmax) = fy(jmax)

           do j=jmax-1,js,-1
              fy(j) = (fy(j)-cy(j)*fy(j+1))/by(j)
              px(i,k,j) = fy(j)
           enddo

        enddo   !do indx

     enddo   ! do icolor

  enddo   !END   do  icount = 1,nuu

!
  rcolor=2
  call color_swap(current_mesh,rcolor)

  ib = merge(1, 0, ALL(type_bc(4,1,:) == 0) .AND. is_proc_bnd(1) )

  if(ib == 1)then  !do not compute the first point
     wrk(0,:,:) = zero
  endif

  IF(resid_flag) then
!
!EVAL residual for transfering to lower levels
!
     myres = zero

     do i=ib,xe
        do k=0,ze


           if(jb == 1) then
              j=0
              wrk(i,k,0) = fx(i,k,0) + cj(i,k,1)*px(i+1,k,0) + cj(i,k,2)*px(i+1,k,1) &
                   + cj(i,k,3)*px(i,k,1)+ cj(i,k,4)*px(i-1,k,1)+ cj(i,k,5)*px(i-1,k,0) &
                   - c9(i,k,j)*px(i,k,0)

              if(abs(wrk(i,k,0)) > myres)then
                 myres = abs(wrk(i,k,0))
                 imax = i;kmax=k;jmax =j
              endif


           endif

           do j=jb,ye

              wrk(i,k,j) =  fx(i,k,j) + c1*(px(i+1,k,j) + px(i-1,k,j)) + ey(j) * ( &
                   cd(i,k,1) *(px(i+1,k,j+1)- px(i+1,k,j-1)) &
                   +cd(i,k,2) *(px(i-1,k,j+1)- px(i-1,k,j-1)) &
                   +c3(i,k,j) * px(i,k,j+1) + c7(i,k,j)*px(i,k,j-1) )&
                   -c9(i,k,j)*px(i,k,j)
              myres = max (myres, abs(wrk(i,k,j)))

           enddo
        enddo
     enddo
     if(current_mesh%mesh_num == 1) then
        CALL MPI_ALLREDUCE(myres,resmax,1,MPI_DOUBLE_PRECISION,MPI_MAX,tmpcomm,ierr)
        maxres(3) = resmax 
!!>        if(myid == 0) then
!!>           print*,icycle,current_mesh%mesh_num,'RESID_MG',resmax
!!>        endif
     endif

  ENDIF   !if (resid_flag)


!---------------------------------------------------------------------------------------------------
  RETURN
END SUBROUTINE gauss_relax
!***************************************************************************************************


!***************************************************************************************************
 SUBROUTINE restrict ( coarse_mesh )

!fine_mesh -->> coarse mesh     T_trans_r_calc (ALI)

  USE data_types
  USE GLOBAL_DATA 
  USE MYMPI

  IMPLICIT NONE

  include "parallel.h"

!---------------------------------------------------------------------------------------------------
  ! Dummy variables:

    TYPE(mesh_pointers) :: coarse_mesh
   
  ! Local Variables

    TYPE(mesh_pointers), POINTER :: fine_mesh
 
    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: px

    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: fx

    INTEGER :: i,j,k,id,kd,jb,itemp
    INTEGER :: ic,kc,jc,mn, col(6)
    INTEGER :: ncx,ncy,ncz,xe,ze,ye

    REAL(KIND=double) :: rm,rj,rp
!---------------------------------------------------------------------------------------------------
!
! This implementation uses full weghted residual restriction which is (little) 
! different from what is used in rocsolid
!
    fine_mesh => coarse_mesh%fine_mesh
    if (fine_mesh%blank) RETURN

    px =>  coarse_mesh%x ; fx =>  coarse_mesh%f
    ncx = coarse_mesh%xe + 1; ncy = coarse_mesh%ye + 1 ;ncz = coarse_mesh%ze + 1;
    xe = fine_mesh%xe ; ye = fine_mesh%ye; ze = fine_mesh%ze;
    mn = fine_mesh%mesh_num; id = id_mg(mn); kd = kd_mg(mn)


    col = 0;col(2) = xe;col(4) = ze; col(6) = ye;

    if ( iscoarse_xz(mn) ) &
      call parallel_swap(wrk,col, fine_mesh)
 
    itemp = 5   !Gross
    CALL DIRICHLET_BC(itemp,fine_mesh)
    if(is_periodic_y) then
       jb = 0
    else
       jb = 1
    endif

    if(coarse_mesh%blank) RETURN


    if ( iscoarse_xz(mn) .AND. iscoarse_y(mn) ) then
!
!   coarsening in x,y,z
!
      do kc=0,ncz-1
        k = kc+kc+kd
      do ic=0,ncx-1
        i = ic+ic+id

!
!       weight on j-1,j,j+1 z planes in rm,rk,rp
!
        DO j=jb-1,ye
           fy(j) = QUARTER*(wrk(i-1,k,j)+ wrk(i+1,k,j) + 2.* wrk(i,k,j))
        ENDDO

        if(jb == 1) &
             fx(ic,kc,0) = third*(fy(1) + two*fy(0) )

        DO jc = jb,ncy-1
           j = jc+jc
           fx(ic,kc,jc) = 0.25*(fy(j-1) + 2.0d0*fy(j) + fy(j+1))
        ENDDO

      end do
      end do       


    elseif ( iscoarse_xz(mn) ) then
     
!
!   coarsening in x,z but not y, 2D coasening
!
      do kc=0,ncz-1
        k = kc+kc+kd
      do jc=0,ncy-1
        j = jc
      do ic=0,ncx-1
        i = ic+ic+id

        rj = QUARTER*(wrk(i-1,k,j)+ wrk(i+1,k,j) + 2.* wrk(i,k,j))

        fx(ic,kc,jc) = rj
      
      end do
      end do
      end do

      elseif ( iscoarse_y(mn) ) then
!
!   coarsening in y but not xy, 1D coasening
!
      do kc=0,ncz-1
        k = kc
      do ic=0,ncx-1
        i = ic
      do jc=jb,ncy-1
        j = jc+jc
  
        fx(ic,kc,jc) = quarter*(wrk(i,k,j+1)+two*wrk(i,k,j)+wrk(i,k,j-1)) 
 
      end do
        if(jb == 1) &
             fx(ic,kc,0) =  third*(wrk(i,k,1)+two*wrk(i,k,0))
      end do
      end do


      else
    
      write(*,*) 'NO COARSENING, CHECK THE CALL!!'

      endif

!----------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE restrict
!**********************************************************************************


!**********************************************************************************
  SUBROUTINE prolong (coarse_mesh ,fine_mesh)

  USE data_types
  USE GLOBAL_DATA 
  USE MYMPI

  IMPLICIT NONE

  include "parallel.h"
!---------------------------------------------------------------------------------------------------
  ! Dummy variables:

    TYPE(mesh_pointers) :: coarse_mesh,fine_mesh

   
  ! Local Variables

    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: cx, fx    !fx is the solution on the fine mesh

    INTEGER :: i,j,k
    INTEGER :: ic,kc,jc,mn,id,kd
    INTEGER :: ncx,ncy,ncz
    INTEGER :: nfx,nfy,nfz

    REAL(KIND=double) :: rm,rj,rp
!---------------------------------------------------------------------------------------------------

    if(fine_mesh%blank) RETURN
!
!   initialize pointers and LOOP limits
!
    cx => coarse_mesh%x ; fx => fine_mesh%x  
    ncx = coarse_mesh%xe + 1; ncy = coarse_mesh%ye + 1 ;ncz = coarse_mesh%ze + 1;
    nfx = fine_mesh%xe + 1; nfy = fine_mesh%ye + 1 ;nfz = fine_mesh%ze + 1;
    mn = fine_mesh%mesh_num; id = id_mg(mn); kd = kd_mg(mn)

!
! ASSUMING pdims(2) > pdims(1) and nx=nz these two are equivalent
!
!    if(any(coarse_mesh%pskip) .AND. iscoarse_xz(mn)) &
    if(coarse_mesh%pskip(2) .AND. iscoarse_xz(mn)) &
       CALL fill_voids(coarse_mesh,fine_mesh)

    if ( iscoarse_xz(mn) .AND. iscoarse_y(mn) ) then
!
!   Interpolating in x,y,z
!   It is done in three step, 1-D interp, 2-D interp and finally 3-D interp
!
      do kc=-kd,ncz
        k = kc+kc+kd
      do jc=0,ncy
        j = jc+jc
      do ic=-id,ncx
        i = ic+ic+id
   
        wrk(i,k,j) = cx(ic,kc,jc)
        if(ic.lt.ncx)&
        wrk(i+1,k,j) = 0.5*(cx(ic,kc,jc)+cx(ic+1,kc,jc))
        if(kc.lt.ncz)&
        wrk(i,k+1,j) = 0.5*(cx(ic,kc,jc)+cx(ic,kc+1,jc))
        if(jc.lt.ncy)&
        wrk(i,k,j+1) = 0.5*(cx(ic,kc,jc)+cx(ic,kc,jc+1))
 
      enddo
      enddo
      enddo
     

      do kc=-kd,ncz
        k = kc+kc+kd
      do jc=0,ncy
        j = jc+jc
      do ic=-id,ncx
        i = ic+ic+id
    
        if(ic.lt.ncx.and.kc.lt.ncz)&
        wrk(i+1,k+1,j) = 0.5*(wrk(i+1,k,j)+wrk(i+1,k+2,j))
        if(ic.lt.ncx.and.jc.lt.ncy)&
        wrk(i+1,k,j+1) = 0.5*(wrk(i+1,k,j)+wrk(i+1,k,j+2))
        if(kc.lt.ncz.and.jc.lt.ncy)&
        wrk(i,k+1,j+1) = 0.5*(wrk(i,k+1,j)+wrk(i,k+1,j+2))

      enddo
      enddo
      enddo

      do kc=-kd,ncz-1
        k = kc+kc+kd
      do jc=0,ncy-1
        j = jc+jc
      do ic=-id,ncx-1
        i = ic+ic+id
    
        wrk(i+1,k+1,j+1) = 0.5*(wrk(i+1,k,j+1)+wrk(i+1,k+2,j+1))
      
      enddo
      enddo
      enddo

      elseif( iscoarse_xz(mn)) then
!
!   interpolating in x,z
!   It is done in two steps, 1-D interp, 2-D interp 
!
      do kc=-kd,ncz
        k = kc+kc+kd
      do jc=0,ncy
        j = jc
      do ic=-id,ncx
        i = ic+ic+id
   
        wrk(i,k,j) = cx(ic,kc,jc)
        if(ic.lt.ncx)&
        wrk(i+1,k,j) = 0.5*(cx(ic,kc,jc)+cx(ic+1,kc,jc))
        if(kc.lt.ncz)&
        wrk(i,k+1,j) = 0.5*(cx(ic,kc,jc)+cx(ic,kc+1,jc))
! 
      enddo
      enddo
      enddo
     

      do jc=0,ncy
        j = jc
      do kc=-kd,ncz-1
        k = kc+kc+kd
      do ic=-id,ncx-1
        i = ic+ic+id
    
        wrk(i+1,k+1,j) = 0.5*(wrk(i+1,k,j)+wrk(i+1,k+2,j))

      enddo
      enddo
      enddo
      !wrk(2*ncx-1,2*ncz-1,:)= 0.5*(cx(ncx,ncz-1,:)+cx(ncx-1,ncz,:))
         
      ELSEIF ( iscoarse_y(mn) ) then
!
!   interpolating in y only
!
      do kc=0,ncz
        k = kc
      do ic=0,ncx
        i = ic
      do jc=0,ncy-1
        j = jc+jc

        wrk(i,k,j) = cx(ic,kc,jc)
        wrk(i,k,j+1) = 0.5*(cx(ic,kc,jc)+cx(ic,kc,jc+1))
 
      enddo
      enddo
      enddo

      ELSE
    
      write(*,*) 'NO COARSENING, CHECK THE CALL!!'

      ENDIF

!
!   Inject the coarse grid correction into the fine mesh approximation
!
      do i=0,nfx-1
      do k=0,nfz-1
      do j=0,nfy-1

         fx(i,k,j) = fx(i,k,j)+wrk(i,k,j)

      enddo
      enddo
      enddo

!---------------------------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE prolong
!***************************************************************************************************


!***************************************************************************************************
  SUBROUTINE transfer ( coarse_mesh )

!This subroutine is the same as restrict (more or less). It is separate (for now)
! because at the start of the algorithm we just transfer the known term down

  USE data_types
  USE GLOBAL_DATA 

  IMPLICIT NONE
 
!---------------------------------------------------------------------------------------------------
  ! Dummy variables:

    TYPE(mesh_pointers) :: coarse_mesh
   
  ! Local Variables

    TYPE(mesh_pointers), POINTER :: fine_mesh

    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: cx
    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: cf
    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: fx

    INTEGER :: i,j,k,id,kd
    INTEGER :: ic,kc,jc,mn
    INTEGER :: ncx,ncy,ncz,incy,incxz
    INTEGER :: xe,ze,ye
!---------------------------------------------------------------------------------------------------
!
!   RETURN IMMEDIATELY IF THE PROC HAS NO POINTS
!
    if(coarse_mesh%blank) RETURN

    fine_mesh => coarse_mesh%fine_mesh
    cx => coarse_mesh%x ; cf =>  coarse_mesh%f
    fx =>   fine_mesh%x
    ncx = coarse_mesh%xe + 1; ncy = coarse_mesh%ye + 1 ;ncz = coarse_mesh%ze + 1;
    mn = fine_mesh%mesh_num; id = id_mg(mn); kd = kd_mg(mn)

!
    if ( iscoarse_xz(mn) .AND. iscoarse_y(mn) ) then
       incxz=2
       incy=2
    else if (iscoarse_xz(mn)) then
       incxz=2
       incy=1
    else if (iscoarse_y(mn)) then
       incxz=1
       id=0;kd=0
       incy=2
    else
       write(*,*) 'ERROR in TRANSFER, mesh # = ',mn
    endif
!
!   coarsening in x,y,z
!
    do kc=0,ncz-1
      k = kc*incxz+kd
    do ic=0,ncx-1
      i = ic*incxz+id
    do jc=0,ncy-1
      j = jc*incy
!
! JUST inject
!
      cf(ic,kc,jc) = wrk(i,k,j)
      cx(ic,kc,jc) =  fx(i,k,j)

    end do
    end do
    end do       

!
!   This extra communication is necessary only on the corasest grid
!   when iex is odd
!
    if (mn == num_mg_meshes-1 .AND. isoddcoarse)&
       call color_swap(coarse_mesh,1)

!!!    cx = zero
      
!---------------------------------------------------------------------------------------------------
    RETURN
 END SUBROUTINE transfer
!***************************************************************************************************


!***************************************************************************************************
 SUBROUTINE output_results( current_mesh )

  USE data_types
  
  USE mg_solver

  USE GLOBAL_DATA 

  IMPLICIT NONE

!---------------------------------------------------------------------------------------------------
  ! Dummy variables:

    TYPE(mesh_pointers) :: current_mesh

   
  ! Local Variables

    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: px

    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: fx

    INTEGER :: i,j,k
    INTEGER :: ic,kc,jc,mn
    INTEGER :: ncx,ncy,ncz,incy
!---------------------------------------------------------------------------------------------------

    px =>  current_mesh%x ; fx =>  current_mesh%f
    ncx = current_mesh%xe + 1; ncy = current_mesh%ye + 1 ;ncz = current_mesh%ze + 1;
    mn = current_mesh%mesh_num

    do jc=0,ncy-1
    do kc=0,ncz-1
    do ic=0,ncx-1

       if(myid == 0) &
          write(*,*)myid,ic,kc,jc,'OUT',px(ic,kc,jc),fx(ic,kc,jc),ncy

    end do
    end do
    end do       
   
!---------------------------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE output_results
!***************************************************************************************************


!***************************************************************************************************
  SUBROUTINE SOLVE_DIRECT( current_mesh )

    USE data_types  
    USE mg_solver
    USE GLOBAL_DATA 
    USE LU_VARS

    IMPLICIT NONE

!---------------------------------------------------------------------------------------------------
! Dummy variables:

    TYPE(mesh_pointers) :: current_mesh


! Local Variables

    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: px

    REAL(KIND=double), DIMENSION(:,:,:), POINTER :: fx

    INTEGER :: i,j,k, ib, jb,itemp
    INTEGER :: j_LU, mn, elcount
    INTEGER :: ncx,ncy,ncz, istart, iend
    REAL(KIND=double) :: maxres_direct
!---------------------------------------------------------------------------------------------------

    px =>  current_mesh%x ; fx =>  current_mesh%f
    ncx = current_mesh%xe + 1; ncy = current_mesh%ye + 1 ;ncz = current_mesh%ze + 1;
    mn = current_mesh%mesh_num

!--- ALL gather the RHS on the coarsest grid
    ib = merge(1, 0, ALL(type_bc(4,1,:) == 0) .AND. is_proc_bnd(1))
    k = 0
    elcount = 0 
    do i=0,ncx-1
       do j=0,ncy-1
          elcount = elcount+1  !starts from 1
          if(i< ib) then
             buff3(elcount) = zero
          else
             buff3(elcount) =  fx(i,k,j)
          endif
       enddo
    enddo

    call MPI_ALLGATHERV(buff3(1),elcount,MPI_DOUBLE_PRECISION, &
       &     LU_F(1),LU_size,LU_disp,MPI_DOUBLE_PRECISION, &
       &     comm3d,ierr)

!solves and returns the solution vector in the same buffer the rhs was provided in
    call LU_SPARSE_SOLVE 

    k = 0
    j_LU = max( (ib_mg(mn)-1)*LU_ny , 0 )
    istart = merge(0,-1,is_proc_bnd(1))
    iend = merge(ncx-1,ncx,is_proc_bnd(2))
    do i = istart, iend
       do j = 0, ncy-1
          j_LU = j_LU+1
          px(i,k,j) = LU_F(j_LU)
       enddo
    enddo

!--- update dirichlet BC for pressure
    itemp = 4  !Gross
    call DIRICHLET_BC(itemp,current_mesh)
    
    if(issymmetric) then
       IF(is_proc_bnd(1)) px(-1,:,:)   = px(+1,:,:)
       IF(is_proc_bnd(2)) px(current_mesh%xe+1,:,:) = px(current_mesh%xe,:,:)
    else
       IF(is_proc_bnd(1)) px(-1,:,:)   = px(current_mesh%xe,:,:)
       IF(is_proc_bnd(2)) px(current_mesh%xe+1,:,:) = px(0,:,:)
    endif
    
!
!evaluate residual this is just a check
!
    jb = 1
    maxres_direct = zero
    residual_direct = zero
    do i=ib,ncx-1
       if(jb == 1) then
          j=0
          wrk(i,k,0) = fx(i,k,0) + current_mesh%cj(i,k,1)*px(i+1,k,0) &
               + current_mesh%cj(i,k,2)*px(i+1,k,1) &
               + current_mesh%cj(i,k,3)*px(i,k,1)   &
               + current_mesh%cj(i,k,4)*px(i-1,k,1) &
               + current_mesh%cj(i,k,5)*px(i-1,k,0) &
               - current_mesh%c9(i,k,j)*px(i,k,0)

          maxres_direct = max(maxres_direct,abs( wrk(i,k,j)))
       endif
       do j=jb,ncy-1
          wrk(i,k,j) =  fx(i,k,j) + current_mesh%c1*(px(i+1,k,j) + px(i-1,k,j)) + &
               current_mesh%ey(j) * ( &
                   current_mesh%cd(i,k,1) *(px(i+1,k,j+1)- px(i+1,k,j-1)) &
                  +current_mesh%cd(i,k,2) *(px(i-1,k,j+1)- px(i-1,k,j-1)) &
                  +current_mesh%c3(i,k,j) * px(i,k,j+1) &
                  +current_mesh%c7(i,k,j) * px(i,k,j-1))&
                  -current_mesh%c9(i,k,j)*px(i,k,j)
          maxres_direct = max(maxres_direct,abs( wrk(i,k,j)))
       enddo
    enddo
    CALL MPI_ALLREDUCE(maxres_direct,residual_direct,1,&
         MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
!---------------------------------------------------------------------------------------------------
    RETURN
  END SUBROUTINE SOLVE_DIRECT
!***************************************************************************************************

