!*********************************************************************
 SUBROUTINE PROJECTION

   USE GLOBAL_DATA
   USE mg_solver
   USE MYMPI
   
   IMPLICIT NONE

   include "parallel.h"
   
!-----------------------------------------------------------------------
!  Local variables
   INTEGER :: i, k, j, eqn, col(7), error
   INTEGER :: xs,zs,xe,ze,ys,ye
   REAL*8  :: gx,gy,gz,dxa,dya,dza
!-----------------------------------------------------------------------

   xs = 0
   xe = finest_mesh%xe
   zs = 0
   ze = 0
   ys = 1
   ye = finest_mesh%ye + 1
!
!  the array col contains information about the color of the pints to be updated
!  col(1:2), the dimension to be updated, col(3:4), and how many point in the
!  j direction need to be updated col(5:6)
!
   col(1)=1; col(2)=2; col(3)=1; col(4) = ndim; 
   col(5)=drange(5); col(6)=drange(6); col(7) = 1

   dxa = 1.0d0/dx
   dza = 1.0d0/dz

   CALL VELOCITY(4)

   CALL UPDATE_CONDITIONS_BC  !uses vvel,uvel
!
!  evaluate the rhs of the pressure Poisson equation (PPE)
!
   do j=ys,ye
      dya = detady(j)/dy
   do k=zs,ze
   do i=xs,xe
      gx = (uvel(i,k,j)-uvel(i-1,k,j))*dxa
      gy = (vvel(i,k,j)-vvel(i,k,j-1))*dya
      divt(i,k,j)= gx+gy  + rate(i,k,j,4)
    enddo
    enddo
    enddo
    p = zero
!
    divt(:,:,ny) = divt(:,:,ny-1)  !This value should never be used
    divt(:,:,0) = zero

    CALL get_coeff_ppe

!
!  SHIFT the vectors to be used in the MGsolver (easy to coarsen)
!
    divt=cshift(divt,1,3)
    p=cshift(p,1,3)

    finest_mesh%x => p
    finest_mesh%f => divt

    CALL fmg_solve   !CALL THE MULTIGRID SOLVER

    divt=cshift(divt,-1,3)
    p=cshift(p,-1,3)

    CALL UPDATE_Q

    CALL parallel_swap4(col,finest_mesh)

    do i = 1,ndim
       CALL  DIRICHLET_BC(i,finest_mesh)
    enddo
!-----------------------------------------------------------------------
   RETURN
 END SUBROUTINE PROJECTION
!*************************************************************************


!*************************************************************************
 SUBROUTINE UPDATE_Q

   USE GLOBAL_DATA
   
   IMPLICIT NONE

!-------------------------------------------------------------------------   
! Local variables
   INTEGER :: i, k, j, xs,zs,xe,ze,ys,ye
   REAL*8  :: gx,gy,gz,dxa,dya,dza,dya4
   REAL*8  :: pressure,poo,coe_rho,rhogas
!-------------------------------------------------------------------------   
   
   xs = drange(1)
   xe = drange(2)
   zs = drange(3)
   ze = drange(4)
   ys = drange(5)+1
   ye = drange(6)-1

   dxa = 1.0d0/dx
   dza = 1.0d0/dz    
 
   poo = PRESSURE(tcyc)
   coe_rho = dim_fact*poo

   p(:,:,ny) = zero
   if(.NOT. is_periodic_y) then
      p(:,:,0) = two*p(:,:,1)-p(:,:,2)
   endif
 
   do j =0,ny
   do i=-1,finest_mesh%xe+1
   do k=-1,finest_mesh%ze+1
      pold(i,k,j) = pold(i,k,j) + p(i,k,j)
   enddo
   enddo
   enddo

!--- UPDATE VELOCITY AND PRESSURE GRADIENT
   do j=ys,nyv(1)
      dya4 = quarter*detady(j)/dy
   do i=xs,nxv(1)
   do k=zs,ze
      rhogas = two*coe_rho/(f(i,k,j,1)+f(i+1,k,j,1))
      gy = (p(i,k,j+1)+p(i+1,k,j+1)-p(i,k,j-1)-p(i+1,k,j-1))*dya4
      gx = (p(i+1,k,j)-p(i,k,j))*dxa
      q(i,k,j,1) = q(i,k,j,1) + (dphidxa(i,k)*gy - gx)/rhogas
      rate(i,k,j,1) = rate(i,k,j,1) + (dphidxa(i,k)*gy - gx)/rhogas
   enddo
   enddo
   enddo

   do j=ys,nyv(3)
      dya = detadya(j)/dy
   do i=xs,nxv(3)
   do k=zs,ze
      rhogas = two*coe_rho/(f(i,k,j+1,1)+f(i,k,j,1))
      gy = (p(i,k,j+1)-p(i,k,j))*dya
      q(i,k,j,3) = q(i,k,j,3) - gy/rhogas
      rate(i,k,j,3) = rate(i,k,j,3) - gy/rhogas
   enddo
   enddo
   enddo


   CALL MGTEST
   if(.FALSE.) &
   CALL RHSPTEST

!-------------------------------------------------------------------------
   RETURN
 END SUBROUTINE UPDATE_Q
!*************************************************************************


!*************************************************************************
 SUBROUTINE get_coeff_ppe

   USE data_types
   USE global_data
   USE LUSOLVER

   IMPLICIT NONE

!-------------------------------------------------------------------------
! Local variables:
   TYPE(mesh_pointers), POINTER :: current_mesh
   REAL(KIND=double), DIMENSION(:,:,:), POINTER :: cd,c3,c7,cj
   REAL(KIND=double), DIMENSION(:), POINTER :: ey,eya

   INTEGER :: xe,ze,ye,i,k,j,jb,mn,mnx,mnz,m,itemp
   INTEGER :: i_LU,j_LU

   REAL(KIND=double) :: dxf,dzf,dyf,dxc,dzc,dyc,dxa,dza,dya
   REAL(KIND=double) :: dxa2,dza2,dy2a,dxdy4,dzdy4
   REAL(KIND=double) :: term,c379
!--------------------------------------------------------------------------

  
   current_mesh => finest_mesh

   xe = current_mesh%xe ; ze = current_mesh%ze; ye = current_mesh%ye
   dxf = current_mesh%dx ; dzf = current_mesh%dz; dyf = current_mesh%dy
   dxa = one/current_mesh%dx ; dza = one/current_mesh%dz; dya = one/current_mesh%dy

   cd => current_mesh%cd  
   c3 => current_mesh%c3  
   c7 => current_mesh%c7 
   cj => current_mesh%cj
   ey => current_mesh%ey 
   eya => current_mesh%eya

   dxa2=dxa*dxa;dza2=dza*dza;dy2a = dya*dya
   dxdy4 = dxa*dya/4.0d0;dzdy4 = dza*dya/4.0d0

   if(is_periodic_y) then
      jb = 0
   else
      jb = 1
   endif

!
!   b1=c1
!
   current_mesh%c1 = -dxa2

   do i=0,xe
      do k=0,ze
         cd(i,k,1) =  (dphidxa(i,k)+dphidx(i,k))*dxdy4
         cd(i,k,2) = -(dphidxa(i-1,k)+dphidx(i,k))*dxdy4
         cd(i,k,3) =  (dphidza(i,k)+dphidz(i,k))*dzdy4
         cd(i,k,4) = -(dphidza(i,k-1)+dphidz(i,k))*dzdy4
         term = (one+dphidx(i,k)**2+dphidz(i,k)**2)*dy2a
         c379 = (dphidxa(i,k)-dphidxa(i-1,k))*dxdy4 + (dphidza(i,k)-dphidza(i,k-1))*dzdy4

         cj(i,k,1) = -dxa2+(dphidx(i,k)-two*dphidxa(i,k))*dxdy4*ey(0)    
         cj(i,k,2) = (dphidx(i,k)+two*dphidxa(i,k))*dxdy4*ey(0)    
         cj(i,k,3) = ey(0) * ( two*c379-term*eya(0))
         cj(i,k,4) = -(dphidx(i,k)+two*dphidxa(i-1,k))*dxdy4*ey(0)
         cj(i,k,5) = -dxa2-(dphidx(i,k)-two*dphidxa(i-1,k))*dxdy4*ey(0)

         cj(i,k,6) = -dza2+(dphidz(i,k)-two*dphidza(i,k))*dzdy4*ey(0)
         cj(i,k,7) = (dphidz(i,k)+two*dphidza(i,k))*dxdy4*ey(0)
         cj(i,k,8) = -(dphidz(i,k)+two*dphidza(i,k-1))*dzdy4*ey(0)
         cj(i,k,9) = -dza2-(dphidz(i,k)-two*dphidza(i,k-1))*dzdy4*ey(0)


         if(jb ==1) then
            c3(i,k,0) = ( two*c379-term*eya(0))
            c7(i,k,0) = zero
         endif

         do j=jb,ye
            c3(i,k,j)=  ( c379-term*eya(j))
            c7(i,k,j)=  (-c379-term*eya(j-1))
         enddo
      enddo
   enddo

!
!  Discretize the PPE on each mesh level
!
   DO

      IF ( .NOT. ASSOCIATED(current_mesh%coarse_mesh) ) EXIT
      current_mesh => current_mesh%coarse_mesh

      xe = current_mesh%xe ; ze = current_mesh%ze; ye = current_mesh%ye
      dxc = current_mesh%dx ; dzc = current_mesh%dz; dyc = current_mesh%dy
      mnx = anint(dxc/dxf); mnz = anint(dzc/dzf);
      m= anint( log10 (dble(mnx)) / log10 (2.0))

      cd => current_mesh%cd ; c3 => current_mesh%c3 ; c7 => current_mesh%c7 ;
      ey => current_mesh%ey ; cj => current_mesh%cj; eya => current_mesh%eya

      dxa2 = 1.0d0/(dxc*dxc); dza2 = 1.0d0/(dzc*dzc); dy2a = 1.0d0/(dyc*dyc)
      dxdy4 = 1.0d0/(dxc*dyc*4.0d0); dzdy4 = 1.0d0/(dzc*dyc*4.0d0)

!
!     b1=c1
!
      current_mesh%c1 = -dxa2

      do i=0,xe
         do k=0,ze


            if( mnx <= 1 ) then               ! didnt coarse in xz direction 

               term = (1.0+dphidx(i,k)**2+dphidz(i,k)**2)*dy2a
               cd(i,k,1) = (dphidxa(i,k)+dphidx(i,k))*dxdy4
               cd(i,k,2) =-(dphidxa(i-1,k)+dphidx(i,k))*dxdy4
               cd(i,k,3) = (dphidza(i,k)+dphidz(i,k))*dzdy4
               cd(i,k,4) =-(dphidza(i,k-1)+dphidz(i,k))*dzdy4
               c379 = (dphidxa(i,k)-dphidxa(i-1,k))*dxdy4 + (dphidza(i,k)-dphidza(i,k-1))*dzdy4

               cj(i,k,1) = -dxa2+(dphidx(i,k)-two*dphidxa(i,k))*dxdy4*ey(0)
               cj(i,k,2) = (dphidx(i,k)+two*dphidxa(i,k))*dxdy4*ey(0)    
               cj(i,k,3) = ey(0) * ( two*c379-term*eya(0))
               cj(i,k,4) = -(dphidx(i,k)+two*dphidxa(i-1,k))*dxdy4*ey(0)    
               cj(i,k,5) = -dxa2-(dphidx(i,k)-two*dphidxa(i-1,k))*dxdy4*ey(0)    

               cj(i,k,6) = -dza2+(dphidz(i,k)-two*dphidza(i,k))*dzdy4*ey(0)    
               cj(i,k,7) =  (dphidz(i,k)+two*dphidza(i,k))*dxdy4*ey(0)    
               cj(i,k,8) = -(dphidz(i,k)+two*dphidza(i,k-1))*dzdy4*ey(0)    
               cj(i,k,9) = -dza2-(dphidz(i,k)-two*dphidza(i,k-1))*dzdy4*ey(0)

            else

               term = (1.0+dphidx_mg(i,k,m)**2+dphidz_mg(i,k,m)**2)*dy2a
               cd(i,k,1) = (dphidxa_mg(i,k,m)+dphidx_mg(i,k,m))*dxdy4
               cd(i,k,2) =-(dphidxa_mg(i-1,k,m)+dphidx_mg(i,k,m))*dxdy4
               cd(i,k,3) = (dphidza_mg(i,k,m)+dphidz_mg(i,k,m))*dzdy4
               cd(i,k,4) =-(dphidza_mg(i,k-1,m)+dphidz_mg(i,k,m))*dzdy4
               c379 = (dphidxa_mg(i,k,m)-dphidxa_mg(i-1,k,m))*dxdy4 +&
                    (dphidza_mg(i,k,m)-dphidza_mg(i,k-1,m))*dzdy4

               cj(i,k,1) = -dxa2+(dphidx_mg(i,k,m)-two*dphidxa_mg(i,k,m))*dxdy4*ey(0)
               cj(i,k,2) =  (dphidx_mg(i,k,m)+two*dphidxa_mg(i,k,m))*dxdy4*ey(0)
               cj(i,k,3) = ey(0) * ( two*c379-term*eya(0))
               cj(i,k,4) = -(dphidx_mg(i,k,m)+two*dphidxa_mg(i-1,k,m))*dxdy4*ey(0)
               cj(i,k,5) = -dxa2-(dphidx_mg(i,k,m)-two*dphidxa_mg(i-1,k,m))*dxdy4*ey(0)

               cj(i,k,6) = -dza2+(dphidz_mg(i,k,m)-two*dphidza_mg(i,k,m))*dzdy4*ey(0)
               cj(i,k,7) =  (dphidz_mg(i,k,m)+two*dphidza_mg(i,k,m))*dxdy4*ey(0)
               cj(i,k,8) = -(dphidz_mg(i,k,m)+two*dphidza_mg(i,k-1,m))*dzdy4*ey(0)
               cj(i,k,9) = -dza2-(dphidz_mg(i,k,m)-two*dphidza_mg(i,k-1,m))*dzdy4*ey(0)

          

            endif
            if(jb ==1) then
               c3(i,k,0) = (two*c379-term*eya(0))
               c7(i,k,0) = zero
            endif

            do j=jb,ye
               c3(i,k,j)= ( c379-term*eya(j))
               c7(i,k,j)= (-c379-term*eya(j-1))
            enddo

         enddo
      enddo

   END DO   


   current_mesh => finest_mesh
   MESHES : DO

      mn = current_mesh%mesh_num
      xe = current_mesh%xe ; ze = current_mesh%ze; ye = current_mesh%ye
      
      do i=0,xe
         do k=0,ze
            do j = 0,ye
               current_mesh%c9(i,k,j) = two*current_mesh%c1 + &
                    current_mesh%ey(j)*(current_mesh%c3(i,k,j)+current_mesh%c7(i,k,j))
            enddo
         enddo
      enddo

!      itemp = -4  !Gross
      CALL DIRICHLET_BC(itemp,current_mesh)  !4 stands for pressure

      IF ( .NOT. ASSOCIATED(current_mesh%coarse_mesh) ) EXIT
      current_mesh => current_mesh%coarse_mesh


   END DO MESHES

   
!
!NOW Perform a sparse LU decomposition on the coarsest mesh 
!
   IF((mod(ncyc,LU_frequency) == 0 .OR. .NOT. done_LU_decomposition &
        .OR. residual_direct > 10d0) .AND. do_solve_direct ) THEN
      
      print*,ncyc,'PERFORMING LU DECOMPOSITION',LU_frequency,done_LU_decomposition,&
           'NUMBER OF DIAGONALS',ndiag
      
      done_LU_decomposition = .TRUE.

      if(num_mg_meshes == 1) then
         current_mesh => finest_mesh
      else
         current_mesh => coarsest_mesh
      endif

      mn = current_mesh%mesh_num
      xe = current_mesh%xe ; ze = current_mesh%ze; ye = current_mesh%ye

      cm(:,:,ndiag-1) = - current_mesh%c1
      cm(:,:,ndiag) = - current_mesh%c1
      do j=0,ye
         do i = 0, xe;  k=0
            i_LU= i+1;j_LU=j+1

            cm(i_LU,j_LU,1) = current_mesh%c9(i,k,j)

            cm(i_LU,j_LU,2) =  - current_mesh%c7(i,k,j) * current_mesh%ey(j)
            cm(i_LU,j_LU,3) =  - current_mesh%c3(i,k,j) * current_mesh%ey(j)

            if(ndiag == 9) then
               if(j==0)then
                  cm(i_LU,j_LU,4) = zero
                  cm(i_LU,j_LU,5) = -current_mesh%cj(i,k,4)
                  cm(i_LU,j_LU,6) = zero
                  cm(i_LU,j_LU,7) = -current_mesh%cj(i,k,2)
                  cm(i_LU,j_LU,8) = -current_mesh%cj(i,k,5)
                  cm(i_LU,j_LU,9) = -current_mesh%cj(i,k,1)
               else
                  cm(i_LU,j_LU,4) = current_mesh%ey(j)*current_mesh%cd(i,k,2)
                  cm(i_LU,j_LU,5) = -current_mesh%ey(j)*current_mesh%cd(i,k,2)
                  cm(i_LU,j_LU,6) = current_mesh%ey(j)*current_mesh%cd(i,k,1)
                  cm(i_LU,j_LU,7) = -current_mesh%ey(j)*current_mesh%cd(i,k,1)
               endif
            endif                 
                  
         enddo
      enddo


      CALL LU_GATHER(xe+1,ye+1)  !->proc 0

!
! this because I would fix the solution at i=0
!
      if(ALL(type_bc(1,1,:) == 1)) then
         cm(1,:,2:ndiag) = zero
         cm(1,:,1) = one
      endif
      
       
      CALL LU_DECOMPOSE

   ENDIF

!---------------------------------------------------------------------
   RETURN
 END SUBROUTINE get_coeff_ppe
!**********************************************************************


!**********************************************************************
SUBROUTINE RHSPTEST

   USE GLOBAL_DATA
   USE data_types
   USE MYMPI
   
   IMPLICIT NONE

   include "parallel.h"
   
!-----------------------------------------------------------------------
!  Local variables
   INTEGER :: i, k, j, eqn, col(7)
   INTEGER :: xs,zs,xe,ze,ys,ye
   REAL*8  :: gx,gy,gz,dxa,dya,dza
   REAL*8  :: pressure,poo,coe,maxq(2),allmaxq(2),rhogas
!-----------------------------------------------------------------------
   
   xs = drange(1)
   xe = drange(2)
   zs = drange(3)
   ze = drange(4)
   ys = drange(5)+1
   ye = drange(6)-1
   col(1)=1; col(2)=2; col(3)=1; col(4) = ndim; 
   col(5)=drange(5); col(6)=drange(6); col(7) = 0

   dxa = 1.0d0/dx
   dza = 1.0d0/dz

   poo = PRESSURE(tcyc)
   coe = two*dim_fact*poo

   do j=ys,ye
      dya = detadya(j)/dy
   do i=xs,xe
   do k=zs,ze
      rhogas = coe/(f(i,k,j+1,1)+f(i,k,j,1))
      gy = (p(i,k,j+1)-p(i,k,j))*dya
      gx = 0.25*(p(i+1,k,j+1)-p(i-1,k,j+1)+p(i+1,k,j)-p(i-1,k,j))*dxa
      gz = 0.25*(p(i,k+1,j+1)-p(i,k-1,j+1)+p(i,k+1,j)-p(i,k-1,j))*dza
      uvel(i,k,j) = uvel(i,k,j) + (dphidx(i,k)*gy - gx)/rhogas
      wvel(i,k,j) = wvel(i,k,j) + (dphidz(i,k)*gy - gz)/rhogas
   enddo
   enddo
   enddo


   CALL parallel_swap4(col,finest_mesh)

   CALL VELOCITY(6)

!
!  TEST the projection
!
   wrk = zero
   do j=ys,ye
      dya = detady(j)/dy
   do k=zs,ze
   do i=xs,xe
      gx = (uvel(i,k,j)-uvel(i-1,k,j))*dxa
      gz = (wvel(i,k,j)-wvel(i,k-1,j))*dza
      gy = (vvel(i,k,j)-vvel(i,k,j-1))*dya
      wrk(i,k,j)= gx+gy+gz 
   enddo
   enddo
   enddo
!
   slambda = zero
   do j=0,ny
   do k=zs,ze
   do i=xs,xe
     slambda(j) = slambda(j) + vvel(i,k,j)
   enddo
   enddo
   enddo

   maxq(1) = maxval(abs(vvel(0:xe,0:ze,0:ye))); maxq(2) = maxval(abs(wrk))
   CALL MPI_ALLREDUCE(maxq,allmaxq,2,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
   if(myid == 0)&
   write(*,*) ncyc,'RHPTEST,maxvvel,maxresid', allmaxq
   CALL MPI_ALLREDUCE(slambda(0),rlambda(0),50,MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
   if(myid == 0) then
     do j=0,ny
        write(*,*) j,'VELSUM',rlambda(j)
     enddo
   endif
!
!-----------------------------------------------------------------------
  RETURN
END SUBROUTINE RHSPTEST
!**********************************************************************


!**********************************************************************
SUBROUTINE MGTEST

   USE data_types
   USE GLOBAL_DATA
   USE MYMPI
   
   IMPLICIT NONE

!-----------------------------------------------------------------------
!  Local variables
   INTEGER :: i, k, j, eqn
   INTEGER :: xs,zs,xe,ze,ys,ye
   REAL*8  :: gx,gy,gz,dxa,dya,dza, myerr_mg,err_mg
   REAL(KIND=double) :: dxa2,dza2,dy2a,dxdy4,dzdy4
   REAL(KIND=double) :: cd(4),cj(9),c3,c7,ey,c1,c9,term,c379
   REAL*8  :: pressure,poo,coe,rhogas
!-----------------------------------------------------------------------
   
   xs = drange(1)
   xe = drange(2)
   zs = drange(3)
   ze = drange(4)
   ys = drange(5)+1
   ye = drange(6)-1

   dxa = 1.0d0/dx
   dza = 1.0d0/dz
   dya = 1.0d0/dy
   dxa2=dxa*dxa; dza2=dza*dza; dy2a = dya*dya
   dxdy4 = dxa*dya/4.0d0; dzdy4 = dza*dya/4.0d0

   poo = PRESSURE(tcyc)
   coe = two*dim_fact*poo

   c1 = -dxa2

!
    wrk= zero

    do i=0,xe
    do k=0,ze
       j=1

      term = (1.0+dphidx(i,k)**2+dphidz(i,k)**2)*dy2a
      c379 = (dphidxa(i,k)-dphidxa(i-1,k))*dxdy4 + (dphidza(i,k)-dphidza(i,k-1))*dzdy4
      cd(1) = (dphidxa(i,k)+dphidx(i,k))*dxdy4
      cd(2) =-(dphidxa(i-1,k)+dphidx(i,k))*dxdy4
      cd(3) = (dphidza(i,k)+dphidz(i,k))*dzdy4
      cd(4) =-(dphidza(i,k-1)+dphidz(i,k))*dzdy4

      cj(1) = -dxa2+(dphidx(i,k)-two*dphidxa(i,k))*dxdy4*detady(1)
      cj(2) = (dphidx(i,k)+two*dphidxa(i,k))*dxdy4*detady(1)
      cj(3) = detady(1) * ( two*c379-term*detadya(1))
      cj(4) = -(dphidx(i,k)+two*dphidxa(i-1,k))*dxdy4*detady(1)
      cj(5) = -dxa2-(dphidx(i,k)-two*dphidxa(i-1,k))*dxdy4*detady(1)

      cj(6) = -dza2+(dphidz(i,k)-two*dphidza(i,k))*dzdy4*detady(1)
      cj(7) = (dphidz(i,k)+two*dphidza(i,k))*dxdy4*detady(1)
      cj(8) = -(dphidz(i,k)+two*dphidza(i,k-1))*dzdy4*detady(1)
      cj(9) = -dza2-(dphidz(i,k)-two*dphidza(i,k-1))*dzdy4*detady(1)

       c9 = four*c1+ cj(3)
       wrk(i,k,1) = (divt(i,k,1) + cj(1) * p(i+1,k,1) + cj(2) * p(i+1,k,2)+&
                   cj(3)*p(i,k,2)+ cj(4)*p(i-1,k,2)+ cj(5)*p(i-1,k,1)+&
                   cj(6)*p(i,k+1,1)+ cj(7)*p(i,k+1,2)+ cj(8)*p(i,k-1,2)+&
                   cj(9)*p(i,k-1,1)) - c9*p(i,k,1)

    do j=2,ye

      c3 =  ( c379-term*detadya(j))
      c7 =  (-c379-term*detadya(j-1))
      ey = detady(j)
     
      c9 = four*c1+ ey*(c3+c7)


      wrk(i,k,j)= (divt(i,k,j) + c1*(p(i+1,k,j) + p(i-1,k,j)+&
                    p(i,k+1,j) + p(i,k-1,j))+ ey * ( &
                    cd(1) *(p(i+1,k,j+1)- p(i+1,k,j-1))&
                   +cd(2) *(p(i-1,k,j+1)- p(i-1,k,j-1))&
                   +cd(3) *(p(i,k+1,j+1)- p(i,k+1,j-1))&
                   +cd(4) *(p(i,k-1,j+1)- p(i,k-1,j-1))&
                   +c3    * p(i,k,j+1) + c7*p(i,k,j-1) ))&
                    -c9 * p(i,k,j)

    enddo
    enddo
    enddo

    myerr_mg= maxval(abs(wrk))
    CALL MPI_ALLREDUCE(myerr_mg,err_mg,1,MPI_DOUBLE_PRECISION,MPI_MAX,comm3d,ierr)
!
!-----------------------------------------------------------------------
  RETURN
END SUBROUTINE MGTEST
!**********************************************************************

!**********************************************************************
 SUBROUTINE INIT_PRESS

   USE GLOBAL_DATA
   
   IMPLICIT NONE

!-------------------------------------------------------------------------   
   ! Local variables
   INTEGER :: i, k, j, xs,zs,xe,ze,ys,ye
   REAL*8  :: dya,pressure,poo,coe,rhogas,rhogas0
!-------------------------------------------------------------------------   
   
   xs = drange(1)
   xe = drange(2)
   zs = drange(3)
   ze = drange(4)
   ys = drange(5)+1
   ye = drange(6)-1

   p(:,:,0) = 2.0*p(:,:,1)-p(:,:,2)
   p(:,:,ny) = zero

   poo = PRESSURE(tcyc)
   coe = two*dim_fact*poo

   do j=8,1,-1
      dya = detadya(j)/dy
   do i=xs-1,xe+1
   do k=zs-1,ze+1
      rhogas0 = coe/(f(i,k,1,1)+f(i,k,0,1))
      rhogas = coe/(f(i,k,j+1,1)+f(i,k,j,1))
      p(i,k,j) = p(i,k,j+1) +  (rhogas0*q(i,k,0,3) - rhogas*q(i,k,j,3))/dya
!      p(i,k,j) = p(i,k,j+1) +  (vvel(i,k,0) - vvel(i,k,j))/dya 
   enddo
   enddo
   enddo

!   CALL MGTEST
!-------------------------------------------------------------------------
   RETURN
 END SUBROUTINE INIT_PRESS
!******************************************************************************

! ********************************************************************
 SUBROUTINE PROJECTION_SOURCE

   USE GLOBAL_DATA
   USE MYMPI

   IMPLICIT NONE

!--------------------------------------------------------------------------------
! Local variables
   INTEGER ::  i, j, jp, k, eqn, ys,ye,xs,xe,zs,ze
   REAL*8  :: dx1, dx2, dy1, dy2, dz1, dz2, dx1a, dx2a, dy1a, dy2a, dz1a, dz2a,term
   REAL*8  :: dxy1, dzy1,gx,gy,gxx,gyy,gxy,gz,gzz,gzy
   REAL*8  :: conv,diff
!----------------------------------------------------------------------------------

!
! Initialize variables
!
   dx1 = dx
   dx2 = dx*dx
   dy1 = dy
   dy2 = dy*dy
   dz1 = dz
   dz2 = dz*dz

   dx1a = 1.0d0/(2.0d0*dx1)
   dx2a = 1.0d0/(1.0d0*dx2)
   dy1a = 1.0d0/(2.0d0*dy1)
   dy2a = 1.0d0/(1.0d0*dy2)
   dz1a = 1.0d0/(2.0d0*dz1)
   dz2a = 1.0d0/(1.0d0*dz2)
   dxy1 = 1.0d0/(4.0d0*dx1*dy1)
   dzy1 = 1.0d0/(4.0d0*dz1*dy1)

   rate = zero


   ys = drange(5)+1
   ye = nyv(0)
   zs = drange(3)
   ze = drange(4)
   xs = drange(1)
   xe = nxv(0)

! evaluate diffusion and reaction terms for the first equation
!
!I am not reupdating the advection velocity of the lambda
!seeking a fist order update

   eqn = 1
   if(neqmax > 1)  CALL RATE_EXP

   do j = ys, ye
      do k = zs, ze
         do i = xs, xe

            gx = (f(i+1,k,j,eqn)-f(i-1,k,j,eqn))*dx1a
            gy = (f(i,k,j+1,eqn)-f(i,k,j-1,eqn))*dy1a
            gz = (f(i,k+1,j,eqn)-f(i,k-1,j,eqn))*dz1a

            gxx = (f(i+1,k,j,eqn)-2.0d0*f(i,k,j,eqn)&
                 +f(i-1,k,j,eqn))*dx2a
            gyy = (f(i,k,j+1,eqn)-2.0d0*f(i,k,j,eqn)&
                 +f(i,k,j-1,eqn))*dy2a
            gzz = (f(i,k+1,j,eqn)-2.0d0*f(i,k,j,eqn)&
                 +f(i,k-1,j,eqn))*dz2a

            gxy = (f(i+1,k,j+1,eqn)-f(i+1,k,j-1,eqn)&
                 -f(i-1,k,j+1,eqn)+f(i-1,k,j-1,eqn))*dxy1
            gzy = (f(i,k+1,j+1,eqn)-f(i,k+1,j-1,eqn)&
                 -f(i,k-1,j+1,eqn)+f(i,k-1,j-1,eqn))*dzy1

            term=1.0d0+dphidx(i,k)**2+dphidz(i,k)**2

            conv = - vvel(i,k,j)*detady(j)*gy &
                 - uvel(i,k,j)*gx-wvel(i,k,j)*gz
            diff =  lambdag(i,k,j)*(gxx + gzz +&
                 (term)*(&
                 (detady(j)**2)*gyy+deta2dy2(j)*gy)&
                 -2.0d0*dphidx(i,k)*detady(j)*gxy&
                 -2.0d0*dphidz(i,k)*detady(j)*gzy&
                 -(dphi2dx2(i,k)+dphi2dz2(i,k))*detady(j)*gy)/cp&
                 +((dlgdx(i,k,j)-dphidx(i,k)*detady(j)*dlgdy(i,k,j))*&
                 (gx-dphidx(i,k)*detady(j)*gy)&
                 +(dlgdz(i,k,j)-dphidz(i,k)*detady(j)*dlgdy(i,k,j))*&
                 (gz-dphidz(i,k)*detady(j)*gy)&
                 +(detady(j)**2)*dlgdy(i,k,j)*gy)/cp
            rate(i,k,j,4) = -(conv+diff+rate(i,k,j,eqn))/f(i,k,j,1)
            rate(i,k,j,1:3) =  zero

         end do
      end do
   end do

!--------------------------------------------------------------------------
   RETURN
 END SUBROUTINE PROJECTION_SOURCE
!*********************************************************************
