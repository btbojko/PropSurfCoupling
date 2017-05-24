!*******************************************************
subroutine reinitia(n_iter,space_ord,time_ord)
!*******************************************************
!     
!     reinitialization of level set fun!tion by solving 
!     Hamilton equation
!     
!-----------------------------------------
USE GLOBAL_DATA
IMPLICIT NONE

integer space_ord,time_ord
integer i,j,ire,n_iter,i_rk
!     
REAL*8 ::  dt_re,cour_cond
real*8 ::  a_20,a_21,a_30,a_32,b_2,b_3
real*8 ::  time_reinitia

time_reinitia = zero

cour_cond = 0.11
dt_re = min(dx,dy/detady(0))*cour_cond

!     
do i = -1 , drange(2)+1
    do j = -3 , ny+3   !Gross 
        if(j==-3)then
            ff(i,j) = sign(one, psi(i,0,j+1))
        else    
            ff(i,j) = sign(one, psi(i,0,j))
        endif
    end do
end do

DO ire = 1, n_iter 
    time_reinitia = time_reinitia + dt_re
    call pde_function
    do i =0,drange(2)
        do j = 0,ny_l_re
            ff(i,j) = ff(i,j) + dt_re*L_phi(i,j)
        enddo
    enddo  
    bc_ptr => ff
    do i =1,2
        call bc_p(i)
    enddo
!
ENDDO

return
end subroutine reinitia
!     
!=======================================================!     
!     
subroutine pde_function

USE GLOBAL_DATA
IMPLICIT NONE    
!     
!-----------------------------------------
!        
integer :: i,j
real*8 ::  sij,dp,dm,prd,sm
real*8 ::  df_x,df_y,dyavg
logical cond1,cond2,cond3
!
bc_ptr => ff
do i =1,2
    call bc_p(i)
enddo
!
call eno_2d
!
do  i = drange(1), drange(2)
    do  j=0,ny_l_re
        sij = sign(one,psi(i,0,j))
        dyavg = 0.5d0*(dyp(i,j)+dym(i,j))
        dp = dxp(i,j) - dphidx(i,0)*dyavg
        dm = dxm(i,j) - dphidx(i,0)*dyavg
        prd = dp*dm
        sm = sij*(dp+dm)
        cond1 = dp*sij .lt. zero .and. sm .lt. zero
        cond2 = dm*sij .gt. zero .and. sm .gt. zero
        cond3 = prd .lt. zero .and. dm*sij .lt. zero
        if(cond1) then
            df_x = dp
        elseif(cond2) then
            df_x = dm
        elseif(cond3) then
            df_x = zero
        else
            df_x = zero
        endif
        dp = dyp(i,j)
        dm = dym(i,j)
        prd = dp*dm
        sm = sij*(dp+dm)
        cond1 = dp*sij .lt. zero .and. sm .lt. zero
        cond2 = dm*sij .gt. zero .and. sm .gt. zero
        cond3 = prd .lt. zero .and. dm*sij .lt. zero
        if(cond1) then
            df_y = dp
        elseif(cond2) then
            df_y = dm
        elseif(cond3) then
            df_y = zero
        else
            df_y = zero
        endif

!            df_x = df_x - dphidx(i,0)*df_y

        L_phi(i,j) = sij*(one-sqrt(df_x*df_x + df_y*df_y))

    enddo
enddo


return
end subroutine pde_function
!
!========================================================
!

    subroutine minabs(a,b,c,kmin)
      real*8 ::  a,b,c
      integer kmin

      if(abs(a) .le. abs(b)) then
         c = a
         kmin = kmin - 1
      else
         c = b
      endif


      return
    end subroutine minabs

!
!========================================================
!
    subroutine div_diff_table (d1,d2,ijst,ijnd)

      USE GLOBAL_DATA
      IMPLICIT NONE    



      real*8 ::  dx_i,dy_j
      integer i,j,d1,d2,ijst,ijnd

      bc_ptr => tab


      if(abs(d1).eq.1) then
         do i = drange(1),drange(2)
            do j = -3,ny_l_re+3
               dx_i = d1 / (xy(i+ijnd,j) - xy(i+ijst,j))
               tab(i,j) = (ff(i+d1,j) - ff(i,j))*dx_i
            enddo
         enddo
         call bc_p(1)
      elseif(abs(d2).eq.1) then
         do j = 0,ny_l_re
            do i = -1,drange(2)+1
               dy_j = d2 / (xy(i,j+ijnd) - xy(i,j+ijst))
               tab(i,j) = (ff(i,j+d2) - ff(i,j))*dy_j
            enddo
         enddo
         call bc_p(2)
      endif

      return
    end subroutine div_diff_table


!    
!========================================================
!     
    subroutine eno_2d

      USE GLOBAL_DATA
      USE MYMPI
      IMPLICIT NONE    

      integer i,j,l,m,kmin,klm1,inegpos,n_ord,iorder
      integer limited_order
      real*8 ::  tmp_deriv(2)

!---------xdirection
      do i =  drange(1)-1, drange(2) +1
         xy(i,-1:ny+1) = x(i)
      enddo

      call div_diff_table(1,0,0,1)

      do i = 0,drange(2)           !--- xdirection
         do j = 0,ny_l_re
            do inegpos = 1,2

               kmin = i + inegpos - 2
               tmp_deriv(inegpos) = tab(kmin,j)


            enddo

            dxm(i,j) = tmp_deriv(1)
            dxp(i,j) = tmp_deriv(2)
         enddo
      enddo

!-----------------YDRIECTION
      do j =  0, ny_l_re+3
         do i = -1,drange(2)+1
            xy(i,j) = y(j)
         enddo
      enddo

      call div_diff_table(0,1,0,1)

      do i = 0,drange(2)  
         do j = 0,ny_l_re
            do inegpos = 1,2
               kmin = j + inegpos - 2
               tmp_deriv(inegpos) = tab(i,kmin)

            enddo

            dym(i,j) = tmp_deriv(1)
            dyp(i,j) = tmp_deriv(2)
         enddo
      enddo

      return
    end subroutine eno_2d

!
!========================================================
!

    subroutine bc_p (flag)


      USE GLOBAL_DATA
      USE MYMPI
      IMPLICIT NONE 
!----------------------------------------- 

      integer ::  xs,xe
      integer i,j,elcount,flag
      INTEGER ::  status(MPI_STATUS_SIZE)

      xs = 0
      xe = drange(2)

      if( flag == 1 )then

         elcount = 0
         DO j = 0, ny_l_re
            slambda(elcount) =  bc_ptr(xe,j)  
            elcount = elcount + 1
         END DO
         CALL MPI_SENDRECV(slambda,elcount,&
              MPI_DOUBLE_PRECISION,dnbr(2),0,rlambda,&
              elcount,MPI_DOUBLE_PRECISION,dnbr(1),0,&
              comm3d,status,ierr)
         elcount = 0
         DO j = 0, ny_l_re
            bc_ptr(xs-1,j)  =  rlambda(elcount)   
            slambda(elcount) = bc_ptr(xs,j)
            elcount = elcount + 1
         END DO

         CALL MPI_SENDRECV(slambda,elcount,&
              MPI_DOUBLE_PRECISION,dnbr(1),1,rlambda,&
              elcount,MPI_DOUBLE_PRECISION,dnbr(2),1,&
              comm3d,status,ierr)


         elcount = 0
         DO j = 0, ny_l_re
            bc_ptr(xe+1,j) = rlambda(elcount)   
            elcount = elcount + 1
         END DO

      elseif( flag == 2 ) then


         do i = drange(1)-1,drange(2)+1
            do j = 1,2
               bc_ptr(i,-j) = bc_ptr(i,0)
               bc_ptr(i,ny_l_re+j) = bc_ptr(i,ny_l_re)
            enddo
         enddo

      endif
!

      return
    end subroutine bc_p
!
!=======================================================
!

    real*8 function minmod(x0,y0)
      real*8 ::  x0,y0
!
      if (x0*y0 .ge. 0) then
         if(dabs(x0) .le. dabs(y0)) then
            minmod=x0
         else
            minmod=y0
         end if
      else
         minmod=0.d0
      end if
!
      return
    end function minmod

!
!========================================================
!

    real*8 function Heaviside(val1,val2,phi0,limit)
      real*8 :: val1,val2,phi0,limit,dval
!
      if (phi0 <= -limit) then
         Heaviside = val1
      elseif(phi0 >= limit) then
         Heaviside = val2
      else
         dval = val2-val1
         Heaviside = val1 + dval*(phi0+limit)/(2.0d0*limit)
      endif
!
      return
    end function Heaviside

!
!========================================================
!
