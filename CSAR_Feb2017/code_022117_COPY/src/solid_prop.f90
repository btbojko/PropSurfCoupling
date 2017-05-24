! ********************************************************************
! ********************************************************************
! Updated: TLJ; 12/19/2016
! Filename: solid_prop.f90
! ********************************************************************
! ********************************************************************
!     ________________________________________________________________
!
!     BRIEF CODE DESCRIPTION:
!     =======================
!
!     This set of subroutines sets the mass fractions at the
!     surface from the solid side, as well as compute the
!     solid volume properties Qs, rho, lambda, Da, and theta.
!
!     subroutines: solid_prop
!                  blend
!     functions:   lambda_3d
!                  lamroot3
!
! ********************************************************************

subroutine solid_prop

  use global_data

  implicit none

!---------------------------------------------------------------------
! local variables
  integer :: i,j,k,m
  real*8  :: tmp_x, tmp_h, tmp_t,tmp_f,factorf_pyr,qheat_p1
!---------------------------------------------------------------------

  do k=drange(3)-1,drange(4)+1
  do i=drange(1)-1,drange(2)+1
!
!    AP surface
!
     xsurf(i,k,2) = 1.0d0            ! X; oxidizer/AP
     xsurf(i,k,3) = 0.0d0            ! Z; decomposition
     xsurf(i,k,4) = 0.0d0            ! Y; fuel/binder
     xsurf(i,k,5) = 0.0d0            ! Y_al; aluminum
     xsurf(i,k,6) = 0.0d0            ! Y_oxide; oxide

     psurf(i,k,1) = 0d0
     psurf(i,k,2) = 0d0
     psurf(i,k,3) = 1d0

     if (psi(i,k,0) < zero) then
!
!    binder surface
!
        xsurf(i,k,2) = alp_w(1)
        xsurf(i,k,3) = 0d0
        xsurf(i,k,4) = alp_w(3)
        xsurf(i,k,5) = alp_w(2)
        xsurf(i,k,6) = 0d0

        psurf(i,k,1) = 0d0
        psurf(i,k,2) = 0d0
        psurf(i,k,3) = alp_w(1)
     endif
  enddo
  enddo
!
! surface and volume properties
!
  do k=drange(3)-1,drange(4)+1
  do i=drange(1)-1,drange(2)+1
  do j=drange(5),drange(6)

     qheats(i,k,j)  = qheat_ox
     rhos(i,k,j)    = rho_ap
     lambdas(i,k,j) = lambda_ap
     da_s(i,k,j)    = da_ap*press**nus_ap
     theta_s(i,k,j) = theta_ap

     ! binder
     if ( psi(i,k,j) < zero ) then
        qheats(i,k,j)  = alp_w(1)*qheat_ox + alp_w(3)*qheat_binder
        rhos(i,k,j)    = rho_eff
        lambdas(i,k,j) = lambda_eff
        da_s(i,k,j)    = (da_ap*press**nus_ap)**alp_v(1)&
                         & * (da_binder*press**nus_binder)**alp_v(3)
        theta_s(i,k,j) = alp_v(1)*theta_ap+alp_v(3)*theta_binder
     endif
  enddo
  enddo
  enddo
!---------------------------------------------------------------------
  return
end subroutine solid_prop
! ********************************************************************
! ********************************************************************

subroutine blend(t3,x3,y3,f,r,nele)
  USE MYMPI
  implicit none
  integer,intent(in) :: nele
  real*8,intent(in)  :: t3(nele),x3(nele),y3(nele)
  real*8,intent(out)  :: f,r
  integer :: ierr
  real*8 :: t,x,lamroot3
  real*8 :: a,b,c,d,p,q,p1,al,pp1,zterm1,r1,det3
  integer :: iter,maxi

  if(nele == 1) then
     f = x3(1)
     r = y3(1)
  else
     t = t3(1)/sum(t3(1:2))
     x = x3(1)/x3(2)
     if(lamroot3(x,t) < -1d0) then
        print*,'Multiple solution lamroot3_2',x3,t3
        call mpi_finalize(ierr)
        stop 'blend'
     end if
     f = x3(2)*lamroot3(x,t)
     r = sum(t3(1:2)*y3(1:2))

     if(nele == 3 .and. t3(3) > 1d-10) then
        t = t3(3)/sum(t3(1:3))
        x = x3(3)/f
        if(lamroot3(x,t) < -1d0) then
           print*,'Multiple solution lamroot3_2',x3,t3
           call mpi_finalize(ierr)
           stop 'blend'
        end if
        f = f*lamroot3(x,t)
        r = r + t3(3)*y3(3)
     end if
  end if


  return
end subroutine blend
! ********************************************************************
! ********************************************************************

real *8 function lambda_3d(x,t,f)
  implicit none
  real*8,intent(in)  :: t, x
  real*8,intent(inout)  ::f

  lambda_3d = ((1.0d0-t)*(1.0d0-x)/(f-x))**3 * f - 1.0d0;

  return
end function lambda_3d

real*8 function lamroot3(x,t)

  implicit none
  integer :: ierr
  real*8,intent(in) :: x,t
  real*8 :: a,b,c,d,p,q,det3,al,p1,pp1,r1,zterm

  d = -x**3
  c = 3d0*x+x**3+3d0*t-1d0+9d0*t*x**2-3d0*t*x**3-3d0*t**2+9d0*t**2*x&
       &-9d0*t**2*x**2+3d0*t**2*x**3+t**3-3d0*t**3*x+ &
       &3d0*t**3*x**2-t**3*x**3-9d0*t*x
  b = -3d0*x;
  a= 1d0
  p = (3d0*a*c-b**2)/(9d0*a**2);
  q = (2d0*b**3-9d0*a*b*c+27d0*a**2*d)/(27d0*a**3)
  det3 = 4d0*p**3+q**2
  if(det3 < 0d0) then
     write(*,*) 'multiple solution in sub blend','x,t',x,t
      lamroot3 = - 99d0
      return
  end if
  al = (-q+sqrt(q**2+4d0*p**3))/2d0;
  p1 = 8d0*al;
  pp1 = sign(1d0,p1)*abs(p1)**(1d0/3d0)
  r1 = 5d-1*pp1 - 2d0*p/pp1
  lamroot3 = r1-b/3d0/a;
end function lamroot3
