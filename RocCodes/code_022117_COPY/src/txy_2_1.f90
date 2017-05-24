!********************************************************************
SUBROUTINE UPDATE_TXY

  USE GLOBAL_DATA
  USE implicit
  USE MYMPI
  use timing

  IMPLICIT NONE

!---------------------------------------------------------------------
! Local Variables
  INTEGER :: i, j, k, n, ys,ye, eqn,sub_iter,newton_max,gauss_cycles
  REAL*8 :: dminus1zero,pressure,poo,coe_rho
  LOGICAL :: outcome,restarted
!---------------------------------------------------------------------

!---  Initialize variables
  ys = drange(5)
  ye = drange(6)
  if(dSerialRange(5).eq.0) ys = ys + 1
  if(dSerialRange(6).eq.ny) ye = ye - 1

  outcome = .TRUE.

  if (skipTemperature) then
     finest_mesh%w => f
     return
  endif

!.....................................................................
!.....................................................................

!--- Set iterative method; this over writes implicit.txt
  iterativemethod = 1     ! use gmres, only for radiation
  if (NORADIATION) iterativemethod = 2  ! use newton otherwise


!--- SET NUMBER of cycles for temp-species iterations
  gauss_cycles = 1
  sub_iter=merge(5,5,ncyc < ncyc_init)

!--- Setup the pointers
  finest_mesh%w => dff
  finest_mesh%r => dfdt

  !
  ! Store old solution
  !
  DO i = -1, drange(2)+1
  DO k = -1, drange(4)+1
  DO j = 0, drange(6)
     DO n = 1, neqmax
        oldsoln(i, k, j, n) = f(i, k, j, n)
     END DO
  END DO
  END DO
  END DO


!--- Start main temperature and species solver
     Call PYROLYSIS
     CALL LAMBDA(0)
     CALL VELOCITY(0)
     CALL COEFFMATRIX   ! fill --> coeff; clsd

     ! fills dfdt in exp_rhs (crank Nicholson)

     CALL RATESETUP     ! fill --> Bm
     CALL EXP_RHS       ! fill --> dfdt

     CALL ADJUST_LHS    ! This adds coeff + Bm

     dff = 0d0

     converged = .FALSE.

     if(iterativemethod==1)then
        Do icycle = 1,GMRES_restarts
           call GMres_solve(iterLINEARsys,gauss_GMRES,icycle>1)
        end Do
     else
        do icycle=1,gauss_cycles
           CALL gauss_relax4 (finest_mesh,dfdt,sub_iter,.TRUE.)
           IF(converged) EXIT
        enddo
     end if

  do i = -1,drange(2)+1
  do k = -1,drange(4)+1
  do j = 0,drange(6)
     do n = 1,neqmax
        f(i,k,j,n) = f(i,k,j,n) + dff(i,k,j,n)
     enddo
  enddo
  enddo
  enddo

  finest_mesh%w => f

END SUBROUTINE UPDATE_TXY
!*********************************************************************


!*********************************************************************
SUBROUTINE COEFFMATRIX

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!---------------------------------------------------------------------
! Local variables
  INTEGER ::  i, j, k, eqn, ys,ye,xs,xe,zs,ze,L
  REAL*8  :: dx1, dx2, dy1, dy2, dz1, dz2, dx1a, dx2a, dy1a, dy2a, dz1a, dz2a,term
  REAL*8  :: dxy1, dzy1,dlx1,dlx2,dlz1,dlz2,dly1,dly2,term1,dya,dyya,dyyya
  REAL*8  :: dtcpp,dxcp1a,dzcp1a,dycp1a,dtloc,alldtmin
  REAL*8  :: invcp,cpr, convx,convz,convy,dlcsi,dlzet,vterm
  REAL*8  :: pressure,solid_dt
  REAL*8  :: c2,c3,c4,c5,c6,c7,c8,c9,c11,sumc
  REAL*8  :: termxzsq,termxxzz,dlambdady,diffus
  REAL*8  :: t1,t31,t32,t33,DNM,meshincrements(8)
  LOGICAL :: drop2first(3)
  logical  ::outcome
!----------------------------------------------------------------------------------------

  dx1 = dx
  dx2 = dx*dx
  dy1 = dy
  dy2 = dy*dy
  dz1 = dz
  dz2 = dz*dz

  ys = drange(5)
  ye = drange(6)
  if(dSerialRange(5).eq.0) ys = ys + 1
  if(dSerialRange(6).eq.ny) ye = ye - 1
  zs = drange(3)
  ze = drange(4)
  xs = drange(1)
  xe = drange(2)

!.................................................................................
! THE FIRST PART IS FOR GAS ONLY THE SECOND IS FOR SOLID ONLY
!.................................................................................

  eqn=1
  invcp = one/cp

  dx1a = one/(2.0d0*dx1) 
  dxcp1a = cp/(2.0d0*dx1) 
  dx2a = one/(1.0d0*dx2)
  dy1a = one/(2.0d0*dy1)
  dy2a = one/(1.0d0*dy2)
  dz1a = one/(2.0d0*dz1)
  dzcp1a = cp/(2.0d0*dz1) 
  dz2a = one/(1.0d0*dz2)
  dxy1 = one/(4.0d0*dx1*dy1)
  dzy1 = one/(4.0d0*dz1*dy1)

! coeff(1,:,:,:) is ON the LHS of the equation.

  meshincrements = (/dx2a,dz2a,dy2a,dxy1,dzy1,dx1a,dz1a,dy1a/)
  mcomp = 1
  !!$  coeff(eqn,10,:,:,:) = coeff(eqn,1,:,:,:)

  DO eqn = 1,neqgas
     do k = zs,ze
        do i = xs,xe

!independent of j
           termxzsq = dphidx(i,k)**2 + dphidz(i,k)**2
           termxxzz = dphi2dx2(i,k)  + dphi2dz2(i,k)

           do j = ys,ye

!
              diffus = lambdag(i,k,j)/lewis(eqn)
              dlambdady = dlgdy(i,k,j)*detady(j)/lewis(eqn)
              t1 = MQchi(j,mcomp)*termxzsq
              t31 = (one +  MQchipr(j,mcomp)*(MQphi(i,k)))
              t32 = t31*t31
              t33 = t31*t32
              DNM = one/t33
              MQnonC(1:2) = diffus
              MQnonC(3) = diffus/t32*(one+t1*MQchi(j,mcomp))
              MQnonC(4) = -two*MQchi(j,mcomp)*dphidx(i,k)*diffus/t31
              MQnonC(5) = -two*MQchi(j,mcomp)*dphidz(i,k)*diffus/t31
              MQnonC(6) = dlgdx(i,k,j)/lewis(eqn) - MQchi(j,mcomp)*dphidx(i,k)*dlambdady/t31
              MQnonC(7) = dlgdz(i,k,j)/lewis(eqn) - MQchi(j,mcomp)*dphidz(i,k)*dlambdady/t31
              MQnonC(8) = DNM *( diffus*( t1* ( two*MQchipr(j,mcomp) + &
                   MQphi(i,k)*(two*MQchipr(j,mcomp)**2-MQchi(j,mcomp)*MQchisc(j,mcomp))) &
                   -MQchi(j,mcomp)*termxxzz*t32 - MQphi(i,k)*MQchisc(j,mcomp)) &
                   +dlambdady*t31*(one+MQchi(j,mcomp)*t1)&
                   -MQchi(j,mcomp)*t32*(dlgdx(i,k,j)*dphidx(i,k)+dlgdz(i,k,j)*dphidz(i,k))/lewis(eqn)  )

!change terms because of the stretch mapping
              MQnonC(8) = MQnonC(8)*detady(j) + MQnonC(3)*deta2dy2(j)
              MQnonC(3) = MQnonC(3) *detady(j)**2
              MQnonC(4:5) = MQnonC(4:5)*detady(j)
!divide by the specific heat
              MQnonC = MQnonC/cp
!add convective terms
              MQnonC(6) = MQnonC(6) - uvel(i,k,j)
              MQnonC(7) = MQnonC(7) - wvel(i,k,j)
              MQnonC(8) = MQnonC(8) - vvel(i,k,j)*detady(j)
!multiply by coe_exp*dt/rho
              MQnonC = MQnonC*dtx(i,k,j)

!multiply by mesh increments
              MQnonC = MQnonC* meshincrements

!these are needed to determine the upwinding
              uwvterm(1) = (uvel(i,k,j))*dx1a*dtx(i,k,j)
              uwvterm(2) = (wvel(i,k,j))*dz1a*dtx(i,k,j)
              uwvterm(3) = (vvel(i,k,j))*detady(j)*dy1a*dtx(i,k,j)

!assemble matrix
              coeff(eqn,2,i,k,j) = MQnonC(1) - MQnonC(6) 
              coeff(eqn,3,i,k,j) = MQnonC(1) + MQnonC(6) 
              coeff(eqn,4,i,k,j) = MQnonC(2) - MQnonC(7) 
              coeff(eqn,5,i,k,j) = MQnonC(2) + MQnonC(7) 
              coeff(eqn,6,i,k,j) = MQnonC(3) - MQnonC(8) 
              coeff(eqn,7,i,k,j) = MQnonC(3) + MQnonC(8) 
              coeff(eqn,8:9,i,k,j) = MQnonC(4:5)

              iconv(i,k,j,1:ndim) = merge( (/ (1,L=1,ndim)/),(/ (0,L=1,ndim)/),uwvterm(1:ndim)<=0)
              dconv(i,k,j,1:ndim) = zero
              drop2first(1) = (any(coeff(eqn,2:3,i,k,j)< zero) .or. lewis(eqn)>1d1)
              drop2first(2) = (any(coeff(eqn,4:5,i,k,j)< zero) .or. lewis(eqn)>1d1)
              drop2first(3) = (any(coeff(eqn,6:7,i,k,j)< zero) .or. lewis(eqn)>1d1)

              if(drop2first(1) ) then   !drop to first order
                 dconv(i,k,j,1) = abs(uwvterm(1))   
                 coeff(eqn,2:3,i,k,j) =  coeff(eqn,2:3,i,k,j) + dconv(i,k,j,1)
              endif
              if(drop2first(2) ) then   !drop to first order 
                 dconv(i,k,j,2) = abs(uwvterm(2))
                 coeff(eqn,4:5,i,k,j) =  coeff(eqn,4:5,i,k,j) + dconv(i,k,j,2)
              endif
              if(drop2first(3) ) then   !drop to first order  
                 dconv(i,k,j,3) = abs(uwvterm(3))
                 coeff(eqn,6:7,i,k,j) =  coeff(eqn,6:7,i,k,j) + dconv(i,k,j,3)
              endif

              coeff(eqn,1,i,k,j) =  one + sum(coeff(eqn,2:7,i,k,j))
              !!$           coeff(eqn,10,i,k,j) =  one + sum(abs(coeff(eqn,2:9,i,k,j)))

           ENDDO
        ENDDO
     ENDDO
  ENDDO
!
!--- SOLID JACOBIAN MATRIX. NOTE, for solid only coeff(eqn,1) = 1/coeff(eqn,1)
!
  solid_dt = coe_dt*dt

  dx1a = solid_dt/(2.0d0*dx1)
  dx2a = solid_dt/(1.0d0*dx2)
  dy1a = solid_dt/(2.0d0*dy1)
  dy2a = solid_dt/(1.0d0*dy2)
  dz1a = solid_dt/(2.0d0*dz1)
  dz2a = solid_dt/(1.0d0*dz2)
  dxy1 = solid_dt/(4.0d0*dx1*dy1)
  dzy1 = solid_dt/(4.0d0*dz1*dy1)

  do j = ys, ye
     do k = zs, ze
        do i = xs, xe


           term=one+dphidx(i,k)**2+dphidz(i,k)**2
           term1=(dphi2dx2(i,k)+dphi2dz2(i,k))*detady(j)*dy1a

           dlx1 = 2.0d0*lambdas(i+1,k,j)*lambdas(i,k,j)/ &
                (lambdas(i+1,k,j) + lambdas(i,k,j))
           dlx2 = 2.0d0*lambdas(i,k,j)*lambdas(i-1,k,j)/ &
                (lambdas(i,k,j) + lambdas(i-1,k,j))
           dlz1 = 2.0d0*lambdas(i,k+1,j)*lambdas(i,k,j)/ &
                (lambdas(i,k+1,j) + lambdas(i,k,j))
           dlz2 = 2.0d0*lambdas(i,k,j)*lambdas(i,k-1,j)/ &
                (lambdas(i,k,j) + lambdas(i,k-1,j))
           dly1 = 2.0d0*lambdas(i,k,j+1)*lambdas(i,k,j)/ &
                (lambdas(i,k,j+1) + lambdas(i,k,j))
           dly2 = 2.0d0*lambdas(i,k,j)*lambdas(i,k,j-1)/ &
                (lambdas(i,k,j) + lambdas(i,k,j-1))

           cpr=cp*rhos(i,k,j)
           vterm = ft0(i,k)*detady(j)*dy1a
           d_vterm(i,k) = d_ft0(i,k)*dy1a

           csld(2,i,k,j) = dlx2*dx2a/cpr
           csld(3,i,k,j) = dlx1*dx2a/cpr
           csld(4,i,k,j) = dlz2*dz2a/cpr
           csld(5,i,k,j) = dlz1*dz2a/cpr
           csld(7,i,k,j) = vterm+&
                (term*(detady(j)**2*dly1*dy2a+deta2dy2(j)*lambdas(i,k,j)*dy1a)+term1&
                *lambdas(i,k,j))/cpr
           csld(6,i,k,j) = term*detady(j)**2*dy2a*(dly1 + dly2)/cpr-csld(7,i,k,j)
           csld(8,i,k,j) = dphidx(i,k)*detady(j)*dxy1/cpr
           csld(9,i,k,j) = dphidz(i,k)*detady(j)*dzy1/cpr
!         csld(1,i,k,j) = one/(one + sum(csld (2:7,i,k,j)))   !for point relaxation
           csld(1,i,k,j) = coe_d + sum(csld (2:7,i,k,j))

        enddo
     enddo
  enddo

!------------------------------------------------------------------------------------
  RETURN 
END SUBROUTINE COEFFMATRIX
!****************************************************************************************

!*************************************************
SUBROUTINE CURB_TXY

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!-----------------------------------------------------------
!!!   Local Variables
  integer i, k, j, l,n
  REAL*8 :: f_old
  integer ybegin,ytop,ncurb,mynclip(maxsize+2)
!-----------------------------------------------------------


  max_f(1) = 4.0d3
  min_f(1) = 500.0d0

  do n=2,neqgas
     min_f(n) = -5.d-2
     min_f(n) = 0.d0
  enddo
  if(neqgas > 4) then
     min_f(5:6) = -1d-4
     min_f(5:6) = 0d0
  end if
  max_f(2:neqgas) = 1.d0 - min_f(2:neqgas)
  max_f(2:neqgas) = 1.d0

  max_f(neqmax) = 2550.0d0
  min_f(neqmax) = tcold - 100.0d0

  ybegin = drange(5)
  ytop = drange(6)
!
  mynclip = 0
  do j=ybegin,ytop
     do k=drange(3)-1,drange(4)+1
        do i=drange(1)-1,drange(2)+1
           ncurb = 0
           do n=1,neqmax
              f_old = f(i,k,j,n)
              f(i,k,j,n) = max(f(i,k,j,n), min_f(n))
              f(i,k,j,n) = min(f(i,k,j,n), max_f(n))
              if(f(i,k,j,n) /=  f_old) then
                 mynclip(n) = mynclip(n)+1
                 if(n == neqmax .and. k >= drange(3).and. k <= drange(4)) then
!!>                    print'(3i3,a,1p123e12.4)',i,k,j,'Curb T_solid',f(i,k,j,n),f_old
                 elseif(n > 4 .and. k >= drange(3).and. k <= drange(4)) then
!!>                    print'(4i3,a,1p123e12.4)',n,i,k,j,'CURB n>4',f(i,k,j,n),f_old  !cancel (uncommnet)
                 end if
              end if
           enddo
        enddo
     enddo
  enddo
  f(:,:,0,neqmax) = f(:,:,0,1)

!clobber small data every 2*writemod iterations
!!>  if (mod(ncyc_run-1,2*writemod) == 0) then
  if (mod(ncyc_run-1,1) == 0) then
     do l = lbound(f,4),ubound(f,4)
        do j = lbound(f,3),ubound(f,3)
           do k = lbound(f,2),ubound(f,2)
              do i = lbound(f,1),ubound(f,1)
                 if(abs(f(i,k,j,l))<= 1d-35)  f(i,k,j,l) = 0d0
              end do
           end do
        end do
     end do
  end if


  call MPI_ALLREDUCE(mynclip,nclip,size(nclip),MPI_INTEGER,  &
       MPI_SUM, comm3d,ierr)
  if(maxval(nclip) > 1) then
     nclipwrite = nclipwrite+1
     if(nclipwrite <= 2) then
        if( myid == 0) write(*,*)ncyc,'CLIPPED, sorted by equation',nclip
     end if
  end if
  if(mod(ncyc,max(10*writemod,100)) == 0) nclipwrite = 0
    

!--------------------------------------------------------
  RETURN
END SUBROUTINE CURB_TXY
!**********************************************************************
!


!**************************************************************
SUBROUTINE GET_TIME_COEFF

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!-----------------------------------------------------------
!     Local Variables
  real*8 :: time_coeff_min,time_coeff_max
!-----------------------------------------------------------

!GAS
  time_coeff = merge(one,half,is_firstorder)
  time_coeff_min = time_coeff/dble(Tfrac)
  time_coeff_max = time_coeff

  if(Tjumps*dble(Tncyc) > 1d-14) then
     time_coeff = min(time_coeff_min + (time_coeff_max-time_coeff_min)*&
          &floor(dble(ncyc)/dble(Tncyc))/Tjumps,time_coeff_max)
  endif

!-- SOLID
  if(ipack == 0) coe_d = zero   !always steady
  coe_dt = half 
  if(coe_d == 0.0d0) coe_dt = one/dt
  if(ncyc < 2*NCYCRAD)  coe_dt = 1d2

!PROJECTION
  doprojection = (ncyc > ncyc_oseen .AND. nx > 10)

  if(mod(ncyc,500) == 0 .AND. myid == 0) write(*,*)&
       'time_coeff,coe_dt,doprojection',&
       time_coeff,coe_dt,doprojection,dt

  if(coe_d == 0.0d0) then
     flag_comp = ' steady'
  else
     flag_comp = ' un_steady'
  endif
  if(doprojection) then
     flag_comp = TRIM(flag_comp) // ' Navier Stokes'
  else
     flag_comp = TRIM(flag_comp) // ' Oseen'
  endif

!----------------------------------------------------
  RETURN
END SUBROUTINE GET_TIME_COEFF
!**********************************************************************
