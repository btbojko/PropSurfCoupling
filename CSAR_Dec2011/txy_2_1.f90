!********************************************************************
SUBROUTINE UPDATE_TXY

USE GLOBAL_DATA
USE MYMPI

IMPLICIT NONE

!---------------------------------------------------------------------
! Local Variables
INTEGER :: i, j, k, n, ys,ye, eqn,sub_iter,newton_max
REAL*8 :: dminus1zero,pressure,poo,coe_rho
LOGICAL :: outcome
!---------------------------------------------------------------------

!---  Initialize variables
ys = drange(5)
ye = drange(6)
if(ddrange(5).eq.0) ys = ys + 1
if(ddrange(6).eq.ny) ye = ye - 1

outcome = .TRUE.

!--- Setup the pointers
finest_mesh%w => f
var_updated = 0
finest_mesh%cm => coeff
finest_mesh%r => dfdt

!.....................................................................
! START THE NEWTON--GAUSS LOOP
!.....................................................................

!--- SET NUMBER of cycles for temp-species iterations
gauss_cycles = 2
sub_iter=merge(8,8,ncyc < ncyc_init)
newton_max = merge(1,1,ncyc < 100); 

!----STORE OLD SOLUTION
do i = -1,nxv(0)+1
    do k = 0,drange(4)
        do j = 0,nyv(0)+1
            do n = 1,neqmax
                oldsoln(i,k,j,n) = f(i,k,j,n)
            enddo
        enddo
    enddo
enddo

dt = timestep

CALL VELOCITY(0)

CALL LAMBDA(0) 

!-- coefficient matrix
CALL COEFFMATRIX
!
! update the rhs of the equations
!
CALL EXP_RHS

DO it_new=1,newton_max
    CALL RATESETUP
    CALL ADJUST_LHS
!    save the lu decomp; the matrix changes every newton iteration
    saved_LU_vars = .FALSE.     
    converged = .FALSE.

    itime_sos = 0
    do icycle=1,gauss_cycles
        if(itime_sos == 0) THEN
            CALL gauss_relax4 (finest_mesh,sub_iter,.TRUE.)
        else
            CALL gauss_relax5 (finest_mesh,sub_iter,.TRUE.)
        endif
        IF(converged) EXIT
    enddo
ENDDO     ! DO it_new=1,max_new
!

call NAN_FIELD(outcome,.FALSE.)


if(ipack <= 0 .AND. coe_d == 0 .and. &
       mod(ncyc,merge(10,500,ipack == 0)) == 0) call RHS_test
!-------------------------------------------------------------------------
RETURN
END SUBROUTINE UPDATE_TXY
!*************************************************************************


!****************************************************************************************
SUBROUTINE COEFFMATRIX

USE GLOBAL_DATA
USE MYMPI

IMPLICIT NONE

!----------------------------------------------------------------------------------------
! Local variables
INTEGER ::  i, j, k, eqn, ys,ye,xs,xe,zs,ze
REAL*8  :: dx1, dx2, dy1, dy2, dz1, dz2, dx1a, dx2a, dy1a, dy2a, dz1a, dz2a,term
REAL*8  :: dxy1, dzy1,dlx1,dlx2,dlz1,dlz2,dly1,dly2,term1,dya,dyya,dyyya
REAL*8  :: dtcp,dxcp,dzcp,dycp,dtloc,alldtmin
REAL*8  :: invcp,cpr, convx,convz,convy,dlcsi,dlzet,vterm
! Fluids vars
REAL*8  :: pressure,poo,coe_invrho,solid_dt
REAL*8  :: c2,c3,c4,c5,c6,c7,c8,c9,c11,sumc
LOGICAL :: drop2first
!----------------------------------------------------------------------------------------

poo = PRESSURE(tcyc)
coe_invrho = one / (dim_fact*poo)

dx1 = dx
dx2 = dx*dx
dy1 = dy
dy2 = dy*dy
dz1 = dz
dz2 = dz*dz

ys = drange(5)
ye = drange(6)
if(ddrange(5).eq.0) ys = ys + 1
if(ddrange(6).eq.ny) ye = ye - 1
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
dxcp = cp/(2.0d0*dx1) 
dx2a = one/(1.0d0*dx2)
dy1a = one/(2.0d0*dy1)
dy2a = one/(1.0d0*dy2)
dz1a = one/(2.0d0*dz1)
dzcp = cp/(2.0d0*dz1) 
dz2a = one/(1.0d0*dz2)
dxy1 = one/(4.0d0*dx1*dy1)
dzy1 = one/(4.0d0*dz1*dy1)

! RECALL coeff(1,:,:,:) is the LHS of the equation.

  dtmin = dt
DO j = ys, nyv(0)
    dya =   detady(j)*dy1a
    dycp =  dya*cp
    dyya =  detady(j)**2*dy2a
    dyyya = deta2dy2(j)*dy1a
    DO k = zs, ze
        DO i = xs, nxv(0)
            dlcsi = dlgdx(i,k,j)-dphidx(i,k)*detady(j)*dlgdy(i,k,j)
            dlzet = dlgdz(i,k,j)-dphidz(i,k)*detady(j)*dlgdy(i,k,j)

            convx = (uvel(i,k,j) )*dxcp - dx1a*dlcsi
            convz = (wvel(i,k,j) )*dzcp - dz1a*dlzet
            convy = (vvel(i,k,j) )*dycp &
                + dya*dlcsi*dphidx(i,k) +&
                dya*dlzet*dphidz(i,k) - dya*dlgdy(i,k,j)*detady(j)

            term = one+dphidx(i,k)**2+dphidz(i,k)**2
            term1 = (dphi2dx2(i,k)+dphi2dz2(i,k))*dya

            coeff(2,i,k,j) = lambdag(i,k,j)*dx2a + convx
            coeff(3,i,k,j) = lambdag(i,k,j)*dx2a - convx

            coeff(4,i,k,j) = convy + lambdag(i,k,j)*(&
                term*(dyya-dyyya)+term1)
            coeff(5,i,k,j) = two*lambdag(i,k,j)*term*dyya - coeff(4,i,k,j)

            dconv(i,k,j) = zero
            drop2first = coeff(5,i,k,j) < zero .OR. coeff(4,i,k,j) < zero
            if(is_TVD .AND. drop2first ) then   !drop to first order
                dconv(i,k,j) = vvel(i,k,j)*dycp
                coeff(4,i,k,j) =  coeff(4,i,k,j) + abs(dconv(i,k,j))
                coeff(5,i,k,j) =  coeff(5,i,k,j) + abs(dconv(i,k,j))
            endif
            coeff(6,i,k,j) = 0.0d0
            coeff(7,i,k,j) = 0.0d0
            coeff(8,i,k,j) = -two*lambdag(i,k,j)*dphidx(i,k)*detady(j)*dxy1
            coeff(9,i,k,j) = 0.0d0 
        ENDDO
    ENDDO
ENDDO

!
!--- SOLID JACOBIAN MATRIX. NOTE, for solid only coeff(1) = 1/coeff(1)
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

do j = ys, nyv(0)
    do k = zs, ze
        do i = xs, nxv(0)
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

            csld(2,i,k,j) = dlx2*dx2a/cpr                   !(lambda*Tx)x
            csld(3,i,k,j) = dlx1*dx2a/cpr                   !(lambda*Tx)x
            csld(5,i,k,j) = vterm+&
                (term*(detady(j)**2*dly1*dy2a+deta2dy2(j)*lambdas(i,k,j)*dy1a)+term1&
                *lambdas(i,k,j))/cpr                        !ftTeta+(1+fx^2)*(lambdaTeta)eta-fxx*lambda*Teta
            csld(4,i,k,j) = term*detady(j)**2*dy2a*(dly1 + dly2)/cpr-csld(5,i,k,j)
            csld(8,i,k,j) = dphidx(i,k)*detady(j)*dxy1/cpr
!           csld(1,i,k,j) = one/(one + sum(csld (2:7,i,k,j)))
        enddo
    enddo
enddo

!  CALL DIRICHLET_BC(-10, finest_mesh)

DO j = ys, nyv(0)
    DO k = zs, ze
        DO i = xs, nxv(0)
            dtcp = dtx(i,k,j)*invcp

            coeff(2,i,k,j) = coeff(2,i,k,j)*dtcp
            coeff(3,i,k,j) = coeff(3,i,k,j)*dtcp
            coeff(4,i,k,j) = coeff(4,i,k,j)*dtcp
            coeff(5,i,k,j) = coeff(5,i,k,j)*dtcp
            coeff(8,i,k,j) = coeff(8,i,k,j)*dtcp

!           coeff(1,i,k,j) = one + sum(coeff(2:7,i,k,j))   !for point relaxation
           
            coeff(1,i,k,j) = one + coeff(2,i,k,j) + coeff(3,i,k,j) &
                + coeff(4,i,k,j) + coeff(5,i,k,j)
!           csld(1,i,k,j) = coe_d + sum(csld (2:7,i,k,j))
            csld(1,i,k,j) = coe_d + csld(2,i,k,j) + csld(3,i,k,j) &
                + csld(4,i,k,j) + csld(5,i,k,j)
            if(is_TVD ) then 
                dconv(i,k,j) = dconv(i,k,j)*dtcp
            endif
        ENDDO
    ENDDO
ENDDO

!------------------------------------------------------------------------------------
RETURN 
END SUBROUTINE COEFFMATRIX
!****************************************************************************************

!****************************************************************************************
SUBROUTINE ADJUST_LHS

  USE GLOBAL_DATA

  IMPLICIT NONE

  include 'parallel.h'

!----------------------------------------------------------------------------------------
! Local variables
  TYPE(mesh_pointers), POINTER :: current_mesh
  INTEGER ::  i,j,k,ii,jj,kk, mn, ys,ye,xs,xe,zs,ze
  INTEGER ::  m,n, neqn, ip,kp,jp, col(4)
  REAL*8  ::  corrflux,s1,s2,s3,sig1,sig2,nuy
!----------------------------------------------------------------------------------------

  neqn = neqmax-1

  current_mesh => finest_mesh

  ys = 1
  ye = drange(6)-1
  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)

!
! N-S Boundary cond for B
!
  Bm(1:neqgas,1:neqgas,:,:,0) = 0.0d0
  Bm(1:neqgas,1:neqgas,:,:,ny) = 0.0d0
  do m = 1,neqgas
     Bm(m,m,:,:,0) = 1.0d0
     Bm(m,m,:,:,ny) = 1.0d0
  enddo

  DO j = ys, nyv(0)
     do k = zs, ze
        do i = xs, nxv(0)


           do m=1,neqn

!           Bm(m,m,i,k,j) = one/(coeff(1,i,k,j) - Bm(m,m,i,k,j))   !for point relaxation
              Bm(m,m,i,k,j) = coeff(1,i,k,j) + Bm(m,m,i,k,j)          !for line relaxation

           enddo

        enddo
     enddo
  ENDDO

!------TVD stuff in this version TVD operators only in y direction

  if(.NOT.is_TVD) RETURN

  DO j = ys, nyv(0)
     do k = zs, ze
        do i = xs, nxv(0)

           nuy = abs(dconv(i,k,j))

           if(nuy < 1.d-12) CYCLE

           if(dconv(i,k,j) > 0) then
              jp = j
           else
              jp = j+1
           endif

           do m=1,neqn

              s1 = (f(i,k,jp+1,m) - f(i,k,jp,m))
              s2 = (f(i,k,jp,m) - f(i,k,jp-1,m))
              s3 = (f(i,k,jp-1,m) - f(i,k,jp-2,m))
              if(jp==1) s3=s2;
              if(jp==ny) s1=s2

              sig1 = (sign(half,s1)+sign(half,s2))*min(abs(s1),abs(s2))
              sig2 = (sign(half,s2)+sign(half,s3))*min(abs(s2),abs(s3))

              corrflux = - nuy * (sig1-sig2)

              dfdt(i,k,j,m) = dfdt(i,k,j,m) + corrflux

           enddo

        enddo
     enddo
  ENDDO


!------------------------------------------------------------------------------------
  RETURN 
END SUBROUTINE ADJUST_LHS
!****************************************************************************************


!**************************************************************
SUBROUTINE GET_TIME_COEFF

USE GLOBAL_DATA
USE MYMPI

IMPLICIT NONE

!-----------------------------------------------------------
!     Local Variables
!-----------------------------------------------------------

!GAS
time_coeff = merge(one,half,is_firstorder)

!SOLID
coe_dt = half 
if(ipack == 0) then
    if(ncyc < ncyc_steady)then
        coe_d = zero
    else
	    coe_d = one
    endif
endif    

if(coe_d == 0.0d0) coe_dt = one/dt

!PROJECTION
doprojection = (ncyc > ncyc_oseen .AND. nx > 10) &
                    .AND. tcyc > time_oseen

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


!*************************************************
SUBROUTINE CURB_TXY
        
  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!-----------------------------------------------------------
!!!   Local Variables
      integer i, k, j, n
      double precision max_f(neqmax),min_f(neqmax), f_old
      integer ybegin,ytop,ncurb, mynclip
!-----------------------------------------------------------


      max_f(1) = 3500.0d0
      min_f(1) = 300.0d0
      do n=2,neqmax-1
         min_f(n) = -5.d-2
         max_f(n) = 1.d0 - min_f(n)
      enddo
      max_f(neqmax) = 2400.0d0
      min_f(neqmax) = tcold - 10.0d0

      ybegin = drange(5)
      ytop = drange(6)
 
      mynclip = 0
      do j=ybegin,ytop
         do k=drange(3)-1,drange(4)+1
            do i=drange(1)-1,drange(2)+1
               ncurb = 0
               do n=1,neqmax
                  f_old = f(i,k,j,n)
                  f(i,k,j,n) = max(f(i,k,j,n), min_f(n))
                  f(i,k,j,n) = min(f(i,k,j,n), max_f(n))
                  if(f(i,k,j,n) /=  f_old) ncurb = ncurb + 1
               enddo
               if(ncurb >  0) mynclip = mynclip+1
            enddo
         enddo
      enddo

      if(mod(ncyc,500) == 0) then
         call MPI_ALLREDUCE(mynclip,nclip,1,MPI_INTEGER,  &
              MPI_SUM, comm3d,ierr)
         if( myid == 0) write(*,*)'CLIPPED',nclip
      endif

 
!--------------------------------------------------------
      RETURN
   END SUBROUTINE CURB_TXY
!**********************************************************************
!

