!
!========================================================
!
subroutine space_avg_dump(t,tout)
  USE GLOBAL_DATA
  USE MYMPI

!     Local variables
  integer j, n, m
  integer xbeg,xstp,xrange
  integer yrange
  REAL*8 :: dy1a,dy2a
  REAL*8 :: gy,gyy,g2eta2
  REAL*8 :: t,tout,per_test,dt_step
  REAL*8 :: pressure, poo, tdump
!-------------------------------------------------


  xbeg = drange(1)
  xstp = drange(2)

  xrange = xstp - xbeg + 1
  yrange = drange(6) - drange(5) + 1

  CALL averagePHIT_space(t,tout)

!    
!     CHECK FOR NEW PERIOD
!

  per_test = dint(phi0/period)

!>      if(per_test .lt. nperiod-1) return



  dt_step = tout - t      
  dt_integ = dt_integ + dt_step

!>      CALL averageT_space(drange,t,tout)
!>
!>      CALL prturb_space(drange)
!>      CALL termseval_space(drange,t,tout)
  CALL surf_termseval_space(t,tout)

  if(mod(ncyc,n_print_terms) .ne. 0) return

  ndump_s_terms = ndump_s_terms + 1 

!      Eval space  avg. terms and surface terms
!
  do n = 1,n_surf_terms
     surf_terms0(n) = surf_terms(n)/dt_integ
     surf_terms(n) = 0.0d0
  enddo

  tdump = tout - dt_integ/2.0d0
  poo = pressure(tdump)
  if(myid == 0) then
     write(*,*)ncyc,tdump,'DUMP',ndump_s_terms,'  TIME INTERV',dt_integ
     write(*,*) ndump_s_terms,tdump,poo, (surf_terms0(m),m=1,n_surf_terms)
  endif


  dt_integ = 0.0d0

  if(myid .eq. 0) then
     write(28,*) ndump_s_terms,tdump,poo, (surf_terms0(m),m=1,n_surf_terms)
  endif



  return
end subroutine space_avg_dump
!
!========================================================
!
subroutine ICaverage_space
  USE GLOBAL_DATA
  USE MYMPI

!     Local variables
!    
  integer j,n
  REAL*8 :: tmp_x,tmp_t,tmp_h
!     

  iper = 0
  pos_period = 0.0d0
  time_period = 0.0d0
  do n = 1,n_surf_terms
     surf_terms(n) = 0.0d0
     surf_terms0(n) = 0.0d0
  enddo

  rb00 = 0.0d0
  rb0save = 0.0d0   !correspond to rb0
  ndump_s_terms = 0
  dt_integ = 0.0d0

  tmp_x=lambda_ap/lambda_binder
  tmp_t=0.8201292467014901
  rblend = tmp_t*rho_ap   +(1.0d0-tmp_t)*rho_binder

  tmp_h=tmp_x+0.5d0*(1.0d0-tmp_x)**2.0d0 *(1.0d0-tmp_t)**2.0d0
  if(tmp_x.lt.1.0d0) then
     lblend=lambda_binder*     (tmp_h+sqrt(tmp_h*tmp_h-tmp_x*tmp_x))
  else
     lblend=lambda_binder* (tmp_h-sqrt(tmp_h*tmp_h-tmp_x*tmp_x))
  endif
  alfa_blend = lblend*fact_blend

  open(28,file='spacesurfterms.dat')

  return
end subroutine ICaverage_space
!
!========================================================

!
!========================================================
!
subroutine averagePHIT_space(t,tout)
  USE GLOBAL_DATA
  USE MYMPI

!     Local variables
  integer i
  integer xbeg,xstp,xrange
  integer yrange,totrange
  REAL*8 :: my_phi,all_phi
  REAL*8 :: t,tout,dt_step
  REAL*8 :: ftmp,per_test,per_old
!
  dt_step = tout - t
  time_period = time_period + dt_step

  xbeg = drange(1)
  xstp = drange(2)

  xrange = xstp - xbeg + 1
  yrange = drange(6) - drange(5) + 1


  ftmp = 0.0d0
  do i =  xbeg,xstp
     ftmp = ftmp + phit(i,0)
  enddo
  my_phi = ftmp*dt_step

  call MPI_ALLREDUCE (my_phi,all_phi,1,MPI_DOUBLE_PRECISION,&
       MPI_SUM,comm3d,ierr) 
  call MPI_ALLREDUCE (xrange,totrange,1,MPI_DOUBLE_PRECISION,&
       MPI_SUM,comm3d,ierr) 

  rb00 = rb00 + all_phi/nx   !space-time average 


  per_test = mod(-phi(0,0),-period)
  per_old = mod(-oldphi(0,0),-period)

  call MPI_BCAST(per_test,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)
  call MPI_BCAST(per_old,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)

  if(per_old .gt. per_test) then
     rb0save = abs(rb00/time_period)
     rb00 = 0.0d0
     time_period = 0.0d0
     ndump_s_terms = 0
  endif

  return
end subroutine averagePHIT_space

!
!======================================================
!   
!
!========================================================
!
!
!the next subroutine is identical to the previous one,but
!considers surface source terms which are still 
!integrated in space and time
!========================================================
!
subroutine surf_termseval_space(t,tout)
  USE GLOBAL_DATA
  USE MYMPI

!     Local variables
  integer i, j, n, e
  integer xbeg,xstp,xrange,totrange
  REAL*8 :: termtmp(n_surf_terms)
  REAL*8 :: t,tout,dt_step,const
  REAL*8 :: myterm(n_surf_terms)
  REAL*8 :: allterm(n_surf_terms)
  REAL*8 :: hx,hy,g1,g2,g3,g4,tx,ty
!

!>      if(dint(phi0/period) .lt. 1.0) return

  dt_step = tout - t

  xbeg = drange(1)
  xstp = drange(2)



  xrange = xstp - xbeg + 1


  hx = 1.0d0/(2.0d0*dx)
  hy = - detady(0)/(2.0d0*dy)

  e = neqmax
  j=0

  do n = 1,n_surf_terms
     termtmp(n) = 0
  enddo
  do i =  xbeg,xstp  

     const = sqrt(1.0+dphidx(i,0)**2)

     Tx = (f(i+1,0,0,e)-f(i-1,0,0,e))*hx
     Ty = (-3.0d0*f(i,0,0,e)+4.0d0*f(i,0,1,e) -f(i,0,2,e))*hy

     g1 = phit(i,0)
     g2 = (const*Ty - dphidx(i,0) * Tx)
     g3 = rhos(i,0,j)*phit(i,0)
     g4 = f(i,0,0,e)


     termtmp(1) = termtmp(1) + g1    !a
     termtmp(2) = termtmp(2) + g2    !b
     termtmp(3) = termtmp(3) + g3    !c
     termtmp(4) = termtmp(4) + g4    !c
  enddo

  do n = 1,n_surf_terms
     myterm(n) = termtmp(n)
  enddo


  call MPI_ALLREDUCE (myterm(1),allterm(1),n_surf_terms,&
       MPI_DOUBLE_PRECISION, MPI_SUM, comm3d,ierr)

  call MPI_ALLREDUCE (xrange,totrange,1,MPI_DOUBLE_PRECISION,&
       MPI_SUM,comm3d,ierr) 


  do n = 1,n_surf_terms
     surf_terms(n) = surf_terms(n) + allterm(n)*dt_step/nx
  enddo

!

  if(myid .eq. 0 .AND. mod(ncyc,2*n_print_terms) .eq. 0) then
     write(*,*)'SURFACE TERMS',(allterm(n)/nx,n=1,n_surf_terms)
     write(*,*)'SURFACE TERMS_SPACE',&
          (surf_terms(n),n=1,n_surf_terms)
  endif
!

  return
end subroutine surf_termseval_space
!
!========================================================
