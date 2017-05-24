!**********************************************************************
subroutine global_parameter
  use global_data 

  implicit none
!---------------------------------------------------------------------
! local variables
  integer :: k
  real*8  ::  th_st_c,th_ap_c
  real*8  ::  sigma_ap, sigma_st, cp_st
  real*8  ::  rblimit
  logical ::  unity_mks
!--------------------------------------------------------------------

  cp_ap = cp
  cp_binder = cp_ap
  cp_gas = cp_ap

  beta = rho_m_st/(1.0d0-rho_m_st)

  rblimit = 25d0
  Tlimit = (/theta_ap/log(da_ap/rblimit),theta_binder/log(da_binder/rblimit)/)

!--cp_st is used to find the heats of reaction from the adiabatic flame temps.
!  cp_st = rho_m_st*cp_ap + (1.0d0-rho_m_st)*cp_binder   !solid cp stoich 
!
! qg(1) = cp_gas*th_ap - cp_ap*tcold - qheat_ox
! qg(3) = (cp_gas*th_st - cp_st*tcold) &  
!      / (1.0d0-rho_m_st)             &  
!      - rho_m_st*(qg(1)+qheat_ox)      & 
!      / (1.0d0-rho_m_st)             &
!      - qheat_binder
! qg(2) = (cp_gas*th_st - cp_st*tcold) & 
!      / (1.0d0-rho_m_st)             &
!      - rho_m_st*(qheat_ox)          &
!      / (1.0d0-rho_m_st)             &
!      - qheat_binder

!  th_st_c = cp_st*tcold + rho_m_st*(qg(1)+qheat_ox)  &
!       +  (1.0d0-rho_m_st)*(qg(3)+qheat_binder)
!  th_ap_c = cp_ap*tcold + (qg(1)+qheat_ox)
!  th_st_c = th_st_c / cp_gas
!  th_ap_c = th_ap_c / cp_gas

! IF (myid == 0) THEN
!    WRITE (6, *) 'T_H = ', th_st_c, th_ap_c
!    WRITE (6, *) 'beta = ', beta
!    WRITE (6, *)
! END IF

!     
!     set coefficients of reaction terms 
!     
!    yf = zero

     yf(1,1:nsteps) = qg(1:nsteps)/cp_gas

!    hardwired for 3 stepp
!    yf(2,1) = -1.0d0
!    yf(2,2) = -beta
!    yf(2,3) =  0.0d0

!    yf(3,1) = +1.0d0
!    yf(3,2) = 0d0
!    yf(3,3) = -beta

!    yf(4,1) = 0d0
!    yf(4,2) = -1d0
!    yf(4,3) = -1d0
 
!    This needs to be cleaned up
     if (NORADIATION.eqv..FALSE.) then
        qg(irxn) = 4500.0d0      ! AL + OX -> Al2O3
        yf(1,irxn) = qg(irxn)/cp_gas
        yf(2:irxn,irxn) =  0.0d0
        yf(neqgas-1,irxn) = -1d0
        yf(neqgas,irxn) = beta_al
        do k = 1,3
           mu(k,5:ubound(mu,2)) = 0d0
        end do
        mu(irxn,1:irxn) = 0d0
        mu(irxn,irxn+1) = 1d0
        mu(irxn,irxn+2) = 0d0
     endif


  if(myid == 0) then
     write(*,'(80("#")/,a,/,20("-"))')'output from  global parameters:'
     write(*,*)'output from  global parameters:'
     write(6,*)'equivalence_ratio = ',equivalence_ratio
     write(6,*)'beta = ',beta
     write(6,*)'alpha_pack = ',alpha_pack
     write(6,*)'da    = ',da(1:nsteps)
     write(6,*)'qg1,qg2,qg3    = ',qg(1:nsteps)
     write(6,*)'thetay = ',thetay(1:nsteps)
     write(6,*)'qheat_solid  = ',qheat_ox,qheat_binder
     write(6,*)'lambda_solid  = ',lambda_ap,lambda_binder
     write(6,*)'nx,ny  = ',nx,ny
     write(6,*)'xend,yend = ',xend,yend
     write(6,*)'c1y = ',  c1y
     write(6,*)'theta_ap,theta_binder = ',theta_ap,theta_binder
     write(6,*)'da_ap,da_binder = ',da_ap,da_binder
     write(6,*)'Limit Temperatures',Tlimit
     write(*,'(20("-"),/a/,80("#"))')'end global parameters'
  endif


  return
end subroutine global_parameter
!******************************************************************************
!******************************************************************************
subroutine read_chem_info

  use global_data

  implicit none
!---------------------------------------------------------------------
! local variables
  integer :: n
!--------------------------------------------------------------------
  !open(unit=16,file="chem_info.txt")
  open (unit=16,file=chem_inputfile_name)

!-set maximum size

  read(16,*) nsteps
  read(16,*) rho_m_st

  read(16,*) tflame1
  read(16,*) tflame2
  read(16,*) Tcold

  do n=1,nsteps
     read(16,*) da(n)
     read(16,*) thetay(n)
     read(16,*) np(n)
     read(16,*) zetay(n)
     read(16,*) qg(n)
  enddo

  do n=1,nsteps
     read(16,*) yf(n+1,1:nsteps)
  enddo

  read(16,*) lewis(1:4)

  do n = 1,nsteps
     read(16,*) mu(n,1:neqgas)
  enddo

  read(16,*) a_cpg,b_cpg,e_cpg,t_cpg
  read(16,*)  a_lamg,b_lamg,e_lamg

  read(16,*) rho_ap
  read(16,*) cp_ap
  read(16,*) lambda_ap
  read(16,*) qheat_ox

  read(16,*) rho_binder
  read(16,*) cp_binder
  read(16,*) lambda_binder
  read(16,*) qheat_binder

  read(16,*) da_ap
  read(16,*) theta_ap
  read(16,*) nus_ap
  read(16,*) da_binder
  read(16,*) theta_binder
  read(16,*) nus_binder

  qsolidpar = 0d0
  qsolidpar(1) = qheat_ox

!additional read to nullify defaults
  f_pyr(1:4) = 0d0
  f_pyr(2) = 1d0


  rho_al  = 2.7d0
  rho_ala = 3.97
  beta_rhoal = rho_ala/rho_al
  lambda_al  = 0.56606477d0
  lambda_ala = 0.0430


100 continue

  close(16)

  if(myid /= 0) goto 200

  write(6,*)'irxn=',irxn
  write(6,*)'neqgas=',neqgas

  open(unit=16,file='chem_info.chk')

  write(16,'(1pe14.6,"   ........   rho_m_st")') rho_m_st

  write(16,'(1pe14.6,"   ........   tflame1")') tflame1
  write(16,'(1pe14.6,"   ........   tflame2")') tflame2

  write(16,'(1pe14.6,"   ........   da")') da(1:nsteps)
  write(16,'(1pe14.6,"   ........   theta1")') thetay(1:nsteps)
  write(16,'(1pe14.6,"   ........   n1")') np(1:nsteps)
  write(16,'(1pe14.6,"   ........   z1")') zetay(1:nsteps)

  write(16,'(1p4e14.6,"   ........   lewis(1:4)")') lewis(1:4)

  do n = 1,irxn
     write(16,'(1p4e14.6,"   ........   nu(n,1:neqgas)")') mu(n,1:neqgas)
  enddo

  write(16,'(1p4e14.6,"   ........   a_cpg,b_cpg,e_cpg,t_cpg")') a_cpg,b_cpg,e_cpg,t_cpg
  write(16,'(1p3e14.6,"   ........   a_lamg,b_lamg,e_lamg")') a_lamg,b_lamg,e_lamg

  write(16,'(1pe14.6,"   ........   rho_ap")') rho_ap
  write(16,'(1pe14.6,"   ........   cp_ap")') cp_ap
  write(16,'(1pe14.6,"   ........   lambda_ap")') lambda_ap
  write(16,'(1pe14.6,"   ........   qheat_ox")') qheat_ox

  write(16,'(1pe14.6,"   ........   rho_binder")') rho_binder
  write(16,'(1pe14.6,"   ........   cp_binder")') cp_binder
  write(16,'(1pe14.6,"   ........   lambda_binder")') lambda_binder
  write(16,'(1pe14.6,"   ........   qheat_binder")') qheat_binder

  write(16,'(1pe14.6,"   ........   da_ap")') da_ap
  write(16,'(1pe14.6,"   ........   theta_ap")') theta_ap
  write(16,'(1pe14.6,"   ........   nus_ap")') nus_ap
  write(16,'(1pe14.6,"   ........   da_binder")') da_binder
  write(16,'(1pe14.6,"   ........   theta_binder")') theta_binder
  write(16,'(1pe14.6,"   ........   nus_binder")') nus_binder

  write(16,'(1p3e14.6,"   ........   qsolidpar(3vals)")') qsolidpar
  write(16,'(1p4e14.6,"   ........   f_pyr(4vals)")')f_pyr(1:4)

  close(16)

200 continue

  lewis(5:ubound(lewis,1)) = 1.2d1

!-----------------------------------------------------------------------------
  return
end subroutine read_chem_info
!******************************************************************************
