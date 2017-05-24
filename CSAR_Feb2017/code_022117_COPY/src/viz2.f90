
!**********************************************************
subroutine gather_xz

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

!  Local variables
  integer i, mystart, myend, myxsize
  integer mysize,error
!--------------------------------------------------------


!note change to diplacemtent1 and allsize 1 to keep info about the other values
  mystart = drange(1)
  myend = drange(2)
  if(rowid(1) == rowsize(1)-1) myend = drange(2)+1
  myxsize = (myend - mystart) + 1
  mysize = myxsize

  allsize1 = 0
  call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allsize1,1,MPI_INTEGER,&
       rowcomm(1),ierr)

  displacements1(1) = 0
  do i = 2,nproc
     displacements1(i) = allsize1(i-1) + displacements1(i-1)
  enddo
  do i = mystart,myend
     myphi1(i-mystart) = x(i)
  enddo

  call MPI_ALLGATHERV(myphi1(0),mysize,MPI_DOUBLE_PRECISION, &
       allphi1(0),allsize1,displacements1,MPI_DOUBLE_PRECISION, &
       rowcomm(1),ierr)
!!
  do i = 0,nx
     G_x(i) = allphi1(i)
  enddo

  if(is2d) return

!
!=============================
!

  mystart = drange(3)
  myend = drange(4)
  if(rowid(1) == rowsize(2)-1) myend = drange(4)+1
  myxsize = (myend - mystart) + 1
  mysize = myxsize

  allsize1 = 0
  call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allsize1,1,MPI_INTEGER,&
       rowcomm(2),ierr)

  displacements1(1) = 0
  do i = 2,nproc
     displacements1(i) = allsize1(i-1) + displacements1(i-1)
  enddo
  do i = mystart,myend
     myphi1(i-mystart) = z(i)
  enddo

  call MPI_ALLGATHERV(myphi1(0),mysize,MPI_DOUBLE_PRECISION, &
       allphi1(0),allsize1,displacements1,MPI_DOUBLE_PRECISION, &
       rowcomm(2),ierr)
  do i = 0,nx
     G_z(i) = allphi1(i)
  enddo

!____________________________________________________
  return
end subroutine gather_xz
!**************************************************************


!*****************************************************************************
subroutine gatherPsiSurf

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

!     Local variables
  integer i, j,k,l,m, n,pcs,flag,nIN
  integer elcount,myzsize,el,ixz, myxsize
  integer mysize,elm,elsize, lastd
  integer mydrange(4),qplet(4)
!---------------------------------------------------------------

  mydrange(1) = drange(1)
  mydrange(2) = drange(2)
  if(rowid(1) == rowsize(1)-1) mydrange(2) = drange(2)+1

  mydrange(3) = drange(3)
  mydrange(4) = drange(4)
  if(rowid(2) == rowsize(2)-1) mydrange(4) = drange(4)+1
!
!  myxsize = (mydrange(2) - mydrange(1)) + 1
!  myzsize = (mydrange(4) - mydrange(3)) + 1
!  mysize = myxsize * myzsize
!  call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allsize,1,MPI_INTEGER,&
!       comm3d,ierr)
!
!  displacements(1) = 0
!  do i = 2,nproc
!     displacements(i) = allsize(i-1) + displacements(i-1)
!  enddo
!
!  lastd =  allsize(nproc) + displacements(nproc)

  if(any(alldrange < 0)) then
     print*,myid,'gatherPsiSurf NEGATIVE DRANGE',alldrange,rowid
     stop 'gatherPsiSurf'
  endif

  elcount = 0
  do i = mydrange(1),mydrange(2)
     do k = mydrange(3),mydrange(4)
        myphi2(elcount) = sign(1.0d0,psi(i,k,0))
        elcount = elcount + 1
     enddo
  enddo
  call MPI_ALLGATHERV(myphi2(0),mysizeB,MPI_DOUBLE_PRECISION,&
       allphi2(0),allsize,displacements,MPI_DOUBLE_PRECISION,&
       comm3d, ierr)

  elcount = 0
  do pcs = 0, nproc-1
     el = 1 + (pcs*6)
     do i = alldrange(el),alldrange(el+1)
        do k = alldrange(el+2),alldrange(el+3)
           G_psi(i,k) = allphi2(elcount)
           elcount = elcount + 1 
        enddo
     enddo
  enddo

!-----------------------------------------------
  return
end subroutine gatherPsiSurf
!**********************************************************

!*****************************************************************************
subroutine gatherPhiSurf

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

!     Local variables
  integer i, j,k,l,m, n,pcs,flag,nIN
  integer elcount,myzsize,el,ixz, myxsize
  integer mysize,elm,elsize, lastd
  integer mydrange(4),qplet(4),mySdrange(4)
!---------------------------------------------------------------

  if(.not. allocated(displacements)) allocate(displacements(nproc+1))
  if(.not. allocated(alldrange)) allocate(alldrange(4*nproc))


  mydrange(1) = drange(1)
  mydrange(2) = drange(2)
  mydrange(3) = drange(3)
  mydrange(4) = drange(4)

  if(mydrange(1).eq.0) mydrange(1) = 0
  if(mydrange(2)+ib_mg(1).eq.nx-1) mydrange(2) = mydrange(2)+1
  if(mydrange(3).eq.0) mydrange(3) = 0
  if(mydrange(4)+kb_mg(1).eq.nz-1) mydrange(4) = mydrange(4)+1

  mySdrange(1:2) = mydrange(1:2) + ib_mg(1)
  mySdrange(3:4) = mydrange(3:4) + kb_mg(1)
!
  myxsize = (mydrange(2) - mydrange(1)) + 1
  myzsize = (mydrange(4) - mydrange(3)) + 1
  mysize = myxsize * myzsize
  call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allsize,1,MPI_INTEGER,comm3d,ierr)
  call MPI_ALLGATHER(mySdrange,4,MPI_INTEGER,alldrange,4,MPI_INTEGER,comm3d,ierr)

!
  displacements(1) = 0
  do i = 2,nproc
     displacements(i) = allsize(i-1) + displacements(i-1)
  enddo

  lastd =  allsize(nproc) + displacements(nproc)

  if(any(alldrange < 0)) then
     print*,myid,'gatherPhiSurf NEGATIVE DRANGE',alldrange,rowid
     stop 'gatherPhiSurf '
  endif

  elcount = 0
  do i = mydrange(1),mydrange(2)
     do k = mydrange(3),mydrange(4)
        myphi2(elcount) = phi(i,k)
        elcount = elcount + 1
     enddo
  enddo
  call MPI_ALLGATHERV(myphi2(0),mysizeB,MPI_DOUBLE_PRECISION,&
       allphi2(0),allsize,displacements,MPI_DOUBLE_PRECISION,&
       comm3d, ierr)

  elcount = 0
  do pcs = 0, nproc-1
     el = 1 + (pcs*4)
     do i = alldrange(el),alldrange(el+1)
        do k = alldrange(el+2),alldrange(el+3)
           G_phi(i,k) = allphi2(elcount)
           elcount = elcount + 1 
        enddo
     enddo
  enddo

!-----------------------------------------------
  return
end subroutine gatherPhiSurf
!**********************************************************

!*****************************************************************************
subroutine gatherPhitSurf

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

  !     Local variables
  integer i, j,k,l,m, n,pcs,flag,nIN
  integer elcount,myzsize,el,ixz, myxsize
  integer mysize,elm,elsize, lastd
  integer mydrange(4),qplet(4)
  !---------------------------------------------------------------

  mydrange(1) = drange(1)
  mydrange(2) = drange(2)
  if(rowid(1) == rowsize(1)-1) mydrange(2) = drange(2)+1

  mydrange(3) = drange(3)
  mydrange(4) = drange(4)
  if(rowid(2) == rowsize(2)-1) mydrange(4) = drange(4)+1
  !
  !  myxsize = (mydrange(2) - mydrange(1)) + 1
  !  myzsize = (mydrange(4) - mydrange(3)) + 1
  !  mysize = myxsize * myzsize
  !  call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allsize,1,MPI_INTEGER,&
  !       comm3d,ierr)
  !
  !  displacements(1) = 0
  !  do i = 2,nproc
  !     displacements(i) = allsize(i-1) + displacements(i-1)
  !  enddo
  !
  !  lastd =  allsize(nproc) + displacements(nproc)

  if(any(alldrange < 0)) then
     print*,myid,'gatherPhitSurf NEGATIVE DRANGE',alldrange,rowid
     stop 'gatherPhitSurf'
  endif

  do m = 1,2
     elcount = 0
     do i = mydrange(1),mydrange(2)
        do k = mydrange(3),mydrange(4)
           myphi2(elcount) = merge(phit(i,k),rb(i,k),m == 1)
           elcount = elcount + 1
        enddo
     enddo
     call MPI_ALLGATHERV(myphi2(0),mysizeB,MPI_DOUBLE_PRECISION,&
          allphi2(0),allsize,displacements,MPI_DOUBLE_PRECISION,&
          comm3d, ierr)

     elcount = 0
     do pcs = 0, nproc-1
        el = 1 + (pcs*6)
        do i = alldrange(el),alldrange(el+1)
           do k = alldrange(el+2),alldrange(el+3)
              if(m == 1) then
                 G_phi(i,k) = allphi2(elcount)
              else
                 G_psi(i,k) = allphi2(elcount)
              end if
              elcount = elcount + 1 
           enddo
        enddo
     enddo
  end do

  !-----------------------------------------------
  return
end subroutine gatherPhitSurf
!**********************************************************



!**********************************************************!
SUBROUTINE PRN_GLOB_VARS

  USE GLOBAL_DATA

  IMPLICIT NONE

  INTEGER :: w_unit
!-------------------------------------------------

  if(myid /= nproc-1) RETURN

  w_unit = 112

  if(w_unit /= 6) open(UNIT=w_unit,FILE='run_parameters.dat',form='formatted') 

  write(w_unit,*) 'NUMBER OF POINTS'
  write(w_unit,*) nx
  write(w_unit,*) nz
  write(w_unit,*) ny
  write(w_unit,*) 'INCREMENTS'
  write(w_unit,*) dx
  write(w_unit,*) dz
  write(w_unit,*) dy
  write(w_unit,*) 'Smallest dy',y(1)-y(0)
  write(w_unit,*) '-------------'
  write(w_unit,*) 'DRANGE IS',drange,'PROC--',myid
  write(w_unit,*) 'XEND,ZEND,YEND'
  write(w_unit,*) xend,myid
  write(w_unit,*) zend,myid
  write(w_unit,*) yend,myid
  write(w_unit,*) 'C1Y'
  write(w_unit,*) c1y,myid

  write(w_unit,*)'MG parameters',size(dx_mg)
  write(w_unit,*)'DX',dx_mg
  write(w_unit,*)'DZ',dz_mg
  write(w_unit,*)'DY',dy_mg
  write(w_unit,*) '-------------'
  write(w_unit,*) 
  write(w_unit,*)'SANDWICH PARAMETERS'
  write(w_unit,*)'XLOC',xloc


  write(w_unit,*) '-------------'
  write(w_unit,*)
  write(w_unit,*)'SOLID PPROPERTIES'
  write(w_unit,*)'RHO_AP,RHO_BIND',rho_ap,rho_binder
  write(w_unit,*)'LAMBDA_AP,LAMBDA_BIND',lambda_ap,lambda_binder
  write(w_unit,*) 'HOMOGENIZATION LAMBDA',lambda_eff
  write(w_unit,*) 'HOMOGENIZATION RHO',rho_eff
  write(w_unit,*) 'HOMOGENIZATION ALPHA',alp_V

  write(w_unit,*) '-------------'
  write(w_unit,*) ''
  write(w_unit,*) 'SURFACE REACTIONS'
  write(w_unit,*) 'Da_ap = ',da_ap, 'Da_bind = ',da_binder
  write(w_unit,*) 'Theta_ap = ',theta_ap, 'Theta_bind = ',theta_binder
  write(w_unit,*) ''
  write(w_unit,*) 'Limiting teperatures on surf for OX and binder',Tlimit(1:2)
  write(w_unit,*) 'Factor RB',factorRB
  write(w_unit,*) '-------------'
  write(w_unit,*) 
  write(w_unit,*) 'GAS PHASE REACTIONS'
  write(w_unit,*) 'STOICHIOMETRIC PARAMETERS:'
  write(w_unit,*) 'Stoichimetric mass fraction Map/Mtot',rho_m_st
  write(w_unit,*) 'Stoichimetric mass ratio Map/Mbinder (beta) ',beta
  write(w_unit,*) 'PACK mass ratio Map/Mbinder (beta_pack) ',beta_pack
  write(w_unit,*) 'HOMOGENIZED VOLUME FRACTION (alpha) ',alpha
  write(w_unit,*) '-------------'
  write(w_unit,*) 
  write(w_unit,*) 'HEAT OF REACTION:'
  write(w_unit,*) 'USING FLAME TEMPERATURES Tap,Tblend(rho_m_st)',th_ap,th_st
  write(w_unit,*) 'FIRST  REACTION, Qg1',qg1
  write(w_unit,*) 'SECOND REACTION, Qg2',qg2
  write(w_unit,*) 'THIRD  REACTION, Qg3',qg3
  write(w_unit,*) '................'
  write(w_unit,*) 'KINETIC PARAMETERS (Calibrated)'
  write(w_unit,*) da1,'        da1'
  write(w_unit,*) theta1,'      theta1'
  write(w_unit,*) n1,'       n1'
  write(w_unit,*) da2,'       da2'
  write(w_unit,*) theta2,'      theta2'
  write(w_unit,*) n2,'        n2'
  write(w_unit,*) da3,'        da3'
  write(w_unit,*) theta3,'      theta3'
  write(w_unit,*) n3,'         n3'
  write(w_unit,*) lewis(1:4),'         Lewis numbers' 
  write(w_unit,*) '-------------'
  write(w_unit,*) 

  
  write(w_unit,*) 'PRESSURE PERTURBATION PARAMETERS'  
  write(w_unit,*) 'BASE pressure',press  
  write(w_unit,*) 'Amplitude',epi_p 
  write(w_unit,*) 'Frequency',omg_p
  write(w_unit,*) 'PHASE',phase_p
  write(w_unit,*) '-------------'
  write(w_unit,*) 


  write(w_unit,*) 'RESTART FLAG ...', irestart
  write(w_unit,*) 'DUMP FILE IN ',TRIM(outputfile)
  if(irestart .eq. 1) write(w_unit,*) 'RESTART FILE IN ',TRIM(restartfile)
  write(w_unit,*) '-------------'
  write(w_unit,*) 
 


  write(w_unit,*) 'TIME DELAYS STEADY-- OSEEN',time_steady_pressure,time_oseen
  write(w_unit,*) 'CYCLES DELAYS STEADY-- OSEEN',ncyc_steady,ncyc_oseen
  write(w_unit,*) '-------------'
  write(w_unit,*) 
  
  write(w_unit,*) 'Local HOMOGENIZATION INFO'  
  write(w_unit,*) 'alp_V',alp_V
  write(w_unit,*) 'alp_W',alp_W
  write(w_unit,*) '-------------'
  write(w_unit,*) 

  write(w_unit,*) 'PACK Global HOMOGENIZATION INFO'  
  write(w_unit,*) 'alpha_pack',alpha_pack  
  write(w_unit,*) 'rho_pack',rho_pack   
  write(w_unit,*) 'alphaH_pack (w/ homogenized)',alphaH_pack  
  write(w_unit,*) 'rhoH_pack (w/ homogenized)',rhoH_pack  
  write(w_unit,*) 'lambda_pack',lambda_pack  
  write(w_unit,*) 'Q Heat Binder',alp_W*qheat_ox+(1.0-alp_W)*qheat_binder
  write(w_unit,*) '-------------'
  write(w_unit,*) 

  write(w_unit,*) 'READ FROM pack3.txt'
  write(w_unit,*) 'rhoV_RocPack, by Volume',rhoV_RocPack
  write(w_unit,*) 'alphaW_Rocpack, by weight',alphaW_Rocpack
  write(w_unit,*) 'TheoreticalPackingDensity',TheoreticalPackingDensity     !    Set periodicity
  write(w_unit,*) '-------------'
  write(w_unit,*) 



  write(w_unit,*) ''
  write(w_unit,*) 'END of FILE'
  if(w_unit /= 6)  close(w_unit)

!-----------------------------------------------------------
  RETURN
END SUBROUTINE PRN_GLOB_VARS
!***********************************************************


!******************************************************************
SUBROUTINE PRINT_1D

  USE GLOBAL_DATA
  IMPLICIT NONE

!---------------------------------------------------------------------
! Dummy Variables
  REAL *8 ::  t

! Local Variables
  INTEGER ::  P_mpierr, w_unit,w_unit2, j,n,i,k
  CHARACTER*80 filnam,filnam2
  logical :: condw
  real*8,allocatable :: yfprint(:,:)
!---------------------------------------------------------------------


  i = 0;k=0
  w_unit = 88
  w_unit2 = 95
  IF (ipack < 0) THEN
     WRITE(filnam,'("ic",I2.2,".dat")') int(press)
     WRITE(filnam2,'("rat",I2.2,".dat")') int(press)
     write(*,*)myid,'PRINT_1D ',ncyc,trim(filnam),'ALPHA',alpha
     write(*,*)ncyc,'BURN RATE',rb(0,0)!,'ALP_W',alp_w
     OPEN (w_unit,file=TRIM(filnam),form='formatted',action="write")
     OPEN (w_unit2,file=TRIM(filnam2),form='formatted',action="write")
     do j = ny,0,-1
        WRITE(w_unit,'(i3,1p123e14.6)') j,-y(j), f(1,0,j,neqmax),0,0,0
     enddo
     DO j=0,ny
        WRITE(w_unit,'(i3,1p123e14.6)') j,y(j),f(1,0,j,1:4)
     ENDDO
     DO j=0,ny
        WRITE(w_unit2,*)'RRR'
        WRITE(w_unit2,101) rate(0,0,j,1:4),vvel(i,k,j),rb(i,k)
        WRITE(w_unit2,*)'BMM'
        do n = 1,4
           WRITE(w_unit2,101) Bm(n,1:4,0,0,j)
        enddo
     ENDDO
  ENDIF

  CLOSE(w_unit)
  CLOSE(w_unit2)
  conv_error = abs(conv_error - rb(0,0))
  write(*,*)'CONVERGENCE ERROR',conv_error
  if(ncyc > frame .AND. conv_error < 1.d-9 .and. all(nclip == 0))then
     steadyconvergence = .true.
  endif
  conv_error = rb(0,0)

  condw = .true.

  if(.not. allocated(yfPrint)) allocate(yfPrint(ubound(yf,1),ubound(yf,2)))

  yfPrint = yf
  yfPrint(1,1:3) = yf(1,1:3)*cp_gas
  
  

  If(condw) Write(*,*)' '
  If(condw) Write(*,'(75("="))')
  If(condw) Write(*,*)'Intermediate output from PRINT_1D'
  if(condw) write(*,*)'Rb,Ts',rb(0,0),f(0,0,0,1),'PRESS',press,'RHO_M',alpha
  if(condw) write(*,*)'da1,da2,da3,theta1,theta2,theta3',da1,da2,da3,theta1,theta2,theta3
  if(condw) write(*,*)'CP_END',cp_gas,cp_ap,cp_binder
  if(condw) write(*,*)'Tflame',f(0,0,ny-2:ny,1)
  if(condw) write(*,*)'MAX,MIN temp',maxval(f(:,:,1:ny,1)),minval(f(:,:,1:ny,1))
  if(condw) write(*,*)'Qg',Qg1,Qg2,Qg3
  if(condw) write(*,*)'F_Pyr,qheats_cp',f_pyr
  if(condw) write(*,*)'Residual @ surf>>'
  if(condw) write(*,*)'PRE_123',da1*press**n1,da2*press**n2,da3*press**n3
  if(condw) write(*,*)'THETA_123',theta1,theta2,theta3
  if(condw) write(*,*)'Yf',yfprint
  if(condw) write(*,*)'Steadyconvergence ',Steadyconvergence
  If(condw) Write(*,'(75("*")///)')
  
  if(Steadyconvergence.or. .true.) then
     
     close(76)
     WRITE(filnam,'(a,f4.2,"_",f4.1,a)') 'solfile',alpha,press,'.DAT'
     open(76,file = filnam ,form = 'unformatted') 
     do j=0,ny 
        write(76)j,y(j),f(0,0,j,1:neqgas),press,alpha,rate(0,0,j,1:neqgas)
     end do
     close(76)
  end if



101 FORMAT(50e14.4)
102 FORMAT(1i4,1p50e20.11)

!---------------------------------------------------------------------
  RETURN 
END SUBROUTINE PRINT_1D
!*************************************************************

!*************************************************************
subroutine VizGrid 

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

  INTEGER :: i,j,k
  REAL*8 :: xx,yy, yys
  !------------------------------------------------------


  if(myid == 0) then
     open(10,file = 'grid.dat')
     do i = drange(1), drange(2)
        do k = drange(3), drange(4)
           do j = 0, ny
              yy = y(j) + fmax + MQchi(j,1)*MQphi(i,k)
              write(10,'(1p5e12.4)') x(i),z(k),yy
           end do
        end do
     enddo
     close(10)
  endif
  call MPI_finalize(ierr)
  stop 'VizGrid'
end subroutine VizGrid
