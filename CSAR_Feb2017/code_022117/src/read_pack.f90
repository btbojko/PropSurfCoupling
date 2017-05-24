! *****************************************************************
Subroutine SAND_PACK (packtime)

  USE GLOBAL_DATA
  USE radiation
  USE MYMPI

  IMPLICIT NONE

!-------------------------------------------------------------------------
! Dummy variables
  REAL*8  ::  packtime

! Local variables
  INTEGER ::  nn, i, i1,i2,i3, mm, n
  INTEGER ::  ii,imax, w_unit,k
  INTEGER ::  sum1, sum2, sum_K
  INTEGER ::  nmode,Ncomments,myploc
  REAL*8  ::  vol,vol_tot,vol_ap,vol_bind,vol_pack
  REAL*8  ::  ddlx,ddly,ddlz
  REAL*8  ::  rvtot_K,rvexp_K
  REAL*8  ::  ymax,zmax,xmax,radmax
  REAL*8  ::  th,beta_inp,diff_betas
  REAL*8  ::  dummy1,dummy2,dummy3
  REAL*8  ::  prdx,prdy,prdz
  character(5) :: dummychar
  character(90) comments(500)
!-------------------------------------------------------------------------

  ipack = 0
  dlength = xend
  dly = dlength 

  !call RADIATION_hardwire
  !CALL BLEND((/ts_AP_H,ts_BD_H,ts_alum_H/),(/lambda_ap,lambda_binder,lambda_al/),&
  !     &(/rho_ap,rho_binder,rho_al/),lambda_eff,rho_eff,3)
  
  alpha_pack = beta/(1.0d0+beta)

! xloc = equivalence_ratio/   &
!      ((1.0-alp_V)*(1.0d0+beta*rho_binder/rho_ap))
!  xloc = xloc*dlength
! if (NORADIATION.eqv..FALSE.) then
!    xloc = (1d0-ts_AP)*dlength
! else
!    xloc = 130.0d0/2.0d0/1.d4
!    xloc = 50.0/2.0/1.e4
! endif
! BinderThickness = two*xloc
  if(equivalence_ratio < 10) then
     xloc = (1d0-ts_AP)*dlength
     ipack = 0
  else
     xloc = equivalence_ratio*1.d-4
     ipack = 0
  endif
!    xloc = 130.0d0/2.0d0/1.d4
    xloc = 50.0/2.0/1.e4
  BinderThickness = two*xloc
  write(6,*)'Sandwich Properties'
  write(6,*) xloc,BinderThickness

  packtime = MPI_WTIME() - packtime

       
!----------------------------------------------------------------------------
  RETURN
END SUBROUTINE SAND_PACK


!****************************************************************************
! *****************************************************************
Subroutine READ_DISCS

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!-------------------------------------------------------------------------

! Local variables
  INTEGER ::  i,n,i5,mm, w_unit,npacks,ii,imax
  REAL*8  ::  area,area_tot,area_bind,dummys
  REAL*8  ::  area_ap,u3,v3,area_pack,area_al,totmass
  REAL*8  ::  print_a_v,print_a_w
  REAL*8  ::  beta00,percent,alpha_test,aspect_ratio
  REAL*8  ::  scalefactor,radmax,ddlx,ddly
  REAL*8  ::  sum1,sum2,area_d,area_h
  REAL*8  ::  ddlz,prdx,prdy
  REAL*8  ::  rmax,xmax,ymax
!-------------------------------------------------------------------------

  if (myid ==0) then
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'READING 2D PACK'
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'****************************************************'
  endif

  is2D = .true.

! open(14,file='pack2D.txt')
  open(14,file='pack200.txt')
! open(14,file='pack2_0.70_200.txt')
!
!*****************************************************************
!*****************************************************************
!
  read(14,*) ddlx,aspect_ratio
  read(14,*) rho1
  read(14,*) npacks,ncircles

  ddly = ddlx*aspect_ratio
!
! Scale factor; scale domain by this value
  scalefactor = 1.0d0
  ddlx = ddlx*scalefactor
  ddly = ddly*scalefactor
  ddlz = 1.0

  lengthscale = (ddlx/2.0d0) * 1.0d-4

  dlx = lengthscale
  dly = lengthscale*dly
  dlz = 1.0
  dlength = dlx
  period = -2.0*dly

  xend = dlx
  zend = dlz
  write(6,*)'xend=',xend

! this is for the SUBGRID CODE  ! ! compute alpha, the precent of blend that is ap
  alpha = (0.79 - rho1)/(1.0d0-rho1)
  alpha = 0.78 - rho1
  alpha_2dpack = alpha
  alpha=0.0
  if (alpha.le.0.0d0) alpha=0.0d0

  ts_AP = rho1
  ts_alum = 0.0
  ts_BD = 1d0 - ts_alum - ts_AP

  ts_AP_H = alpha
  if (alpha.le.0.0) ts_AP_H = 0.0
  ts_alum_H = 0.0
  ts_BD_H = 1d0 - ts_alum_H - ts_AP_H

  totmass = rho_ap*ts_AP_H+rho_binder*ts_BD_H+rho_al*ts_alum_H
  ws_AP_H = ts_AP_H*rho_ap/totmass
  ws_alum_H = ts_alum_H*rho_al/totmass
  ws_BD_H = 1d0 - ws_alum_H-ws_AP_H


  alp_V(1:3) = (/ts_AP_H,ts_alum_H,ts_BD_H/)
  alp_W(1:3) = (/ws_AP_H,ws_alum_H,ws_BD_H/)
  

!-----------------------------
!allocation for the pack variables
  if(ncircles <= 0) then
     write(*,*)'ERROR in reading PACK information'
     write(*,*)'Ncircles < 0,',ncircles
  else
     CALL ALLOCATE_PACK
  endif
!-----------------------------
!
!
!*****************************************************************
!*****************************************************************
!
!           here I am rescaling domain
   dlx = lengthscale
   dly = lengthscale*aspect_ratio
   period = -2.0*dly
   dlength = dlx
!
!     divide dlength by 2 to get half length so that x = [-dlength,dlength]
!     scale radius by dlength to get dimensional values
!     the cube is mapped into -dl < x < dl; 0 < y < -2*dl
  pi = acos(-1.0d0)
  area = 0.0d0
  sum1 = 0
  sum2 = 0
  radmax = 0.0d0

   alpha_test = 1.0  
!   alpha_test = 0.8

  print *,'before the loop ...'
  do i=1,ncircles
     read(14,*) mm,mm,mm,ycenter(i),xcenter(i),rad(i)
     rad(i)     = rad(i)    *dlength
     xcenter(i) = xcenter(i)*dlength*alpha_test
     ycenter(i) = (ycenter(i)*dlength + period/2.0)/alpha_test

     write(19,101) i,xcenter(i),xcenter(i),ycenter(i),rad(i)
  end do
101 format(2x,i6,2x,4f12.6)
  close(14)
  close(15)

  if(myid == 0)then
     write(6,*)'READING THE PACK ...'
     write(6,*)
     write(6,*)'dlength = ',dlength
     write(6,*)'period = ',period
     write(6,*)'aspect_ratio = ',aspect_ratio
     write(6,*)'Lx, Ly = ',dlx,dly
     write(6,*)'Sum = ',sum1,sum2
     write(6,*)'Max diameter (microns,----) = ',radmax*2*10000.0d0
     write(6,*)'Packing fraction = ',area/(4.0*dlx*dly)
     write(6,*) dlength,yend
  
     write(6,*) 'alpha_test = ', alpha_test

  endif
113 format(2x,3f12.6)
!
!*****************************************************************
!*****************************************************************
!
  area_tot = 4.0d0*dlx*dly
  area_d = area
  area_h = area_tot - area_d
  area_ap = area_d + alpha*area_h
  area_bind = (1.0d0-alpha)*area_h
  area_pack = area_ap/area_tot
  beta_pack = rho_ap*area_ap/(rho_binder*area_bind)
!
!*****************************************************************
!*****************************************************************
!
!
!          
  print_a_v = alpha_pack
  print_a_w = print_a_v/(print_a_v+(1.d0-print_a_v)*rho_binder/rho_ap)

  alpha_pack = print_a_w

  !alp_V=0d0
  !alp_W=0d0
  !if(iwflag == 0) then
  !   alp_V(1) = alpha
  !   alp_W(1) = alp_V(1)/(alp_V(1)+(1.d0-alp_V(1))*rho_binder/rho_ap)
  !else 
  !   alp_W(1) = alpha
  !   alp_V(1) = alp_W(1)/(rho_ap/rho_binder*(1.0d0-alp_W(1))+alp_W(1))
  !endif
  !alp_V(2)=0.d0
  !alp_V(3)=1.d0
  !alp_W(1)=0.d0
  !alp_W(2)=0.d0
  !alp_W(3)=1.d0
  !write(6,*)'alpha=',alpha
  !write(6,*)'iwflag=',iwflag
  !write(6,*)'here', alp_V(1:3)
  !write(6,*)'here', alp_W(1:3)

  ii=0
  n = 0
  do i=1,ncircles
     ymax = -10.0
     n = n+1
     do ii=n,ncircles
        if(ycenter(ii).ge.ymax) then
          ymax = ycenter(ii)
          xmax = xcenter(ii)
          rmax = rad(ii)
          imax = ii
        endif
     enddo
     ycenter(imax)=ycenter(i)
     xcenter(imax)=xcenter(i)
     rad(imax)=rad(i)
     ycenter(i)=ymax
     xcenter(i)=xmax
     rad(i)=rmax
  enddo
  write(6,*)'ncircles=',ncircles
  do i=1,ncircles
     write(26,101) i,xcenter(i),ycenter(i)
  enddo
  close(26)


!----------------------------------------------------------------------------
  RETURN
END SUBROUTINE READ_DISCS
!****************************************************************************

!****************************************************************************
!****************************************************************************
Subroutine READ_3DPACK (packtime)

  USE GLOBAL_DATA
  USE MYMPI
  USE radiation

  IMPLICIT NONE

  !-------------------------------------------------------------------------
  ! Dummy variables
  REAL*8  ::  packtime

  ! Local variables
  INTEGER ::  nn, i, i1,i2,i3, mm, n
  INTEGER ::  ii,imax, w_unit,k
  INTEGER ::  sum1, sum2, sum_K
  INTEGER ::  nmode,Ncomments
  INTEGER ::  inunit
  REAL*8  ::  vol,vol_tot,vol_ap,vol_bind,vol_pack
  REAL*8  ::  ddlx,ddly,ddlz
  REAL*8  ::  rvtot_K,rvexp_K
  REAL*8  ::  ymax,zmax,xmax,radmax
  REAL*8  ::  th,beta_inp,diff_betas
  REAL*8  ::  dummy1,dummy2,dummy3
  REAL*8  ::  prdx,prdy,prdz,alp_V_ele,ts_Mt_HT
  REAL*8  ::  ss,rhoth,totmass
  character(5 ) :: dummychar
  character(90) :: comments(500),cctrim
  character(8) :: elstr,elstrL
  !-------------------------------------------------------------------------

  issymmetric = .FALSE.
  is2D = .FALSE.

  if (myid ==0) then
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'READING 3D PACK'
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'****************************************************'
  endif

  packtime = MPI_WTIME()

  pi = ACOS(-1.0d0)

! Read in pack
  inunit = 14
  open (inunit,file=pack_inputfile_name)

  !--- READ PACK HEADER
  if(myid == 0) write(*,*)'Reading pack'

  commloop: do i =1,500
     read(inunit,'(a)')comments(i)
     cctrim = adjustL(comments(i))
     if(myid == 0 .and. len_trim(comments(i)) > 5) write(*,*) trim(comments(i)),k
     if(len_trim(comments(i)) > 5 .and. cctrim(1:3) /= '***') then
        read(cctrim,*) elstr,alp_V_ele
        elstrL = adjustL(elstr)
        if(elstrL(1:2) == 'al' .or.  elstrL(1:2) == 'Al' .or.  elstrL(1:2) == 'AL') then
           ts_alum_HT = max(alp_V_ele,0d0)
        elseif(elstrL(1:2) == 'ap' .or.  elstrL(1:2) == 'Ap' .or.  elstrL(1:2) == 'AP') then
           ts_AP_HT = max(alp_V_ele,0d0)
        elseif(index(elstr,'Binder') > 0) then
           ts_BD_HT = max(alp_V_ele,0d0)
           ts_Mt_HT = ts_BD_HT+ts_AP_HT+ts_alum_HT
           ts_alum_H = ts_alum_HT/ts_Mt_HT
           ts_AP_H   = ts_AP_HT/ts_Mt_HT
           ts_BD_H   = 1.0 - ts_AP_H - ts_alum_H
           ts_alum = ts_alum_HT
           ts_BD   = ts_BD_HT
           Ncomments = i
           sCASE = 'PACK'
           if(myid == 0)print*,'EXITING commloop'
           exit commloop
        endif
     endif
  enddo commloop

  if (NORADIATION.eqv..FALSE.) call RADIATION_hardwire
  CALL BLEND((/ts_AP_H,ts_BD_H,ts_alum_H/),(/lambda_ap,lambda_binder,lambda_al/),&
       &(/rho_ap,rho_binder,rho_al/),lambda_eff,rho_eff,3)

  totmass = rho_ap*ts_AP_H+rho_binder*ts_BD_H+rho_al*ts_alum_H
  ws_AP_H = ts_AP_H*rho_ap/totmass
  ws_alum_H = ts_alum_H*rho_al/totmass
  ws_BD_H = 1d0 - ws_alum_H-ws_AP_H

  alp_V(1:3) = (/ts_AP_H,ts_alum_H,ts_BD_H/)
  alp_W(1:3) = (/ws_AP_H,ws_alum_H,ws_BD_H/)

  if (myid==0) then
   write(6,*) 'Hom fractions=',ts_AP_H,ts_BD_H,ts_alum_H
   write(6,*) 'lambdas=',lambda_ap,lambda_binder,lambda_al
   write(6,*) 'densities=',rho_ap,rho_binder,rho_al
   write(6,*) 'lambda_eff=',lambda_eff
   write(6,*) 'rho_eff=',rho_eff
   write(6,*) 'alp_V=',alp_V(1:3)
   write(6,*) 'alp_W=',alp_W(1:3)
  endif

  
  
  read(inunit,*) ddlx
  read(inunit,*) ddlz
  read(inunit,*) ddly
  read(inunit,*) ncircles,ncirclesIN
  read(inunit,*) rhoV_Rocpack
  read(inunit,*) TheoreticalPackingDensity !theoretical value
  read(inunit,*) nmode
  do i = 1, nmode
     read(inunit,*) nn,dummychar,mm,dummy1
  enddo


  ! note that ddlx = ddlz always
  ! here I am fixing the depth of solid to be 1000 microns
  ! and then converting to cm
  ! lengthscale = (500.0d0/2.0d0) * 1.0d-4
  ! here I am rescaling domain

  ddlx = ddlx*1d-4/2.d0
  ddlz = ddlz*1d-4/2.d0
  ddly = ddly*1d-4/2.d0
  lengthscale = ddlx

  ! divide dlength by 2 to get half length so that x = [-dlength,dlength]
  ! scale radius by dlength to get dimensional values
  ! Caution: the 3d packing code assumes burning in the z-direction,
  ! and so prints out x,y,z; however, this 3d combustion
  ! code burns in the y-direction, and so reads in x,z,y

  dlx = ddlx
  dly = ddly
  dlz = ddlz
  dlength = lengthscale
  prdx = ddlx/dlength
  prdy = ddly/dlength
  prdz = ddlz/dlength



  alpha = alp_V(1)
  xend = dlx
  zend = dlz
  period = -2.0d0*dly
  periodXZ = 2.0d0*dlx   !dlx and dlz should be the same

  CALL ALLOCATE_PACK

  !--- READ PACK SPHERES

  sum1 = 0;sum2 = 0;vol =0
  do i=1,ncircles
     read(inunit,*) i1,i2,OxidizerType(i),xcenter(i),&
                    zcenter(i),ycenter(i),rad(i)


     !       Mapping II
     !       y lies in -2*prdy < y < 0
     !
     ycenter(i) = ycenter(i)-prdy

     !       now check for proper number of spheres inside domain
     !
     if (abs(xcenter(i)).le.prdx .and.&
          abs(zcenter(i)).le.prdz .and.&
          ycenter(i).le.0.0d0 .and.&
          ycenter(i).ge.-2.0d0*prdy) then
        sum1 = sum1 + 1
     endif

     !       Mapping III
     !       scale everything by lengthscale
     !
     xcenter(i) = xcenter(i)*lengthscale
     ycenter(i) = ycenter(i)*lengthscale
     zcenter(i) = zcenter(i)*lengthscale
     rad(i)     = rad(i)    *lengthscale
     if(myid == 0) write(19,100) i,xcenter(i),zcenter(i),ycenter(i),rad(i)

     !       now check for proper number of spheres inside domain
     !
     if (abs(xcenter(i)).le.dlx .and.&
          abs(zcenter(i)).le.dlz .and.&
          ycenter(i).le.0.0d0 .and.&
          ycenter(i).ge.-2.0d0*dly) then
        sum2 = sum2 + 1
        vol = vol + 4.0d0*pi*rad(i)*rad(i)*rad(i)/3.0d0
     endif
  end do

  close(inunit)

  vol_tot  = 8.0d0*dlx*dly*dlz
  vol_ap   = vol
  vol_bind = vol_tot - vol_ap
  beta_pack = (vol_ap/vol_bind)
  rhoH_pack = beta_pack/(1.0d0+beta_pack)
  vol_ap  = vol_ap + alp_V(1)*vol_bind
  vol_bind = ts_BD_H*vol_bind
  vol_pack = vol_ap/vol_tot
  beta_pack = (vol_ap/vol_bind)
  rhoH_pack = beta_pack/(1.0d0+beta_pack)

  if (myid == 0) then
    write(6,*)'vol_ap,vol_bind=',vol_ap*1e12,vol_bind*1e12
    write(6,*)'packing fraction=',vol_ap/vol_tot
    write(6,*)'beta_pack=', beta_pack
    write(6,*)'alp_V(1)=', alp_V(1)
    write(6,*)'vol_ap,vol_bind=',vol_ap*1e12,vol_bind*1e12
  endif

  !
  !    SORT the circles against ycenter
  !
  ii=0
  n = 0
  do i=1,ncircles
     ymax=-10.0
     n=n+1
     do ii=n,ncircles
        if(ycenter(ii).ge.ymax) then
           ymax=ycenter(ii)
           zmax=zcenter(ii)
           xmax=xcenter(ii)
           radmax=rad(ii)
           imax=ii
        endif
     enddo
     ycenter(imax)=ycenter(i)
     zcenter(imax)=zcenter(i)
     xcenter(imax)=xcenter(i)
     rad(imax)=rad(i)
     ycenter(i)=ymax
     zcenter(i)=zmax
     xcenter(i)=xmax
     rad(i)=radmax
  enddo


  if(myid == 0) then
     write(*,'(80("*")/,a,/,20("-"))')'output from READPACK:'
     write(6,*)'SUMS',sum1,sum2,ncirclesIN
     write(6,*)'Total volume     = ',vol_tot*1.0d12
     write(6,*)'AP volume        = ',vol_ap*1.0d12
     write(6,*)'Binder volume    = ',vol_bind*1.0d12
     write(6,*)'AP vol frac      = ',vol_ap/vol_tot
     write(6,*)'Bind vol frac    = ',vol_bind/vol_tot
     write(6,*)'Packing fraction = ',beta_pack
     write(6,*)'Pack density     = ',rhoH_pack
     write(6,*)'Calculated beta  = ',(vol_ap/vol_bind)*rho_ap/rho_binder
     write(6,*)'Homogenization alphas:alp_V = ',alp_V(1:3)
     write(6,*)'Homogenization alphas:alp_W = ',alp_W(1:3)
     write(*,*)'PERIOD == ',period
     write(*,'(20("-"),/a/,80("*"))')'END READPACK:'
  endif
100 format(2x,i5,2x,4f15.8)

  packtime = MPI_WTIME() - packtime


  !----------------------------------------------------------------------------
  RETURN
END SUBROUTINE READ_3DPACK
!****************************************************************************


