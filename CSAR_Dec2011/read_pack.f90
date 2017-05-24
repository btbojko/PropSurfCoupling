! *****************************************************************
Subroutine READ_DISCS

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

!-------------------------------------------------------------------------

! Local variables
  INTEGER ::  i,n,i5,mm, w_unit,npacks
  REAL*8  ::  area,area_tot,area_bind
  REAL*8  ::  area_ap,u3,v3,area_al
  REAL*8  ::  print_a_v,print_a_w
  REAL*8  ::  beta00,percent
  REAL*8  ::  scalefactor,radmax,ddlx,ddly,lengthscale
  REAL*8  ::  sum1,sum2,area_d,area_h,area_fdf
!-------------------------------------------------------------------------

  if (myid ==0) then
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'READING PACK'
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'****************************************************'
     print*,'****************************************************'
  endif

  open(14,file='pack2.txt')
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

  lengthscale = (ddlx/2.0d0) * 1.0d-4

  if (alpha.le.0.0d0) alpha=0.0d0

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


!-modi0

   alpha_test = 1.0  
!   alpha_test = 0.8

  do i=1,ncircles
     read(14,*) i5,mm,mm,xcenter(i),ycenter(i),u3,v3,rad(i)
     rad(i)     = rad(i)    *dlength
     xcenter(i) = xcenter(i)*dlength*alpha_test
     ycenter(i) = (ycenter(i)*dlength + period/2.0)/alpha_test
  end do
101 format(2x,i6,2x,4f12.6)
  close(15)

  dly = dlength/alpha_test 
  period = -2.0*dly

  dlength = dlength*alpha_test  
  dlx = dlength

!-------end modi0


!           now check for proper number of spheres inside domain
  do i=1,ncircles
     if (abs(xcenter(i)).le.dlx .and.   &
        &             ycenter(i).le.0.0d0 .and.&
        &             ycenter(i).ge.-2.0d0*dly) then
        sum2 = sum2 + 1
        area = area + pi*rad(i)*rad(i)
     endif
     radmax=dmax1(radmax,rad(i))
  end do
  close(14)

!
!*****************************************************************
!*****************************************************************
!
  area_tot = 4.0d0*dlx*dly
  area_d = area
  area_h = area_tot - area_d
  area_ap = area_d + alp_V*area_h          
  area_bind = (1.0d0-alp_V)*area_h
  alpha_fdf = rho_ap*area_d/(rho_premixed*area_h+rho_ap*area_d)
  alpha_pack = rho_ap*area_ap/(rho_binder*area_bind+rho_ap*area_ap)
!
!*****************************************************************
!*****************************************************************
!
  if(myid == 0)then
     write(6,*)'READING THE PACK ...'
     write(6,*)
     write(6,*)'dlength = ',dlength
     write(6,*)'period = ',period
     write(6,*)'aspect_ratio = ',aspect_ratio
     write(6,*)'Lx, Ly = ',dlx,dly
     write(6,*)'Sum = ',sum1,sum2
     write(6,*)'Max diameter (microns,----) = ',radmax*2*10000.0d0
     write(6,*) dlength,yend
     write(6,*) 'alpha_test = ', alpha_test
     write(6,*) 'alpha_pack = ', alpha_pack
     write(6,*) 'alpha_FDF = ', alpha_fdf    

  endif
113 format(2x,3f12.6)
!----------------------------------------------------------------------------
  RETURN
END SUBROUTINE READ_DISCS
!****************************************************************************

