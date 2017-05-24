Module radiation
  USE data_types
  use global_data
  USE MYMPI

  implicit none
  SAVE
  real*8 :: etastar,C_GAMMA,C_K1,C_K2,C_A1,ctvlims(2)
  real*8 :: C_rho_oxide,C_rho_alum,C_D_oxide,C_D_alum,C_D_oxideT,C_D_alumT,C_alpha_oxide,C_alpha_alum,C_eps_oxide,C_eps_alum
  real*8 :: C_rho_AP,C_rho_BD,C_alpha_AP,C_alpha_BD,C_eps_BD,C_eps_AP
  real*8 :: C_sigma,C_sigopi,C_fourpi,C_pi
  real*8 :: G_inf
  logical :: set_constant = .false.
  real*8,parameter :: R2_limit=1d-3,mic2cm = 1d-4
  character(LEN=1) :: BCTYPE_RAD,TVAR

contains
  subroutine S_etastar
    implicit none
    integer :: i,k,j,jj
    real*8,allocatable :: sumterm(:),allterm(:)

!    allocate(sumterm(ny+1),allterm(lbound(sumterm,1):ubound(sumterm,1)))
    allocate(sumterm(ny+1))
	allocate(allterm(lbound(sumterm,1):ubound(sumterm,1)))

    sumterm=0d0
    do j = 1,ny-1
       sumterm(j) = sum(rate(drange(1):drange(2),drange(3):drange(4),j,2))
    end do
    sumterm(ny) = (drange(2)-drange(1)+1)*(drange(4)-drange(3)+1)
    sumterm(ny+1) = sum(f(drange(1):drange(2),drange(3):drange(4),ny,1))

    CALL MPI_ALLREDUCE(sumterm,allterm,size(sumterm),MPI_DOUBLE_PRECISION,MPI_SUM,comm3d,ierr)
    allterm = allterm/allterm(ny)


    G_inf = merge(F_Ginf(finest_mesh%c2f(0,0,0,ny),allterm(ny+1)),0d0,BCTYPE_RAD == 'D')

    jj = 0d0
    do j = 1,ny-1
       if(allterm(j)>R2_limit) jj = j
    end do
    etastar = y(jj)
    call RADIATION_hardwire

    etastar = .1d0/4d0!cancel!
!!>    etastar = 1d99

    return
  end subroutine S_etastar
!
  subroutine constants
    implicit none
    real*8,parameter :: joule2cal = 0.238845897
    if(set_constant) return

    C_K1 = 1d0
    C_K2 = 0.015d0*press/13d0  !assuming the pressure in atmospheres
!
    C_alpha_alum  = 1d-1
    C_alpha_oxide = 0.45d0
    C_rho_alum    = 0.90d0
    C_rho_oxide   = 0.55d0
    C_eps_alum    = 1d0
    C_eps_oxide   = .45d0
    C_D_alumT      = 1d0  *mic2cm
    C_D_oxideT     = 0.3d0*mic2cm
    C_D_alum      = 1d0 !for aluminum found from particle model
    C_D_oxide     = 1d0 
!
    C_rho_AP = 6d-1
    C_rho_BD = 3d-1
    C_alpha_AP = 1d0-C_rho_AP
    C_alpha_BD = 1d0-C_rho_BD
    C_eps_BD = C_alpha_AP
    C_eps_AP = C_alpha_BD
!
    C_GAMMA = 1d2!/mic2cm
    C_A1 = 0d0
    C_sigma = 5.67d-12*joule2cal
    C_sigopi = C_sigma / acos(-1d0)
    C_fourpi = 4d0* acos(-1d0)
    C_pi = acos(-1d0)

    ctvlims = (/0.7d0,1.3d0/)

    BCTYPE_RAD = 'N'
    TVAR = 'T'
!!>    sCASE = 'M28'
!!>    sCASE = 'none'

    set_constant = .true.
  end subroutine constants
!---------------------------------

  subroutine RADIATION_hardwire
    real*8 :: totmass

    call constants

    if(sCASE == 'M28') then
       ts_AP  = 52.29d0/1d2
       ts_alum = 4.79d0/1d2
    elseif(sCASE == 'M29') then
       ts_AP   = 57.8d0/1d2
       ts_alum = 4.23d0/1d2
    elseif(sCASE == 'baseline') then
       ts_AP  = 62.25d0/1d2
       ts_alum = 13.32d0/1d2
    elseif(sCASE == 'none') then
       ts_AP  = 69d0/1d2
       ts_alum = 0d0
    end if
    

    if(sCase /= 'PACK') then
       ts_BD = 1d0 - ts_alum-ts_AP
       if(ipack == 0) then
          ts_AP_H = 0d0
          ts_AP = equivalence_ratio*(1d0-ts_alum)/(1.0d0+beta*rho_binder/rho_ap)
          ts_alum_H = ts_alum/(1d0-ts_AP)
          ts_BD = 1d0 - ts_alum-ts_AP
       elseif(ipack == 1) then
          ts_AP_H = 0d0
          ts_alum_H = ts_alum/(ts_BD+ts_alum)*(1d0-ts_AP_H)
       else
          ts_AP_H =   ts_AP
          ts_alum_H = ts_alum
       end if
    end if
    ts_BD_H = 1d0 - ts_alum_H-ts_AP_H
       

    totmass = rho_ap*ts_AP_H+rho_binder*ts_BD_H+rho_al*ts_alum_H
    ws_AP_h = ts_AP_H*rho_ap/totmass
    ws_alum_H = ts_alum_H*rho_al/totmass
    ws_BD_H = 1d0 - ws_alum_H-ws_AP_H

    alp_V(1:3) = (/ts_AP_H,ts_alum_H,ts_BD_H/)
    alp_W(1:3) = (/ws_AP_H,ws_alum_H,ws_BD_H/)


    beta_al = 27d0/102d0

    Initial_Diameter = 3d0   !this is in micron

    return
  end subroutine RADIATION_hardwire
!---------------------------------
  double precision function t_alum(t)
    type (coarse_to_fine) :: t
    integer :: i,k,j,m
    i = t%indx(1)
    k = t%indx(2)
    j = t%indx(3)
    m=1
    if(is2d) then
       k = 0
       t_alum = Vod(i,k,j,m)*t%alpha(1)*t%alpha(3)  &
            & + Vod(i+1,k,j,m)*(1d0-t%alpha(1))*t%alpha(3)&
            & + Vod(i,k,j+1,m)*t%alpha(1)*(1d0-t%alpha(3))&
            & + Vod(i+1,k,j+1,m)*(1d0-t%alpha(1))*(1d0-t%alpha(3))
    else
       t_alum = Vod(i,k,j,m)*t%alpha(1)*t%alpha(2)*t%alpha(3)  &
            & + Vod(i+1,k,j,m)*(1d0-t%alpha(1))*t%alpha(2)*t%alpha(3)&
            & + Vod(i,k,j+1,m)*t%alpha(1)*t%alpha(2)*(1d0-t%alpha(3))&
            & + Vod(i+1,k,j+1,m)*(1d0-t%alpha(1))*t%alpha(2)*(1d0-t%alpha(3))&
            & + Vod(i,k+1,j,m)*t%alpha(1)*(1d0-t%alpha(2))*t%alpha(3)  &
            & + Vod(i+1,k+1,j,m)*(1d0-t%alpha(1))*(1d0-t%alpha(2))*t%alpha(3)&
            & + Vod(i,k+1,j+1,m)*t%alpha(1)*(1d0-t%alpha(2))*(1d0-t%alpha(3))&
            & + Vod(i+1,k+1,j+1,m)*(1d0-t%alpha(1))*(1d0-t%alpha(2))*(1d0-t%alpha(3))
    end if
    t_alum = max(t_alum,1d-12/C_D_alumT)
    
  end function t_alum
!--------------------------------
  double precision function t_oxide(t)
    type (coarse_to_fine) :: t
    integer :: i,k,j,m
    i = t%indx(1)
    k = t%indx(2)
    j = t%indx(3)
    m=2
    if(is2d) then
       k = 0
       t_oxide = Vod(i,k,j,m)*t%alpha(1)*t%alpha(3)  &
            & +  Vod(i+1,k,j,m)*(1d0-t%alpha(1))*t%alpha(3)&
            & +  Vod(i,k,j+1,m)*t%alpha(1)*(1d0-t%alpha(3))&
            & +  Vod(i+1,k,j+1,m)*(1d0-t%alpha(1))*(1d0-t%alpha(3))
    else
       t_oxide = Vod(i,k,j,m)*t%alpha(1)*t%alpha(2)*t%alpha(3)  &
            & +  Vod(i+1,k,j,m)*(1d0-t%alpha(1))*t%alpha(2)*t%alpha(3)&
            & +  Vod(i,k,j+1,m)*t%alpha(1)*t%alpha(2)*(1d0-t%alpha(3))&
            & +  Vod(i+1,k,j+1,m)*(1d0-t%alpha(1))*t%alpha(2)*(1d0-t%alpha(3))&
            & +  Vod(i,k+1,j,m)*t%alpha(1)*(1d0-t%alpha(2))*t%alpha(3)  &
            & +  Vod(i+1,k+1,j,m)*(1d0-t%alpha(1))*(1d0-t%alpha(2))*t%alpha(3)&
            & +  Vod(i,k+1,j+1,m)*t%alpha(1)*(1d0-t%alpha(2))*(1d0-t%alpha(3))&
            & +  Vod(i+1,k+1,j+1,m)*(1d0-t%alpha(1))*(1d0-t%alpha(2))*(1d0-t%alpha(3))
    end if
    t_oxide = max(t_oxide,1d-12/C_D_oxideT)

  end function t_oxide
!--------------------------------
  double precision function F_ka(t)
    type (coarse_to_fine) :: t
    F_Ka = 3d0*(t_oxide(t)*C_alpha_oxide/C_D_oxide + C_K1*t_alum(t)*C_alpha_alum/C_D_alum)
    F_Ka = F_Ka+F_kp(t)
  end function F_ka
!--------------------------------
  double precision function F_ks(t)
    type (coarse_to_fine) :: t
    F_Ks = 3d0*(t_oxide(t)*C_rho_oxide/C_D_oxide + C_K1*t_alum(t)*C_rho_alum/C_D_alum)
  end function F_ks
!--------------------------------
  double precision function F_kem(t)
    type (coarse_to_fine) :: t
    F_Kem = 3d0*(t_oxide(t)*C_eps_oxide/C_D_oxide + C_K1*t_alum(t)*C_eps_alum/C_D_alum)
!!>    F_Kem = max(F_Kem,F_kp(t))
    F_Kem = F_Kem+F_kp(t)
  end function F_kem
!--------------------------------
  double precision function F_kp(t)
    type (coarse_to_fine) :: t
    F_Kp = 0.3d0*press
  end function F_kp
!--------------------------------
  double precision function F_omega(t)
    type (coarse_to_fine) :: t
    if(F_ks(t)+F_ka(t) > 1d-12) then
       F_omega = F_ks(t)/(F_ks(t)+F_ka(t))
    else
       F_omega = 0d0
    end if
  end function F_omega
!--------------------------------
  double precision function F_Ib(t)
    real*8 :: t
    F_Ib = C_sigopi * t**4
  end function F_Ib
!--------------------------------
  double precision function F_BC(eta,psi,temp,flag)
    type (coarse_to_fine) :: eta
    real*8 :: psi,temp
    integer :: flag
    real*8 :: aa,ee,term1,term2

    aa = merge(C_alpha_AP,C_alpha_BD,psi>0)
    ee = merge(C_eps_AP,C_eps_BD,psi>0)

    term1 = (F_ks(eta)+F_ka(eta))*(3d0-C_A1*F_omega(eta))/(aa-2d0)/2d0
    if(flag == 1) then
       term2 = C_fourpi*ee*F_Ib(temp)
    elseif(flag == 2) then
       term2 = -aa
    else 
       stop 'function F_BC'
    end if
    F_BC = term1*term2
    return
  end function F_BC
!--------------------------------
  double precision function F_HBC(eta,psi,temp,flag)
    type (coarse_to_fine) :: eta
    real*8 :: psi,temp
    integer :: flag
    real*8 :: aa,ee,term1,term2

    aa = merge(C_alpha_AP,C_alpha_BD,psi>0)
    ee = merge(C_eps_AP,C_eps_BD,psi>0)

    term1 = .5d0/(2d0-aa)
    if(flag == 1) then
       term2 = C_fourpi*ee*F_Ib(temp)
    elseif(flag == 2) then
       term2 = -aa
    else 
       stop 'function F_BC'
    end if
    F_HBC = term1*term2
    return
  end function F_HBC
!--------------------------------
  double precision function F_LAP_RHS_PRE(t)
    type (coarse_to_fine) :: t
    F_LAP_RHS_PRE = (F_ks(t) + F_ka(t)) * (3d0-C_A1*F_omega(t))
    return
  end function F_LAP_RHS_PRE
!--------------------------------
  double precision function F_Ginf(eta,temp)
    type (coarse_to_fine) :: eta
    real*8 :: psi,temp
    integer :: flag
    real*8 :: aa,ee,term1,term2
!
    if(F_ka(eta) > 1d-20) then
       F_Ginf = C_fourpi*F_kem(eta)*F_Ib(temp)/F_ka(eta)
    else
       F_Ginf = 0d0
    end if
    return
  end function F_GINF
!--------------------------------
  double precision function F_P(t)
    type (coarse_to_fine) :: t
    F_P = -(F_ks(t) + F_ka(t)) * (3d0-C_A1*F_omega(t))
    return
  end function F_P
!--------------------------------
!*************************************************************************

!Laplacian coeffs


!*************************************************************************
  SUBROUTINE RADIATION_SOURCE
    USE mg_solver


    IMPLICIT NONE

!-------------------------------------------------------------------------
! Local variables:
    TYPE(mesh_pointers), POINTER :: current_mesh
    TYPE(cppevec), DIMENSION(:,:), POINTER :: CPPE
    REAL(KIND=double), DIMENSION(:,:), POINTER :: ccmg
    REAL(KIND=double), DIMENSION(:), POINTER :: ey,eya,eta,etaa
    REAL(KIND=double), DIMENSION(:,:), POINTER :: chi_mg,chipr_mg

    INTEGER :: xs,zs,ys,lbb(4)
    INTEGER :: xe,ze,ye,i,k,j,mnx,mnz,m,mm,mmm,ipm,kpm

    REAL(KIND=double) :: dxf,dzf,dyf,dxc,dzc,dyc,dxa,dza,dya
    REAL(KIND=double) :: dxa2,dza2,dy2a,dxdy4,dzdy4
    REAL(KIND=double) :: term,term_m1,c379,sqterm,Aterm,Bterm
    REAL(KIND=double) :: dya4,dyxz4,dytermP(2),chterm(2),chtermP(2),zzz(3)
    REAL(KIND=double) :: MQt(5),fxza(5),coetv(7)
    type (coarse_to_fine) :: t,tv(7)
!--------------------------------------------------------------------------

    call constants
    call S_etastar
    if(max(ts_alum,ts_alum_H) <= 1d-12 .or. NORADIATION) then
       radheat = 0d0
       BCradheat = 0d0
       return
    end if

    divt = 0d0
!
    xs = drange(1)
    xe = drange(2)
    zs = drange(3)
    ze = drange(4)
    ys = drange(5)+1
    ye = drange(6)-1
!
    dya = 1d0/dy
    dxa = 1d0/dx
    dza = 1d0/dz
!
    m = 0
    do k=zs,ze
       do i=xs,xe
          do j=ys,ye
             t = finest_mesh%c2f(0,i,k,j)
             divt(i,k,j) =  (F_Ka(t)*G_inf-C_fourpi*F_kem(t)*F_Ib(f(i,k,j,1)))*F_LAP_RHS_PRE(t)
             if(.not.abs(divt(i,k,j)) < 1d88) print*,i,k,j,myid, divt(i,k,j),&
                  &F_Ka(t),G_inf,C_fourpi,F_kem(t),F_Ib(f(i,k,j,1)),F_LAP_RHS_PRE(t),f(i,k,j,1),t%indx,t%alpha,'VODs',&
                  &Vod(t%indx(1),0,t%indx(3),1:2),Vod(t%indx(1)+1,0,t%indx(3),1:2),&
                  &Vod(t%indx(1),0,t%indx(3)+1,1:2),Vod(t%indx(1)+1,0,t%indx(3)+1,1:2),'VFRACS',&
                  &Vfrac(t%indx(1),0,t%indx(3),1:2),Vfrac(t%indx(1)+1,0,t%indx(3),1:2),&
                  &Vfrac(t%indx(1),0,t%indx(3)+1,1:2),Vfrac(t%indx(1)+1,0,t%indx(3)+1,1:2),'DIAMPs',&
                  &DiamP(t%indx(1),0,t%indx(3),1:2),DiamP(t%indx(1)+1,0,t%indx(3),1:2),&
                  &DiamP(t%indx(1),0,t%indx(3)+1,1:2),DiamP(t%indx(1)+1,0,t%indx(3)+1,1:2),'F',&
                  &F(t%indx(1),0,t%indx(3),5:6),F(t%indx(1)+1,0,t%indx(3),5:6),&
                  &F(t%indx(1),0,t%indx(3)+1,5:6),F(t%indx(1)+1,0,t%indx(3)+1,5:6)
          end do
          j = ys
          t = finest_mesh%c2f(0,i,k,j)
          term_m1 = dphidx_mg(i,k,m)**2+dphidz_mg(i,k,m)**2
          term = 1d0+term_m1
          sqterm = sqrt(term)
          Bterm = F_BC(t,psi_mg(i,k,m),temp_mg(i,k,m),1) &
               &+ F_BC(t,psi_mg(i,k,m),temp_mg(i,k,m),2)*G_inf
          divt(i,k,j) = divt(i,k,j) + sqterm*Bterm*detady(j)*dya
       end do
    end do

    if( .not. all(abs(divt) < 1d88)) then
       do k=zs,ze
          do i=xs,xe
             do j=ys,ye
                if(.not.abs(divt(i,k,j)) < 1d88) print*,i,k,j,myid, divt(i,k,j)
             end do
          end do
       end do
    end if


    divt(:,:,ny) = divt(:,:,ny-1)  
    divt(:,:,0) = zero


    divt=cshift(divt,1,3)
    p = Gfield

    finest_mesh%x => p
    finest_mesh%f => divt


    current_mesh => finest_mesh
    dxf = finest_mesh%dx ; dzf = finest_mesh%dz; dyf = finest_mesh%dy

    DO

       xe = current_mesh%xe ; ze = current_mesh%ze; ye = current_mesh%ye
       dxc = current_mesh%dx ; dzc = current_mesh%dz; dyc = current_mesh%dy
       mnx = anint(dxc/dxf); mnz = anint(dzc/dzf);
       m = anint( log10 (dble(mnx)) / log10 (2.0) )

       dxa = 1.0/dxc ; dza = 1.0/dzc; dya = 1d0/dyc
       dxa2 = dxa*dxa; dza2 = dza*dza; dy2a = dya*dya
       dxdy4 = dxa*dya/4.0d0;dzdy4 = dza*dya/4.0d0
       dya4 = dya/4.0d0;


       mcomp   = 1
       mgradjM = 0
       mgradj  = 2

       eta      => current_mesh%y
       etaa     => current_mesh%ya
       ey       => current_mesh%ey 
       eya      => current_mesh%eya
       chi_mg   => current_mesh%chi
       chipr_mg => current_mesh%chipr
       CPPE     => current_mesh%cppe
!
!   b1=c1
!
       current_mesh%c1 = -dxa2

       iloop: do i=0,xe
          kloop: do k=0,ze


             ccmg => cPPE(i,k)%ccmg

             jloop: DO j = 0,ye

                do ipm=-1,1
                   MQterm(i+ipm,k,mcomp)   = (one +  chipr_mg(j,1)  *  phi_mg(i+ipm,k,m) )
                   MQterm(i+ipm,k,mgradjM) = (one +  chipr_mg(j-1,2)*  phi_mg(i+ipm,k,m) )
                   MQterm(i+ipm,k,mgradj)  = (one +  chipr_mg(j,2)  *  phi_mg(i+ipm,k,m) )
                end do
                do kpm=-1,1,2
                   MQterm(i,k+kpm,mcomp)   = (one +  chipr_mg(j,1)  *  phi_mg(i,k+kpm,m) )
                   MQterm(i,k+kpm,mgradjM) = (one +  chipr_mg(j-1,2)*  phi_mg(i,k+kpm,m) )
                   MQterm(i,k+kpm,mgradj)  = (one +  chipr_mg(j,2)  *  phi_mg(i,k+kpm,m) )
                enddo


                dyxz4 = ey(j)*dya4*dxa

                tv(1) = finest_mesh%c2f(0,i,k,j)
                tv(2) = finest_mesh%c2f(1,i-1,k,j)
                tv(3) = finest_mesh%c2f(1,i,k,j)
                tv(4) = finest_mesh%c2f(2,i,k-1,j)
                tv(5) = finest_mesh%c2f(2,i,k,j)
                tv(6) = finest_mesh%c2f(3,i,k,j-1)
                tv(7) = finest_mesh%c2f(3,i,k,j)

                do mm = 2,7
                   coetv(mm) = (2d0-F_P(tv(mm))/F_P(tv(1)))
                   coetv(mm) = max(min(coetv(mm),ctvlims(2)),ctvlims(1))
                end do

                MQt(2) = MQterm(i-1,k,mcomp)*coetv(2)
                MQt(3) = MQterm(i+1,k,mcomp)*coetv(3)
                MQt(4) = MQterm(i,k-1,mcomp)*coetv(4)
                MQt(5) = MQterm(i,k+1,mcomp)*coetv(5)

                fxza(2) = dphidxa_mg(i-1,k,m)*coetv(2)
                fxza(3) = dphidxa_mg(i,k,m)*coetv(3)
                fxza(4) = dphidza_mg(i,k-1,m)*coetv(4)
                fxza(5) = dphidza_mg(i,k,m)*coetv(5)

                if(j == 0) then

                   term_m1 = dphidx_mg(i,k,m)**2+dphidz_mg(i,k,m)**2
                   term = 1d0+term_m1
                   sqterm = sqrt(term)
                   Aterm = F_BC(tv(1),psi_mg(i,k,m),temp_mg(i,k,m),2)


                   dyterm = ey(j)*  (/ eya(j-1)*coetv(6), eya(j)*coetv(7) /)  *dya**2
                   chterm = (/ chi_mg(j-1,mgradj)*coetv(6), chi_mg(j,mgradj)*coetv(7) /)


                   ccmg(2,j) = current_mesh%c1*MQterm(i-1,k,mcomp) + &
                        dyxz4*(2d0*dphidxa_mg(i-1,k,m)*chi_mg(j,mcomp) - dphidx_mg(i,k,m) )

                   ccmg(3,j) = current_mesh%c1*MQterm(i+1,k,mcomp) +  &
                        dyxz4*(-2d0*dphidxa_mg(i,k,m)*chi_mg(j,mcomp) + dphidx_mg(i,k,m) )

                   ccmg(4,j) = current_mesh%c1*MQterm(i,k-1,mcomp) + &
                        dyxz4*(2d0*dphidza_mg(i,k-1,m)*chi_mg(j,mcomp) - dphidz_mg(i,k,m) )

                   ccmg(5,j) = current_mesh%c1*MQterm(i,k+1,mcomp) + &
                        dyxz4*(-2d0*dphidza_mg(i,k,m)*chi_mg(j,mcomp) + dphidz_mg(i,k,m) )



                   ccmg(6,j) = zero
                   ccmg(7,j) = ( -term_m1/MQterm(i,k,mgradj)*chi_mg(j,mgradj)*chi_mg(j+1,mcomp) - &
                        one/MQterm(i,k,mgradj)  )*dyterm(2) + 2d0*(dphidxa_mg(i,k,m)-dphidxa_mg(i-1,k,m) +&
                        & dphidza_mg(i,k,m)-dphidza_mg(i,k-1,m) ) *chi_mg(j+1,mcomp)*dyxz4 -&
                        0.5d0*sqterm*Aterm*ey(j)*dya


                   ccmg(8,j) = zero
                   ccmg(9,j) = zero
                   ccmg(10,j) =( 2d0*dphidxa_mg(i,k,m)*chi_mg(j+1,mcomp)  +(2d0*chi_mg(j,mgradj)/&
                        &MQterm(i,k,mgradj)*MQterm(i+1,k,mgradj)-1d0)*dphidx_mg(i,k,m))*dyxz4
                   ccmg(11,j) =(-2d0*dphidxa_mg(i-1,k,m)*chi_mg(j+1,mcomp)-(2d0*chi_mg(j,mgradj)/&
                        &MQterm(i,k,mgradj)*MQterm(i-1,k,mgradj)-1d0)*dphidx_mg(i,k,m))*dyxz4

                   ccmg(12,j) = zero
                   ccmg(13,j) = zero
                   ccmg(14,j) =( 2d0*dphidza_mg(i,k,m)*chi_mg(j+1,mcomp) +(2d0*chi_mg(j,mgradj)/&
                        MQterm(i,k,mgradj)*MQterm(i,k+1,mgradj)-1d0)*dphidz_mg(i,k,m))*dyxz4
                   ccmg(15,j) =(-2d0*dphidza_mg(i,k-1,m)*chi_mg(j+1,mcomp) -(2d0*chi_mg(j,mgradj)/&
                        MQterm(i,k,mgradj)*MQterm(i,k-1,mgradj)-1d0)*dphidz_mg(i,k,m))*dyxz4

                   ccmg(1,j) = ( (-dphidx_mg(i,k,m)**2-dphidz_mg(i,k,m)**2)*chi_mg(j,mgradj)*&
                        &chi_mg(j,mcomp) - one  )*dyterm(2)/MQterm(i,k,mgradj) &
                        &- 4d0* MQterm(i, k, mcomp)*dxa**2 &
                        &-2d0*(-dphidxa_mg(i,k,m)+dphidxa_mg(i-1,k,m)-dphidza_mg(i,k,m)+dphidza_mg(i,k-1,m))*&
                        &chi_mg(j,mcomp)*dyxz4  -&
                        1.5d0*sqterm*Aterm*ey(j)*dya


                else

                   dyterm = ey(j)*  (/ eya(j-1)*coetv(6), eya(j)*coetv(7) /)  *dya**2
                   chterm = (/ chi_mg(j-1,mgradj)*coetv(6), chi_mg(j,mgradj)*coetv(7) /)

                   ccmg(2,j) = current_mesh%c1*MQt(2) - &
                        &(chterm(2)  *MQterm(i-1,k,mgradj) /MQterm(i,k,mgradj) &
                        &- chterm(1)*MQterm(i-1,k,mgradjM)/MQterm(i,k,mgradjM))&
                        &*dphidx_mg(i,k,m)*dyxz4
                   ccmg(3,j) = current_mesh%c1*MQt(3) + &
                        &( chterm(2)  *MQterm(i+1,k,mgradj) /MQterm(i,k,mgradj) &
                        & -chterm(1)*  MQterm(i+1,k,mgradjM)/MQterm(i,k,mgradjM))&
                        &*dphidx_mg(i,k,m)*dyxz4

                   ccmg(4,j) = current_mesh%c1*MQt(4) - &
                        &(chterm(2)  *MQterm(i,k-1,mgradj) /MQterm(i,k,mgradj) &
                        & -chterm(1)*MQterm(i,k-1,mgradjM)/MQterm(i,k,mgradjM))&
                        &*dphidz_mg(i,k,m)*dyxz4
                   ccmg(5,j) = current_mesh%c1*MQt(5) + &
                        &(chterm(2)  *MQterm(i,k+1,mgradj) /MQterm(i,k,mgradj) -&
                        &  chterm(1)*MQterm(i,k+1,mgradjM)/MQterm(i,k,mgradjM))&
                        &*dphidz_mg(i,k,m)*dyxz4

                   ccmg(6,j) = ( (-dphidx_mg(i,k,m)**2-dphidz_mg(i,k,m)**2)/MQterm(i,k&
                        &,mgradjM)*chi_mg(j-1,mgradj)*chi_mg(j-1,mcomp) - one/MQterm(i&
                        &,k,mgradjM)  )*dyterm(1) + (-fxza(3)+fxza(2) -&
                        & fxza(5)+fxza(4))*chi_mg(j-1,mcomp)*dyxz4
                   ccmg(7,j) = ( (-dphidx_mg(i,k,m)**2-dphidz_mg(i,k,m)**2)/MQterm(i,k&
                        &,mgradj)*chi_mg(j,mgradj)*chi_mg(j+1,mcomp) - one/MQterm(i&
                        &,k,mgradj)  )*dyterm(2) + (fxza(3)-fxza(2) +&
                        & fxza(5)-fxza(4) ) *chi_mg(j+1,mcomp)*dyxz4
                   ccmg(8,j) = (fxza(2)*chi_mg(j-1,mcomp) + chterm(1)&
                        &/MQterm(i,k,mgradjM)*MQterm(i-1,k,mgradjM)*dphidx_mg(i,k,m))*dyxz4
                   ccmg(9,j) = (-fxza(3)*chi_mg(j-1,mcomp)-chterm(1)&
                        &/MQterm(i,k,mgradjM)*MQterm(i+1,k,mgradjM)*dphidx_mg(i,k,m))*dyxz4
                   ccmg(10,j) =( fxza(3)*chi_mg(j+1,mcomp)+chterm(2)/&
                        &MQterm(i,k,mgradj)*MQterm(i+1,k,mgradj)*dphidx_mg(i,k,m))*dyxz4
                   ccmg(11,j) =(-fxza(2)*chi_mg(j+1,mcomp)-chterm(2)&
                        &/MQterm(i,k,mgradj)*MQterm(i-1,k,mgradj)*dphidx_mg(i,k,m))*dyxz4
                   ccmg(12,j) =(chterm(1)/MQterm(i,k,mgradjM)*&
                        &MQterm(i,k-1,mgradjM)*dphidz_mg(i,k,m)+fxza(4)*chi_mg(j-1,mcomp))*dyxz4
                   ccmg(13,j) =(-chterm(1)/MQterm(i,k,mgradjM)*MQterm(i,k+1,mgradjM)&
                        &*dphidz_mg(i,k,m)-fxza(5)*chi_mg(j-1,mcomp))*dyxz4
                   ccmg(14,j) =(chterm(2)/MQterm(i,k,mgradj)*MQterm(i,k+1,mgradj)&
                        &*dphidz_mg(i,k,m)+fxza(5)*chi_mg(j+1,mcomp))*dyxz4
                   ccmg(15,j) =(-chterm(2)/MQterm(i,k,mgradj)*MQterm(i,k-1,mgradj)&
                        &*dphidz_mg(i,k,m)-fxza(4)*chi_mg(j+1,mcomp))*dyxz4
                   ccmg(1,j) = sum(ccmg(2:15,j))
                endif
             ENDDO jloop
             do j = 0,ye
                t = current_mesh%c2f(0,i,k,j)
                ccmg(1,j) = ccmg(1,j) - F_LAP_RHS_PRE(t)*F_Ka(t)
                ccmg(6:7,j) = - ccmg(6:7,j)
             enddo
          enddo kloop
       enddo iloop

       IF ( .NOT. ASSOCIATED(current_mesh%coarse_mesh) ) EXIT
       current_mesh => current_mesh%coarse_mesh
    END DO



    CALL fmg_solve(BCTYPE_RAD)   !CALL THE MULTIGRID SOLVER


!

    Gfield = p
    divt=cshift(divt,-1,3)
    p=cshift(p,-1,3)

    p = p+G_inf

!
!p should conatin the radiation intensity solution
!
    m = 0
    do j = 1,ny-1
       do k = drange(3),drange(4)
          do i = drange(1),drange(2)
             t = finest_mesh%c2f(0,i,k,j-1)
             radheat(i,k,j) = -C_fourpi*F_kem(t)*F_Ib(f(i,k,j,1)) + F_ka(t)*p(i,k,j)
          enddo
       enddo
    enddo
!!>
!!>    if(myid == 0) print'(1p5234e16.8)',y(1:ny-1)
!!>    if(myid == 0) print'(1p5234e16.8)',radheat(0,0,1:ny-1)
!!>    if(myid == 0) print'(1p5234e16.8)',p(0,0,1:ny-1)
!!>    if(myid == 0) print'(1p5234e16.8)',f(0,0,1:ny-1,1)
!!>    if(myid == 0) print'(1p5234e16.8)',divt(0,0,1:ny-1)
!!>    if(myid == 0) print'(a,1p1234e16.8)','G_inf',G_inf
!!>    stop



    j=1;
    m = 0
    do i = drange(1),drange(2)
       do k = drange(3),drange(4)
          t = finest_mesh%c2f(0,i,k,0)
          BCradheat(i,k) = -(F_HBC(t,psi_mg(i,k,m),temp_mg(i,k,m),1)+&
               & F_HBC(t,psi_mg(i,k,m),temp_mg(i,k,m),2) * p(i,k,j))*&
               &sqrt(1.0d0+dphidx_mg(i,k,m)**2+dphidz_mg(i,k,m)**2)
       end do
    end do

!---------------------------------------------------------------------
    RETURN
  END SUBROUTINE RADIATION_SOURCE
!**********************************************************************
end Module radiation
