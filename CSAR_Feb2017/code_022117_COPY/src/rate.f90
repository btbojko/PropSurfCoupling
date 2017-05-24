! ********************************************************************
! ********************************************************************
! Updated: TLJ; 1/17/2017
! Filename: rate.f90
! ********************************************************************
! ********************************************************************
!     ________________________________________________________________
!
!     BRIEF CODE DESCRIPTION:
!     =======================
!
!     This set of subroutines computes the rates for kinetics
!
!     subroutines: RATESETUP
!
! ********************************************************************
!
!     Set up rate for three step kinetics with aluminum
!        f(i,k,j,1)  = T    (gas phase temperature)
!        f(i,k,j,2)  = Y_O  (gas phase oxidizer X)
!        f(i,k,j,3)  = Y_Z  (gas phase intermediate Z)
!        f(i,k,j,4)  = Y_F  (gas phase fuel Y)
!        f(i,k,j,5)  = Y_al (gas phase aluminum Y_al)
!        f(i,k,j,6)  = Y_ox (gas phase oxide Y_oxide)
!        f(i,k,j,7)  = T    (solid phase temperature)
!
!        maxsize  = 7    number of gas+solid phase variables
!        neqgas   = 6    number of gas phase variables
!        irxn     = 4    number of reaction rates
!
!     Set up rate for three step kinetics without aluminum
!        f(i,k,j,1)  = T    (gas phase temperature)
!        f(i,k,j,2)  = Y_O  (gas phase oxidizer X)
!        f(i,k,j,3)  = Y_Z  (gas phase intermediate Z)
!        f(i,k,j,4)  = Y_F  (gas phase fuel Y)
!        f(i,k,j,5)  = T    (solid phase temperature)
!
!        maxsize  = 5    number of gas+solid phase variables
!        neqgas   = 4    number of gas phase variables
!        irxn     = 3    number of reaction rates
!
! ********************************************************************

SUBROUTINE RATESETUP

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE
!
!---------------------------------------------------------------------
! Local Variables
  INTEGER :: i, j, k,ys,ye,xs,xe,zs,ze
  INTEGER :: l,m,n,inn1,inn2
  INTEGER :: main_oxidizer
  REAL*8  :: pressure, poo, tsqinv
  REAL*8  :: ratey(maxsize), pre(maxsize),pre_spec(maxsize) 
  REAL*8  :: abf(maxsize,maxsize)
  REAL*8  :: dratey(maxsize,maxsize)
  REAL*8  :: rmax, allrmax,gy
  REAL*8  :: eps,fact_sign
  REAL*8  :: expHermsen, expDiamL, expTempG, expPresG 
  REAL*8  :: expXiEffG, xiCO2, xiO2, xiH2O, xiH2
  REAL*8  :: xiEffG, dOxH2RdOxnoH2, diffRel,expMF,fff
  REAL*8  :: sss(maxsize),TLIMIT_ALflame(2)
!---------------------------------------------------------------------------

! get the value of pressure
  poo = PRESSURE(tcyc(2))

! set the main oxidizer for aluminum product
  main_oxidizer = 2

  eps=1.0d-60

! For ultrafine aluminum in the binder
  if (NORADIATION.eqv..FALSE.) then
     expHermsen = 1.9d0
     expDiamL   = 3.0d0-expHermsen   
     expMF      = expDiamL/3d0 - 1d0

     expTempG  = 1.57d0   
     expPresG  = 0.20d0   
     expXiEffG = 0.39d0   

     xiCO2 = 0.20d0   
     xiO2  = 0.02d0   
     xiH2O = 0.20d0   
     xiH2  = 0.20d0

     xiEffG = ( xiO2 + 0.58d0 * xiH2O + 0.22d0 * xiCO2 )   

     dOxH2RdOxnoH2 = 3.5d0   

     diffRel = 1.0d0 + xiH2 * ( dOxH2RdOxnoH2 - 1.0d0 )

     da(irxn) = 2.636d0/expHermsen*Initial_Diameter**(-expHermsen)&
           &*xiEffG**(expXiEffG)*diffRel*beta_al**expMF
     da(irxn) = da(irxn)*KrierFACTOR  !kreir

     if (NOQ4) da(irxn) = 0.0d0

     thetay(irxn) = 0d0
     np(irxn) = expPresG
     zetay(irxn) = expTempG

     mu(irxn,1:maxsize) = 0d0 ! initialize to zero
     mu(irxn,main_oxidizer) = expXiEffG   ! use the final product
                                          ! main oxydating species Y_O
     mu(irxn,neqgas-1) = 1d0+expMF
     mu(irxn,neqgas) = -expMF

  endif

! for aluminum, the combustion is turned on when
! the temperature reaches the melt temperature
!   T_mely_oxide = 2350 K
! and the combustion is turned off when
!   T > T_boil = 4100 K
! See equation (64) of JPP 2008 paper

  if(ncyc > NCYCRAD) then
     TLIMIT_ALflame = (/2350d0,4.0d3/)   !cancel
  else
     TLIMIT_ALflame = (/2050d0,4.0d3/)
  end if

! Determine range for each CPU
  ys = drange(5) + 1
  ye = drange(6) - 1
  xs = drange(1)
  xe = drange(2)
  zs = drange(3)
  ze = drange(4)

! Determine pre-exponential for each reaction
  pre(1:irxn) = da(1:irxn)*poo**np(1:irxn)

  dratey = 0d0
  ratey = 0d0

  do j = ys, ye
  do k = zs, ze
  do i = xs, xe

     tsqinv = one/(f(i,k,j,1)*f(i,k,j,1))

     abf = 1d0
     do m=1,irxn     ! Do for reaction steps
        pre_spec(m) = 1d0
        do n=2,neqgas     ! Do for species

           if (f(i,k,j,n) < 0.0d0) f(i,k,j,n) = 0.0d0
           if (f(i,k,j,n) > 1.0d0) f(i,k,j,n) = 1.0d0

           if( abs(mu(m,n)) <= 1d-12 ) cycle

           abf(m,n) = f(i,k,j,n)**mu(m,n)

           ! This is for aluminum combustion
           if (NORADIATION.eqv..FALSE.) then
           if(m == irxn) then 
               if( n == neqgas) then
                   fff = f(i,k,j,neqgas-1)*beta_al + f(i,k,j,neqgas)
               elseif(n == main_oxidizer) then
                   fff = min(max(1d0 -sum( f(i,k,j,2:neqgas)),0d0),1d0)
               else
                   fff = f(i,k,j,n)
               end if
               if(abs(fff) > 1d-12) then
                   abf(m,n) =  sign(1.0d0,fff)  * abs(fff)**mu(m,n)
               else
                   abf(m,n) = 0d0
               end if
           end if
           end if
           pre_spec(m) = pre_spec(m)*abf(m,n)
        enddo
     enddo
!
     do m=1,irxn

!---    calculate reaction rates save the pressure dependent term
        ratey(m) = pre_spec(m)* exp(-thetay(m)/f(i,k,j,1)) * f(i,k,j,1)**zetay(m) 

!---    needed for output only to calculate total heat release
        rate_out(i,k,j,m) =  pre(m) * ratey(m) 

!---    multiply by dt/rho and by the pressure dependent term
        ratey(m) =  dtx(i,k,j) * pre(m) * ratey(m) 

     end do

!--- Take derivatives for Jacobain, derivative of rxns wrt T
     dratey = 0d0
     dratey(1:irxn,1) = thetay(1:irxn)*ratey(1:irxn)&
                *tsqinv + zetay(1:irxn)*ratey(1:irxn)/f(i,k,j,1)

!--- Take derivatives for Jacobain, derivative of rxns wrt Y
     dratey(1,2) =  mu(1,2) *ratey(1)/merge(sign(eps,f(i,k,j,2)),&
                f(i,k,j,2),f(i,k,j,2)==0.0d0)&
                   *merge(2d0,1d0,f(i,k,j,2)<0d0 .and. mu(1,2) < 1d0)

     if (NORADIATION.eqv..FALSE.) then
     dratey(1,4) =  -mu(1,4) *ratey(1) &
                /merge(eps,one-f(i,k,j,4),f(i,k,j,4)==one)
     do m = 2,irxn
        do n = 2,neqgas
           mugt0if: if(abs(mu(m,n)) <= 1d-12 ) then
           else
              sss = 0d0
              if(m==irxn .and. n==main_oxidizer) then
                          fff = min(max(1d0 -sum( f(i,k,j,2:neqgas)),0d0),1d0)
                          sss =-1
                          inn1 = 2
                          inn2 = neqgas
              elseif(m==irxn .and. n==neqgas) then
                          fff = f(i,k,j,neqgas-1)*beta_al + f(i,k,j,neqgas)
                          sss(neqgas-1:neqgas) = (/beta_al,1d0/)
                          inn1 = neqgas-1
                          inn2 = neqgas
              else
                          fff = f(i,k,j,n)
                          sss = 1
                          inn1 = n
                          inn2 = n
              end if
              dratey(m,inn1:inn2) = dratey(m,inn1:inn2) + &
                   &sss(inn1:inn2) * mu(m,n) * ratey(m) /&
                   &merge(sign(eps,fff),fff,fff==0.0d0)&
                   &*merge(2d0,1d0,fff<0d0 .and. mu(m,n) < 1d0)
           end if mugt0if
        enddo
     enddo
     endif
     if (NORADIATION) then
     do m = 2,irxn
        do n = 2,neqgas
           if(abs(mu(m,n)) <= 1d-12 ) then
           else
              fff = f(i,k,j,n)
              sss = 1
              inn1 = n
              inn2 = n
              dratey(m,inn1:inn2) = dratey(m,inn1:inn2) + &
                   &sss(inn1:inn2) * mu(m,n) * ratey(m) /&
                   &merge(sign(eps,fff),fff,fff==0.0d0)&
                   &*merge(2d0,1d0,fff<0d0 .and. mu(m,n) < 1d0)
           end if
        enddo
     enddo
     endif

     if (NORADIATION.eqv..FALSE.) then
        ratey(irxn) = ratey(irxn) * rhog(i,k,j)
        dratey(irxn,:) = dratey(irxn,:) * rhog(i,k,j)
        if(f(i,k,j,1)<TLIMIT_ALflame(1) .or. &
                &f(i,k,j,1)>TLIMIT_ALflame(2)) then
            ratey(irxn) = 0d0
            dratey(irxn,:) = 0d0
        end if
     end if
!
     do n = 1,neqgas
        do m = 1,neqgas
           Bm(m,n,i,k,j) = - sum(yf(m,1:irxn) * dratey(1:irxn,n))
        enddo
        Bm(n,n,i,k,j) = max(Bm(n,n,i,k,j),-one)
        !Bm(n,n,i,k,j) = max(Bm(n,n,i,k,j),-two)

        rate(i,k,j,n) = sum(yf(n,1:irxn)*ratey(1:irxn))
     enddo

  end do
  end do
  end do
!
  do j = drange(5),drange(6) 
  do k =drange(3),drange(4)
  do i =drange(1),drange(2)   
     gy = (f(i,k,j+1,neqmax)-f(i,k,j-1,neqmax))
     csld(10,i,k,j) = d_vterm(i,k)*gy*detady(j)
     ! Add radiation contribution to energy equation
     rate(i,k,j,1) = rate(i,k,j,1) + dtorhocp(i,k,j)*radheat(i,k,j)
  enddo
  enddo
  enddo

!
!  Boundary Eqauations Jacobians is using implicit_bc
!  CALL BC_JAC
  if(allocated(dfdt))then
!--- Boundary Eqauations Jacobians
     CALL BC_JAC
  end if



!---------------------------------------------------------------------
  RETURN
END SUBROUTINE RATESETUP
!*********************************************************************
