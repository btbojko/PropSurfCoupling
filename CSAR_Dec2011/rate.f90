!     Set up rate for four step kinetics
!     f(i,k,j,1)  = T    (gas phase temperature)
!     f(i,k,j,2)  = Y_F  (gas phase fuel Y)
!     f(i,k,j,3)  = Y_O  (gas phase oxidizer X)
!     f(i,k,j,4)  = Y_Z1 (gas phase intermediate Z1 from AP)
!     f(i,k,j,5)  = Y_Z2 (gas phase intermediate Z2 from Binder)
!     f(i,k,j,6)  = T    (solid phase temperature)
! ********************************************************************
SUBROUTINE RATESETUP
  
USE GLOBAL_DATA
USE MYMPI
IMPLICIT NONE

!---------------------------------------------------------------------------
! Local Variables
INTEGER :: i, j, k,ys,ye,xs,xe,zs,ze
INTEGER :: m,n
REAL*8  :: pressure, poo, tsqinv
REAL*8  :: addrate(neqgas)
REAL*8  :: ratey(irxn), pre(irxn), abf(irxn,neqgas),dratey(irxn,neqgas)
REAL*8  :: rmax, allrmax
  
!---------------------------------------------------------------------------

! get the value of the pressure
poo = PRESSURE(tcyc)
  
! Determine range for each CPU  
ys = drange(5) + 1
ye = nyv(0)
xs = drange(1)
xe = nxv(0)
zs = drange(3)
ze = drange(4)

! Determine pre-exponential for each reaction
pre(1:irxn) = da(1:irxn)*poo**np(1:irxn)

do j = ys, ye
    do k = zs, ze
        do i = xs, xe
            tsqinv = one/(f(i,k,j,1)*f(i,k,j,1))   ! 1/Tgas^2
            do m=1,irxn       ! Do for reaction
                do n=2,neqgas      !Do for species
                    if(f(i,k,j,n) < 0.0d0) f(i,k,j,n) = 0.0d0
                    if(f(i,k,j,n) > 1.0d0) f(i,k,j,n) = 1.0d0
                    abf(m,n) = sign(1.0d0,f(i,k,j,n))* abs(f(i,k,j,n))**mu(m,n)
                enddo
            enddo
        
! Concentration and exponent terms of reaction rate
            ratey(1:irxn) = exp(-thetay(1:irxn)/f(i,k,j,1))
            do m = 1,irxn
                do n = 2,neqgas
                    ratey(m) = ratey(m)*abf(m,n)
                enddo
            enddo

!---  multiply by dt/rho and by the pressure dependent term
            ratey(1:irxn) =  dtx(i,k,j) * pre(1:irxn) * ratey(1:irxn) 

! Take derivatives for jacobain, derivative of 4 rxns with respect to T
            dratey(1:irxn,1) = thetay(1:irxn)*ratey(1:irxn)*tsqinv

            dratey(1:irxn,2:neqgas) = 0.0d0
  
            do m=1,irxn
                do n=2,neqgas
                    if(yf(n,m) < 0.0d0)then
                        if(f(i,k,j,n) == 0) then
                            dratey(m,n) = 0.0d0
                        else
                            dratey(m,n) = mu(m,n) * ratey(m) / f(i,k,j,n)
                        endif
                    endif
                enddo
            enddo
                                   
            do m = 1,neqgas
                do n = 1,neqgas
                    Bm(m,n,i,k,j) = - sum(yf(m,1:irxn) * dratey(1:irxn,n))  !yf(1,rxn) = Q/cp, yf(2:4,rxn) are rxn moles of each species
                enddo
                Bm(m,m,i,k,j) =  max(Bm(m,m,i,k,j),-one)
                addrate(m) = sum(Bm(m,1:neqgas,i,k,j)*f(i,k,j,1:neqgas))   !sum(Q*rate*species)
                dfdt(i,k,j,m) = dqdt(i,k,j,m) + sum(yf(m,1:irxn)*ratey(1:irxn)) + addrate(m)
            enddo

        end do
    end do
end do

if(itime_sos ==1) then
    do i = drange(1),nxv(0)
        do k = drange(3),drange(4)
            do j = drange(5)+1,drange(6)-1
                csld(10,i,k,j) = d_vterm(i,k) * (f(i,k,j+1,neqmax) - f(i,k,j-1,neqmax))*detady(j)
                dfdt(i,k,j,neqmax) = dqdt(i,k,j,neqmax) - csld(10,i,k,j)*f(i,k,0,neqmax)
            enddo
        enddo
    enddo
else
    dfdt(:,:,:,neqmax) = dqdt(:,:,:,neqmax)
endif


!...................................................................................
! EVALUATE THE JACOBIAN TERMS THAT APPEAR IN THE IMPLICIT BOUNDARY CONDITION
!...................................................................................
CALL BC_JAC


!--------------------------------------------------------------------
RETURN
END SUBROUTINE RATESETUP
!************************************************************************************


!*********************************************************************
SUBROUTINE RATE_EXP  !  [explicit rate]
USE GLOBAL_DATA
IMPLICIT NONE

!---------------------------------------------------------------------
! Local Variables
INTEGER :: i, j, k, m, n, ys,ye,xs,xe,zs,ze
REAL*8  :: pressure, poo
REAL*8  :: ratey(irxn), pre(irxn), abf(irxn,neqgas)
!---------------------------------------------------------------------

! get the value of the pressure
poo = PRESSURE(tcyc)

ys = drange(5) + 1
ye = nyv(0)
xs = drange(1)
xe = nxv(0)
zs = drange(3)
ze = drange(4)
  
!     Set up rate for two step kinetics
!     f(i,k,j,1)  = T    (gas phase temperature)
!     f(i,k,j,2)  = Y_F  (gas phase fuel Y)
!     f(i,k,j,3)  = Y_O  (gas phase oxidizer X)
!     f(i,k,j,4)  = Y_Z1 (gas phase intermediate Z1 from AP)
!     f(i,k,j,5)  = Y_Z2 (gas phase intermediate Z2 from Binder)
!     f(i,k,j,6)  = T    (solid phase temperature)

pre(1:irxn) = da(1:irxn)*poo**np(1:irxn)

do j = ys, ye
    do k = zs, ze
        do i = xs, xe
        
            do m=1,irxn
                do n=2,neqgas
                    if(f(i,k,j,n) < 0.0d0) f(i,k,j,n) = 0.0d0
                    if(f(i,k,j,n) > 1.0d0) f(i,k,j,n) = 1.0d0
                    abf(m,n) = sign(1.0d0,f(i,k,j,n))* abs(f(i,k,j,n))**mu(m,n)
                enddo
            enddo

            ratey(1:irxn) = pre(1:irxn)*exp(-thetay(1:irxn)/f(i,k,j,1))
            
            do m = 1,irxn
                do n = 2,neqgas
                    ratey(m) = ratey(m)*abf(m,n)
                enddo
            enddo
           
            do m = 1,neqgas-1
                rate(i,k,j,m) = sum(yf(m,1:irxn)*ratey(1:irxn))
            enddo
          
        end do
    end do
end do
!----------------------------------------------------------     
RETURN
END SUBROUTINE RATE_EXP
!**************************************************************************


!*********************************************************************
SUBROUTINE RATE_PRINT  !  [print rate]
USE GLOBAL_DATA
IMPLICIT NONE

!---------------------------------------------------------------------
! Local Variables
INTEGER :: i, j, m, n, ys,ye,xs,xe
REAL*8  ::  pressure, poo
REAL*8  :: ratey(irxn), pre(irxn), abf(irxn,neqgas)
!---------------------------------------------------------------------

! get the value of the pressure
poo = PRESSURE(tcyc)

ys = drange(5) + 1
ye = nyv(0)
xs = drange(1)
xe = nxv(0)

!     Set up rate for two step kinetics
!     f(i,k,j,1)  = T    (gas phase temperature)
!     f(i,k,j,2)  = Y_F  (gas phase fuel Y)
!     f(i,k,j,3)  = Y_O  (gas phase oxidizer X)
!     f(i,k,j,4)  = Y_Z  (gas phase intermediate Z)
!     f(i,k,j,5)  = T    (solid phase temperature)

pre(1:irxn) = da(1:irxn)*poo**np(1:irxn)

do j = 0, ny
    do i = 0, nx

        do m=1,irxn
            do n=2,neqgas
                if(G_f(i,j,n) < 0.0d0) G_f(i,j,n) = 0.0d0
                if(G_f(i,j,n) > 1.0d0) G_f(i,j,n) = 1.0d0
                abf(m,n) = sign(1.0d0,G_f(i,j,n))* abs(G_f(i,j,n))**mu(m,n)
            enddo
        enddo

        ratey(1:irxn) = pre(1:irxn)*exp(-thetay(1:irxn)/G_f(i,j,1))
            
        do m = 1,irxn
            do n = 2,neqgas
                ratey(m) = ratey(m)*abf(m,n)
            enddo
        enddo
        
        do m = 1,irxn
            G_rate(i,j,m) = yf(1,m)*ratey(m) * 1.d-3
        enddo

    end do
end do
!----------------------------------------------------------     
RETURN
END SUBROUTINE RATE_PRINT
!**************************************************************************


! ********************************************************************
!maybe it would be useful to use the old suroutine with just a flag
SUBROUTINE RATE_SPLIT

USE GLOBAL_DATA
USE MYMPI
IMPLICIT NONE

!---------------------------------------------------------------------------
! Local Variables
INTEGER :: i, j, k,ys,ye,xs,xe,zs,ze
INTEGER :: m,n
REAL*8  :: pressure, poo, tsqinv
REAL*8  :: addrate(neqgas)
REAL*8  :: ratey(irxn), pre(irxn), abf(irxn,neqgas),dratey(irxn,neqgas)
REAL*8  :: rmax, allrmax
!---------------------------------------------------------------------------

!-- get pressure
!  poo = PRESSURE(press,epi_p,omg_p,period_p,pi,tcyc)
poo = PRESSURE(tcyc)  !gross


ys = drange(5) + 1
ye = nyv(0)
xs = drange(1)
xe = nxv(0)
zs = drange(3)
ze = drange(4)

pre(1:irxn) = da(1:irxn)*poo**np(1:irxn)

do j = ys, ye
    do k = zs, ze
        do i = xs, xe

            if(iconv(i,k,j) ==1) CYCLE

            tsqinv = one/(f(i,k,j,1)*f(i,k,j,1))

            do m=1,irxn
                do n=2,neqgas
                    if(f(i,k,j,n) < 0.0d0) f(i,k,j,n) = 0.0d0
                    if(f(i,k,j,n) > 1.0d0) f(i,k,j,n) = 1.0d0
                    abf(m,n) = sign(1.0d0,f(i,k,j,n))* abs(f(i,k,j,n))**mu(m,n)
                enddo
            enddo

            ratey(1:irxn) = pre(1:irxn)*exp(-thetay(1:irxn)/f(i,k,j,1))
            
            do m = 1,irxn
                do n = 2,neqgas
                    ratey(m) = ratey(m)*abf(m,n)
                enddo
            enddo

            ratey(2) = max(ratey(2),0.0d0)

! Take derivatives for jacobain, derivative of 4 rxns with respect to T
            dratey(1:irxn,1) = thetay(1:irxn)*ratey(1:irxn)*tsqinv

            dratey(1:irxn,2:neqgas) = 0.0d0
  
            do m=1,irxn
                do n=2,neqgas
                    if(yf(n,m) < 0.0d0)then
                        if(f(i,k,j,n) == 0) then
                            dratey(m,n) = 0.0d0
                        else
                            dratey(m,n) = mu(m,n) * ratey(m) / f(i,k,j,n)
                        endif
                    endif
                enddo
            enddo

            do m = 1,neqgas
                do n = 1,neqgas
                    Bm(m,n,i,k,j) = sum(yf(m,1:irxn) * dratey(1:irxn,n))
                enddo
                rate(i,k,j,m) = sum(yf(m,1:irxn)*ratey(1:irxn))        
            enddo
        end do
    end do
end do
!--------------------------------------------------------------------
RETURN
END SUBROUTINE RATE_SPLIT
!************************************************************************************ 
! ********************************************************************
SUBROUTINE RATE_CANC
  
USE GLOBAL_DATA
USE MYMPI
IMPLICIT NONE

!---------------------------------------------------------------------------
! Local Variables
INTEGER :: i, j, k,ys,ye,xs,xe,zs,ze
INTEGER :: m,n
REAL*8  :: pressure, poo, tsqinv
REAL*8  :: addrate(neqgas)
REAL*8  :: ratey(irxn), pre(irxn), abf(irxn,neqgas),dratey(irxn,neqgas)
REAL*8  :: rmax, allrmax
REAL*8  :: rhogas,coe_rho
  
!---------------------------------------------------------------------------

! get the value of the pressure
poo = PRESSURE(tcyc)
coe_rho = dim_fact*poo
 
ys = drange(5) + 1
ye = nyv(0)
xs = drange(1)
xe = nxv(0)
zs = drange(3)
ze = drange(4)

pre(1:irxn) = da(1:irxn)*poo**np(1:irxn)

do j = ys, ye
    do k = zs, ze
        do i = xs, xe

            tsqinv = one/(f(i,k,j,1)*f(i,k,j,1))

            rhogas = coe_rho/f(i,k,j,1)	
            dtx(i,k,j) = 1.d-6/rhogas
            do m=1,irxn
                do n=2,neqgas
                    if(f(i,k,j,n) < 0.0d0) f(i,k,j,n) = 0.0d0
                    if(f(i,k,j,n) > 1.0d0) f(i,k,j,n) = 1.0d0
                    abf(m,n) = sign(1.0d0,f(i,k,j,n))* abs(f(i,k,j,n))**mu(m,n)
                enddo
            enddo
            
! Concentration and exponent terms of reaction rate
            ratey(1:irxn) = exp(-thetay(1:irxn)/f(i,k,j,1))
            do m = 1,irxn
                do n = 2,neqgas
                    ratey(m) = ratey(m)*abf(m,n)
                enddo
            enddo

!---  multiply by dt/rho and by the pressure dependent term
            ratey(1:irxn) =  dtx(i,k,j) * pre(1:irxn) * ratey(1:irxn) 

! Take derivatives for jacobain, derivative of 4 rxns with respect to T
            dratey(1:irxn,1) = thetay(1:irxn)*ratey(1:irxn)*tsqinv

            dratey(1:irxn,2:neqgas) = 0.0d0
  
            do m=1,irxn
                do n=2,neqgas
                    if(yf(n,m) < 0.0d0)then
                        if(f(i,k,j,n) == 0) then
                            dratey(m,n) = 0.0d0
                        else
                            dratey(m,n) = mu(m,n) * ratey(m) / f(i,k,j,n)
                        endif
                    endif
                enddo
            enddo
            
            do m = 1,neqgas
                do n = 1,neqgas
                    Bm(m,n,i,k,j) = - sum(yf(m,1:irxn) * dratey(1:irxn,n))
                enddo
                Bm(m,m,i,k,j) =  max(Bm(m,m,i,k,j),-one)
                addrate(m) = sum(Bm(m,1:neqgas,i,k,j)*f(i,k,j,1:neqgas))
                rate(i,k,j,m) = sum(yf(m,1:irxn)*ratey(1:irxn))
            enddo

        end do
    end do
end do
!--------------------------------------------------------------------
RETURN
END SUBROUTINE RATE_CANC
!************************************************************************************

