c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c   
c
      subroutine w2v
      include 'rocprepack_3.h'
c
c     local variables
c
      integer i,j,indx,n,np
      double precision a_local,b_local,d_local
      dimension a_local(500,500),b_local(500),indx(6)
c
c     This subroutine converts weight fraction to volume
c        fraction using the appropriate densities
c
c     initialize arrays to zero
      indx = 0
      a_local = 0
      b_local = 0
c
      do i=1,iblend_max
         do j=1,iblend_max
            a_local(i,j) = w_frac(i)*(density(j)-density_binder)
         enddo
         a_local(i,i) = a_local(i,i) - density(i)
         b_local(i) = -w_frac(i)*density_binder
      enddo
c
      n = iblend_max
      np = n
      call ludcmp(a_local,n,np,indx,d_local)
      call lubksb(a_local,n,np,indx,b_local)
c
      rho_v = sum(b_local)
      rho_w = sum(w_frac)
      write(606,*)
      do i=1,iblend_max
         v_frac(i) = b_local(i)
         write(606,101) i,v_frac(i),w_frac(i)
      enddo
      write(606,101) iblend_max+1,1.0-sum(v_frac),1.0-sum(w_frac)
 101  format(2x,'Cut, Volume, Weight = ',i3,2x,2f12.8)
      write(606,*)
      write(606,*)' Total Weight Solids Fraction = ',rho_w
      write(606,*)' Total Volume Solids Fraction = ',rho_v
      write(606,*)' Total Volume Binder Fraction = ',1.0-rho_v
      write(6,*)' Total Weight Solids Fraction = ',rho_w
      write(6,*)' Total Volume Solids Fraction = ',rho_v
      write(6,*)' Total Volume Binder Fraction = ',1.0-rho_v
      write(606,*)
c
      return
      end
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c   
      subroutine lognormal(iblend)
      include 'rocprepack_3.h'
c
c     local variables
c
      integer i,nBin,j,iblend
      double precision diamsLog,diams,weights,
     &   diamMin,diamMax,diamPeak,power,
     &   diamMinLog,diamMaxLog,diamDeltaLog,ddd1,ddd2,
     &   acoeff,weightsSum

      dimension diamsLog(500), diams(500), weights(500)
c
c     This subroutine computes the log-normal for an initial
c     distribution of AP or aluminum given the range
c     [D_min,D_max], peak D_peak, number of bins
c
c     Print Input Data
c
      diamMin = diam_Min(iblend)
      diamMax = diam_Max(iblend)
      diamPeak = diam_Peak(iblend)
      nBin = 80
      nBin = 12
      nBin = 20
c     nBin = 200
      power = 3.0d0

      WRITE(606,*) 'Minimum Diameter, diamMin  = ', diamMin
      WRITE(606,*) 'Maximum Diameter, diamMax  = ', diamMax
      WRITE(606,*) 'Peak    Diameter, diamPeak = ', diamPeak
      WRITE(606,*) 'B-Coefficient,    bcoeff   = ', bcoeff(iblend)
      WRITE(606,*) 'Number of Bins,   nBin     = ', nBin
      WRITE(606,*) 'Power,            power    = ' ,power

c     Determine equidistant distribution in Log space

      diamMinLog = LOG(diamMin)
      diamMaxLog = LOG(diamMax)

      diamDeltaLog =  (diamMaxLog-diamMinLog)/float(nBin-1)

      diamsLog(1) = diamMinLog
      do i = 2, nBin
         diamsLog(i) = diamsLog(i-1) + diamDeltaLog
      enddo
c
c     Get diameter in physical space
c
      do i=1,nBin
         diams(i) = EXP(diamsLog(i))
      enddo
c
c     Compute aCoeff based on diamPeak and diamMax
c       aCoeff is found by setting df/dD = at D = diamPeak
c
      ddd1 = (diamPeak/diamMin)**power
      ddd2 = (diamPeak/diamMax)**power

      aCoeff = LOG(diamPeak) + bcoeff(iblend)**2 * power *
     &    ( ddd1*(1.0-ddd2) + ddd2*(1.0-ddd1) ) / 
     &    ( (1.0-ddd1)*(1.0-ddd2) )

      WRITE(606,*) 'aCoeff = ', aCoeff
c
c     Determine weights
c
      weightsSum = 0.0d0
      do i = 1, nBin
         weights(i) = ( 1.0 - (diams(i)/diamMax)**power )
     &       *        ( 1.0 - (diams(i)/diamMin)**power )
     &       * EXP( -( LOG( diams(i) ) -aCoeff )**2 
     &       /( 2.0*bcoeff(iblend)**2))
         weightsSum = weightsSum + weights(i)
      enddo 
c
c     Scale Weights
c
c     weightsSum = SUM(weights)
      weights = weights/weightsSum
c
c     Sum Weights
c
      write(606,*) 'Sum of Weights: weightsScaled = ', weightsSum
c 
c     Write output
c
      WRITE(606,*) '#Log Normal Distribution'
      j = 0
      do i = nBin,1,-1
         j = j + 1
         di(iblend,j) = diams(i)
         percent(iblend,j) = abs(weights(i)*100)
      enddo
      nmode(iblend) = nBin
      do i=1,nmode(iblend)
         write(606,105) di(iblend,i),percent(iblend,i)
      enddo
      do i=1,nmode(iblend)-1
         di(iblend,i) = di(iblend,i+1)
         percent(iblend,i) = percent(iblend,i+1)
      enddo
      nmode(iblend) = nmode(iblend)-2
c
 105  format(3X,F8.2,4X,1PE12.5)
c
      return
      end
