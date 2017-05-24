c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c   
      subroutine cuts(diam,dnj,Vfrackept,dN)
c
c
c     Compute the number of particles within each cut
c        Also, store total volume V_total of solids 
c        (excluding binder region)
c     This subroutine takes into account the
c        volume fraction discarded
c
c     On input, diam, dnj, Vfrackept are needed
c     On output, dN is determined, as well as V_total
c
      include 'rocprepack_3.h'
c
      pi = acos(-1.0)
c
c
      if (iblend_max.ge.2) then

         do iblend=1,iblend_max
            sum = 0
            do i=1,nmode(iblend)
               sum = sum + dnj(iblend,i)*diam(iblend,i)**idim
            enddo
            if (Vfrackept(iblend).lt.1.0d-12) Vfrackept(iblend)=1.0d-12
            b(iblend) = pi*sum/(2.0*idim)/Vfrackept(iblend)
         enddo

         den = 0.0
         do i=1,iblend_max
            c(i) = v_frac(1)/v_frac(i)
            if (b(i).gt.1.0d-12) den = den + 1.0/b(i)/c(i)
            write(6,*) c(i)
         enddo
     
         V(1) = ntot/den
         do i=2,iblend_max
            V(i) = V(1)/c(i)
         enddo

         do i=1,iblend_max
            if (b(i).gt.1.0d-12) dN(i) = V(i)/b(i)
         enddo
  
         sum = 0.0
         do i=1,iblend_max
           sum = sum + V(i)
         enddo
         V_total = sum
   
      else
         dN(1) = ntot
         sum = 0
         do i=1,nmode(1)
            sum = sum + dnj(1,i)*diam(1,i)**idim
         enddo
         b(1) = pi*sum/(2.0*idim)/Vfrackept(1)
         V_total = b(1)*dN(1)
      endif

      sum = 0.0
      write(6,*)
      do iblend=1,iblend_max
         sum = sum + dN(iblend)
         write(6,*)'Total particles for cut',iblend,'  = ',dN(iblend)
      enddo
      write(6,*)'Total particles for all cuts = ',sum
      write(6,*)
c
      return
      end
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c   
      subroutine consolidate(diam,dN,dnj,j,iblend)
c
c
c     Consolidation is done BEFORE the various cuts have been
c        combined. This subroutine treats aluminum and AP
c        separately.
c
c
      include 'rocprepack_3.h'
c
c     This subroutine consolidates particles within each cut,
c        then combines all cuts into one distribution; it does
c        not combine like diameters from different cuts.
c
c     On output:
c        itotal      = total number of particles
c        isolid(i)   = solid type for particle i
c        dN_final(i) = Number of particles for particle i
c        D_final(i)  = diameter of particle i
c
c     Consolidate aluminum and AP separately
c
      pi = acos(-1.0)
      j = 1
      do iblend=1,iblend_max

         helpsum = 0
         summing = 0
         do k=1,nmode(iblend)

            if (summing.eq.0) then
               dN_final(j) = dnj(iblend,k)*dN(iblend)
               helpsum = dN_final(j)*(diam(iblend,k)**idim)
               icut(j) = iblend
            elseif (summing.eq.1) then
               dN_final(j) = dN_final(j)
     &               +dnj(iblend,k)*dN(iblend)
               helpsum = helpsum
     &               +dnj(iblend,k)*dN(iblend)*(diam(iblend,k)**idim)
            endif
   
            if (dN_final(j).lt.1.0.and.k.lt.nmode(iblend)) then
               summing = 1
            else
               isolid(j) = isolid_type(iblend)
               D_final(j) = (helpsum/dN_final(j))**(1.0/idim)
               summing = 0
               j = j + 1
            endif

         enddo

      enddo
      itotal = j-1
c
c
      return
      end
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c   
      subroutine output
      include 'rocprepack_3.h'
c
c     This subroutine prints out relevant information
c       to run the packing code
c
      pi = acos(-1.0)
      dntot = 0
      Ntot = 0
      vol_ap0  = 0.0  ! volume based on integer number of particles
      vol_al0  = 0.0  ! volume based on integer number of particles
      vol_rdx0 = 0.0  ! volume based on integer number of particles
      vol_ap1  = 0.0  ! volume based on fractional number of particles
      vol_al1  = 0.0  ! volume based on fractional number of particles
      vol_rdx1 = 0.0  ! volume based on fractional number of particles
      do i = 1,itotal
         dntot = dntot + dN_final(i)
         Ni = dN_final(i) + 0.5
         Ntot = Ntot + Ni
         if (isolid(i).eq.1) then
            vol_ap0 = vol_ap0 + Ni*D_final(i)**idim
            vol_ap1 = vol_ap1 + 
     &           dN_final(i)*D_final(i)**idim
         endif
         if (isolid(i).eq.2) then
            vol_al0 = vol_al0 + Ni*D_final(i)**idim
            vol_al1 = vol_al1 + 
     &           dN_final(i)*D_final(i)**idim
         endif
         if (isolid(i).eq.3) then
            vol_rdx0 = vol_rdx0 + Ni*D_final(i)**idim
            vol_rdx1 = vol_rdx1 + 
     &           dN_final(i)*D_final(i)**idim
         endif
      enddo
      write(6,*)
      write(6,*)'Total number of particles before roundoff = ',dntot
      write(6,*)'Total number of particles after  roundoff = ',Ntot
      write(6,*)
 100  format(2x,4i8,2x,f20.12)
c
c     Compute total volume fraction discarded, and
c       the value 'alpha' of the blend that is AP
c
      Vfrac_kept0 = (pi*vol_ap0/(2.0*idim))/V_total
      Vfrac_kept1 = (pi*vol_ap1/(2.0*idim))/V_total
      Vfrac_kept2 = (pi*vol_rdx1/(2.0*idim))/V_total
c
      x = rho/(1.0 + (1.0-Vfrac_kept1)/Vfrackept_total)
      alpha = (rho-x)/(1.0-x)
      w = Vfrac_kept1
      write(6,*)'Values based on fractional number of particles'
      write(6,*)'   Fraction of solids kept      = ',w
      write(6,*)'   Fraction of solids discarded = ',1.0-w
      write(6,*)'   Theoretical packing fraction       = ',rho
      write(6,*)'   Updated packing fraction           = ',x
      write(6,*)'   Volume fraction AP in binder       = ',alpha
      x = rho/(1.0 + (1.0-Vfrac_kept0)/Vfrackept_total)
      alpha = (rho-x)/(1.0-x)
      w = Vfrac_kept0
      write(6,*)
      write(6,*)'Values based on integer number of particles'
      write(6,*)'   Fraction of solids kept      = ',w
      write(6,*)'   Fraction of solids discarded = ',1.0-w
      write(6,*)'   Theoretical packing fraction       = ',rho
      write(6,*)'   Updated packing fraction           = ',x
      write(6,*)'   Volume fraction AP in binder       = ',alpha
c
      if (alpha.lt.0.0) then
         write(6,*)
         write(6,*)'**********************************'
         write(6,*)'CAUTION: ALPHA NEGATIVE; REDO PACK'
         write(6,*)'**********************************'
         write(6,*)
      endif
c
c     Compute dimensions of rectangle (2D) or cube (3D)
c
      vol_ap0 = vol_ap0 * pi/(2.0*idim)
      vol_ap1 = vol_ap1 * pi/(2.0*idim)
      vol_al0 = vol_al0 * pi/(2.0*idim)
      vol_al1 = vol_al1 * pi/(2.0*idim)
      vol_rdx0 = vol_rdx0 * pi/(2.0*idim)
      vol_rdx1 = vol_rdx1 * pi/(2.0*idim)
      vol_tot0 = (vol_ap0 + vol_al0 + vol_rdx0)/x
      vol_tot1 = (vol_ap1 + vol_al1 + vol_rdx1)/rho
      vol_binder = vol_tot-(vol_ap0 + vol_al0 + vol_rdx0)
      write(6,*)
      if (idim.eq.2) then
         if (height.gt.0) then
            dly = height
            dlx = vol_tot0/height
         else
            dly = sqrt(vol_tot0)
            dlx = dly
         endif
         write(6,*)'Values based on fractional number of particles'
         write(6,*)'   Area fraction of AP  = ',100.0*vol_ap1/vol_tot1
         write(6,*)'   Area fraction of Al  = ',100.0*vol_al1/vol_tot1
         write(6,*)'   Area fraction of RDX = ',100.0*vol_rdx1/vol_tot1
         write(6,*)'Values based on integer number of particles'
         write(6,*)'   Area fraction of AP  = ',100.0*vol_ap0/vol_tot0
         write(6,*)'   Area fraction of Al  = ',100.0*vol_al0/vol_tot0
         write(6,*)'   Area fraction of RDX = ',100.0*vol_rdx0/vol_tot0
         write(6,*)
         write(6,*)'Number of cylinders  = ',ntot
         write(6,*)'Total Overall Length = ',dlx
         write(6,*)'Total Overall Height = ',dly
      else
         if (height.gt.0) then
            dlz = height
            dlx = sqrt(vol_tot0/height)
            dly = dlx
         else
            dlx = (vol_tot0)**(1.0/idim)
            dly = dlx
            dlz = dlx
         endif
         write(6,*)'Values based on fractional number of particles'
         write(6,*)'  Volume fraction of AP  = ',100.0*vol_ap1/vol_tot1
         write(6,*)'  Volume fraction of Al  = ',100.0*vol_al1/vol_tot1
         write(6,*)'  Volume fraction of RDX = ',100.0*vol_rdx1/vol_tot1
         write(6,*)'Values based on integer number of particles'
         write(6,*)'  Volume fraction of AP  = ',100.0*vol_ap0/vol_tot0
         write(6,*)'  Volume fraction of Al  = ',100.0*vol_al0/vol_tot0
         write(6,*)'  Volume fraction of RDX = ',100.0*vol_rdx0/vol_tot0
         write(6,*)
         write(6,*)'Number of spheres = ',ntot
         write(6,*)'Length dlx = ',dlx
         write(6,*)'Width  dly = ',dly
         write(6,*)'Height dlz = ',dlz
      endif
c
      open(UNIT=15,FILE='rocprepack_out')
      if (idim.eq.2) then
         write(15,*) dlx,'	! Length'
         write(15,*) dly,'	! Height'
      endif
      if (idim.eq.3) then
         write(15,*) dlx,'	! Length'
         write(15,*) dly,'	! Width'
         write(15,*) dlz,'	! Height'
      endif
      write(15,*) rho,'	! Theoretical packing density'
      write(15,*) x,'	! Packing density'
      write(15,*) alpha,'	! Fraction of binder that is Solids'
      write(15,*) ntot
      write(15,*) itotal
      isum1 = 0
      do i= 1,itotal
         Ni = dN_final(i) + 0.5
         if (Ni.eq.0) then
            Ni=1
            write(6,*)'*** increasing Ni',Ni
         endif
         write(15,110) icut(i),isolid(i),Ni,D_final(i)
         if (isolid(i).eq.2) 
     &      write(20,110) icut(i),
     &                       isolid(i),Ni,D_final(i)
         if (isolid(i).eq.2.and.icut(i).eq.1) 
     &      write(30,110) icut(i),
     &                       isolid(i),Ni,D_final(i)
         if (isolid(i).eq.2.and.icut(i).eq.2) 
     &      write(31,110) icut(i),
     &                       isolid(i),Ni,D_final(i)
         if (isolid(i).eq.2.and.icut(i).eq.3) 
     &      write(32,110) icut(i),
     &                       isolid(i),Ni,D_final(i)
         if (isolid(i).eq.2.and.icut(i).eq.4) 
     &      write(33,110) icut(i),
     &                       isolid(i),Ni,D_final(i)
         write(22,110) icut(i),isolid(i),Ni,D_final(i)
         if (isolid(i).eq.2) isum1 = isum1 + Ni
      enddo
      write(6,*)'Total number of aluminum particles = ',isum1
 110  format(2x,3i8,2x,f20.12)
      close(15)
c
      return
      end
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
      subroutine mix_and_rnd
      include 'rocprepack_3.h'
c
      dimension Dtemp(200000)
      dimension dntemp(200000)
      dimension Ditemp(200000)
      dimension DNi(200000)
      dimension Vi(200000)
c
c
c     This subroutine combines like particles, and is only
c       used as an output device.
c
c     This subroutine should not be called if there is aluminum
c       in the propellant, as it will combine like diameters of
c       AP and Al. This routine would need some work to
c       avoid this.
c
c
      pi = acos(-1.0)
      dntot = 0
      do i= 1,itotal
         dntot = dntot + dN_final(i)
         dntemp(i) = dN_final(i)
         Dtemp(i) = D_final(i)
      enddo
c
c     Sort's particles by size; not really needed
c       if there is only one cut, they are already sorted
c
      ipass=1
  10  continue
      if (ipass.le.itotal) then
         do i=1,itotal-1
            if (Dtemp(i).lt.Dtemp(i+1)) then
               temp = Dtemp(i+1)
               Dtemp(i+1) = Dtemp(i)
               Dtemp(i)=temp
               temp = dntemp(i+1)
               dntemp(i+1) = dntemp(i)
               dntemp(i) = temp
            endif
         enddo
         ipass = ipass + 1
         go to 10
      end if
c
c
      tol=1.0e-5
      Ditemp(1)=Dtemp(1)
      DNi(1)=dntemp(1)

      ii=2

      hlp_eps=1.0e-12
      do k=2,itotal
         if (abs(Dtemp(k)-Dtemp(k-1)).lt.tol) then
            DNi(ii-1)=DNi(ii-1)+dntemp(k)
         else
            if (dntemp(k).gt.hlp_eps) then
               Ditemp(ii)=Dtemp(k)
               DNi(ii)=dntemp(k)
               ii=ii+1
            endif
         endif
      enddo

      nmax2 = ii-1

      Ntot = 0
      do k=1,nmax2
         Ni = Dni(k)+0.5
         Ntot = Ntot + Ni
      enddo
c
      write(16,*) dlx,dly,dlz
      write(16,*) alpha
      write(16,*) ntot
      write(16,*) nmax2
      Ntot = 0
      do k=1,nmax2
         Ni = Dni(k)+0.5
         write(16,30) k,Ditemp(k),DNi(k),Ni
      enddo
  30  format(2x,i10,2f12.4,3x,i8)
c
c
      return
      end
c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c
      subroutine num_fraction(iflag,iblend,nmode,
     &                        diam,Vfract,dnj,idim)
      dimension a(500,500),b(500),indx(500)
      dimension diam(6,500),Vfract(6,500),dnj(6,500)
      dimension nmode(6)
c
c     This subroutine computes the number fractions for each cut
c       given Nlocal; since we are only interested in the number
c       fraction, the value of Nlocal is irrelevant
c
      nmode0 = nmode(iblend)
      Nlocal = 10000
c
      pi = acos(-1.0)
      do i=1,200
         indx(i)=0
         b(i) = 0.0
         do j=1,200
            a(i,j) = 0.0
         enddo
      enddo
c
      do i=1,nmode0
         a(i,i) = 1.0
         a(i,i+nmode0) = -pi*(diam(iblend,i)**idim)/(2.0*idim)
         if (iflag.eq.1) a(i,i+nmode0) = -diam(iblend,i)
      enddo
      do i = nmode0+1,2*nmode0-1
         a(i,1) = 1.0
         a(i,i-nmode0+1) = -Vfract(iblend,1)/Vfract(iblend,i-nmode0+1)
      enddo
      do j=nmode0+1,2*nmode0
         a(2*nmode0,j) = 1.0
      enddo
c
      do i=1,2*nmode0-1
         b(i) = 0.0
      enddo
      b(2*nmode0) = Nlocal
c
      n = 2*nmode0
      np = n
      call ludcmp(a,n,np,indx,d)
      call lubksb(a,n,np,indx,b)
c
      do i=1,nmode0
         dnj(iblend,i) = abs(b(nmode0+i))/Nlocal
      enddo
c
      return
      end
c
c
c
c
c
      SUBROUTINE ludcmp(a,n,np,indx,d)
      INTEGER indx(500)
      dimension a(500,500)
      PARAMETER (NMAX=500,TINY=1.0e-20)
      dimension vv(NMAX)
      d=1.
      do 12 i=1,n
        aamax=0.
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
        if (aamax.eq.0.) stop ! 'singular matrix in ludcmp'
        vv(i)=1./aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.
        do 16 i=j,n
          sum=a(i,j)
          do 15 k=1,j-1
            sum=sum-a(i,k)*a(k,j)
15        continue
          a(i,j)=sum
          dum=vv(i)*abs(sum)
          if (dum.ge.aamax) then
            imax=i
            aamax=dum
          endif
16      continue
        if (j.ne.imax)then
          do 17 k=1,n
            dum=a(imax,k)
            a(imax,k)=a(j,k)
            a(j,k)=dum
17        continue
          d=-d
          vv(imax)=vv(j)
        endif
        indx(j)=imax
        if(a(j,j).eq.0.)a(j,j)=TINY
        if(j.ne.n)then
          dum=1./a(j,j)
          do 18 i=j+1,n
            a(i,j)=a(i,j)*dum
18        continue
        endif
19    continue
      return
      END
      SUBROUTINE lubksb(a,n,np,indx,b)
      INTEGER indx(500)
      dimension a(500,500),b(500)
      ii=0
      do 12 i=1,n
        ll=indx(i)
        sum=b(ll)
        b(ll)=b(i)
        if (ii.ne.0)then
          do 11 j=ii,i-1
            sum=sum-a(i,j)*b(j)
11        continue
        else if (sum.ne.0.) then
          ii=i
        endif
        b(i)=sum
12    continue
      do 14 i=n,1,-1
        sum=b(i)
        do 13 j=i+1,n
          sum=sum-a(i,j)*b(j)
13      continue
        b(i)=sum/a(i,i)
14    continue
      return
      END
