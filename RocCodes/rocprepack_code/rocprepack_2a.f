c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c   
      subroutine cuts
      include 'rocprepack_3.h'
c
c     local variables
c
      integer iblend,i
      double precision V,b,c,sum1,den,Dn_sum
      dimension V(max_cut)
      dimension b(max_cut)
      dimension c(max_cut)
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
c
c
      if (iblend_max.ge.2) then

         do iblend=1,iblend_max
            sum1 = 0
            do i=1,nmode(iblend)
               sum1 = sum1 + dnj(iblend,i)*diam(iblend,i)**idim
            enddo
            if (Vfrackept(iblend).lt.1.0d-12) Vfrackept(iblend)=1.0d-12
            b(iblend) = pi*sum1/(2.0*idim)/Vfrackept(iblend)
         enddo

         den = 0.0
         do i=1,iblend_max
            c(i) = v_frac(i)/v_frac(1)
            if (b(i).gt.1.0d-12) den = den + c(i)/b(i)
         enddo

         V(1) = ntot/den
         do i=2,iblend_max
            V(i) = V(1)*c(i)
         enddo

         do i=1,iblend_max
            if (b(i).gt.1.0d-12) dN(i) = V(i)/b(i)
         enddo

         V_solids = sum(V(1:iblend_max))
         V_binder = V_solids*(1.0-rho_v)/rho_v
         V_total = V_solids + V_binder

      else
         dN(1) = ntot
         sum1 = 0
         do i=1,nmode(1)
            sum1 = sum1 + dnj(1,i)*diam(1,i)**idim
         enddo
         b(1) = pi*sum1/(2.0*idim)/Vfrackept(1)
         V_solids = b(1)*dN(1)
         V_binder = V_solids*(1.0-rho_v)/rho_v
         V_total = V_solids + V_binder
      endif

      write(606,*)
      dN_sum = sum(dN(1:iblend_max))
      do iblend=1,iblend_max
         write(6,*)'      Total particles for cut',iblend,
     &               '(',TRIM(cut_type(iblend)),')',
     &               '	= ',dN(iblend)
         write(606,*)'      Total particles for cut',iblend,
     &               '(',TRIM(cut_type(iblend)),')',
     &               '	= ',dN(iblend)
      enddo
      write(606,*)'      Total particles for all cuts 	= ',dN_sum
      write(606,*)
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
      subroutine consolidate
      include 'rocprepack_3.h'
c
c     local variables
c
      integer iblend,j,k,Ni,i
      double precision helpsum,summing
c
c
c     Consolidation is done BEFORE the various cuts have been
c        combined. This subroutine treats aluminum and AP
c        separately.
c
c     This subroutine consolidates particles within each cut,
c        then combines all cuts into one distribution; it does
c        not combine like diameters from different cuts.
c
c     The consolidation process preserves volume
c
c     On output:
c        itotal         = total number of particles
c        solidname(i)   = solid name for particle i
c        dN_final(i)    = Number of particles for particle i
c        D_final(i)     = diameter of particle i
c
c     Consolidate different cuts separately
c
      j = 1
      do iblend=1,iblend_max

         helpsum = 0
         summing = 0
         do k=1,nmode(iblend)

            write(19,11) dnj(iblend,k),dN(iblend),diam(iblend,k)
  11        format(2x,3f15.8)
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
               solidname(j) = cut_type(iblend)
               D_final(j) = (helpsum/dN_final(j))**(1.0/idim)
               summing = 0
               j = j + 1
            endif

         enddo

      enddo
      itotal = j-1
c
c     reset dN_final so that Ni can never be zero
c
c     do i = 1,itotal
c        Ni = dN_final(i) + 0.5
c        if (Ni.eq.0) dN_final(i) = 0.9
c     enddo
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
c     local variables
c
      integer i,ii,Ni,iblend,j,icase
      integer isum(max_solid_type)
      double precision vol_solids_0,vol_solids_1,alp_0,alp_1,
     &   vol_0,vol_1,rho_v_0,rho_v_1,vol_tot0,vol_tot1,
     &   v_blend_0,v_blend_1,alpha_0,alpha_1,sum0,sum1,
     &   dntot,V_kept0,V_kept1

      dimension vol_solids_0(20),vol_solids_1(20)
      dimension alp_1(20),alp_0(20)
      character*(200) comment_line
      character*(200) color
c
c     This subroutine prints out relevant information
c       to run the packing code
c
c
c
c     v_solids = theoretical volume of particles
c     v_binder = theoretical volume of pure binder
c     v_total  = theoretical total volume
c
c     vol_solids_0(k) = vol of integer    num. particles for solid_type k
c     vol_solids_1(k) = vol of fractional num. particles for solid_type k
c     vol_0 = total solids (only) vol of integer    number of particles
c     vol_1 = total solids (only) vol of fractional number of particles
c     rho_v_0  = packing fraction for integer    number of particles
c     rho_v_1  = packing fraction for fractional number of particles
c     vol_tot0 = total vol of integer    number of particles
c     vol_tot1 = total vol of fractional number of particles
c
c     note: vol_0,vol_1 <= v_solids
c
c     v_solids-vol_0 = discard volume to be added to binder
c     v_solids-vol_1 = discard volume to be added to binder
c
c     dntot = total number of spheres before roundoff
c     Ntot  = total number of spheres after roundoff
c
c
c
c     compute total number of particles
c
      Ntot = sum(int(dN_final + 0.5))
      dntot = sum(dN_final)
c
c     compute volume for each propellant ingredient
c
      vol_solids_0 = 0.0d0
      vol_solids_1 = 0.0d0
      do i=1,Ntype
         do ii=1,itotal
            if (solidname(ii).eq.solidtype(i)) then
               Ni = dN_final(ii) + 0.5
               vol_solids_0(i) = vol_solids_0(i) + 
     &            Ni*D_final(ii)**idim
               vol_solids_1(i) = vol_solids_1(i) +
     &            dN_final(ii)*D_final(ii)**idim
            endif
         enddo
      enddo
c
      vol_solids_0 = vol_solids_0 * pi/(2.0*idim)
      vol_solids_1 = vol_solids_1 * pi/(2.0*idim)
      vol_0 = sum(vol_solids_0)
      vol_1 = sum(vol_solids_1)
      rho_v_0 = vol_0/v_total
      rho_v_1 = vol_1/v_total
      vol_tot0 = vol_0/rho_v_0
      vol_tot1 = vol_1/rho_v_1
c
      v_blend_0 = v_binder + (v_solids-vol_0)
      v_blend_1 = v_binder + (v_solids-vol_1)
      alpha_0 = (v_solids-vol_0)/v_blend_0
      alpha_1 = (v_solids-vol_1)/v_blend_1
c
c     compute fraction of solids kept
c
      V_kept0 = vol_0/V_solids
      V_kept1 = vol_1/V_solids
c
c     compute fraction of solids that is to go into the binder
c        to form the blend
c
      do i=1,Ntype
         sum0 = 0.0d0
         sum1 = 0.0d0
         do iblend=1,iblend_max
            if (cut_type(iblend).eq.solidtype(i)) then
               sum0 = sum0 + v_frac(iblend)*v_solids/rho_v
               sum1 = sum1 + v_frac(iblend)*v_solids/rho_v
            endif
         enddo
         alp_0(i) = (sum0 - vol_solids_0(i))/v_total
         alp_1(i) = (sum1 - vol_solids_1(i))/v_total
      enddo
c
c
      write(606,*)
      write(606,*)'Total number of particles before roundoff = ',dntot
      write(606,*)'Total number of particles after  roundoff = ',Ntot
      write(606,*)
c
c
      write(6,*)
      write(6,*)'   Theoretical packing fraction     = ',rho_v
      write(6,*)'   Updated packing fraction         = ',rho_v_0
      write(6,*)
      write(606,*)'Values based on fractional number of particles'
      write(606,*)'   Fraction of solids kept          = ',V_kept1
      write(606,*)'   Fraction of solids discarded     = ',1.0-V_kept1
      write(606,*)'   Theoretical packing fraction     = ',rho_v
      write(606,*)'   Updated packing fraction         = ',rho_v_1
      write(606,*)'   Vol fraction of solids in binder = ',alpha_1
      write(606,*)
      write(606,*)'Values based on integer number of particles'
      write(606,*)'   Fraction of solids kept          = ',V_kept0
      write(606,*)'   Fraction of solids discarded     = ',1.0-V_kept0
      write(606,*)'   Theoretical packing fraction     = ',rho_v
      write(606,*)'   Updated packing fraction         = ',rho_v_0
      write(606,*)'   Vol fraction of solids in binder = ',alpha_0
c
      if (alpha_0.lt.0.0) then
         write(606,*)
         write(606,*)'**********************************'
         write(606,*)'CAUTION: ALPHA NEGATIVE; REDO PACK'
         write(606,*)'         INCREASE NTOT'
         write(606,*)'**********************************'
         write(606,*)

         write(*,*)
         write(*,*)'**********************************'
         write(*,*)'CAUTION: ALPHA NEGATIVE; REDO PACK'
         write(*,*)'         INCREASE NTOT'
         write(*,*)'**********************************'
         write(*,*)

c        pause
      endif
c
c     Compute dimensions 
c
      write(606,*)
         if (TRIM(container).eq.'Cuboid' .or.
     &       TRIM(container).eq.'Cube' .or.
     &       TRIM(container).eq.'cuboid' .or.
     &       TRIM(container).eq.'cube') then
            if (height.gt.0) then ! use actual height, in microns
               dly = height
               dlx = sqrt(v_total/height)
               dlz = dlx
            else ! use aspect ratio
               dlx = (v_total/abs(height))**(1.0/3.0)
               dly = dlx*abs(height)
               !dlz = 0.24*dlx
               dlz = dlx
               write(6,*) dlx,dlz,dly,height
            endif
            write(606,*)'Values based on fractional number of particles'
            do i=1,Ntype
               sum1 = vol_solids_1(i)/vol_tot1
               write(606,55) TRIM(solidtype(i)),100.0*sum1
            enddo
            write(606,*)'Values based on integer number of particles'
            do i=1,Ntype
               sum0 = vol_solids_0(i)/vol_tot0
               write(606,55) TRIM(solidtype(i)),100.0*sum0
            enddo
            write(606,*)
            write(606,*)'Number of spheres = ',ntot
            write(606,*)'Length dlx = ',dlx
            write(606,*)'Width  dlz = ',dlz
            write(606,*)'Height dly = ',dly
            write(6,*)
            write(6,*)'Number of spheres = ',ntot
            write(6,*)'Length dlx = ',dlx
            write(6,*)'Width  dlz = ',dlz
            write(6,*)'Height dly = ',dly
            icase = 1
         endif
         if (TRIM(container).eq.'Cylinder' .or.
     &       TRIM(container).eq.'cylinder') then
            if (height.gt.0) then ! use actual height, in microns
               dlx = sqrt(v_total/pi/height)
               dly = height
               dlz = dlx
            else ! use aspect ratio
               height = abs(height)
               dlx = (0.5d0*v_total/pi/height)**(1.0d0/3.0d0)
               dly = 2.0d0*height*dlx
               dlz = dlx
            endif
            write(606,*)'Values based on fractional number of particles'
            do i=1,Ntype
               sum1 = vol_solids_1(i)/vol_tot1
               write(606,55) TRIM(solidtype(i)),100.0*sum1
            enddo
            write(606,*)'Values based on integer number of particles'
            do i=1,Ntype
               sum0 = vol_solids_0(i)/vol_tot0
               write(606,55) TRIM(solidtype(i)),100.0*sum0
            enddo
            write(606,*)
            write(606,*)'Number of spheres = ',ntot
            write(606,*)'Radius  dlx = ',dlx
            write(606,*)'Height  dly = ',dly
            write(6,*)
            write(6,*)'Number of spheres = ',ntot
            write(6,*)'Radius  dlx = ',dlx
            write(6,*)'Height  dly = ',dly
            icase = 2
         endif
         if (TRIM(container).eq.'Annulus' .or.
     &       TRIM(container).eq.'annulus') then
            dlx = rad_inner
            dly = height
            dlz = sqrt(dlx*dlx+v_total/(pi*dly))
            write(606,*)'Values based on fractional number of particles'
            do i=1,Ntype
               sum1 = vol_solids_1(i)/vol_tot1
               write(606,55) TRIM(solidtype(i)),100.0*sum1
            enddo
            write(606,*)'Values based on integer number of particles'
            do i=1,Ntype
               sum0 = vol_solids_0(i)/vol_tot0
               write(606,55) TRIM(solidtype(i)),100.0*sum0
            enddo
            write(606,*)
            write(606,*)'Number of spheres = ',ntot
            write(606,*)'Inner radius dlx = ',dlx
            write(606,*)'Outer radius dlz = ',dlz
            write(606,*)'Height       dly = ',dly
            write(6,*)
            write(6,*)'Number of spheres = ',ntot
            write(6,*)'Inner radius dlx = ',dlx
            write(6,*)'Outer radius dlz = ',dlz
            write(6,*)'Height       dly = ',dly
            icase = 3
         endif
  54  format(2x,'Area fraction of ',A3,' = ',f12.6)
  55  format(2x,'Vol. fraction of ',A3,' = ',f12.6)
c
c     Compute sample mass, result in gm, and density gm/cm^3
c
      sample_mass = sum(density(1:iblend_max) * v_frac(1:iblend_max))
     &     * V_solids/rho_v + density_binder*V_binder
      sample_density = sample_mass/V_total
      sample_mass = sample_mass*1d-12
c
c--------------------------------------------------------------
c
c     now print out info to rocprepack_out
c         specifically for packing code
c
      open(UNIT=15,FILE='rocprepack_out')
      open(UNIT=35,FILE='newpackformat_out')

      write(35,*)'# pack -vD -f',rho_v_0
      write(35,*)
      write(35,*)'# create a pack'
      write(35,*)'import "shapes/adn"'
      write(35,*)'import "shapes/cube"'
      write(35,*)'import "shapes/hmx"'
      write(35,*)'import "shapes/petn"'
      !write(35,*)'import "shapes/cl20"'
      !write(35,*)'import "shapes/rdx"'
      write(35,*)

      write(35,*)
      write(35,*)'set packing_fraction = ',rho_v_0,';'
      write(35,*)'set seed = 13;'
      write(35,*)'set growth_rate = 1.0;'
      if (icase.eq.1) then
         write(35,*)'boundary {box',dlx/dlx,dly/dlx,dlz/dlx,' periodic}'
      endif
      if (icase.eq.2) then
         write(35,*)'boundary {cylinder ',dlx/dlx,dly/dlx,'}'
      endif
      if (icase.eq.3) then
         write(35,*)'boundary {annulus ',dlx/dlz,dlz/dlz,dly/dlz,'}'
      endif
      write(35,*)

      do i= 1,itotal
         Ni = dN_final(i) + 0.5
         if (icut(i).eq.1) color = '0.6 0.6 0.6'
         !if (icut(i).eq.2) color = '0.6 0.8 0.99'
         !if (icut(i).eq.2) color = '0.0 0.4 0.8'
         if (icut(i).eq.2) color = '1  .7 .4'
         if (icut(i).eq.3) color = '0 .5 1'
         !if (icut(i).eq.3) color = '.4 .7 1.'
         !if (icut(i).eq.3) color = '.99 .6 .2'
         if (TRIM(solidname(i)).eq.'Al') color='1 0 0'
         !if (TRIM(solidname(i)).eq.'Al') color='.2 .4 1.'
         if (TRIM(solidname(i)).eq.'AP'
     &        .or.TRIM(solidname(i)).eq.'Al') then
            write(35,*)'create ',Ni,' sphere ',
     &       ' size ',D_final(i),'tag',i,'color ',TRIM(color),' ;'
         else
            write(35,*)'create ',Ni,TRIM(solidname(i)),
     &       ' size ',D_final(i),'tag',i,'color ',TRIM(color),' ;'
         endif
      enddo
      write(35,*)
      write(35,*)
      write(35,*)


      write(15,*) '***   ',TRIM(prop_name)
      write(15,*) '***   ',TRIM(prop_description)
      write(35,*) '***   ',TRIM(prop_name)
      write(35,*) '***   ',TRIM(prop_description)
c     write(15,*) '***   ',TRIM(container)
      if (val(7).ge.10) then
         write(15,'(1x,a,i5,2x,a,i3,2x,a,i3,2x,3(a,i2))')
     &         '***   Yr:',val(1),'Mo:',val(2),'Day:',val(3),
     &         'Time: ',val(5),':',val(6),':',val(7)
         write(35,'(1x,a,i5,2x,a,i3,2x,a,i3,2x,3(a,i2))')
     &         '***   Yr:',val(1),'Mo:',val(2),'Day:',val(3),
     &         'Time: ',val(5),':',val(6),':',val(7)
      else
         write(15,'(1x,a,i5,2x,a,i3,2x,a,i3,2x,2(a,i2),a,i1)')
     &         '***   Yr:',val(1),'Mo:',val(2),'Day:',val(3),
     &         'Time: ',val(5),':',val(6),':0',val(7)
         write(35,'(1x,a,i5,2x,a,i3,2x,a,i3,2x,2(a,i2),a,i1)')
     &         '***   Yr:',val(1),'Mo:',val(2),'Day:',val(3),
     &         'Time: ',val(5),':',val(6),':0',val(7)
      endif
      write(15,*)
      write(35,*)
      write(15,'(1x,a,3x,a,1x,1pe12.5,2x,a)')'***','Propellant Density'
     &     ,sample_density,'gm/cm^3'
      write(15,'(1x,a,3x,a,4x,1pe12.5,2x,a)')'***','Propellant Mass'
     &     ,sample_mass,'gm'
      write(35,'(1x,a,3x,a,1x,1pe12.5,2x,a)')'***','Propellant Density'
     &     ,sample_density,'gm/cm^3'
      write(35,'(1x,a,3x,a,4x,1pe12.5,2x,a)')'***','Propellant Mass'
     &     ,sample_mass,'gm'
      write(15,*)
      write(35,*)
      do i=1,Ntype
         comment_line = '!(V_{'// TRIM(solidtype(i))//',total} - V_{'//
     &        TRIM(solidtype(i))//',spheres}) / (Length Width Height)'
         write(15,101)TRIM(solidtype(i)),alp_0(i),trim(comment_line)
         write(35,101)TRIM(solidtype(i)),alp_0(i),trim(comment_line)
      enddo
      comment_line = 'Binder'
      comment_line = '! V_{'// TRIM(comment_line) //
     &     '} / (Length Width Height)'
      write(15,101)'Binder',1.d0-sum(v_frac),TRIM(comment_line)
      write(35,101)'Binder',1.d0-sum(v_frac),TRIM(comment_line)
         if (TRIM(container).eq.'Cuboid' .or.
     &       TRIM(container).eq.'Cube' .or.
     &       TRIM(container).eq.'cuboid' .or.
     &       TRIM(container).eq.'cube') then
            write(15,111) dlx
            write(15,112) dlz
            write(15,113) dly
            write(35,111) dlx
            write(35,112) dlz
            write(35,113) dly
         endif
         if (TRIM(container).eq.'Cylinder' .or.
     &       TRIM(container).eq.'cylinder') then
            write(15,114) dlx
            write(15,115) dly
            write(35,114) dlx
            write(35,115) dly
         endif
         if (TRIM(container).eq.'Annulus' .or.
     &       TRIM(container).eq.'annulus') then
            write(15,116) dlx
            write(15,117) dlz
            write(15,118) dly
            write(35,116) dlx
            write(35,117) dlz
            write(35,118) dly
         endif
      write(15,120) rho_v
      write(15,121) rho_v_0
      write(15,122) ntot
      write(15,123) itotal
      write(35,120) rho_v
      write(35,121) rho_v_0
      write(35,122) ntot
      write(35,123) itotal
      isum = 0
      do i= 1,itotal
         Ni = dN_final(i) + 0.5
ccTLJ        write(15,110) icut(i),TRIM(solidname(i)),Ni,D_final(i)
ccTLJ     &        ,coating(icut(i))
             write(15,134) icut(i),Ni,TRIM(solidname(i)),D_final(i)
             write(35,134) icut(i),Ni,TRIM(solidname(i)),D_final(i)
         do j=1,Ntype
            if (solidname(i).eq.solidtype(j)) isum(j)=isum(j)+Ni
         enddo
      enddo
c
c     end print out to rocprepack_out
c
c--------------------------------------------------------------
c     print out integer number of particles for each cut
c
      write(606,*)
      do i=1,Ntype
         write(606,124) TRIM(solidtype(i)),isum(i)
      enddo
c
 101  format(2x,A8,2x,1(2x,f20.8),2x,a)
 102  format(2x,A8,2x,2(2x,f20.8))
 110  format(2x,i8,A8,i12,2x,f24.12,2x,f20.12)
 134  format(2x,i8,i12,2x,A12,2x,f24.12)
 111  format(2x,f24.14,8x,'! Length')
 112  format(2x,f24.14,8x,'! Width')
 113  format(2x,f24.14,8x,'! Height')
 114  format(2x,f24.14,8x,'! Diameter')
 115  format(2x,f24.14,8x,'! Height')
 116  format(2x,f24.14,8x,'! Inner radius')
 117  format(2x,f24.14,8x,'! Outer radius')
 118  format(2x,f24.14,8x,'! Height')
 120  format(2x,f24.14,8x,'! Theoretical packing density')
 121  format(2x,f24.14,8x,'! Packing density')
 122  format(2x,i24,8x,'! Total number of particles')
 123  format(2x,i24,8x,'! Total number of different diameters')
 124  format(2x,'Total number of ',A4,' particles',i24)
      close(15)
      close(35)
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
      subroutine mix_and_rnd
      include 'rocprepack_3.h'
c
c     local variables
c
      integer i,ipass,ii,k,nmax2,Ni
      integer isum(max_solid_type)
      double precision Dtemp,dntemp,Ditemp,Dni,Vi,dntot,
     &   temp,tol,hlp_eps
 
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
      write(16,*) TRIM(prop_name)
      write(16,*) dlx
      write(16,*) dly
      write(16,*) dlz
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
      subroutine num_fraction(iblend)
      include 'rocprepack_3.h'
c
c     local variables
c
      integer iblend,nmode0,Nlocal,i,j,n,np
      integer indx
      double precision a_local,b_local,d_local

      dimension a_local(500,500),b_local(500),indx(500)
c
c     This subroutine computes the number fractions for each cut
c       given Nlocal; since we are only interested in the number
c       fraction, the value of Nlocal is irrelevant
c
      nmode0 = nmode(iblend)
      Nlocal = 10000
      if (nmode0.eq.0) then
         return
      endif
c
      a_local = 0.0d0
      b_local = 0.0d0
c
      do i=1,nmode0
         a_local(i,i) = 1.0d0
         a_local(i,i+nmode0) = -pi*(diam(iblend,i)**idim)/(2.0*idim)
         if (iflag.eq.1) a_local(i,i+nmode0) = -diam(iblend,i)
      enddo
      do i = nmode0+1,2*nmode0-1
         a_local(i,1) = 1.0d0
         a_local(i,i-nmode0+1) = 
     &       -Vfract(iblend,1)/Vfract(iblend,i-nmode0+1)
      enddo
      do j=nmode0+1,2*nmode0
         a_local(2*nmode0,j) = 1.0d0
      enddo
c
      do i=1,2*nmode0-1
         b_local(i) = 0.0d0
      enddo
      b_local(2*nmode0) = Nlocal
c
      n = 2*nmode0
      np = n
      call ludcmp(a_local,n,np,indx,d_local)
      call lubksb(a_local,n,np,indx,b_local)
c
      do i=1,nmode0
         dnj(iblend,i) = abs(b_local(nmode0+i))/Nlocal
      enddo
c
c
      return
      end
c
c
c
c
c
      SUBROUTINE ludcmp(a,n,np,indx,d)
c
c     local variables
c
      implicit none
      integer nmax,i,j,k,n,np,imax
      integer indx(500)
      double precision tiny,d,aamax,sum,dum
      double precision a(500,500)
      PARAMETER (NMAX=500,TINY=1.0e-20)
      double precision vv(NMAX)
c
      d=1.d0
      do 12 i=1,n
        aamax=0.d0
        do 11 j=1,n
          if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
11      continue
c       if (aamax.eq.0.) pause 'singular matrix in ludcmp'
        vv(i)=1.d0/aamax
12    continue
      do 19 j=1,n
        do 14 i=1,j-1
          sum=a(i,j)
          do 13 k=1,i-1
            sum=sum-a(i,k)*a(k,j)
13        continue
          a(i,j)=sum
14      continue
        aamax=0.d0
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
      implicit none
      INTEGER indx(500),n,np,ii,i,j,ll
      double precision a(500,500),b(500)
      double precision sum
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
c
c
      return
      end
