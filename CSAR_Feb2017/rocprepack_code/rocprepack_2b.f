c
c
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
c
c   
c
      subroutine input_pack_setup
      include 'rocprepack_3.h'
c
c     local variables
c
      integer Public_flag,igui,iblend,i,ifill,j
      double precision diam_singlesize
c
c     This subroutine reads in appropriate information
c       to generate packing distribution
c
      open(UNIT=600,file='filenames.list',ACTION='READ')
      read(600,*) inputfile
      do while (scan(inputfile,'!')==1)
         read(600,*) inputfile
      enddo
      write(6,*) TRIM(inputfile)
      prop_name=inputfile

      mix_name=TRIM('Formulations_data/')//TRIM(prop_name)
      write(6,*) TRIM(mix_name)

      open(500,file=mix_name)
      read(500,*) prop_description
      read(500,*) container
      read(500,*) idim
      read(500,*) height
      if (TRIM(container).eq.'Annulus') read(500,*) rad_inner
      read(500,*)
      read(500,*) iblend_max
      read(500,*) ntot
      read(500,*) density_binder

      if (TRIM(container).ne.'Cylinder' .and.
     &    TRIM(container).ne.'cylinder' .and.
     &    TRIM(container).ne.'Annulus' .and.
     &    TRIM(container).ne.'Cuboid' .and.
     &    TRIM(container).ne.'Cube' .and.
     &    TRIM(container).ne.'cuboid' .and.
     &    TRIM(container).ne.'cube') then
         write(6,*)'Incorrect input for container name'
         stop
      endif

      if (TRIM(container).eq.'Cylinder' .or.
     &    TRIM(container).eq.'cylinder') then
         idim = 3
      endif

      write(606,*)
      write(606,*) '***   ',TRIM(prop_name)
      write(606,*) '***   ',TRIM(prop_description)
      write(606,*) '***   ',TRIM(container)
      call date_and_time(VALUES = val)
      if (val(7).ge.10) then
         write(606,'(1x,a,i5,2x,a,i3,2x,a,i3,2x,3(a,i2))')
     &         '***   Yr:',val(1),'Mo:',val(2),'Day:',val(3),
     &         'Time: ',val(5),':',val(6),':',val(7)
      else
         write(606,'(1x,a,i5,2x,a,i3,2x,a,i3,2x,2(a,i2),a,i1)')
     &         '***   Yr:',val(1),'Mo:',val(2),'Day:',val(3),
     &         'Time: ',val(5),':',val(6),':0',val(7)
      endif
      write(606,*)
      write(606,100) idim
      write(606,*)
 100  format(2x,'idim    = ',i8)

      write(606,101)
      write(606,102)
 101  format(2x,'Cut Number',8x,'Name',12x,'Solid',4x,'density',5x,
     &          'w_frac',4x,'dcutoff',4x,'coating')
 102  format(2x,83('-'))
      do iblend=1,iblend_max

         read(500,*) 
         read(500,*) cut_name(iblend)
         read(500,*) cut_type(iblend)
         if (cut_name(iblend).eq.'singlesize') then
            read(500,*) diam_singlesize
            di(iblend,1) = diam_singlesize
            di(iblend,2) = diam_singlesize
            percent(iblend,1) = 100.0
            percent(iblend,2) = 0.0
            nmode(iblend) = 2
         endif
         read(500,*) density(iblend)
         read(500,*) w_frac(iblend)
         read(500,*) dcutoff(iblend)
         read(500,*) coating(iblend)

         write(606,103) iblend,TRIM(cut_name(iblend)),
     &                  TRIM(cut_type(iblend)),
     &                  density(iblend),w_frac(iblend),
     &                  dcutoff(iblend),coating(iblend)
 103     format(2x,'Cut = ',i3,A20,2x,A8,4(2x,f9.5))

         if (TRIM(cut_name(iblend)).ne.'lognormal') then
            cut_name1(iblend)=TRIM('Ingredients_data/')
     &               //TRIM(cut_name(iblend))
         else
            read(500,*) bcoeff(iblend)
            read(500,*) diam_Peak(iblend)
            read(500,*) diam_Min(iblend)
            read(500,*) diam_Max(iblend)
         endif

      enddo
      write(606,104) iblend_max+1,density_binder,1.0-sum(w_frac)
 104  format(2x,'Cut = ',i3,14x,'Binder',6x,'HTPB',2(2x,f9.5))

c
c     compute volume fractions given weight fractions
c       and corresponding densities
c
      iflag = 0
      call w2v
c
      if (v_frac(1).le.1.0d-12) then
         write(606,*)
         write(606,*)'**********************************'
         write(606,*)'ERROR: V_FRAC(1) IS TOO SMALL'
         write(606,*)'       REORDER CUTS'
         write(606,*)'**********************************'
         write(606,*)
         stop
      endif
c
c     generate ingredient table
c
      solidtype = 'Null'
      solidtype(1) = cut_type(1)
      Ntype = 1
      do i=2,iblend_max
         ifill = 1
         do j=1,i-1
            if (cut_type(i).eq.solidtype(j)) ifill = 0
         enddo
         if (ifill.eq.1) then
            Ntype = Ntype + 1
            solidtype(Ntype) = cut_type(i)
         endif
      enddo
c
c
      return
      end
