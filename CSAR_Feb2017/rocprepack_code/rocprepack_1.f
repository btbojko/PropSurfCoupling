c     //////////////////////////////////////////////////
c     Author: Thomas L. Jackson
c     //////////////////////////////////////////////////
c     
c     ***************************************************************
c     PROPRIETARY CODE - DO NOT DISTRIBUTE OR USE THIS CODE WITHOUT
c     EXPLICIT PERMISSION TO DO SO
c     
c     ***************************************************************
c     
c     
c     _______________________________________________________________
c     
c     Thomas L. Jackson
c     
c     16 June 2016
c     _______________________________________________________________
c     
c     Original release date : April, 2005
c     _______________________________________________________________
c     
c     Permission is granted to use this code under the stipulations
c     
c     1. The user will not distributed or sell this code by
c     any means to other parties for any reason whatsoever
c     2. The user agrees to acknowledge the following archival
c     publications in all written manuscripts/publications
c     that are based on the use of this code:
c     
c     [1] G.M. Knott, T.L. Jackson & J. Buckmaster (2001). 
c         "Random packing of heterogeneous propellants". 
c         AIAA Journal, Vol. 39(4), pp. 678-686.
c     [2] S. Kochevets, J. Buckmaster, T.L. Jackson & A. Hegab (2001). 
c         "Random propellant packs and the flames they support."
c         AIAA J. of Propulsion and Power, Vol. 17, pp. 883-891.
c     [3] F. Maggi, S. Stafford, T.L. Jackson & J. Buckmaster (2008).
c         "Nature of packs used in propellant modelling."
c         Physical Review E, Vol. 77(4), 046107.
c     
c     _______________________________________________________________
c     
c     The author of this code assumes no responsibility whatsoever
c     in the use or missuse of this code, or in any results obtained
c     from using this code. The user of this code assumes all risks.
c     _______________________________________________________________
c     
c     BRIEF CODE DESCRIPTION:
c     =======================
c     
c     Code determines distribution of solid particles
c     needed for the in-house packing code, Rocpack
c     
c     Dependencies:
c     1. rocprepack_1.f  
c        - main code
c     2. rocprepack_2a.f
c        - subroutine cuts
c        - subroutine consolidate
c        - subroutine output
c        - subroutine mix_and_rnd
c        - subroutine num_fraction
c        - subroutine ludcmp, lubksb
c     3. rocprepack_2b.f
c        - input_pack_setup
c     4. rocprepack_2c.f
c        - subroutine w2v
c        - subroutine lognormal
c     5. rocprepack_3.h  
c        - variables declaration file
c     6. Ingredients_data 
c        - directory containing predefined
c           propellant distributions
c     7. Formulations_data
c        - directory containing predefined
c          propellant morphologies
c     8. filenames.list
c        - contains names found in Formulations_data
c     
c     ***************************************************************
c     
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     
c     
c     
c     NONMENCLATURE
c     
c     
c     Global Parameters
c     -----------------
c     
c     max_cut = 20; maximum number of cuts that can be combined
c        to generate a propellant matix; most propellants
c        are either bimodal or trimodal AP with aluminum
c        (this results in only 4 cuts maximum)
c     max_particles = 500; maximum number of different diameters
c        within each cut
c     max_solid_type = 20; maximum number of different solid types;
c        e.g., AP, Al, RDX, HMX, ADN
c     
c     
c     User Supplied Parameters
c     ------------------------
c     
c     From input file "filenames.list":
c
c     inputfile = name of predefined propellant morphology
c        in Formulations_data
c
c     container = 'cuboid', 'cylinder', or 'annulus'
c
c     idim = 2 (disks) or (3) spheres
c        Currently idim=2 is disabled since current packing
c        code only computes 3D packs
c
c     height = height of the 3D domain:
c        If (height>0) height given in microns
c        If (height<0) height is the nondimensional aspect ratio
c
c        For a cuboid, the total volume is given by V_T = L_x*L_y*L_z
c           with the domain always a square in the x,z plane.
c        Thus, L_y = height, L_x = L_z width in x,z direction.
c        If (height>0), L_y = height, volume = L_x^2*height
c           so that L_x = L_z = sqrt(volume/height)
c        If (height<0), L_y = |height|*L_x,
c           so that volume = L_x^3*|height|, thus
c           L_x = L_z = (volume/|height|)**(1/3), L_y = height*L_x
c
c        For a cylinder, the total volume is given by
c        V_T = pi*R^2*h, where R is the radius of the
c        cylinder, h is the height. 
c        If (height>0), L_y = height, L_x the radius, then
c           volume = pi*L_x^2*height
c           so that L_x = sqrt(volume/pi/height)
c        If (height<0), let L_x be radius, then L_y = 2*|height|*L_x,
c           so that volume = 2*pi*L_x^3*|height|, thus
c           L_x = (volume/2/pi/|height|)**(1/3), L_y = 2*|height|*L_x
c        Note that L_z is not used.
c
c        For an annulus, the input file has two rows:
c          L_y = height     in microns
c          L_x = rad_inner  in microns
c          Now let L_z = outer radius (microns)
c          The volume is then V = pi*(L_z^2-L_x^2)*L_y
c          so that L_z = sqrt(L_x^2 + V/(pi*L_y))
c     
c     
c        iblend_max = number of cuts to be combined to form a
c           propellant matrix
c     
c        ntot = total number of particles; due to roundoff (we can
c           only have a integer number of particles), the final
c           number of spheres (3D) might not exactly match ntot
c
c        density_binder = density of binder (gm/cc)
c
c
c        cut_name = name of propellant cut; predefined
c          cuts located in directory Ingredinents_data
c
c        cut_type(i) = type of solid
c           Examples: AP, Al, PETB, HMX, ADN
c
c        density(i) = density for cut_type(i) (gm/cc)
c
c        w_frac(i) = weight fraction for cut_type(i)
c          Note: all values are in weight precent; the code then
c          converts these to volume percent using the subroutine
c          w2v, which generates v_frac(i)
c     
c        dcutoff(i) = cut off diameter for cut_type(i);
c           smallest diameter to be considered for the pack 
c
c        coating(i) = coating parameter for cut_type(i);
c           simulates binder deposit over large
c           oxidizer particles. coating is a *radial* displacement: diam/2
c           + coating is the acutal contact radius of the spheres. Note that a
c           too large coating reduces the packing fraction achievable using
c           Rocpack.
c           Set to 0; current code has this disabled
c     
c     
c     Internal Parameters and Variables
c     ---------------------------------
c
c     rho_v = theoretical packing fraction, and assumed to
c           be equal to the total solid volume fraction; due to roundoff this
c           value may or may not be acheived in practice. For single cuts, rho
c           must be specified by the user. For multiple cuts, rho is
c           determined internally using the equation rho = sum (i=1 to
c           iblend_max) v_frac(i)
c     
c     iflag = 0; uses volume to compute number fractions
c           = 1; uses diameter to compute number fraction
c           code hard wires iflag=0
c     
c     iblend = counter, 1 <= iblend <= iblend_max
c     i = counter, 1 <= i <= nmode(iblend)
c     nmode(iblend) = number of different particles for cut iblend
c     diam(iblend,i) = diameter
c     Vfrac(iblend,i) = percent(iblend,i)/100; volume fraction
c     dnj(iblend,i) = number fraction before consolidation
c     
c     dN(iblend)  = number of particles after consolidation for cut
c     iblend
c     Ni = delNj_final(iblend,i) + 0.5; integer value after roundoff
c     
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     
c     
c     Predefined propellant cuts stored in tables
c     (located in the directory Ingredients_data)
c     
c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
c     
c     
c     NOTE: CODE MUST BE COMPILED USING DOUBLE PRECISION
c     
c     
c     
      include 'rocprepack_3.h'
c     
c     local variables
c     
      integer iblend,i,iout
      double precision sum1,w,y,rho_v_1
      double precision Dtemp(5000)
      double precision dNtemp(5000)
      double precision dDtemp(5000)
      double precision delNj(5000)
      double precision Vtemp(5000)
c     
c     Create global output file for relevant user information
c     
      open(606,file='output_info')
c     
      call input_pack_setup     ! subroutine in rocprepack_2b.f
c     
      pi = acos(-1.0d0)
c     
c     
      write(606,102)
      do iblend = 1,iblend_max
         write(606,*)
         write(606,*)'*** starting cut = ',iblend
         write(606,*)
c     
         if (cut_name(iblend).ne.'lognormal') then
            if (cut_name(iblend).ne.'singlesize') then
               open(505,file=cut_name1(iblend))
               do i=1,200
                  read(505,*,end=5) di(iblend,i),percent(iblend,i)
                  if (percent(iblend,i).eq.0.0) goto 5
               enddo
 5             continue
               close(505)
               nmode(iblend) = i-1
            endif
            sum1 = 0.0
            do i=1,nmode(iblend)
               if (ivag.eq.0) then
                  diam(iblend,i) = (di(iblend,i)+di(iblend,i+1))/2.0
                  del_d(iblend,i) = di(iblend,i)-di(iblend,i+1)
               endif

               if (ivag.eq.1) then
                  diam(iblend,i) = di(iblend,i)
                  del_d(iblend,i) = 0.
               endif

               Vfract(iblend,i) = percent(iblend,i)/100.0
               if (diam(iblend,i).le.dcutoff(iblend)) goto 6
               if (percent(iblend,i).eq.0.0) goto 6
               sum1 = sum1 + Vfract(iblend,i)
            enddo
 6          continue
            nmode(iblend) = i-1
         else
            call lognormal(iblend)
            sum1 = 0.0
            do i=1,nmode(iblend)
               diam(iblend,i) = di(iblend,i)
               del_d(iblend,i) = 0.0
               Vfract(iblend,i) = percent(iblend,i)/100.0
               if (diam(iblend,i).le.dcutoff(iblend)) goto 7
               sum1 = sum1 + Vfract(iblend,i)
            enddo
 7          continue
            nmode(iblend) = i-1
         endif
c     
c     Store the volume fraction kept within each cut
c     
         Vfrackept(iblend) = sum1
         if (Vfrackept(iblend).lt.1.0d-12) Vfrackept(iblend)=1.0d-12
         if (Vfrackept(iblend).ge.1.0) Vfrackept(iblend)=1.0
         Vfracdisc(iblend) = 1.0-Vfrackept(iblend)
c     
c     Compute number fractions for all kept
c     particles and for each cut
c     
         call num_fraction(iblend)
c     
c     Write out table
c     
         write(606,10)
 10      format(4x,'Mode',4x,'Mesh Size',2x,'% in channel',
     &        2x,'Mean diam',8x,'nj',14x,'Nj')
         do i=1,nmode(iblend)
            write(606,20) i,di(iblend,i),percent(iblend,i),
     &           diam(iblend,i),dnj(iblend,i),dnj(iblend,i)*ntot
         enddo
         i = nmode(iblend)+1
         write(606,20) i,di(iblend,i),percent(iblend,i)
 20      format(2x,i5,2x,3f12.3,f16.9,f12.3)
         write(606,*)
         write(606,*)'*** completeing cut = ',iblend
         write(606,*)
      enddo
      write(606,102)
 102  format(75('x'))
c     
c     Write out information based on fractional number of particles
c     
      rho_v_1 = sum(v_frac*Vfrackept)
      w = rho_v_1/rho_v
      y = (1.0-w)/w
      alpha = 0.0
      if (rho_v_1.le.0.999) alpha = (rho_v-rho_v_1)/(1.0-rho_v_1)
      write(606,*)
      write(606,*)
      write(606,*)'Values based on fractional number of particles:'
      write(606,*)
      do i=1,iblend_max
         write(606,*)'   Cut = ',i,'   ',TRIM(cut_type(i))
         write(606,*)'      Vfrac kept      = ',Vfrackept(i)
         write(606,*)'      Vfrac discarded = ',Vfracdisc(i)
      enddo
      write(606,*)
      write(606,*)'   Fraction of solids kept         = ',w
      write(606,*)'   Fraction of solids discarded    = ',1.0-w
      write(606,*)'   Mixture fraction discarded      = ',y
      write(606,*)
      write(606,*)'   Theoretical packing fraction    = ',rho_v
      write(606,*)'   Fraction of discrete particles  = ',rho_v_1
      write(606,*)'      This is also the updated vol. fraction'
      write(606,*)'   Fraction of fine particles      = ',rho_v-rho_v_1
      write(606,*)'      Vol. of solids to be added to binder'
      write(606,*)
      write(606,*)'   New volume fraction of blend    = ',1.0-rho_v_1
      write(606,*)'      Vol. of binder plus vol. fine particles'
      write(606,*)'   Fraction of blend that is solid = ',alpha
      write(606,*)
c     
c     Compute number of particles for each cut; the number
c     of particles for each cut adds to ntot
c     
      write(606,*)
      write(606,*)'..... calling cuts'
      call cuts
c     
c     Consolidate large particles and print out final solution
c     
      write(606,*)
      write(606,*)'..... calling consolidate'
      call consolidate
c     
c     Print out final solution; information stored is used for
c     the packing codes
c     
      write(606,*)
      write(606,*)'..... calling output'
      call output
c     
c     Print out table if all AP; information printed out
c     is "user friendly", and is not used as part of the
c     packing codes
c     
      iout=0
      do i=1,itotal
         if (solidname(i).ne.'AP') iout=1
      enddo
      if (iout.eq.0) call mix_and_rnd
c     
c     
      stop
      end
