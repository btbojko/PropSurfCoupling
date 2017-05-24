
      implicit none

c     declarations

      integer max_cut,max_particles,max_solid_type

      parameter (max_cut=20)
      parameter (max_particles=500)
      parameter (max_solid_type=20)

      integer iblend_max,mix
      integer Ntot,Ntype,idim
      integer val(8)
      integer iflag
      integer itotal
      integer nmode,icut
      integer ivag

      double precision rho_w,rho_v
      double precision density,density_binder
      double precision height
      double precision V_total,V_solids,V_binder
      double precision alpha,sample_mass,sample_density
      double precision dlx,dly,dlz
      double precision pi
      double precision di,diam,percent,dnj,Vfrackept,Vfracdisc,
     &     del_d,Vfract,dN_final,D_final,Dj,dn,dcutoff,
     &     Vkept,diam_Min,diam_Max,diam_Peak,bcoeff,
     &     w_frac,v_frac,coating
      double precision rad_inner

c     dimensions

      dimension di(max_cut,max_particles)
      dimension diam(max_cut,max_particles)
      dimension percent(max_cut,max_particles)
      dimension dnj(max_cut,max_particles)
      dimension nmode(max_cut)

      dimension Vfrackept(max_cut)
      dimension Vfracdisc(max_cut)

      dimension del_d(max_cut,max_particles)
      dimension Vfract(max_cut,max_particles)

      dimension dN_final(max_particles)
      dimension D_final(max_particles)

      dimension Dj(max_cut,max_particles)
      dimension dn(max_cut)
      dimension dcutoff(max_cut)
      dimension coating(max_cut)
      dimension Vkept(max_cut)
      dimension icut(max_particles)

      dimension w_frac(max_cut)
      dimension v_frac(max_cut)
      dimension density(max_cut)

      dimension diam_Min(max_cut)
      dimension diam_Max(max_cut)
      dimension diam_Peak(max_cut)
      dimension bcoeff(max_cut)

      character*(200) prop_name
      character*(200) inputfile
      character*(200) prop_description
      character*(200) container
      character*(200) mix_name
      character*(200) cut_name(max_cut)
      character*(200) cut_name1(max_cut)
      character*(200) cut_type(max_cut)
      character*(200) solidtype(max_cut)
      character*(200) solidname(max_particles)

c     common blocks

      common /coeff00/ prop_name,prop_description,container,mix_name,
     &                 cut_name,cut_name1,cut_type,solidtype,
     &                 solidname

      common /coeff01/ iblend_max,mix,Ntot,Ntype,idim,val,iflag,
     &                 itotal,nmode,icut,ivag

      common /coeff02/ rho_w,rho_v,height,V_total,V_solids,V_binder,
     &                 alpha,dlx,dly,dlz,pi,sample_mass,sample_density,
     &                 rad_inner

      common /coeff03/ v_frac,w_frac,density,density_binder

      common /coeff04/ diam_Min,diam_Max,diam_Peak,bcoeff

      common /coeff05/ dN_final,D_final,Dj,dn,dcutoff,Vkept,
     &                 di,diam,percent,dnj,Vfrackept,Vfracdisc,
     &                 del_d,Vfract,coating
