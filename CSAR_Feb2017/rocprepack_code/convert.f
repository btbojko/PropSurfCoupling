c
c     This fortran code converts the output from Rocprepack
c     and the latest version of Rocpack (Guilherme Amadio)
c     to the previous version of Rocpack (Shane Stafford).
c     This is necessary because Rocfire requires the previous
c     version of Rocpack as input.
c
c     Written: TLJ, 12/13/2016
c
c     To run:
c        pack-ls -D newformat.out
c        pack2xyz pack.out pack.xyz
c        gfortran convert.f
c        ./a.out <pack.xyz
c     This creates rocfire_input, which can then be saved
c     and moved to the appropriate Rocfire folder.
c
      implicit double precision (a-h,o-z)
      dimension x(10000),y(10000),z(10000),r(10000)
      character(90) comments(500),cctrim
      character(90) solidtype
      character(8) elstr,elstrL,dummy
c
c     read in pack from Guilhere's code
c     note that the domain used in Rocpack is [-0.5,0.5] in x,y,z
c     Rocfire requires dmain [-1,1]
c     so scale Guilhere's data by 2 to convert [-0.5,0.5] to [-1,1]
      open(11,file='newpack.out')
      do i=1,10000
         read(11,*,end=90) x(i),y(i),z(i),r(i)
         x(i)=x(i)*2
         y(i)=y(i)*2
         z(i)=z(i)*2
         r(i)=r(i)*2
      enddo
  90  continue
      npts=i-1

      open(14,file='rocprepack_out')
      open(16,file='rocfire_input')

      do i =1,500
         read(14,'(a)')comments(i)
         cctrim = adjustL(comments(i))
         write(16,'(a)')TRIM(comments(i))
         if (len_trim(comments(i)) > 5 .and. cctrim(1:3) /= '***') then
            read(cctrim,*) elstr,alp_V_ele
            elstrL = adjustL(elstr)
            if (index(elstr,'Binder') > 0) then
               goto 99
            endif
         endif
      enddo
  99  continue
      read(14,*) ddlx
      read(14,*) ddlz
      read(14,*) ddly
      write(16,101) ddlx
      write(16,102) ddlz
      write(16,103) ddly
 101  format(2x,f19.8,' ! Length')
 102  format(2x,f19.8,' ! Width')
 103  format(2x,f19.8,' ! Height')
c
c     compute number of spheres, including overlaps
c
      isum = 0
      do kx=1,3
        xshift = -2.0d0 + dfloat(kx-1)*2.0d0
       do ky = 1,3
          yshift = -2.0d0 + dfloat(ky-1)*2.0d0
          do kz = 1,3
            zshift = -2.0d0 + dfloat(kz-1)*2.0d0
            do i=1,npts
              xout = x(i) + xshift
              yout = y(i) + yshift
              zout = z(i) + zshift
              rout = r(i)
              if ((abs(xout)).le.1.2d0 .and.
     &            (abs(yout)).le.1.2d0 .and.
     &            (abs(zout)).le.1.2d0) then
                  isum = isum + 1
              end if
            end do
          end do
        end do
      end do

      write(16,104) isum,npts
 104  format(2x,2i6,' ! Number of spheres (Total,Inside)') 
      read(14,*) TheoreticalPackingDensity !theoretical value
      read(14,*) rhoV_Rocpack
      write(16,105) rhoV_Rocpack,rhoV_Rocpack
      write(16,106) TheoreticalPackingDensity
 105  format(2x,2f19.8,' ! packing fraction')
 106  format(2x,f19.8,' ! Theoretical packing density')
      read(14,*) ncircles
      read(14,*) nmode
      write(16,107) nmode
 107  format(2x,i6)
      do i = 1, nmode
         read(14,*) n1,n2,dummy,diam
         write(16,108) n1,TRIM(dummy),n2,diam
      enddo
 108  format(2x,i6,a9,i6,f19.9)

      solidtype='AP'

      isum = 0
      do i=1,npts
      if (i.le.40) write(6,*) i,r(i)
      do kx=1,3
        xshift = -2.0d0 + dfloat(kx-1)*2.0d0
       do ky = 1,3
          yshift = -2.0d0 + dfloat(ky-1)*2.0d0
          do kz = 1,3
            zshift = -2.0d0 + dfloat(kz-1)*2.0d0
              xout = x(i) + xshift
              yout = y(i) + yshift
              zout = z(i) + zshift
              rout = r(i)
              if ((abs(xout)).le.1.2d0 .and.
     &            (abs(yout)).le.1.2d0 .and.
     &            (abs(zout)).le.1.2d0) then
                  isum = isum + 1
                  write(16,100) isum,i,TRIM(solidtype),
     &               xout,yout,zout,rout
              end if
            end do
          end do
      end do
      end do
 100  format(1x,2i8,2x,a5,2x,19f14.6)
      write(6,*)'isum=',isum
c
c     write to file for rocfire input
c
c
      stop
      end
