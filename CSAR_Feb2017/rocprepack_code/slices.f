c
      implicit double precision (a-h,o-z)
      dimension x(50000),y(50000),z(50000),r(50000)
      dimension m(50000)
      character(8) a
c
      write(6,*)'Input starting point'
      read(5,*) iend
      open(unit=98,file='rocfire_input')
      do i=1,iend
         read(98,*)
      enddo
      do i=1,20000
        read(98,*,end=5) i1,i2,a,x(i),z(i),y(i),r(i)
      end do
   5  continue
      close(98)
      npts = i-1
      write(6,*)'npts=', npts
c
      iout = 199
      ncuts = 9
      dzcut = 2.0d0/dfloat(ncuts-1)
      pi = acos(-1.0d0)
c
      do 999 nn=1,ncuts
        j = 0
        iout = iout + 1
        zcut = -1.0d0 + dfloat(nn-1)*dzcut
        area = 0.0d0
        do i=1,npts
          if (abs(z(i)-zcut).lt.r(i)) then
            j = j + 1
            rstar = sqrt(r(i)**2-(z(i)-zcut)**2)
            write(iout,100) j,m(i),x(i),y(i),zcut,rstar,z(i)
            if (abs(x(i)).le.1.0d0.and.abs(y(i)).le.1.0d0) then
              area = area + pi*rstar**2
            end if
          end if
        end do
        total_area = 4.0d0
        ap = area
        binder = total_area - ap
        rho=ap/total_area
        write(300,101) nn,zcut,ap,binder,rho
 999  continue
 100  format(2x,2i8,2x,8f12.6)
 101  format(2x,i8,2x,5f12.6)
c
      ncuts = 1000
      dzcut = 2.0d0/dfloat(ncuts-1)
      pi = acos(-1.0d0)
c
      volume = 0.0d0
      do 997 nn=1,ncuts
        j = 0
        zcut = -1.0d0 + dfloat(nn-1)*dzcut
        area = 0.0d0
        do i=1,npts
          if (abs(z(i)-zcut).lt.r(i)) then
            j = j + 1
            rstar = sqrt(r(i)**2-(z(i)-zcut)**2)
            if (abs(x(i)).le.1.0d0.and.abs(y(i)).le.1.0d0) then
              area = area + pi*rstar**2
            end if
          end if
        end do
        total_area = 4.0d0
        ap = area
        volume = volume + area
        binder = total_area - ap
        rho=ap/total_area
        write(301,101) nn,zcut,ap,binder,rho
 997  continue
      volume = volume*dzcut
      write(6,*)'Volume=',volume
c
c
      stop
      end
