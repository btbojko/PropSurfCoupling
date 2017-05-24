! Heterogeneous Propellant 3D Combustion Code
! Constant Mass Flux
! DOE-ASCI Center for Simulation of Advanced Rockets
! University of Illinois at Urbana-Champaign
! Author: Thomas L. Jackson
! Date last Modified: 16 January 2001
! Parallelized by M. Campbell on May 9, 2001
! //////////////////////////////////////////////////
!
! ********************************************************************
! PROPRIETARY CODE - DO NOT USE WITHOUT PERMISSION OF CSAR!
! ********************************************************************
!
!     ________________________________________________________________
!
!     BRIEF CODE DESCRIPTION:
!     =======================
!
!     Quadrants
!
!     Subroutine breaks the solid up into 4 quadrants
!
! ********************************************************************
SUBROUTINE QUADRANTS

! REVISIT - This whole routine needs some thought.. it cant be done this
!           way now.        
  USE GLOBAL_DATA

  IMPLICIT NONE

!     Local Variables
  INTEGER :: i ,n, j, k, kkk
  INTEGER ::  xs,zs,ys,xe,ze,ye, min_range
  REAL*8  ::   del,rmax, y0, y1,xsc,dxc,eps
  CHARACTER(LEN=21) :: quadfile


  WRITE(quadfile,'("quad",I3.3,".dat")') myid
  if(doquaddump.eq.1) then
     open(UNIT=20,FILE=quadfile,form='formatted')
  endif
!
  xs = drange(1)
  xe = drange(2)-1
  zs = drange(3)
  ze = drange(4)-1
  ys = drange(5)
  ye = drange(6)

  eps=1.d-8

!
! IF the number of cells assigned as input is larger than the minimum
! grid dimension on the particula processor, use the minimum range as 
! number of cells
!
  if(is2D)then
     min_range = max(drange(2),1)
  else
     min_range = max(min(drange(2),drange(4)),1)
  endif
  if(ncells.gt.min_range) then 
     if(myid == 0) write(*,*)'NCELLS CHANGED TO',min_range
     ncells = min_range
  endif
!!!
  rmax=0.0d0
  do n=1,ncircles
     rmax = dmax1(rmax,rad(n))
  enddo

!     Breaks the x-z plane into cells, and defines their
!     looping limits
!
  do n=1,ncells
     xsc=x(xs)
     dxc=x(xe+1)-xsc
     xl(n) = (xsc+dxc*dfloat(n-1)/ncells)
     xm(n) = (xsc+dxc*dfloat(n)/ncells)
     xsc=z(zs)
     dxc=z(ze+1)-xsc
     zl(n) = (xsc+dxc*dfloat(n-1)/ncells)
     zm(n) = (xsc+dxc*dfloat(n)/ncells)
     if(doquaddump.eq.1) then
        write(20,*) n,xl(n),xm(n),zl(n),zm(n)!,z(0),z(1),z(2)
     endif
  enddo
  nxquad(1,1)=0
  nzquad(1,1)=0
  do n=2,ncells
     do i=xs,xe
        if (x(i).le.xl(n)-eps .and.&
             xl(n)-eps.lt.x(i+1)) then
           nxquad(n,1) = i+1
           nxquad(n-1,2) = i
        endif
     enddo
     do i=zs,ze
        if (z(i).le.zl(n)-eps .and.&
             zl(n)-eps.lt.z(i+1)) then
           nzquad(n,1) = i+1
           nzquad(n-1,2) = i
        endif
     enddo
  enddo
  nxquad(ncells,2)=xe+1
  nzquad(ncells,2)=ze+1
  if(doquaddump.eq.1) then
     do n=1,ncells
        write(20,*)n,nxquad(n,1),nxquad(n,2),nzquad(n,1),nzquad(n,2)
     enddo
  endif
!     
!     now put spheres into each cell
!
  del = 1.01*rmax
  if(doquaddump.eq.1) then
     write(20,*)'del = ',del
  endif

  y0=0.0d0
  y1=period           !-2.0d0*dly
  kkk=0
  do i=1,ncells
     do j=1,ncells
        kkk=kkk+1
        k=1
        do n=1,ncircles
           del = 1.01*rad(n)
           if (y1-del.le.ycenter(n) .and.&
                ycenter(n).le.y0+del) then
              if (zl(i)-del.le.zcenter(n) .and.&
                   zcenter(n).le.zm(i)+del .and.&
                   xl(j)-del.le.xcenter(n) .and.&
                   xcenter(n).le.xm(j)+del) then
                 quad(kkk,k)=n
                 nquad(kkk)=k
                 k=k+1
              endif
           endif
        enddo
     enddo
  enddo
!
!
  if(doquaddump.eq.1) then
     write(20,*) ncircles,(nquad(n),n=1,ncells)
  endif
  do j=1,ncells
     do n=1,nquad(j)
        i = quad(j,n)
        if(doquaddump.eq.1) then
           write(20,100) nquad(j),n,i,&
                xcenter(i),zcenter(i),ycenter(i)
        endif
     enddo
  enddo
!     
100 format(2x,3i8,8f14.6)
!
  if(doquaddump.eq.1) then
     close(20)
  endif
!
  return
end SUBROUTINE QUADRANTS
