!*****************************************************************************8
subroutine VISUALIZE_FIELD
USE GLOBAL_DATA
USE MYMPI
IMPLICIT NONE
!------------------------------------------------------

IF(mod(ncyc,frame) /= 0) RETURN

IF(ipack == -1) THEN
    CALL PRINT_1D
ELSE
    CALL PRINT_2D
ENDIF

!------------------------------------------------------
RETURN
END subroutine VISUALIZE_FIELD
!*****************************************************************************8

!*****************************************************************************8
subroutine PRINT_2D

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

  INTEGER :: i,j,eqn,jStar,i1,i2
  REAL*8 :: xx,yy, yys, x_loc, interp, percentRate, oneDrate,finalspec,ave_y,brate,ave_b
!!!
  CHARACTER(LEN=21) :: fieldsfile,solids,surfaces
  CHARACTER(LEN=6) :: title
!!!
!---------------------------------------

  ALLOCATE(G_x(0:nx),G_phi(0:nx),G_speed(0:nx))
  ALLOCATE(G_f(-1:nx+1,-1:ny+1,1:neqmax),G_rate(0:nx,0:ny,1:irxn))
  ALLOCATE(G_q(-1:nx+1,-1:ny+1,1:3),G_psi(-1:nx+5,0:ny+5))


  CALL GATHER_FINAL_SOLUTION

!!!
  WRITE(fieldsfile,'("f",I6.6,".vtk")') ncyc
  WRITE(solids,'("s",I6.6,".vtk")') ncyc
  WRITE(surfaces,'("surf",I6.6,".vtk")') ncyc
!!!

  if(myid.eq.0) then
     close(13)
     open (UNIT=13,FILE=fieldsfile,STATUS='UNKNOWN')
     open (UNIT=25,FILE=solids,STATUS='UNKNOWN')
     open (UNIT=95,FILE=surfaces,STATUS='UNKNOWN')



! -- write fields.dat :: unit = 10
     write(13,1000) nx+1,ny+1,1,(nx+1)*(ny+1)
     write(25,1000) nx+1,ny+1,1,(nx+1)*(ny+1)
     do j = 0, ny
        do i = 0, nx
           xx = G_x(i)*10000.0d0
           yy = (y(j) + G_phi(i))*10000.0d0
           yys = (-y(j) + G_phi(i))*10000.0d0
!           write(13,900) xx, yy, (G_f(i,j,eqn), eqn=1,neqgas),&
!                (G_rate(i,j,eqn), eqn=1,3), (G_q(i,j,eqn), eqn=1,2)
           write(13,1001) xx,yy,zero
!           write(25,900) xx, yys, G_psi(i,j)
           write(25,1001) xx,yys,zero
        end do
     end do

     write(13,1002) (nx+1)*(ny+1)
     write(25,1002) (nx+1)*(ny+1)
     do eqn = 1,neqgas+1
        if(eqn.eq.1) then 
            title = '  Temp'
        elseif(eqn.eq.2) then 
            title = 'HTPB'
        elseif(eqn.eq.3) then 
            title = 'AP'
        elseif(eqn.eq.4)then
            title = 'MONO'
        elseif(eqn.eq.5)then
            title = 'BIND'     
        elseif(eqn.eq.6)then
            title = 'PRIMARY'  
        else
            title = 'FINAL'  
        endif
        if(eqn == 7)then
           write(13,1003) title
            do j = 0, ny
                do i = 0, nx
                    finalspec = 1.0d0 - sum(G_f(i,j,2:neqgas))
                    if(abs(finalspec) < 1.0d-30)then
                        write(13,1004) 0.0d0
                    else
                        write(13,1004) finalspec
                    endif
                enddo
            enddo 
        else
            write(13,1003) title
            do j = 0, ny
                do i = 0, nx
                    if(abs(G_f(i,j,eqn))< 1.0d-30)then
                        write(13,1004) 0.0d0
                    else
                        write(13,1004) G_f(i,j,eqn)
                    endif    
                enddo
            enddo
        endif     
     enddo
     write(25,1003) ' rho_s'
     do j = 0, ny
        do i = 0, nx
           write(25,1004) G_psi(i,j)
        enddo
     enddo
     write(25,1003) 'Temp'
     do j = 0, ny
        do i = 0, nx
           write(25,1004)  G_f(i,j,neqmax)
        enddo
     enddo

! -- write surface.dat :: unit = 95
     write(95,1000) nx+1,1,1,nx+1
     do i = 0, nx
        write(95,1001) G_x(i)*10000.0d0, (y(0) + G_phi(i))*10000.0d0,zero
     enddo
     write(95,1002) (nx+1)
     write(95,1003) 's_Temp'
     do i = 0, nx
        write(95,1004) G_f(i,0,neqmax)
     enddo
     write(95,1003) 'B_Rate'
     do i = 0, nx
        write(95,1004) G_speed(i)
     enddo

! -- write BurningRate.txt, heartbeat output
     !if(ncyc >= ncyc_constant) then
        ave_y = zero
        ave_b = zero
        do i = 0,nx
            ave_y = ave_y + (y(0) + G_phi(i))
            ave_b = ave_b + G_speed(i)
        enddo
        ave_y = ave_y/(nx+1)
        ave_b = ave_b/(nx+1)
        if(ncyc <= ncyc_constant)then
            brate = ave_b
        else
            brate = -ave_y/(tcyc - dt*ncyc_constant)
        endif
        write(hb,1005) ncyc,dt,tcyc,ave_y,brate,ave_b
     !endif
     close(8)
     close(13)
     close(25)
     close(95)

  endif


  DEALLOCATE(G_x,G_phi,G_speed,G_f,G_rate,G_q,G_psi)


899 FORMAT(i6,2f14.5)
900 FORMAT(2x,19f14.5)
902 format(2x,12e15.6)
1000 format('# vtk DataFile Version 3.0'/'CSAR_Trial'/'ASCII'/&
     'DATASET STRUCTURED_GRID'/'DIMENSIONS',1x,3i5/'POINTS',1x,i6,1x,'float')
1001 format(e11.5,1x,e12.5,1x,f3.1)
1002 format(/'POINT_DATA',1x,i6)
1003 format('SCALARS',1x,a6,1x,'float'/&
     'LOOKUP_TABLE default')
1004 format(e14.7)
1005 format(i7,',',2x,es10.2,',',2x,f10.6,',',3x,f12.6,',',3x,f10.4,',',3x,f10.4)
!------------------------------------
  RETURN
END subroutine PRINT_2D
!********************************************************************

!*****************************************************************************8
subroutine gather_final_solution

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

  INTEGER :: i,j,eqn
  REAL*8 :: xx,yy, yys
!---------------------------------------

  call gather_x
  call gather_phi
  call gather_speed
  call gather_f
  call gather_q
  call gather_psi

  call vel_print
  call rate_print



  if(myid == 0) write(*,*)'DONE GATHERING',ncyc

  RETURN
END subroutine gather_final_solution
!********************************************************************

!********************************************************************
subroutine gather_phi
  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

  integer i, mystart, myend, myxsize
  integer elcount
  integer mysize
  INTEGER, ALLOCATABLE  :: allsize(:),displacements(:)
!--------------------------------------------------------

  allocate(myphi(0:nx),allphi(0:nx))
  allocate(allsize(nproc+1),displacements(nproc+1))

  mystart = drange(1)
  myend = drange(2)
  if(myid == nproc-1) myend = drange(2)+1
  myxsize = (myend - mystart) + 1
  mysize = myxsize

  if(mystart.eq.0) then
     mystart = 0
  endif
  if(myend.eq.nx) then 
     myend = nx
  endif
  myxsize = (myend - mystart) + 1
  mysize = myxsize         

  call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allsize,1,MPI_INTEGER,&
       comm3d,ierr)

  displacements(1) = 0
  do i = 2,nproc
     displacements(i) = allsize(i-1) + displacements(i-1)
  enddo

  elcount = 0
  do i = mystart,myend
     myphi(elcount)     = phi(i,0)
     elcount = elcount + 1
  enddo
  call MPI_ALLGATHERV(myphi(0),mysize,MPI_DOUBLE_PRECISION,&
       allphi(0),allsize,displacements,MPI_DOUBLE_PRECISION,&
       comm3d,ierr)
  elcount = 0
  do i = 0,nx
     G_phi(i) = allphi(elcount)
     elcount = elcount + 1
  enddo

  deallocate(myphi,allphi)
  deallocate(allsize,displacements)
  return
end subroutine gather_phi
!*****************************************************************************

!*****************************************************************************
subroutine gather_f

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

!     Local variables
  integer i, j, n, myxsize,pcs
  integer elcount,myysize,el
  integer mysize,elm,elsize,lastd
  integer mydrange(4)
  INTEGER, ALLOCATABLE  :: allsize(:),displacements(:),alldrange(:)

!---------------------------------------------------------------

  allocate(allsize(nproc+1),displacements(nproc+1), alldrange(4*nproc))


  mydrange(1) = drange(1)
  mydrange(2) = drange(2)
  if(myid == nproc-1) mydrange(2) = drange(2)+1

  mydrange(3) = 0
  mydrange(4) = ny
  if(mydrange(1).eq.0) mydrange(1) = 0
  if(mydrange(2).eq.nx) mydrange(2) = nx 
  if(mydrange(3).eq.0) mydrange(3) = 0
  if(mydrange(4).eq.ny) mydrange(4) = ny
  myxsize = (mydrange(2) - mydrange(1)) + 1
  myysize = (mydrange(4) - mydrange(3)) + 1
  mysize = myxsize * myysize
  call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allsize,1,MPI_INTEGER,&
       comm3d,ierr)
  call MPI_ALLGATHER(mydrange,4,MPI_INTEGER,alldrange,4,MPI_INTEGER,&
       comm3d,ierr)
  displacements(1) = 0
  do i = 2,nproc
     displacements(i) = allsize(i-1) + displacements(i-1)
  enddo

  lastd =  allsize(nproc) + displacements(nproc)

!--- ALLOCATION
  allocate(myphi(0:mysize+1), allphi(0:lastd+1))

  do pcs = 1, nproc-1
     elm = 1 + 4 * (pcs-1)
     el =  1 + 4 * pcs
     elsize = alldrange(el+1)-alldrange(el+0)
     alldrange(el+0) = alldrange(elm+1) + 1
     alldrange(el+1) = alldrange(el+0) + elsize
  enddo

  do n = 1,neqmax
     elcount = 0
     do i = mydrange(1),mydrange(2)
        do j = mydrange(3),mydrange(4)
           myphi(elcount) = f(i,0,j,n)
           elcount = elcount + 1
        enddo
     enddo
     call MPI_ALLGATHERV(myphi(0),mysize,MPI_DOUBLE_PRECISION,&
          allphi(0),allsize,displacements,MPI_DOUBLE_PRECISION,&
          comm3d, ierr)

     elcount = 0
     do pcs = 0, nproc-1
        el = 1 + (pcs*4)
        do i = alldrange(el),alldrange(el+1)
           do j = alldrange(el+2),alldrange(el+3)
              G_f(i,j,n) = allphi(elcount)
              elcount = elcount + 1 
           enddo
        enddo
     enddo
  enddo


  deallocate(myphi,allphi)
  deallocate(allsize,displacements)

  return
end subroutine gather_f
!**********************************************************


!**********************************************************
subroutine gather_x

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

!  Local variables
  integer i, mystart, myend, myxsize
  integer mysize
  INTEGER, ALLOCATABLE  :: allsize(:),displacements(:)
!--------------------------------------------------------

  allocate(myphi(0:nx+1),allphi(0:nx+1))
  allocate(allsize(nproc+1),displacements(nproc+1))


  mystart = drange(1)
  myend = drange(2)
  if(myid == nproc-1) myend = drange(2)+1
  myxsize = (myend - mystart) + 1
  mysize = myxsize

  call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allsize,1,MPI_INTEGER,&
       comm3d,ierr)
  displacements(1) = 0
  do i = 2,nproc
     displacements(i) = allsize(i-1) + displacements(i-1)
  enddo
  do i = mystart,myend
     myphi(i-mystart) = x(i)
  enddo

  call MPI_ALLGATHERV(myphi(0),mysize,MPI_DOUBLE_PRECISION, &
       allphi(0),allsize,displacements,MPI_DOUBLE_PRECISION, &
       comm3d,ierr)
  do i = 0,nx
     G_x(i) = allphi(i)
  enddo


  deallocate(myphi,allphi)
  deallocate(allsize,displacements)
!____________________________________________________
  return
end subroutine gather_x
!**************************************************************

!********************************************************************
subroutine gather_speed
  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

  integer i, mystart, myend, myxsize
  integer elcount
  integer mysize
  INTEGER, ALLOCATABLE  :: allsize(:),displacements(:)
!--------------------------------------------------------

  allocate(myphi(0:nx),allphi(0:nx))
  allocate(allsize(nproc+1),displacements(nproc+1))

  mystart = drange(1)
  myend = drange(2)
  if(myid == nproc-1) myend = drange(2)+1
  myxsize = (myend - mystart) + 1
  mysize = myxsize

  if(mystart.eq.0) then
     mystart = 0
  endif
  if(myend.eq.nx) then 
     myend = nx
  endif
  myxsize = (myend - mystart) + 1
  mysize = myxsize         

  call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allsize,1,MPI_INTEGER,&
       comm3d,ierr)

  displacements(1) = 0
  do i = 2,nproc
     displacements(i) = allsize(i-1) + displacements(i-1)
  enddo

  elcount = 0
  do i = mystart,myend
     myphi(elcount)  = phit(i,0)
     elcount = elcount + 1
  enddo
  call MPI_ALLGATHERV(myphi(0),mysize,MPI_DOUBLE_PRECISION,&
       allphi(0),allsize,displacements,MPI_DOUBLE_PRECISION,&
       comm3d,ierr)
  elcount = 0
  do i = 0,nx
     G_speed(i) = allphi(elcount)
     elcount = elcount + 1
  enddo

  deallocate(myphi,allphi)
  deallocate(allsize,displacements)
  return
end subroutine gather_speed
!*****************************************************************************


!*****************************************************************************
subroutine gather_q

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

!     Local variables
  integer i, j, n, myxsize,pcs
  integer elcount,myysize,el
  integer mysize,elm,elsize, lastd
  integer mydrange(4)
  INTEGER, ALLOCATABLE  :: allsize(:),displacements(:),alldrange(:)

!---------------------------------------------------------------

  allocate(allsize(nproc+1),displacements(nproc+1), alldrange(4*nproc))


  mydrange(1) = drange(1)
  mydrange(2) = drange(2)
  if(myid == nproc-1) mydrange(2) = drange(2)+1

  mydrange(3) = 0
  mydrange(4) = ny
  if(mydrange(1).eq.0) mydrange(1) = 0
  if(mydrange(2).eq.nx) mydrange(2) = nx 
  if(mydrange(3).eq.0) mydrange(3) = 0
  if(mydrange(4).eq.ny) mydrange(4) = ny
  myxsize = (mydrange(2) - mydrange(1)) + 1
  myysize = (mydrange(4) - mydrange(3)) + 1
  mysize = myxsize * myysize
  call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allsize,1,MPI_INTEGER,&
       comm3d,ierr)
  call MPI_ALLGATHER(mydrange,4,MPI_INTEGER,alldrange,4,MPI_INTEGER,&
       comm3d,ierr)
  displacements(1) = 0
  do i = 2,nproc
     displacements(i) = allsize(i-1) + displacements(i-1)
  enddo

  lastd =  allsize(nproc) + displacements(nproc)

!--- ALLOCATION
  allocate(myphi(0:mysize+1), allphi(0:lastd+1))

  do pcs = 1, nproc-1
     elm = 1 + 4 * (pcs-1)
     el =  1 + 4 * pcs
     elsize = alldrange(el+1)-alldrange(el+0)
     alldrange(el+0) = alldrange(elm+1) + 1
     alldrange(el+1) = alldrange(el+0) + elsize
  enddo

  do n = 1,ndim
     elcount = 0
     do i = mydrange(1),mydrange(2)
        do j = mydrange(3),mydrange(4)
           myphi(elcount) = q(i,0,j,n)
           elcount = elcount + 1
        enddo
     enddo
     call MPI_ALLGATHERV(myphi(0),mysize,MPI_DOUBLE_PRECISION,&
          allphi(0),allsize,displacements,MPI_DOUBLE_PRECISION,&
          comm3d, ierr)

     elcount = 0
     do pcs = 0, nproc-1
        el = 1 + (pcs*4)
        do i = alldrange(el),alldrange(el+1)
           do j = alldrange(el+2),alldrange(el+3)
              G_q(i,j,n) = allphi(elcount)
              elcount = elcount + 1 
           enddo
        enddo
     enddo
  enddo

  deallocate(myphi,allphi)
  deallocate(allsize,displacements)

!-----------------------------------------------
  return
end subroutine gather_q
!**********************************************************

!*****************************************************************************
subroutine gather_psi

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

!     Local variables
  integer i, j, kplot, n, myxsize,pcs,nghost
  integer elcount,myysize,el
  integer mysize,elm,elsize, lastd
  integer mydrange(4)
  INTEGER, ALLOCATABLE  :: allsize(:),displacements(:),alldrange(:)

!---------------------------------------------------------------

  allocate(allsize(nproc+1),displacements(nproc+1), alldrange(4*nproc))


  kplot = 0
  mydrange(1) = drange(1)
  mydrange(2) = drange(2)
  if(x(mydrange(2)) >= period/two - two*dx) &
       mydrange(2) = drange(2)+1

  mydrange(3) = 0
  mydrange(4) = ny
  if(mydrange(1).eq.0) mydrange(1) = 0
  if(mydrange(2).eq.nx) mydrange(2) = nx 
  if(mydrange(3).eq.0) mydrange(3) = 0
  if(mydrange(4).eq.ny) mydrange(4) = ny
  if(myid == 0)then             !Remove ghost nodes from x direction 
     write(6,*) "*** myid = ***", myid         
     myxsize = (mydrange(2) - mydrange(1)) + 1
    ! myxsize = (mydrange(2) - mydrange(1))
  else
     myxsize = (mydrange(2) - mydrange(1))
  endif
  myysize = (mydrange(4) - mydrange(3)) + 1
  mysize = myxsize * myysize

  call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allsize,1,MPI_INTEGER,&
       comm3d,ierr)
  call MPI_ALLGATHER(mydrange,4,MPI_INTEGER,alldrange,4,MPI_INTEGER,&
       comm3d,ierr)

  displacements(1) = 0
  do i = 2,nproc
     displacements(i) = allsize(i-1) + displacements(i-1)
  enddo
  lastd =  allsize(nproc) + displacements(nproc)
!--- ALLOCATION
  allocate(myphi(0:mysize+1), allphi(0:lastd+1))

  do pcs = 1, nproc-1
     elm = 1 + 4 * (pcs-1)
     el =  1 + 4 * pcs
     elsize = alldrange(el+1)-alldrange(el+0)
     alldrange(el+0) = alldrange(elm+1) + 1
     alldrange(el+1) = alldrange(el+0) + elsize
  enddo

  do n = 1,1
     elcount = 0
     nghost =  mydrange(1)
     if(nproc /= 1 .and. myid /= 0) nghost =  mydrange(1) + 1
     do i = nghost,mydrange(2)
        do j = mydrange(3),mydrange(4)
           myphi(elcount) = rhos(i,0,j)
           elcount = elcount + 1
        enddo
     enddo
     call MPI_ALLGATHERV(myphi(0),mysize,MPI_DOUBLE_PRECISION,&
          allphi(0),allsize,displacements,MPI_DOUBLE_PRECISION,&
          comm3d, ierr)
   
     elcount = 0
     do pcs = 0, nproc-1
        el = 1 + (pcs*4)
        do i = alldrange(el),alldrange(el+1)
           do j = alldrange(el+2),alldrange(el+3)
              G_psi(i,j) = allphi(elcount)
              elcount = elcount + 1
           enddo
        enddo
     enddo
  enddo

  deallocate(myphi,allphi)
  deallocate(allsize,displacements)

!-----------------------------------------------
  return
end subroutine gather_psi
!**********************************************************


!**********************************************************!
SUBROUTINE PRN_GLOB_VARS

  USE GLOBAL_DATA

  IMPLICIT NONE

  INTEGER :: w_unit
!-------------------------------------------------

  if(myid /= nproc-1) RETURN

  w_unit = 12

  if(w_unit /= 6) open(UNIT=w_unit,FILE='run_parameters.dat',form='formatted') 

  write(w_unit,*) 'NUMBER OF POINTS'
  write(w_unit,*) nx
  write(w_unit,*) nz
  write(w_unit,*) ny
  write(w_unit,*) 'INCREMENTS'
  write(w_unit,*) dx
  write(w_unit,*) dz
  write(w_unit,*) dy
  write(w_unit,*) 'Smallest dy',y(1)-y(0)
  write(w_unit,*) '-------------'
  write(w_unit,*) 'DRANGE IS',drange,'PROC--',myid
  write(w_unit,*) 'XEND,ZEND,YEND'
  write(w_unit,*) xend,myid
  write(w_unit,*) zend,myid
  write(w_unit,*) yend,myid
  write(w_unit,*) 'C1Y'
  write(w_unit,*) c1y,myid

  write(w_unit,*)'MG parameters'
  write(w_unit,*)'DX',dx_mg
  write(w_unit,*)'DZ',dz_mg
  write(w_unit,*)'DY',dy_mg
  write(w_unit,*) '-------------'
  write(w_unit,*) 
  write(w_unit,*)'SANDWICH PARAMETERS'
  write(w_unit,*)'XLOC',xloc


  write(w_unit,*) '-------------'
  write(w_unit,*)
  write(w_unit,*)'SOLID PPROPERTIES'
  write(w_unit,*)'RHO_AP,RHO_BIND',rho_ap,rho_binder
  write(w_unit,*)'LAMBDA_AP,LAMBDA_BIND',lambda_ap,lambda_binder

  write(w_unit,*) '-------------'
  write(w_unit,*) ''
  write(w_unit,*) 'SURFACE REACTIONS'
  write(w_unit,*) 'Da_ap = ',da_ap, 'Da_bind = ',da_binder
  write(w_unit,*) 'Theta_ap = ',theta_ap, 'Theta_bind = ',theta_binder
  write(w_unit,*) '-------------'
  write(w_unit,*) 
  write(w_unit,*) 'GAS PHASE REACTIONS'
  write(w_unit,*) 'STOICHIOMETRIC PARAMETERS:'
  write(w_unit,*) 'Stoichimetric mass fraction Map/Mtot',rho_m_st
  write(w_unit,*) 'Stoichimetric mass ratio Map/Mbinder (beta) ',beta
  write(w_unit,*) 'PACK mass ratio Map/Mbinder (beta_pack) ',beta_pack
  write(w_unit,*) 'HOMOGENIZED VOLUME FRACTION (alpha) ',alpha
  write(w_unit,*) 'RHO_1',rho1
  write(w_unit,*) '-------------'
  write(w_unit,*) 
  write(w_unit,*) 'HEAT OF REACTION:'
  write(w_unit,*) 'USING FLAME TEMPERATURES Tap,Tblend(rho_m_st)',th_ap,th_st
  write(w_unit,*) 'FIRST  REACTION, Qg1',qg(1)
  write(w_unit,*) 'SECOND REACTION, Qg2',qg(2)
  write(w_unit,*) 'THIRD  REACTION, Qg3',qg(3)
  write(w_unit,*) '................'
  write(w_unit,*) 'KINETIC PARAMETERS (Calibrated)'
  write(w_unit,*) da(1),'        da1'
  write(w_unit,*) thetay(1),'      theta1'
  write(w_unit,*) np(1),'       n1'
  write(w_unit,*) da(2),'       da2'
  write(w_unit,*) thetay(2),'      theta2'
  write(w_unit,*) np(2),'        n2'
  write(w_unit,*) da(3),'        da3'
  write(w_unit,*) thetay(3),'      theta3'
  write(w_unit,*) np(3),'         n3'
  write(w_unit,*) lewis(1:4),'         Lewis numbers' 
  write(w_unit,*) '-------------'
  write(w_unit,*) 

  
  write(w_unit,*) 'PRESSURE PERTURBATION PARAMETERS'  
  write(w_unit,*) 'BASE pressure',press  
  write(w_unit,*) 'Amplitude',epi_p 
  write(w_unit,*) 'Frequency',omg_p
  write(w_unit,*) 'PHASE',phase_p
  write(w_unit,*) '-------------'
  write(w_unit,*) 


  write(w_unit,*) 'RESTART FLAG ...', irestart
  write(w_unit,*) 'DUMP FILE IN ',TRIM(outputfile)
  if(irestart .eq. 1) write(w_unit,*) 'RESTART FILE IN ',TRIM(restartfile)
  write(w_unit,*) '-------------'
  write(w_unit,*) 
 


  write(w_unit,*) 'TIME DELAYS STEADY-- OSEEN',time_steady_pressure,time_oseen
  write(w_unit,*) 'CYCLES DELAYS STEADY-- OSEEN',ncyc_steady,ncyc_oseen
 


  write(w_unit,*) ''
  write(w_unit,*) 'END of FILE'
  if(w_unit /= 6)  close(w_unit)

!-----------------------------------------------------------
  RETURN
END SUBROUTINE PRN_GLOB_VARS
!***********************************************************

!*****************************************************************************
subroutine TRACK_SURF

  USE GLOBAL_DATA
  USE MYMPI
  IMPLICIT NONE

  INTEGER :: i,j,eqn,iprint
  REAL*8 :: ycut,ycutold
  REAL*8, ALLOCATABLE :: yprint(:)
  CHARACTER*120 filnam
!---------------------------------------

  if(ipack == -1) return

  ycut = phi(0,0) - y(0)
  ycut = mod(ycut,period)
  ycutold=oldphi(0,0) - y(0)
  ycutold=mod(ycutold,period)

  ALLOCATE(yprint(4))
  yprint =(/0.25d0,0.5d0,0.75d0,0.995d0/)
  yprint = yprint*period

  iprint = 0
  do i = 1,ubound(yprint,1)
     if(ycut < yprint(i) .AND. ycutold > yprint(i)) then
        iprint = i
        EXIT
     endif
  enddo


  CALL MPI_BCAST(iprint,1,MPI_DOUBLE_PRECISION,0,comm3d,ierr)

  if(iprint == 0) RETURN


  ALLOCATE(G_x(0:nx),G_phi(0:nx),G_speed(0:nx))
  ALLOCATE(G_f(0:nx,0:ny,1:neqmax),G_rate(0:nx,0:ny,1:irxn),G_q(0:nx,0:ny,1:irxn))
  ALLOCATE(G_psi(0:nx+1,0:ny))

  CALL GATHER_FINAL_SOLUTION 


  if(myid.eq.0) then

     WRITE(filnam,'("surf",I2.2,".dat")') iprint

     OPEN (23,file=TRIM(filnam),form='formatted',action="write")

     write(*,*) TRIM(filnam),ncyc 

     do i = 0, nx
        ycut = y(0) + G_phi(i) - iperiod*period 
        write(23,900) G_x(i), ycut, G_f(i,0,neqmax)
     end do

     close(23)

  endif
  
  DEALLOCATE(G_x,G_phi,G_speed,G_f,G_rate,G_q,G_psi)

  DEALLOCATE(yprint)


900 FORMAT(2x,9f16.7)
902 format(2x,12e15.6)

  if(myid == 0) write(*,*)'DONE TRACK_SURF',ncyc


  RETURN
END subroutine TRACK_SURF
!********************************************************************

!******************************************************************
SUBROUTINE PRINT_1D

  USE GLOBAL_DATA
  IMPLICIT NONE

!---------------------------------------------------------------------
! Dummy Variables
  REAL *8 ::  t,dummy

! Local Variables
  INTEGER ::  P_mpierr, w_unit,w_unit2, j,n,i,k
  CHARACTER*80 filnam,filnam2
!---------------------------------------------------------------------

  call RATE_CANC

  i = 0;k=0
   w_unit = 88
   w_unit2 = 99
   IF (ipack < 0) THEN
      WRITE(filnam,'("ic",I3.3,".dat")') int(press)
      WRITE(filnam2,'("rat",I3.3,".dat")') int(press)
      write(*,*)myid,'PRINT_1D ',ncyc,trim(filnam),'ALPHA',alpha
      write(*,*)ncyc,'BURN RATE',rb(0,0)!,'ALP_W',alp_w
      OPEN (w_unit,file=TRIM(filnam),form='formatted',action="write")
      OPEN (w_unit2,file=TRIM(filnam2),form='formatted',action="write")
      WRITE(w_unit,103) 
      do j = ny,1,-1
         WRITE(w_unit,102) j,-y(j), f(1,0,j,neqmax),zero,zero,zero,zero,zero,zero,rb(i,k)
      enddo
      dummy = 1-sum(f(1,0,j,2:5))

      DO j=0,ny
         WRITE(w_unit,102) j,y(j),f(1,0,j,1:5),dummy,&
              uvel(i,k,j),vcnv(i,k,j)
      ENDDO
      DO j=0,ny
         WRITE(w_unit2,*)'RRR'
         WRITE(w_unit2,101) rate(0,0,j,1:4),vvel(i,k,j),rb(i,k)
         WRITE(w_unit2,*)'BMM'
         do n = 1,4
            WRITE(w_unit2,101) Bm(n,1:4,0,0,j)
         enddo
      ENDDO
   ENDIF

   CLOSE(w_unit)
   CLOSE(w_unit2)
   conv_error = abs(conv_error - rb(0,0))
   write(*,*)'CONVERGENCE ERROR',conv_error
   if(ncyc > 2*frame .AND. conv_error < 1.d-9)then
      write(*,*) '1-D sol is converged, &
           &stopped in file visualization.f90'
      write(*,104)
      DO j=0,ny
         WRITE(*,'(1p50e15.7)') y(j),f(1,0,j,1:neqgas)
      ENDDO
      CALL MPI_FINALIZE(ierr)
      stop
   endif
   conv_error = rb(0,0)
 
 101  FORMAT(50e14.4)
 102  FORMAT(1i4,1p50e20.11)
 103  FORMAT('Node',8x,'y(cm)',15x,'Temp(K)',14x,'HTPB',17x,'AP',16x,'MONO',13x,'BIND',14x,'PRI',15x,'vvel',14x,'y_vel')
 104  FORMAT(6x,'y(cm)',10x,'Temp(K)',9x,'HTPB',11x,'AP',11x,'MONO',10x,'BIND',10x,'PRI')

!---------------------------------------------------------------------
  RETURN 
END SUBROUTINE PRINT_1D
!*************************************************************
