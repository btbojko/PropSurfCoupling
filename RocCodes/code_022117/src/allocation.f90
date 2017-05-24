!***********************************************************
SUBROUTINE ALLOCATE_PACK
  
  USE GLOBAL_DATA
  
  IMPLICIT NONE
  INTEGER :: error

  ALLOCATE(xcenter(ncircles),STAT=error)
  CALL CHECK_ERROR(error,'PACK')
  ALLOCATE(ycenter(ncircles),STAT=error)
  CALL CHECK_ERROR(error,'PACK')
  ALLOCATE(zcenter(ncircles),STAT=error)
  CALL CHECK_ERROR(error,'PACK')
  ALLOCATE(rad(ncircles),STAT=error)
  CALL CHECK_ERROR(error,'PACK')
  allocate(OxidizerType(ncircles),STAT=error)
  CALL CHECK_ERROR(error,'PACK')


  
!-----------------------------------------------------------
  RETURN
END SUBROUTINE ALLOCATE_PACK
!***********************************************************
!***********************************************************
SUBROUTINE ALLOCATE_MATRIX
  
  USE GLOBAL_DATA
  
  IMPLICIT NONE
  INTEGER :: error,n

  n = matrix%nFourier

  ALLOCATE(matrix%M(n),STAT=error)
  CALL CHECK_ERROR(error,'MATRIX')
  ALLOCATE(matrix%F(n),STAT=error)
  CALL CHECK_ERROR(error,'MATRIX')
  ALLOCATE(matrix%A(n),STAT=error)
  CALL CHECK_ERROR(error,'MATRIX')

  
!-----------------------------------------------------------
  RETURN
END SUBROUTINE ALLOCATE_MATRIX
!***********************************************************

!***********************************************************
SUBROUTINE ALLOCATE_MAIN
  
  USE GLOBAL_DATA
  USE BC_DATA

  IMPLICIT NONE
  
  INTEGER :: error,i,j
  error = 0
!-----------------------------------------------------------

  ALLOCATE(x(-1:xr(1)+1),STAT=error)
  CALL CHECK_ERROR(error,'GRID')

  ALLOCATE(y(-1:xr(3)+1),STAT=error)
  ALLOCATE(ya(-1:xr(3)+1),STAT=error)
  CALL CHECK_ERROR(error,'GRID')

  ALLOCATE(z(-1:xr(2)+1),STAT=error)  
  CALL CHECK_ERROR(error,'GRID')

  allocate(x_mg(lbound(x,1):ubound(x,1),num_mg_meshes))
  allocate(z_mg(lbound(x,1):ubound(x,1),num_mg_meshes))

  ALLOCATE(detadya(0:xr(3)+1),STAT=error)
  ALLOCATE(detady(0:xr(3)+1),STAT=error)
  CALL CHECK_ERROR(error,'GRID')

  ALLOCATE(deta2dy2(0:xr(3)+1),STAT=error)
  ALLOCATE(deta2dy2a(0:xr(3)+1),STAT=error)
  CALL CHECK_ERROR(error,'GRID')

  ALLOCATE(slambda(0:maxsize*(maxrange+4)**2),STAT=error)
  CALL CHECK_ERROR(error,'PARALLEL')

  ALLOCATE(rlambda(0:maxsize*(maxrange+4)**2),STAT=error)
  CALL CHECK_ERROR(error,'PARALLEL')

  ALLOCATE(sbuf1(0:maxrange+4),STAT=error)
  CALL CHECK_ERROR(error,'PARALLEL')

  ALLOCATE(sbuf2(0:maxrange+4),STAT=error)
  CALL CHECK_ERROR(error,'PARALLEL')

  ALLOCATE(rbuf1(0:maxrange+4),STAT=error)
  CALL CHECK_ERROR(error,'PARALLEL')

  ALLOCATE(rbuf2(0:maxrange+4),STAT=error)
  CALL CHECK_ERROR(error,'PARALLEL')

  ALLOCATE(myright_ff_sbuffer(0:(maxrange+2)**2),STAT=error)
  CALL CHECK_ERROR(error,'PARALLEL')

  ALLOCATE(myleft_ff_sbuffer(0:(maxrange+2)**2),STAT=error)
  CALL CHECK_ERROR(error,'PARALLEL')

  ALLOCATE(myright_ff_rbuffer(0:(maxrange+2)**2),STAT=error)
  CALL CHECK_ERROR(error,'PARALLEL')

  ALLOCATE(myleft_ff_rbuffer(0:(maxrange+2)**2),STAT=error)
  CALL CHECK_ERROR(error,'PARALLEL')

  ALLOCATE(qheats(-1:xr(1)+1,-1:xr(2)+1,-1:xr(3)+1),STAT=error)
  CALL CHECK_ERROR(error,'PHYSICS')

  ALLOCATE(rhos(-1:xr(1)+1,-1:xr(2)+1,-1:xr(3)+1),STAT=error)
  ALLOCATE(tmdpnt(-1:xr(1)+1,-1:xr(2)+1,-1:xr(3)+1),STAT=error)
  CALL CHECK_ERROR(error,'PHYSICS')

  ALLOCATE(lambdas(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  CALL CHECK_ERROR(error,'PHYSICS')

  ALLOCATE(lambdag(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  ALLOCATE(rhog(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  CALL CHECK_ERROR(error,'PHYSICS_Gas')

  ALLOCATE(dlgdx(0:xr(1),0:xr(2),0:xr(3)),STAT=error)
  CALL CHECK_ERROR(error,'PHYSICS_Gas')

  ALLOCATE(dlgdy(0:xr(1),0:xr(2),0:xr(3)),STAT=error)
  CALL CHECK_ERROR(error,'PHYSICS_Gas')

  ALLOCATE(dlgdz(0:xr(1),0:xr(2),0:xr(3)),STAT=error)
  CALL CHECK_ERROR(error,'PHYSICS_Gas')


  ALLOCATE(lbdavg(-1:xr(1)+1,-1:xr(2)+1,-1:xr(3)+1,ndim),STAT=error)
  CALL CHECK_ERROR(error,'PHYSICS_Gas')
  

  ALLOCATE(da_s(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  CALL CHECK_ERROR(error,'PHYSICS')

  ALLOCATE(theta_s(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  CALL CHECK_ERROR(error,'PHYSICS')

  ALLOCATE(xsurf(-2:xr(1)+2,-2:xr(2)+2,maxsize),STAT=error)
  ALLOCATE(psurf(-2:xr(1)+2,-2:xr(2)+2,maxsize),STAT=error)
  CALL CHECK_ERROR(error,'Surface PHYSICS')


  ALLOCATE(f(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,maxsize),STAT=error)
  CALL CHECK_ERROR(error,'solution')
  ALLOCATE(q(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,ndim),STAT=error)
  ALLOCATE(newqBnd(-2:xr(1)+2,-2:xr(2)+2,ndim),STAT=error)
  ALLOCATE(srfQBnd(-2:xr(1)+2,-2:xr(2)+2,ndim),STAT=error)
  CALL CHECK_ERROR(error,'solution')
  ALLOCATE(oldsoln(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,maxsize+ndim),STAT=error)
  ALLOCATE(dff(lbound(f,1):ubound(f,1),lbound(f,2):ubound(f,2),&
       &lbound(f,3):ubound(f,3),lbound(f,4):ubound(f,4)),STAT=error)
  CALL CHECK_ERROR(error,'solution')

  ALLOCATE(ay(-2:xr(3)+2))
  ALLOCATE(by(-2:xr(3)+2))
  ALLOCATE(cy(-2:xr(3)+2))
  ALLOCATE(zy(-2:xr(3)+2))
  ALLOCATE(fy(-2:xr(3)+2)) 
  ALLOCATE(aq(maxsize,-2:xr(3)+2))
  ALLOCATE(cq(maxsize,-2:xr(3)+2))
  allocate(deltaKron(maxsize,maxsize))
  deltaKron = 0d0
  do i = 1,ubound(deltaKron,1)
     deltaKron(i,i)=1d0
  end do

  ALLOCATE(Bb(maxsize-1, maxsize-1, -2:xr(3)+2))
  ALLOCATE(fyv(maxsize-1, -2:xr(3)+2))
  ALLOCATE(Binv(maxsize-1, maxsize-1))
  ALLOCATE(Vdiag0(maxsize-1))
  ALLOCATE(Vright1(maxsize-1))

  ALLOCATE(iconv(0:xr(1),0:xr(2),0:xr(3),NDIM))
  ALLOCATE(dconv(0:xr(1),0:xr(2),0:xr(3),NDIM))
  

  ALLOCATE(phit(0:xr(1),0:xr(2)),STAT=error)
  CALL CHECK_ERROR(error,'SURFACE')

  ALLOCATE(   phi(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(oldphi(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error,'SURFACE')

  ALLOCATE(dphidx(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error,'SURFACE')

  ALLOCATE(dphi2dx2(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error,'SURFACE')

  ALLOCATE(dphidz(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error,'SURFACE')

  ALLOCATE(dphi2dz2(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error,'SURFACE')


  ALLOCATE(vec(-1:xr(1)+1,-1:xr(2)+1,0:2),STAT=error)
  ALLOCATE(vec0(-1:xr(1)+1,-1:xr(2)+1),STAT=error)
  vec = 1.0d0
  vec0 = 1.0d0
  ALLOCATE(dphidxa(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(dphidza(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error,'SURFACE')

  allocate(taGsmooth(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error,'SURFACE')

  allocate(MQterm(-2:xr(1)+2,-2:xr(2)+2,0:3),STAT=error)
  allocate(MQvelocity(-2:xr(1)+2,-2:xr(2)+2,0:ny,0:3),STAT=error)
  allocate(MQeta(-1:ny,2),MQchi(0:ny,2),MQchipr(0:ny,2),MQchisc(0:ny,2),STAT=error)
  ALLOCATE( MQphi(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error,'MQvariables')

  ALLOCATE(psi(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  CALL CHECK_ERROR(error,'SURFACE')

  ALLOCATE(rb(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(d_rb(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(ft0(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(d_ft0(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(d_vterm(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error,'SURFACE')

  ALLOCATE(nquad((maxcells+1)*(maxcells+1)),STAT=error)
  CALL CHECK_ERROR(error,'QUADRANTS')

  ALLOCATE(quad((maxcells+1)*(maxcells+1),maxpack),STAT=error)
  CALL CHECK_ERROR(error,'QUADRANTS')

  ALLOCATE(nxquad((maxcells+1)*(maxcells+1),2),STAT=error)
  CALL CHECK_ERROR(error,'QUADRANTS')

  ALLOCATE(nzquad((maxcells+1)*(maxcells+1),2),STAT=error)
  CALL CHECK_ERROR(error,'QUADRANTS')

  ALLOCATE(xl((maxcells+1)*(maxcells+1)),STAT=error)
  CALL CHECK_ERROR(error,'QUADRANTS')

  ALLOCATE(xm((maxcells+1)*(maxcells+1)),STAT=error)
  CALL CHECK_ERROR(error,'QUADRANTS')

  ALLOCATE(zl((maxcells+1)*(maxcells+1)),STAT=error)
  CALL CHECK_ERROR(error,'QUADRANTS')

  ALLOCATE(zm((maxcells+1)*(maxcells+1)),STAT=error)
  CALL CHECK_ERROR(error,'QUADRANTS')

  ALLOCATE(phi_movie(0:nz+1),STAT=error)
  CALL CHECK_ERROR(error,'MAIN')

  ALLOCATE(z_movie(0:nz+1),STAT=error)
  CALL CHECK_ERROR(error,'MAIN')
  
  ALLOCATE(p(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  ALLOCATE(pold(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  CALL CHECK_ERROR(error,'Pressure')

  
  ALLOCATE(rate(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,maxsize),STAT=error)
  ALLOCATE(rate_out(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,maxsize),STAT=error)
!!>  ALLOCATE(inirate(lbound(rate,1):ubound(rate,1), lbound(rate,2):ubound(rate,2), lbound(rate,3):ubound(rate,3), lbound(rate,4):ubound(rate,4)),STAT=error)
!!>  ALLOCATE(inif(lbound(f,1):ubound(f,1),lbound(f,2):ubound(f,2),lbound(f,3):ubound(f,3),lbound(f,4):ubound(f,4)),STAT=error)
  CALL CHECK_ERROR(error,'RATE')

  ALLOCATE(radheat(lbound(rate,1):ubound(rate,1),lbound(rate,2):ubound(rate,2),&
       &           lbound(rate,3):ubound(rate,3)),STAT=error)
  ALLOCATE(BCradheat(lbound(radheat,1):ubound(radheat,1),lbound(radheat,2):ubound(radheat,2)),STAT=error)
  allocate(Gfield(lbound(p,1):ubound(p,1),lbound(p,2):ubound(p,2),lbound(p,3):ubound(p,3)),STAT=error)
  allocate(Vfrac(lbound(p,1):ubound(p,1),lbound(p,2):ubound(p,2),lbound(p,3):ubound(p,3),2),STAT=error)
  allocate(DiamP(lbound(p,1):ubound(p,1),lbound(p,2):ubound(p,2),lbound(p,3):ubound(p,3),2),STAT=error)
  allocate(VOD(lbound(p,1):ubound(p,1),lbound(p,2):ubound(p,2),lbound(p,3):ubound(p,3),2),STAT=error)
  radheat = 0d0
  BCradheat = 0d0
  Gfield = 0d0
  CALL CHECK_ERROR(error,'RDIATION SOURCES')

  ALLOCATE(solidProbes(Nprobes),OldSolidProbes(Nprobes))


  ALLOCATE(capa(maxsize, -2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(capb(maxsize, -2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(capc(maxsize, -2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error,'BC')

  ALLOCATE(capd(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(capdd(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(cape(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(capee(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error,'BC')

  ALLOCATE(capr(maxsize-1, -2:xr(1)+2, -2:xr(2)+2),STAT=error)
  ALLOCATE(capT(maxsize-1, -2:xr(1)+2, -2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error,'BC')

  ALLOCATE(ffnew(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error,'BC')

  ALLOCATE(vvel(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  ALLOCATE(d_vvel(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  ALLOCATE(uvel(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  ALLOCATE(vcnv(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  ALLOCATE(wvel(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  CALL CHECK_ERROR(error,'Convection')

  !>>> vic: shear: from juPack
!! aI
  ALLOCATE (shear(-2:xr(3)+2), STAT=error)
  CALL CHECK_ERROR(error,'Shear')
!! aI
  !<<< vic

  ALLOCATE(dtx(-1:xr(1)+1,-1:xr(2)+1,-1:xr(3)+1),STAT=error)
  ALLOCATE(dtorhocp(lbound(dtx,1):ubound(dtx,1),lbound(dtx,2):ubound(dtx,2),&
       &lbound(dtx,3):ubound(dtx,3)),STAT=error)
  CALL CHECK_ERROR(error,'Timesteps')
!------------------------------------------------------------------------
  RETURN
END SUBROUTINE ALLOCATE_MAIN
!***********************************************************************
!***********************************************************
SUBROUTINE ALLOCATE_SOLVER
  
  USE GLOBAL_DATA
  USE IMPLICIT
  USE BC_DATA

  IMPLICIT NONE
  
  INTEGER :: error,sumerr
  error = 0
!-----------------------------------------------------------

  ALLOCATE(dfdt(lbound(dff,1):ubound(dff,1),lbound(dff,2):ubound(dff,2),&
       &lbound(dff,3):ubound(dff,3),lbound(dff,4):ubound(dff,4)),STAT=error)
  ALLOCATE(dqdt(lbound(dfdt,1):ubound(dfdt,1),lbound(dfdt,2):ubound(dfdt,2),&
       &lbound(dfdt,3):ubound(dfdt,3),lbound(dfdt,4):ubound(dfdt,4)),STAT=error)
  ALLOCATE(oldrate(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  ALLOCATE(rhsold(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,maxsize),STAT=error)
  CALL CHECK_ERROR(error,'dfdtdqdtrate')
  dfdt = 0d0
  dqdt = 0d0
  oldrate = 0d0


!******************************************************************************
! MEMORY ALLOCATION FOR JACOBIANS FOR THE IMPLICIT SOLVER

  ALLOCATE(Bm(neqgas,neqgas,0:xr(1),0:xr(2),0:xr(3)),STAT=error)
  CALL CHECK_ERROR(error,'CHemJacobians')
  sumerr = 0
  if(iterativemethod == 1) then
     gmrssize = min(iterLINEARsys,maxsizeGMRES)
     if(associated(coarsest_mesh%gm)) then
        gmrssize2 = max(gmrssize,ubound(finest_mesh%gm,1),ubound(coarsest_mesh%gm,1))
     elseif(associated(finest_mesh%gm)) then
        gmrssize2 = max(gmrssize,ubound(finest_mesh%gm,1))
     else
        gmrssize2 = gmrssize
     end if
     if(myid == 0) Print*,'allocating',gmrssize,gmrssize2
     allocate(vgm(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,maxsize,gmrssize),STAT=error);sumerr = sumerr+error
     allocate(dgm(gmrssize2+1),hgm(gmrssize2+1,gmrssize2+1),c1gm(gmrssize2+1))
     allocate(c2gm(gmrssize2+1),ygm(gmrssize2+1),STAT=error); sumerr = sumerr+error
     ALLOCATE(wrk4(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,maxsize),STAT=error); sumerr = sumerr+error
     CALL CHECK_ERROR(sumerr,'GMres')
     vgm = 0d0
     wrk4 = 0d0
  endif
!option save_LU_VARS increases memory requirement
  ALLOCATE(L_save(ubound(Bm,1),ubound(Bm,2),0:xr(1),0:xr(2),0:xr(3)),STAT=error)
  CALL CHECK_ERROR(error,'JACOBIANS')
  ALLOCATE(Binv_save(ubound(Bm,1),ubound(Bm,2),0:xr(1),0:xr(2),0:xr(3)),STAT=error)
  CALL CHECK_ERROR(error,'JACOBIANS')

  ALLOCATE(coeff(maxsize,10,0:xr(1),0:xr(2),0:xr(3)),STAT=error)
  CALL CHECK_ERROR(error,'JACOBIANS')
  coeff = zero
  coeff(:,1,:,:,:) = one
  ALLOCATE(csld(10,0:xr(1),0:xr(2),0:xr(3)),STAT=error)
  CALL CHECK_ERROR(error,'JACOBIANS')


!******************************************************************************
!MEMORY ALLOCATION FOR FLUIDS VARIABLES



  ALLOCATE(divt(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  ALLOCATE(Postdivt(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  CALL CHECK_ERROR(error,'PROJECTION')


  allocate(Uimposed(Nimposed),Zimposed(Nimposed))
  Uimposed = 0d0
  Zimposed = 0d0

  ALLOCATE(rhogasCont(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,4),STAT=error)
  ALLOCATE(dtCont(2),STAT=error)
  CALL CHECK_ERROR(error,'rhogasCont')



  


!------------------------------------------------------------------------
  RETURN
END SUBROUTINE ALLOCATE_SOLVER
!***********************************************************************


SUBROUTINE CHECK_ERROR(err,msg)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

  INTEGER :: err,allerr
  character(len =*) :: msg
  
  if(err.gt.0) then
     write(*,*) myid, ':: ERROR!  FATAL! ACK! Unable to allocate something.'
     write(*,*) TRIM(msg)
     write(*,*) 'Grid Dimensions',nx,nz,ny
  endif
  
  CALL MPI_ALLREDUCE(err,allerr,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
  
  if(allerr.gt.0) stop
  
  RETURN

END SUBROUTINE CHECK_ERROR
