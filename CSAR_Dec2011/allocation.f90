!***********************************************************
SUBROUTINE ALLOCATE_PACK
  
  USE GLOBAL_DATA
  
  IMPLICIT NONE
  INTEGER :: error

  ALLOCATE(xcenter(ncircles),STAT=error)
  CALL CHECK_ERROR(error)
  ALLOCATE(ycenter(ncircles),STAT=error)
  CALL CHECK_ERROR(error)
  ALLOCATE(itype(ncircles),STAT=error)
  CALL CHECK_ERROR(error)
  ALLOCATE(rad(ncircles),STAT=error)
  CALL CHECK_ERROR(error)

  
!-----------------------------------------------------------
  RETURN
END SUBROUTINE ALLOCATE_PACK
!***********************************************************

!***********************************************************
SUBROUTINE ALLOCATE_INITIAL
  
  USE GLOBAL_DATA
  
  IMPLICIT NONE
  INTEGER :: error

  ALLOCATE(thetay(num_rxn),STAT=error)
  CALL CHECK_ERROR(error)
  ALLOCATE(da(num_rxn),STAT=error)
  CALL CHECK_ERROR(error)
  ALLOCATE(np(num_rxn),STAT=error)
  CALL CHECK_ERROR(error)
  
!-----------------------------------------------------------
  RETURN
END SUBROUTINE ALLOCATE_INITIAL
!***********************************************************

!***********************************************************
SUBROUTINE ALLOCATE_MAIN
  
  USE GLOBAL_DATA
  USE BC_DATA

  IMPLICIT NONE
  
  INTEGER :: error
  error = 0
!-----------------------------------------------------------

  ALLOCATE(x(-1:xr(1)+1),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(y(-1:xr(3)+1),STAT=error)
  ALLOCATE(ya(-1:xr(3)+1),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(z(-1:xr(2)+1),STAT=error)  
  CALL CHECK_ERROR(error)

  ALLOCATE(detadya(0:xr(3)+1),STAT=error)
  ALLOCATE(detady(0:xr(3)+1),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(deta2dy2(0:xr(3)+1),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(slambda(0:maxsize*(maxrange+4)**2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(rlambda(0:maxsize*(maxrange+4)**2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(sbuf1(0:maxrange+4),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(sbuf2(0:maxrange+4),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(rbuf1(0:maxrange+4),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(rbuf2(0:maxrange+4),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(myright_ff_sbuffer(0:(maxrange+2)**2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(myleft_ff_sbuffer(0:(maxrange+2)**2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(myright_ff_rbuffer(0:(maxrange+2)**2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(myleft_ff_rbuffer(0:(maxrange+2)**2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(qheats(-1:xr(1)+1,-1:xr(2)+1,-1:xr(3)+1),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(rhos(-1:xr(1)+1,-1:xr(2)+1,-1:xr(3)+1),STAT=error)
  ALLOCATE(rhos_old(-1:xr(1)+1,-1:xr(2)+1,-1:xr(3)+1),STAT=error)
  ALLOCATE(rhog(-1:xr(1)+1,-1:xr(2)+1,-1:xr(3)+1),STAT=error)
  ALLOCATE(tmdpnt(-1:xr(1)+1,-1:xr(2)+1,-1:xr(3)+1),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(lambdas(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(da_s(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(theta_s(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(xsurf(-2:xr(1)+2,-2:xr(2)+2,maxsize),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(f(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,maxsize),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(q(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,ndim),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(oldsoln(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2,maxsize),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(ay(-2:xr(3)+2))
  ALLOCATE(by(-2:xr(3)+2))
  ALLOCATE(cy(-2:xr(3)+2))
  ALLOCATE(zy(-2:xr(3)+2))
  ALLOCATE(fy(-2:xr(3)+2)) 
  ALLOCATE(aq(-2:xr(3)+2))
  ALLOCATE(cq(-2:xr(3)+2))

  ALLOCATE(Bb(maxsize-1, maxsize-1, -2:xr(3)+2))
  ALLOCATE(fyv(maxsize-1, -2:xr(3)+2))
  ALLOCATE(Binv(maxsize-1, maxsize-1))
  ALLOCATE(Vdiag0(maxsize-1))
  ALLOCATE(Vright1(maxsize-1))
  

!
!  RATE dimension have been augmented to use it as work vector in MG
!
  ALLOCATE(rate(-2:xr(1)+2,-2:xr(2)+2,0:xr(3),maxsize),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(lambdag(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(dlgdx(0:xr(1),0:xr(2),0:xr(3)),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(dlgdy(0:xr(1),0:xr(2),0:xr(3)),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(dlgdz(0:xr(1),0:xr(2),0:xr(3)),STAT=error)
  CALL CHECK_ERROR(error)
    
  ALLOCATE(phit(0:xr(1),0:xr(2)),STAT=error)
  ALLOCATE(  ff (-3:xr(1)+3,-3:xr(3)+3),STAT=error)
  ALLOCATE(  tab(-3:xr(1)+3,-3:xr(3)+3),STAT=error)
  ALLOCATE(  dxp(-3:xr(1)+3,-3:xr(3)+3),STAT=error)
  ALLOCATE(  dxm(-3:xr(1)+3,-3:xr(3)+3),STAT=error)
  ALLOCATE(  dyp(-3:xr(1)+3,-3:xr(3)+3),STAT=error)
  ALLOCATE(  dym(-3:xr(1)+3,-3:xr(3)+3),STAT=error)
  ALLOCATE(L_phi(-3:xr(1)+3,-3:xr(3)+3),STAT=error)
  ALLOCATE(xy   (-3:xr(1)+3,-3:xr(3)+3),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(   phi(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(oldphi(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(dphidx(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(dphi2dx2(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(dphidz(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(dphi2dz2(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error)


  ALLOCATE(psi(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(rb(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(d_rb(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(ft0(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(d_ft0(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(d_vterm(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(dfdt(0:xr(1),0:xr(2),0:xr(3),maxsize),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(capa(2, -2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(capb(2, -2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(capc(2, -2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(capd(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(capdd(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(cape(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(capee(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(capr(maxsize-1, -2:xr(1)+2, -2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(ffnew(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(nquad((maxcells+1)*(maxcells+1)),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(quad((maxcells+1)*(maxcells+1),maxpack),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(nxquad((maxcells+1)*(maxcells+1),2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(nzquad((maxcells+1)*(maxcells+1),2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(xl((maxcells+1)*(maxcells+1)),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(xm((maxcells+1)*(maxcells+1)),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(zl((maxcells+1)*(maxcells+1)),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(zm((maxcells+1)*(maxcells+1)),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(phi_movie(0:nz+1),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(z_movie(0:nz+1),STAT=error)
  CALL CHECK_ERROR(error)

!******************************************************************************
! MEMORY ALLOCATION FOR JACOBIANS FOR THE IMPLICIT SOLVER
  ALLOCATE(oldrate(0:xr(1),0:xr(2),0:xr(3)),STAT=error)
  ALLOCATE(dtx(-1:xr(1)+1,-1:xr(2)+1,-1:xr(3)+1),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(Bm(neqgas,neqgas,0:xr(1),0:xr(2),0:xr(3)),STAT=error)
!option save_LU_VARS increases memory requirement
  ALLOCATE(L_save(neqgas,neqgas,0:xr(1),0:xr(2),0:xr(3)),STAT=error)
  ALLOCATE(Binv_save(neqgas,neqgas,0:xr(1),0:xr(2),0:xr(3)),STAT=error)
!
  CALL CHECK_ERROR(error)

  ALLOCATE(coeff(9,0:xr(1),0:xr(2),0:xr(3)),STAT=error)
  ALLOCATE(csld(10,0:xr(1),0:xr(2),0:xr(3)),STAT=error)
  ALLOCATE(wrk_vec(2,0:xr(1),0:xr(2),0:xr(3)))
  CALL CHECK_ERROR(error)

!******************************************************************************
!MEMORY ALLOCATION FOR FLUIDS VARIABLES


  ALLOCATE(vec(-1:xr(1)+1,-1:xr(2)+1),STAT=error)
  ALLOCATE(dphidxa(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  ALLOCATE(dphidza(-2:xr(1)+2,-2:xr(2)+2),STAT=error)
  CALL CHECK_ERROR(error)

  ALLOCATE(vvel(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  ALLOCATE(vcnv(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  ALLOCATE(uvel(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  ALLOCATE(wvel(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  ALLOCATE(qsurf(-2:xr(1)+2,-2:xr(2)+2,3),STAT=error)
  ALLOCATE(p(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  ALLOCATE(pold(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  ALLOCATE(divt(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  ALLOCATE(v_south(-1:xr(3)+1),STAT=error)
  ALLOCATE(dconv(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  ALLOCATE(iconv(-2:xr(1)+2,-2:xr(2)+2,-2:xr(3)+2),STAT=error)
  CALL CHECK_ERROR(error)

!$4 dimension vector (:,:,:ndim)
  ALLOCATE(lbdavg(-1:xr(1)+1,-1:xr(2)+1,-1:xr(3)+1,ndim),STAT=error)
  ALLOCATE(fdiff(-1:xr(1)+1,-1:xr(2)+1,-1:xr(3)+1,ndim),STAT=error)
  ALLOCATE(dqdt(0:xr(1),0:xr(2),0:xr(3),maxsize),STAT=error)
  CALL CHECK_ERROR(error)

!  WRITE(*,*) mywid,':: Done Allocating Field Data.'

  RETURN

END SUBROUTINE ALLOCATE_MAIN

SUBROUTINE CHECK_ERROR(err)

  USE GLOBAL_DATA
  USE MYMPI

  IMPLICIT NONE

  INTEGER :: err,allerr
  
  if(err.gt.0) then
     write(*,*) mywid, ':: ERROR!  FATAL! ACK! Unable to allocate something.'
  endif
  
  CALL MPI_ALLREDUCE(err,allerr,1,MPI_INTEGER,MPI_MAX,MPI_COMM_WORLD,ierr)
  
  if(allerr.gt.0) stop
  
  RETURN

END SUBROUTINE CHECK_ERROR
