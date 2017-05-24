MODULE LUSOLVER

  USE GLOBAL_DATA, ONLY : Ainv => Binv, comm3d, nproc, ierr, myid
  USE LU_VARS
  USE MYMPI

  IMPLICIT NONE
  REAL *8, DIMENSION (7,7)  :: u,l
  REAL *8, DIMENSION (7)    :: b,y
!----------------------------------------------------------

CONTAINS

  SUBROUTINE LUINV(A,n)

    integer :: i,j,n
    real*8  :: A(n,n)
    real*8  :: den


    l = 1.0
    u = 0.0d0
    u(1,1) = a(1,1)
    u(1,2:n) = a(1,2:n)
    l(2:n,1) = a(2:n,1)/a(1,1)

    do i =2,n-1

       u(i,i) = a(i,i) - sum(l(i,1:i-1)*u(1:i-1,i))
       den = 1.0 / u(i,i)

       do j = i+1,n
          u(i,j) = a(i,j) - sum(l(i,1:i-1)*u(1:i-1,j))
          l(j,i) = (a(j,i) - sum(l(j,1:i-1)*u(1:i-1,i))) * den
       enddo
    enddo

    i=n
    u(i,i) = a(i,i) - sum(l(i,1:i-1)*u(1:i-1,i))


    do j = 1,n

       b(1:n) = 0.0d0
       b(j)   = 1.0d0

       y(1) = b(1)

       do i = 2,n
          y(i) = b(i) - sum(l(i,1:i-1)*y(1:i-1))
       enddo

       b(n) = y(n)/u(n,n)
       do i = n-1,1,-1
          b(i) = (y(i) - sum(u(i,i+1:n)*b(i+1:n)))/u(i,i)
       enddo

       Ainv(1:n,j) = b(1:n)

    enddo


    RETURN 
  END  SUBROUTINE LUINV

!********************************************************
SUBROUTINE LUSOLVE(A,f,n)

  integer :: i,j,n
  real*8  :: A(n,n),f(n)
  real*8  :: d


  l = 1.0
  u = 0.0d0
  u(1,1) = a(1,1)
  u(1,2:n) = a(1,2:n)
  l(2:n,1) = a(2:n,1)/a(1,1)

  do i =2,n-1

     u(i,i) = a(i,i) - sum(l(i,1:i-1)*u(1:i-1,i))
     d = 1.0 / u(i,i)

     do j = i+1,n
        u(i,j) = a(i,j) - sum(l(i,1:i-1)*u(1:i-1,j))
        l(j,i) = (a(j,i) - sum(l(j,1:i-1)*u(1:i-1,i))) * d
     enddo
  enddo

  i=n
  u(i,i) = a(i,i) - sum(l(i,1:i-1)*u(1:i-1,i))


  y(1) = f(1)

  do i = 2,n
     y(i) = f(i) - sum(l(i,1:i-1)*y(1:i-1))
  enddo

  f(n) = y(n)/u(n,n)
  do i = n-1,1,-1
     f(i) = (y(i) - sum(u(i,i+1:n)*f(i+1:n)))/u(i,i)
  enddo


  RETURN 
END SUBROUTINE LUSOLVE
!********************************************************

!********************************************************
SUBROUTINE LU_ALLOC(a,b,c,d,xrange,yrange)

  implicit none

  integer :: a,b,c,d,xrange,yrange
  integer :: i,mysize,eqnrange
!--------------------------------------------------------
  LU_nx = a
  LU_ny = b
  nmaxj = 3*b
  nmaxi = 3*b
  ndiag = c
  sizea = LU_nx*LU_ny
  LU_frequency = d
  done_LU_decomposition = .FALSE.
  
  ALLOCATE(as(sizea,nmaxj),ab(nmaxi,sizea),ls(sizea,nmaxj),&
       ub(nmaxi,sizea),us(sizea,nmaxj))
  ALLOCATE(LU_y(sizea),LU_f(sizea))
  ALLOCATE(cm(LU_nx,LU_ny,ndiag))

  ALLOCATE(ja(sizea,nmaxj),ia(sizea,nmaxi),nja(sizea),nia(sizea),&
       mainD_j(sizea),mainD_i(sizea),njl(sizea),jl(sizea,nmaxj),&
       niu(sizea),iu(sizea,nmaxi),nju(sizea),ju(sizea,nmaxj))

  ALLOCATE(posdiag(ndiag,2))

  ALLOCATE(LU_size(nproc))
  ALLOCATE(LU_disp(nproc))
  ALLOCATE(allcounts(nproc))
  ALLOCATE(displacements(nproc))

  if(ndiag == 5) then
     LU_diag_gather = ndiag - 2
  elseif(ndiag == 9) then
     LU_diag_gather = ndiag 
  else
     print*,'ERROR  wrong number of diags'
     stop 'lusolver.f90'
  endif

  eqnrange = LU_diag_gather

  mysize = xrange*yrange

  call MPI_ALLGATHER(mysize,1,MPI_INTEGER,allcounts,1,MPI_INTEGER, &
       &     comm3d,ierr)

  LU_size(1:nproc) = allcounts(1:nproc)
  allcounts(1:nproc) = allcounts(1:nproc)*eqnrange

  displacements(1) = 0
  LU_disp(1) = 0
  do i = 2, nproc
     displacements(i) = allcounts(i-1) + displacements(i-1)
     LU_disp(i) = LU_size(i-1) +  LU_disp(i-1)
  enddo
  
  lastmem = displacements(nproc) +  allcounts(nproc)
  LU_last =  LU_size(nproc) +  LU_disp(nproc)
  allocate(buff(0:mysize*eqnrange))
  allocate(buff2(0:lastmem+1))
  allocate(buff3(0:LU_last+1))


  residual_direct = 0.0d0
  print*, 'DONE LU_ALLOC',LU_disp,displacements
  
!---------------------------------------------------------
  RETURN 
END SUBROUTINE LU_ALLOC
!*********************************************************


!********************************************************
SUBROUTINE LU_DECOMPOSE

  implicit real*8 (a-h,o-z)
  implicit integer(i-n)
  logical :: inbounds,flag,point_found
!--------------------------------------------------------


  prec = 1.d-7
  posdiag = 0
  posdiag(2,2) = -1
  posdiag(3,2) = +1
  posdiag(4,1) = -1;posdiag(4,2) = -1
  posdiag(5,1) = -1;posdiag(5,2) = 1
  posdiag(6,1) = 1 ;posdiag(6,2) = -1
  posdiag(7,1) = 1 ;posdiag(7,2) = 1
  posdiag(8,1) = -1
  posdiag(9,1) = +1

!
! LU_prd = 0 nullify out_of_bounds  
! LU_prd = 1 set-up periodic BC
! LU_prd = 2 reflect
!
  i=0
  LU_prd(1) = 2;LU_prd(2) = 0
  LU_np(1) = LU_nx;LU_np(2) = LU_ny


  if(LU_prd(1) == 1)then
     write(*,*)'ILLEGAL OPTION, NEED TO POST proocess it after the loop';
     stop 'ASSEMBLING matrix'
  endif

!
! MAIN MATRIX ASSIGNMENT LOOP
  nja(:)= 0
  as = 0.0d0
  do ic = 1,LU_np(1)
     do jc = 1,LU_np(2)
        i=i+1
        do n =1,ndiag
           inbounds = ( jc+posdiag(n,2) >= 1 .AND. jc+posdiag(n,2) <= LU_np(2) ) .OR. LU_prd(2) > 0
           inbounds = inbounds .AND. ( (ic+posdiag(n,1) >= 1 &
                &.AND. ic+posdiag(n,1) <= LU_np(1) )  .OR. LU_prd(1) > 0 )
           inbounds = inbounds .AND. abs(cm(ic,jc,n)) > prec
           if(inbounds) then
              if(LU_prd(1) == 2)then
                 ix = ic+posdiag(n,1)-1
!                 if(ix > LU_np(1)) ix = ix - (ix - LU_np(1))
                 if(ix == -1)ix = +1
                 if(ix == LU_np(1))ix = LU_np(1)-1
              else
                 ix = (ic+posdiag(n,1)-1)
              endif
              ii =  ix*LU_np(2)
              if(LU_prd(2) > 0)then
                 jdiag = mod(jc+posdiag(n,2)-1+LU_prd(2),LU_prd(2)) + 1 + ii
              else
                 jdiag = jc + posdiag(n,2) + ii  !\bc i cannot use mod(n,0) in FORTRAN
              endif

              point_found = .FALSE.
              do jold = 1,nja(i)
                 if(jdiag == ja(i,jold)) then
                    as(i,jold) = as(i,jold) + cm(ic,jc,n)
                    point_found = .TRUE.
                    EXIT
                 endif
              enddo
              if(.NOT. point_found) then
                 nja(i)=nja(i)+1
                 as(i,nja(i)) = cm(ic,jc,n)
                 ja(i,nja(i)) = jdiag
              endif
           endif
        enddo
!SORT the entries against the j,i.e nja(i)
        do n =2,nja(i)
           minj=ja(i,n-1)
           amin = as(i,n-1)
           do m = n,nja(i)
              if(ja(i,m)< minj) then
                 ja(i,n-1)=ja(i,m)
                 as(i,n-1)=as(i,m)
                 ja(i,m)=minj
                 as(i,m)=amin
                 minj=ja(i,n-1)
                 amin = as(i,n-1)
              endif
           enddo
        enddo
!        if(myid == 0) write(23,*)i,ic,jc,as(i,1:nja(i)),ja(i,1:nja(i))
     enddo
  enddo
!  close(23)
   
  print*,'START DECOMPOSITION',' EFFECTIVE NDIAG',maxval(nja(1:sizea))

  nia(:)= 0
  do i =1,sizea
     do ic = 1,nja(i)
        j = ja(i,ic)
        nia(j) = nia(j)+1
        ab(nia(j),j) = as(i,ic)
        ia(j,nia(j)) = i
        if(j == i)then
           mainD_j(i) = ic
           mainD_i(j) = nia(j)
        endif
     enddo
  enddo

  
  njl = 0
  d = 1.0d0/as(1,1)
  do ic = 2,nia(1)
     i = ia(1,ic)
     njl(i) = 1
     ls(i,1) = ab(ic,1)*d
     jl(i,1) = 1
  enddo

  niu = 0
  do ic = 1,nja(1)
     j = ja(1,ic)
     niu(j) = 1
     ub(1,j) = as(1,ic)
     iu(j,1) = 1
  enddo

!
!  MAIN DECOMPOSITION LOOP
!  
  do i = 2,sizea-1


     jj=i;ii=i
     call sparse_sum(niu(i),njl(i),iu(i,:),jl(i,:),ub(:,i),ls(i,:),sumsp,flag)


     niu(jj) = niu(jj)+1
     ub(niu(jj),jj) = as(i,mainD_j(i)) - sumsp
     iu(jj,niu(jj)) = i
     d=1.0d0/ub(niu(jj),jj)


     nstca = mainD_j(i) + 1
     do j = i+1,sizea

        if(niu(j) == 0 .AND. (nia(j) == 0 .OR. ia(j,1) > i) ) cycle

        call sparse_sum(niu(j), njl(i),iu(j,:),jl(i,:),ub(:,j),ls(i,:),sumsp,flag)

        if(ja(i,nstca) == j) then     
           niu(j) = niu(j)+1
           ub(niu(j),j) = as(i,nstca) - sumsp
           iu(j,niu(j)) = i
           nstca = nstca+1
        elseif(flag) then    
           niu(j) = niu(j)+1
           ub(niu(j),j) = - sumsp
           iu(j,niu(j)) = i
        endif
     enddo

     nstca = mainD_i(i) + 1
     do i2 = i+1,sizea

        if(njl(i2) == 0 .AND. (nja(i2) == 0 .OR. ja(i2,1) > jj ) ) cycle

        call sparse_sum(niu(jj), njl(i2),iu(jj,:),jl(i2,:),ub(:,jj),ls(i2,:),sumsp,flag)

        if(ia(jj,nstca) == i2) then     
           njl(i2) = njl(i2)+1
           ls(i2,njl(i2)) = (ab(nstca,jj) - sumsp)*d
           jl(i2,njl(i2)) = jj
           nstca = nstca+1
        elseif(flag) then        
           njl(i2) = njl(i2)+1
           ls(i2,njl(i2)) = (- sumsp)*d
           jl(i2,njl(i2)) = jj
        endif
     enddo

  enddo

  i=sizea;
  call sparse_sum(niu(i), njl(i),iu(i,:),jl(i,:),ub(:,i),ls(i,:),sumsp,flag)
  niu(i) = niu(i)+1
  ub(niu(i),i) = as(i,mainD_j(i)) - sumsp
  iu(i,niu(i)) = sizea

!
! FINDS the matrix us from ub
!
  nju(:)= 0
  do j =1,sizea
     do ic = 1,niu(j)
        i = iu(j,ic)
        nju(i) = nju(i)+1
        us(i,nju(i)) = ub(ic,j)
        ju(i,nju(i)) = j
     enddo
  enddo

  print*,'END DECOMPOSITION'
!---------------------------------------------------------
  RETURN 
END SUBROUTINE LU_DECOMPOSE
!********************************************************

!********************************************************
SUBROUTINE LU_SPARSE_SOLVE

  implicit none
  integer :: i
!--------------------------------------------------------
 
  LU_Y(1)= LU_F(1)
  do i = 2,sizea
     LU_y(i) = LU_F(i) - sum(ls(i,1:njl(i))*LU_y(jl(i,1:njl(i))))
!     if(myid ==0) print*,i,LU_f(i), ls(i,1:njl(i))
  enddo

  LU_f(sizea) = LU_y(sizea) / us( sizea,nju(sizea) )
  do i =sizea-1,1,-1
     LU_f(i) = (LU_y(i) - sum(us(i,2:nju(i))*LU_f(ju(i,2:nju(i))) ))/ us(i,1)
!     print*,i,LU_f(i)
  enddo
!---------------------------------------------------------
  RETURN 
END SUBROUTINE LU_SPARSE_SOLVE
!********************************************************

!********************************************************************
subroutine sparse_sum(nrow,ncol,i,j,vecj,veci,sums,flag)

  implicit none

  integer :: nrow,ncol,i(nrow),j(ncol),icol,irow
  real*8 :: veci(ncol),vecj(nrow),sums
  LOGICAL :: flag

  icol = 1;irow =1

  sums = 0.0d0


  do while(icol <= ncol .AND. irow <= nrow)
     if(i(irow) == j(icol)) then
        sums = sums + veci(icol)*vecj(irow)
        icol = icol+1; irow = irow + 1
     elseif(i(irow) < j(icol))then
        irow = irow+1
     else
        icol = icol+1
     endif
  enddo

  flag = abs(sums) > 1.d-14
  
!----------------------------------------------------
return
end subroutine sparse_sum
!***************************************************

!*****************************************************
subroutine LU_GATHER(xrange,yrange)

  integer :: xrange,yrange
  integer :: el_eqn,i_LU,j_LU,pcs,yoff
  integer :: i,j,n,elcount,mysize



  if(nproc == 1) RETURN

  mysize = xrange*yrange*LU_diag_gather

  elcount = 0
  do i = 0, xrange-1
     i_LU= i+1;
     do j= 0, yrange-1
        j_LU=j+1
        do n = 1, LU_diag_gather
           buff(elcount) = cm(i_LU,j_LU,n)
!           if(myid == 0) then
!              print*,i,j,n,'BUFF',buff(elcount)
!           endif
           elcount = elcount + 1
        enddo
     enddo
  enddo


  if(elcount /= mysize) stop 'ERROR IN GATHER_LU'


  call MPI_ALLGATHERV(buff(0),elcount,MPI_DOUBLE_PRECISION, &
       &     buff2(0),allcounts,displacements,MPI_DOUBLE_PRECISION, &
       &     comm3d,ierr)


  elcount = 0
  do el_eqn = 0, lastmem-1, LU_diag_gather
     j = mod(elcount,LU_ny) + 1
     i = int(elcount/LU_ny) + 1
     cm(i,j,1:LU_diag_gather) = buff2(el_eqn : el_eqn+LU_diag_gather-1)
!     if(myid == 0) print*,i,j,elcount,el_eqn,&
!          'LU_GATHER',cm(i,j,1:LU_diag_gather),buff2(el_eqn:el_eqn+LU_diag_gather-1)
     elcount = elcount+1
  enddo
  

!-----------------------------------------------
  return
end subroutine LU_gather
!*************************************************

!*****************************************************
subroutine LU_scatter(xrange,yrange)

!------------------------------------
  integer :: xrange,yrange
  integer :: mysize
!------------------------------------


  if(nproc == 1) then
     LU_Y = LU_F
  else
     mysize = xrange*yrange
     CALL MPI_SCATTERV(LU_F(1), LU_size, LU_disp, MPI_DOUBLE_PRECISION,&
          LU_Y(1), mysize, MPI_DOUBLE_PRECISION, 0,comm3d , IERR) 
  endif

!-----------------------------------------------
  return
end subroutine LU_scatter
!*************************************************

END MODULE LUSOLVER
