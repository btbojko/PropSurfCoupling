!**********************************************************************
SUBROUTINE global_parameter
USE GLOBAL_DATA 

IMPLICIT NONE
!---------------------------------------------------------------------
! Local variables
INTEGER ::  m
REAL*8  ::  qb_heat,binder,ap_fines,ap_coarse
!--------------------------------------------------------------------
!Premixed binder heat release is the mass average of the 2 ingredients
qb_heat = alp_W*qheat_ap+(1.0-alp_W)*qheat_binder

!Calculate reaction coefficients for final diffusion flame (amount of AP not in binder to reach stoich)
if(alpha == 1.0)then
    yf(4,4) = 0.0d0
elseif(ipack == 0)then
    yf(4,4) = -AP_frac*rho_ap/((1.0d0-AP_frac)*rho_premixed)
else
    binder = 1.0d0 - alpha_pack
    ap_fines = alpha*binder/(1.0d0 - alpha)
    ap_coarse = alpha_pack - ap_fines
    yf(4,4) = -ap_coarse/binder
endif
yf(5,4) = -1.0d0

!Flag if binder is treated as 2 ingredients, not fully implemented
if(premixed == 0) then         
    qb_heat = qheat_binder
    yf(4,4) = alp_W/(1.0d0 - alp_W)
endif
    
!Calculate heat release for each reaction
qg(1) = cp*(th_ap - tcold)* abs(yf(3,1)) - qheat_ap * abs(yf(3,1))                                          !AP monopropellant
qg(2) = cp*(th_bind - tcold)* abs(yf(2,2)) - qb_heat * abs(yf(2,2))                                         !Premixed binder
qg(3) = cp*(th_st - tcold)*(abs(yf(3,3))+abs(yf(2,3))) - qheat_ap * abs(yf(3,3)) - qb_heat * abs(yf(2,3))   !Primary flame
!qg(4) = cp*(th_final - th_ap)* abs(yf(4,4)) + cp*(th_final - th_bind)* abs(yf(5,4))                         !Final diffusion flame
qg(4) = cp*(th_st - th_ap)* abs(yf(4,4)) + cp*(th_st - th_bind)* abs(yf(5,4))                         !Final diffusion flame

if(myid.eq.0) then
    write(6,*)
    write(6,111)'Flame Temps = ',th_final,th_st,th_bind,th_ap
    write(6,*)'beta = ',-yf(3,3)
    write(6,*)
endif
!     
!     SET coefficients of reaction terms 
!     

yf(1,1:irxn) = qg(1:irxn)/cp

if(myid.eq.0) then
    write(6,*) "Input Parameters:"
    write(6,*)
    write(6,*)'pressure = ',press
    write(6,*)'alpha_pack = ',alpha_pack
    write(6,*)'alpha_binder =', alp_W
    write(6,*)'Final Diffusion Flame, AP coeff. = ',-yf(4,4)
    write(6,*)'Reaction Parameters'
    write(6,*) '          A        n       E         Q'
    do m=1,irxn
        write(6,112)'#',m,da(m),np(m),thetay(m),qg(m)
    enddo
    write(6,113)'qheat_solid  = ',qheat_ap,qheat_binder
    write(6,114)'lambda_solid  = ',lambda_ap,lambda_binder
    write(6,*)'nx,ny  = ',nx,ny
    write(6,*)'xend,yend = ',xend,yend
    write(6,*)'c1y = ',  c1y
    write(6,*)'theta_ap,theta_binder = ',theta_ap,theta_binder
    write(6,*)'da_ap,da_binder = ',da_ap,da_binder
    write(6,*)
endif

111 format (1x,a13,4(2x,f7.1))
112 format(1x,a1,i1,3x,e10.3,3x,f3.1,3x,f7.1,3x,f7.1)
113 format(1x,a15,2(2x,f7.1))
114 format(1x,a16,p2(2x,e10.4))
return
end SUBROUTINE global_parameter
!******************************************************************************

!******************************************************************************
subroutine run_parameter

USE GLOBAL_DATA
USE MYMPI

IMPLICIT NONE
!---------------------------------------------------------------------
! Local     Variables
integer ::   i,m,n,flag,num_p,num_d
REAL*8 :: pureap, dummy, p60, t60, mixap, rb_bind, top, bot, tt, temp(5)
REAL*8 :: p_cor(15),d_cor(15),rate_cor(15,15),reduce(5),sum,reduce_tot,new_ap
REAL*8 :: a,b,c,t1,t2          !Bisection variables
!--------------------------------------------------------------------

call read_chem_info 
p60 = 60.0d0      !Reference pressure

sum = 0.0d0
reduce_tot = 0.0d0

if(myid == 0)then
!Calculate the monopropellant flame temperature
    flag = 0
    pureap = 0.9999
    call edwards_input(flag,pureap,temp(1))

!Calculate the premixed binder flame temperature
    if(alp_W == 1.0) then
        dummy = pureap
    else
        dummy = alp_W
    endif
    call edwards_input(flag,dummy,temp(2))

!Calculate primary flame temperature
    call edwards_input(flag,rho_m_st,temp(3))

!Calculate final diffusion flame temperature
    mixap = alpha_pack
    if(alp_W == 1.0) mixap = pureap
    if(ipack == 0) mixap = 0.86d0
    call edwards_input(flag,mixap,temp(4))
    
!Calculate flame temperature for binder rate correlation at 60 atm
    flag = 1
    call edwards_input(flag,dummy,temp(5))
    
!Modify binder values to account for decreased particle size effect with increased pressure
!Decrease the amount of AP homogenized with the binder, based on Diffusion Code calculations
    if(bind_sizes /= 0 .and. alp_W /= 1.0d0)then
        call br_corr(num_p,num_d,p_cor,d_cor,rate_cor)  !Read in correlations
        do i=1,bind_sizes
            if (bind_diam(i)<1.0) bind_diam(i) = 1.0d0
            call polin2(p_cor,d_cor,rate_cor,num_p,num_d,pcyc,bind_diam(i),reduce(i))
            !Sum the amount to reduce the rate based on a simple mass average
            sum = sum + bind_mass(i)
        enddo
        do i=1,bind_sizes
            reduce_tot = reduce_tot + reduce(i)*bind_mass(i)/sum
        enddo
      
        !Calculate reduction in burning rate
        rb_bind = 0.537d0*pcyc**0.833d0*exp(-4056.0d0/temp(2))
        t2 = -4056.0d0/(dlog(reduce_tot*rb_bind/(0.537d0*pcyc**0.833d0)))
        flag = 0
        
        ! Alternate method because of large differences in T at <75% AP content
        if(alp_W < 0.75d0 .and. pcyc > p60) then
            rb_bind = reduce_tot*rb_bind*(p60/pcyc)**0.833d0
            t2 = -4056.0d0/(dlog(rb_bind/(0.537d0*p60**0.833d0)))
            flag = 1
        endif

        !Determine new AP% in binder with bisection method
        a = 0.001
        b = alp_w
        t1 = 0.0d0
        do while (abs(t1-t2) > 1.0)
            c = (a+b)/2.0d0
            call edwards_input(flag,c,t1)
            if(t1-t2 < 0.0d0)then
                a=c
            else
                b=c
            endif
        enddo
        write(*,*) 
        write(*,*) "New AP% in binder = ",c*100.0d0
    
        if(alp_W < 0.75d0 .and. pcyc > p60) then
            call edwards_input(0,c,temp(2))
            temp(5) = t2
        else
            temp(2) = t2
            call edwards_input(1,c,temp(5))
        endif
    endif
endif

call MPI_BCAST(temp,5,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

th_ap = temp(1)
th_bind = temp(2)
th_st = temp(3)
th_final = temp(4)
t60 = temp(5)

!Determine rate correlation for AP (60 atm)
rb_melt = 0.73d0 
da_ap = rb_melt*exp(theta_ap/T_melt)

!Determine rate correlation for AP/HTPB binder
T_bmelt = 976.4827*exp(-105.0361/t60)

rb_bind = 0.537d0*p60**0.833d0*exp(-4056.0d0/t60)   !AP/HTPB reference rates based on Beckstead trend at 60 atm
da_binder = rb_bind*exp(theta_binder/T_bmelt)

!Convert condensed phase lambda
lambda_ap =     lambda_ap / (4.184d0 * 100.0d0)
lambda_binder = lambda_binder / (4.184d0 * 100.0d0)

! add Correlation for the pre-exponential factor (fit using Foster data)   
if(da(2) /= 0.0d0)then      
    da(2) = 917.966615546023d0*exp(3714.51410005178d0/t60)
endif
!-----------------------------------------------------------------------------
return
end subroutine run_parameter
!******************************************************************************


!******************************************************************************
subroutine read_chem_info

USE GLOBAL_DATA
USE STRINGS

IMPLICIT NONE
!---------------------------------------------------------------------
! Local variables
integer :: k,m,n,ios,nargs,nsplit(irxn,2),counter
REAL*8 :: abind,bbind,cbind,dbind,ebind,aa
REAL*8 :: coef_re(irxn,neqgas),coef_pr(irxn,neqgas)
character*8 :: sp1(neqgas),spec(neqgas)
Character*6 :: reacts(irxn,neqgas),prods(irxn,neqgas)
Character*20 :: kinetics(neqgas)
Character*100 :: dummy,reaction,rp(2),params
!--------------------------------------------------------------------

open(UNIT=16,FILE="chem_info.txt")

!Set some parameters
rho_m_st = 0.88856203
lewis(1:neqgas) = 1.0
mu = zero

!Read in first line and format answer
call readline(16,dummy,ios)
call removesp(dummy)
dummy = uppercase(dummy)
if(dummy /= 'SPECIES')then
    write(*,*) "Enter species in input file" ;  Pause;   stop
endif

k=3
! Read in the species
do n=1,neqgas
    call readline(16,sp1(n),ios)
    call removesp(sp1(n))
    sp1(n) = uppercase(sp1(n))
!Put the species in order: HTPB, AP, ...    
    if(sp1(n) /= 'HTPB' .AND. sp1(n) /= 'AP')then
        spec(k) = sp1(n)
        k = k +1
    endif        
enddo
spec(1) = 'HTPB'
spec(2) = 'AP'

! Read in the reactions
call readline(16,dummy,ios)
call removesp(dummy)
dummy = uppercase(dummy)
if(dummy /= 'REACTIONS')then
    write(*,*) "Number of species is incorrect"; Pause; stop
endif

! Set reaction coefficients to 0
yf(2:neqgas,1:irxn) = 0.0d0
coef_re(:,:) = 0.0d0
coef_pr(:,:) = 0.0d0

!Read in each reaction
do m=1,irxn
    call readline(16,reaction,ios)            ! Read in each reaction
    reaction = uppercase(reaction)
    call parse(reaction,'=',rp,k)         ! Parse into reactants and products

!Parse reactants
    call name_split(rp(1),coef_re(m,:),reacts(m,:),nsplit(m,1))

!Parse products     
    call name_split(rp(2),coef_pr(m,:),prods(m,:),nsplit(m,2))
    
!read in kinetic parameters    
    call readline(16,params,ios)
    call compact(params)
    call parse(params,' ',kinetics,k)
    if (k > 3) Then
        write(*,*) "Error in kinetics, parameter error in reaction",n ;pause; stop
    endif
    call value(kinetics(1),da(m),ios)
    call value(kinetics(2),np(m),ios)
    call value(kinetics(3),thetay(m),ios)
enddo        

if(myid.eq.0) then
    do m=1,irxn
        reaction = ' '                           
        call reaction_name(coef_re(m,:),reacts(m,:),coef_pr(m,:),prods(m,:),reaction)
        write(*,'(a60)') reaction
        write(*,311) 'DA(',m,') =',da(m)
        write(*,312) 'np(',m,') =',np(m)
        write(*,313) 'Ea(',m,') =',thetay(m)
    enddo
endif

!Set the reaction coefficients
do m=1,irxn
    counter = 1
    do nargs= 1,nsplit(m,1) !loop for possible number of coefficients
        call removesp(reacts(m,nargs))
        do n = 1,neqgas -1    !loop for each species except final
            if(reacts(m,nargs) == spec(n))then
                yf(n+1,m) = -coef_re(m,nargs)
                mu(m,n+1) = 1.0d0
                counter = counter + 1
            endif
        enddo      
    enddo
    do nargs= 1,nsplit(m,2)
        call removesp(prods(m,nargs))
        do n = 1,neqgas-1     !loop for each species except final
            if(prods(m,nargs) == spec(n))then
                yf(n+1,m) = coef_pr(m,nargs)
                counter = counter + 1
            endif
        enddo
    enddo
enddo

!Read in final properties
call readline(16,dummy,ios)
call removesp(dummy)
dummy = uppercase(dummy)
if(dummy /= 'PROPERTIES')then
    write(*,*) "Input Error, check number of reactions"
    Pause
    stop
endif

read(16,*) qheat_ap
read(16,*) qheat_binder

read(16,*) theta_ap
read(16,*) theta_binder  

read(16,*) T_melt
read(16,*) lambda_ap
read(16,*) lambda_binder
read(16,*) coe_lamb1 !Gas conductivity coe_lamb1*T + coe_lamb2
read(16,*) coe_lamb2 !Gas conductivity
close(16)

311 format (2x,a3,i1,a3,1x,1p,e11.4)
312 format (2x,a3,i1,a3,1x,f7.3)
313 format (2x,a3,i1,a3,1x,f7.1)

!-----------------------------------------------------------------------------
return
end subroutine read_chem_info
!******************************************************************************

subroutine name_split(name,coef,gas,nargs)
USE GLOBAL_DATA
USE STRINGS

IMPLICIT NONE
!---------------------------------------------------------------------
INTEGER, INTENT(OUT) :: nargs
REAL*8, INTENT(OUT) :: coef(neqgas)
CHARACTER*80, INTENT(IN) :: name
CHARACTER*8, INTENT(OUT) :: gas(neqgas)
! Local variables
integer :: m,ios
REAL*8 :: num
CHARACTER*14 :: dummy(neqgas),spliter,before
CHARACTER*26 :: delims
!--------------------------------------------------------------------

! Define the deliminators
delims = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'

! Zero out the output
coef(:) = 0.0d0
gas(:) = ' '

! Parse name
call parse(name,'+',dummy,nargs)     
do m = 1,nargs                            ! Loop through each species looking for numbers
    spliter = ' '
    call split(dummy(m),delims,before,spliter)
    call insertstr(spliter,dummy(m),2)
    call removesp(spliter)
    spliter = adjustl(spliter)
    gas(m) = spliter
    if(before == ' ')then
        num = 1
    else
        call value(before,num,ios)
    endif
    coef(m) = num 
enddo        

!-----------------------------------------------------------------------------
return
end subroutine name_split
!******************************************************************************

subroutine reaction_name(c_r,reac,c_p,prod,name)
USE GLOBAL_DATA
USE STRINGS

IMPLICIT NONE
!---------------------------------------------------------------------
REAL*8, INTENT(IN) :: c_r(neqgas),c_p(neqgas)
CHARACTER*6, INTENT(IN) :: reac(neqgas),prod(neqgas)
CHARACTER*100, INTENT(OUT) :: name
! Local variables
integer :: m,n,plus
CHARACTER*8 :: dummy
!--------------------------------------------------------------------

plus = 0
n = 1
do m=1,neqgas
    if(c_r(m) /= 0.0)then
        if(plus ==1 )  then
            call insertstr(name,'+',n)       
            n = n+6   
        endif
        dummy = ' '
        call writenum(c_r(m),dummy,'f7.3')
        call trimzero(dummy)
        call insertstr(name,dummy,n)
        n = n + 6
        call insertstr(name,reac(m),n)
        n = n + 6
        plus = 1
    endif
enddo

dummy = ' ='
plus = 0    
call insertstr(name,dummy,n)
n = n + 6
do m=1,neqgas
    if(c_p(m) /= 0.0)then
        if(plus ==1 )  then
            call insertstr(name,'+',n)   
            n = n+6       
        endif
        dummy = ' '
        call writenum(c_p(m),dummy,'f7.3')
        call trimzero(dummy)
        call insertstr(name,dummy,n)
        n = n + 6
        call insertstr(name,prod(m),n)
        n = n + 6
        plus = 1
    endif
enddo
call compact(name)

!-----------------------------------------------------------------------------
return
end subroutine reaction_name
!******************************************************************************