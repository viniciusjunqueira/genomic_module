! ===================================================================================
!Copyright (C) 2008, 2012, 2016 Fernando Cardoso, Vinicius Junqueira, Rodrigo Mota
!
!    This program is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with this program.  If not, see <http://www.gnu.org/licenses/>.
! ===================================================================================

PROGRAM intergen
!==========================================================================================!
! Fortran F95 program to fit quantitative genetics models, including
! multiple-trait, multiple-breed and multiple-sires (reduced) animal models, 
! random regression and reaction norms models, using MCMC techniques.
! Requires sparse matrix libraries, modules and subrotines available from I. Misztal at UGA
! Version 1.3 - Feb, 2015
! by Fernando F. Cardoso (FC).
!------------------------------------------------------------------------------------------!
! Historic of development:
! - Metropolis-hastings algorithm included a residual heteroskedasticity using a 
! continuous function on the environmental gradient with unknown covariates (Jan, 2010)
! - Reaction norms with unknown covariates (Dec, 2006)
! Random regression models to study genotype by environment interaction via reaction norms 
! based on Su et al. (2006) were included
! - Robust structural model for residual variances (May, 2003)
! Features to allow heteroscedastic and robust models (Student t & Slash distribution)
! based on Cardoso et al. (2005) were implemented on top of previous developements
! - Multibreed animal model (April, 2002)
! New subroutines were added to allow estimatition of breed specific and segregation 
! genetic variances on crossbred populations (Cardoso & Tempelman, 2004). 
! - Uncertain paternity features (Sep, 2001)
! Modification of Ignacy Misztal GIBBSF90 program to run a reduced animal model (RAM) 
! with uncertain paternity (multiple sires), maternal effects and genetic groups: 
! allows for a FULLY BAYES (Cardoso & Tempelman, 2003) or AVERAGE A (Henderson, 1988) 
! approaches on uncertain paternity 
!==========================================================================================!

!==========================================================================================!
! Historic of modifications:
! - BLUP was implemented using Misztal et al. algorithms (Feb.,2015). It is possible the
! definition of Solve Methods (PCG,FSPAK,SOR), Convergence Criteria and the Stardant Error 
! of Prediction (SEP - performed by the inverse of LHS). Prints were modified too. 
! By Vinicius Junqueira (VJ), Henry Carvalho (HC) and Bruno Teixeira (BT)
! - A subroutine to compute dinv was incorported such that user dont have to put it in pedi-
! gree file for analyses that include multiple sires (MS) (Feb.,2015). By Vinicius Junqueira (VJ)
! - Included an adapted subroutine from Misztal BLUPF90 to compute yhat, adjusted y and residuals. 
! (May.,2015). By Vinicius Junqueira (VJ)
! - Included a subroutine which adjust a normal prior effect for fixed effects . 
! (October 2015). By Vinicius S Junqueir (VJ)
! - Included subroutines from CBLUPTHR multivariate Threshold models. The first trait defined
! on section OBSERVATION(S) must be the categorical trait. (July 2016). By Vinicius Junqueira
!==========================================================================================!
use model; use gibbs; use prob; use pcg; use textop; use ranlib

implicit none

real (r8), parameter:: pi=3.14159265358979

real,allocatable :: y(:),&                  ! observation value
                    weight_y(:),&           ! FC weights for records - trait specific
                    indata(:,:)             ! FC all lines of input data

real (r8)::         weight_r                ! FC weights for residuals

type (sparse_hashm)::xx
type (sparse_ija):: xx_ija
type(sparse_hashm),allocatable::ainv(:)
type(sparse_ija),allocatable::ainv_ija(:)

real (rh), allocatable :: xy(:),sol(:),&		!storage for equations
                          vb(:)             !storage for normal prior effects

real,allocatable :: weight_cov(:,:)

integer,allocatable:: address(:,:),&    ! start and address of each effect
                      df_random(:)

integer :: neq,io,&                   ! number of equations and io-status
           data_len,&                 ! length of data record to read
           i,j,k,l                    ! extra variables          

real (rh) ::  val
real :: dat_eff,second   
real (r8) :: conv
real (r8),allocatable::orig_r(:,:),orig_g(:,:,:)	! keep original r and g 


! FC Model choice vectors and variables for CPO-PBF and DIC
real (r8),allocatable:: logl_obs(:,:)   ! logl_obs has, for each obs, the harmonic 
real (r8):: logl_rnd		            ! mean of l (col=1) the mean logl (col=2) &
										! logl at posterior means of theta (col=3) & e/sigmae (col=4) .
										! logl_rnd has logl for all obs at each round


! FC Multiple-sire (ms) and reduced animal model (ram) variables and arrays
real (rh), allocatable:: sol_n(:,:),&		! ram - storage for nonparent solutions
			 probs(:),&			! ms - probalities of each msire
                         ms(:),&            ! ms info: n of possible sires,
 		         ms_alphas(:)       ! sire ID 1, prob 1, ..., prob n;
                                            ! ms - hyperpar of Dirichlet distn of sire assignments

real,allocatable :: dinv(:)			! the inverse of proportion of var due to mendelian sampling

integer,allocatable:: ped(:,:,:),&  ! pedigree data (could be multiple pedigrees)
 		      p_ms(:), &    ! pointer to ms vector, 
		      ms_stat(:)    ! indicates if animals has multiple-sires by (1)

integer :: irec,naddr,&         ! # of addresses including ms.
	   maxnaddr,sirenew     ! maximum # of addresses, sampled sire

! FC Multiple-breed (mb) variables and arrays
real (rh), allocatable:: breed(:,:,:)       ! animal's breed composition 

integer,allocatable:: mb_dbelief(:,:,:)   ! degrees of belief of proposal densities in MH

! FC Structutal heterogeneous variances
real (r8), allocatable:: r_sol(:),r_nrec(:),& ! solutions for structural effects on r & # of records
			 r_alpha(:), r_alpha_varprop(:),& ! or sum of covariate values for each effect on r_sol
			                                  ! Heterogeneity parameter and variance of proposal density
 		       & h_varprop(:) ! and for unknown covariate heteroskedastic

integer,allocatable:: r_dbeliefprop(:)    ! degrees of belief of proposal densities in r MH

! FC Robust errors
real (rh), allocatable:: weight_robust(:)

real (r8):: w_alpha, sumw_alphaaccep,w_alpha_varprop, & ! robustness parameters, its acceptance rate,
			sumw_alpha,ssw_alpha					    ! variance of proposal density, sum and ss

integer:: w_dbeliefprop    ! degrees of belief of proposal densities in w_alpha MH

! FC save samples mean and variance
real (r8),allocatable:: sumsol(:),sssol(:),&         ! solutions  
			sumsol_n(:,:),sssol_n(:,:),& ! for non-parents,  
			summs(:),ssms(:),&           ! sire assignment proportion & probs   
                        sumG(:,:,:),sumR(:,:),ssG(:,:,:),ssR(:,:),& ! and varcomps
                        sumr_sol(:),ssr_sol(:),&        ! solutions for structural effects on r 
			summb_accep(:,:,:),&          ! acceptance rates for mb Metropolis
			sumh_accep(:),&       ! acceptance rates for unknown covariate
			sumr_accep(:),sumr_alphaaccep(:),&  ! acceptance rates for r and r_alpha MH
			sumr_alpha(:),ssr_alpha(:),&  !Heterogeneity parameters sum and ss
			sumweight_robust(:),ssweight_robust(:)  !weights for robustness

! FC cpu time
real:: c_start,c_now,c_last  ! keeps time

! VJ Threshold
integer, allocatable :: place(:)
real (rh)::xthr,lt
real (rh), allocatable :: xrhs(:),xlhs(:,:),xmi2(:),bwb(:,:),&
                          rci(:,:),rc(:,:),b(:)


call copyright

call read_parameters

call print_parameters('INTERGEN 1.3')

call set_seed(seed)

maxnaddr=neff+maxnsires
neq=ntrait*sum(nlev)
data_len=max(maxval(pos_weight),maxval(pos_y),maxval(pos_eff),maxval(nestedcov))
print*,'Data record length = ',data_len

allocate (xy(neq+nthr),sol(neq+nthr),vb(neq+nthr),address(maxnaddr,ntrait),&
          weight_cov(maxnaddr,ntrait),y(ntrait),weight_y(ntrait),indata(nrec,data_len),&
          ainv(neff),ainv_ija(neff),df_random(neff),ped(neff,maxval(lenped),3),&
          p_ms(maxval(lenped)),ms_stat(maxval(lenped)),probs(maxnsires+2),&
          orig_r(ntrait,ntrait),orig_g(size(g,dim=1),size(g,dim=2),size(g,dim=2)),&
          ginv(size(g,dim=1),size(g,dim=2),size(g,dim=2)),&
          sumG(neff,maxcorr,maxcorr),sumR(ntrait,ntrait),&
          ssG(neff,maxcorr,maxcorr), ssR(ntrait,ntrait),&
          sumSOL(neq+nthr),ssSOL(neq+nthr),logl_obs(nrec,3+ntrait),&
          place(ntrait),xmi2(ntrait),bwb(ntrait,ntrait),XLHS(NTHR,NTHR),&
          XRHS(NTHR), b(ntrait+1))

call zerom(xx,neq+nthr); xy=0; sol=0; vb=0; orig_g=g; orig_r=r; ms_stat=0; probs=1.
sumG=0;sumR=0;ssG=0;ssR=0;sumSOL=0;ssSOL=0;logl_obs=0
xlhs=0;xrhs=0

! FC Structural heterogeneous residual
if (neffr > 0) then
  allocate (r_sol(neq),r_nrec(neq),sumr_sol(neq),ssr_sol(neq),sumr_accep(neffr),&
            r_dbeliefprop(neffr),r_alpha(neffr),sumr_alphaaccep(neffr),sumr_alpha(neffr),&
            ssr_alpha(neffr),r_alpha_varprop(neffr),sumh_accep(neq),h_varprop(neff))

  r_sol=1.;r_nrec=0.;sumr_sol=0;ssr_sol=0;sumr_accep=0;r_dbeliefprop=2000;
  r_alpha=3;sumr_alphaaccep=0.;sumr_alpha=0.;ssr_alpha=0.;r_alpha_varprop=.1;
  sumh_accep=0; h_varprop=100
endif

! FC Robustness
if (residualtype == r_slash .or. residualtype == r_struct_slash .or. &
    residualtype == r_studentt .or. residualtype == r_struct_studentt) then
  allocate (weight_robust(nrec),sumweight_robust(nrec),ssweight_robust(nrec))
  weight_robust=1.; sumweight_robust=0.; ssweight_robust=0.
  w_alpha=1.5;sumw_alphaaccep=0.;sumw_alpha=0.;ssw_alpha=0.;w_alpha_varprop=.1;
endif

! FC read pedigree and ms (if uncertain paternity present) files
call read_files

if (rstart == 'y') call restart(start) !get # info previous chain

open(io_d,file=datafile)         !data file

! FC time control
call cpu_time(c_start)
c_last=c_start

!VJ - HG - Splits the program in MCMC or BLUP analisys
! Variable nround is initialize assuming zero. If MCMC, nround will be greater than zero
if (nround > 0) then
   do round = start,nround
     call setup_g                    ! invert G matrices
     xy=0; logl_rnd=0;
     !FC Uses average paternity in round 1 so that all possible positions in 
     !xx are initialized in round 1 avoiding problems with the link_hash_aij
     xx%x(3,:)=0	! zero only the data part
     !Sets up mixed model equations
     if (any(effecttype == effrnorm)) then
       call setup_eq_rnorm 
       call gibbsvar
     else 
       call setup_eq
       call gibbssol
       call gibbsvar
     endif
   end do

   round=nround !FC to write the correct number of rounds

   if (round > nburn) call deviance

   call store_solutions

   call cpu_time(c_now)
   print*,' Time elapsed ',(c_now-c_start)/60,' min',&
	        '(',(nround-start)*60/(c_now-c_start), 'rounds/min )' 

! VJ - HCG If nround equals zero, blup will be performed
else
   ! Linear BLUP analysis
   if(cat_value.eq.-1) then        
     call setup_g                    ! invert G matrices
     call setup_eq                   ! setup equations
     ! Solution methods
     if (solv_method == 'FSPAK') then      
        xx_ija=xx
        call fspak90('factorize',xx_ija)
        call fspak90('solve',xx_ija,rs8=xy,sol8=sol)
     elseif(solv_method == 'SOR') then
        call default_iter(conv=conv_crit,maxround=maxrounds,relax=r_factor,zerosol=lzerosol)
        call solve_iterm(xx,xy,sol)
     elseif(solv_method == 'PCG') then      
        call default_iter(conv=conv_crit,maxround=maxrounds,zerosol=lzerosol)
        call solve_pcg(xx,xy,sol,blksize)
     else
        print*,'Wrong solve_method ',solv_method
        stop
     endif
   ! Categorical BLUP analysis
   else
     call setup_eq_threshold
   endif

   if (neq <= 10) print  '( '' solution:'' ,10(/10f8.2))', sol
   !Store Blup solutions
   print*,'' 
   call store_solutions_blup
   !Print elapsed time to solve MME
   call cpu_time(c_now)
   !print*,'' 
   print ('(1x,A,1x,F5.2,1x,A)'), 'Elapsed time to solve MME',(c_now-c_start)/60,' minutes'
   !print*,''

   ! VJ Run function to compute yhat, adjusted y and residuals
   if(cat_value.eq.-1) then        
     call write_adjusted_data
   endif

endif



contains


 subroutine find_rinv
! From Blupf90
! calculates inv(Q R Q), where Q is an identity matrix zeroed for 
! elements corresponding to y(i)=miss
  interface
    subroutine ginv1(a,n,nmax,tol,rank)
     use kinds
     integer::n,nmax,rank
     real (rh)::a(nmax,nmax),tol
    end subroutine
  end interface

  integer :: i,irank
!
  rinv=r
  do i=1,ntrait
     if (y(i) == miss) then
        rinv(i,:)=0; rinv(:,i)=0
     endif
  enddo
  call ginv1(rinv,ntrait,ntrait,real(1d-8,rh),irank)
  end subroutine
 

subroutine setup_eq
    integer :: i,j,k,l
    real (rh) ::  val
  ! setup xx and xy
	do irec=1,nrec
	   ! FC keep data in memory 
	   if (round == start) then 
		 read(io_d,*,iostat=io) indata(irec,:)
		 if (io.ne.0) then
		   print*,'total number of record in data set: ', irec-1, ' is not equal &
			 &  to the value you have entered: ', nrec
		   stop
		 endif
	   endif
	   !FC

! if (round == start .and. irec==1) print*,'before decode_record'

	   call decode_record
	   call find_addresses
	   !FC if necessary changes r for ram or structural variances

     ! VJ Perform a function selection depending if MCMC or Blup analysis
	 if(nround>0) then
        call predict_missing_y
	    call decode_r
     else
        call find_rinv
     endif

	   !FC
	   if (round == start .and. irec<10) then
		 print*,'    trait    eff    addr  weight_cov  weight_y' 
		 write (*,102) ((j, i, address(i,j), weight_cov(i,j), weight_y(j), j=1, ntrait), i = 1, naddr) 
          ! VJ Changes spaces format
		 102 format (10X, 3I8, 2F10.3)

	   endif
	   !FC
	   do i=1,naddr
       do j=1,naddr
		     do k=1,ntrait
		       do l=1,ntrait
		         val=weight_cov(i,k)*weight_cov(j,l)*sqrt(weight_y(k))*sqrt(weight_y(l))*rinv(k,l)
	!                 if (val /= 0) &
			       call addm(val,address(i,k),address(j,l),xx)
		       enddo
		     enddo     
	     enddo
	     do k=1,ntrait
		     do l=1,ntrait
		        xy(address(i,k))=xy(address(i,k))+sqrt(weight_y(k))*sqrt(weight_y(l))*rinv(k,l)*y(l)*weight_cov(i,k)
		     enddo 
	     enddo     
	   enddo
	enddo

	if (round == start) then
	  call cpu_time(c_now)
	  print*,'read ',nrec,' records in ',c_now-c_start,' s,  ',xx%filled,' nonzeroes'
	endif

	! Random effects' contributions
	do i=1,neff
	   select case (randomtype(i))
	     case (g_fixed)
		      continue                ! fixed effect, do nothing
	     case (g_diag)
		      call add_g_diag(i)
	     case (g_A, g_As, g_A_UPG,g_A_UPG_INB)
		      rewind i+io_off
		      call add_g_add(randomtype(i),i)
	     case (g_A_MS)		! FC add average A to parent equations
		      call add_g_aav(i,maxnsires)
	     case (g_A_MB, g_diag_MB)  ! FC add mbreed Ginv to xx
		      if(max_msires==0) then
			    if(nbreed(i)==1) then
                            rewind i+io_off
                            call add_g_add(g_A,i)
			    else
			      call add_g_mb(i,maxnsires,randomnumb(i)*ntrait)
	            endif
		      else
			    call add_g_mbs(i,maxnsires,randomnumb(i)*ntrait)
		      endif
	     case (g_usr)
        	rewind i+io_off
         	call add_g_usr(i,'')
             case (g_usr_inv)
         	rewind i+io_off
                call add_g_usr(i,'inv')
             case(g_normal)                ! VJ add contributions due to normal prior effects
                rewind i+io_off
                call add_g_norm_prior1(i)
	     case default
	       print*,'unimplemented random type',randomtype(i)
	   endselect
	enddo

	if (round == start) then
	  call cpu_time(c_now)
	  print*,'finished peds in ',c_now-c_start,' s,  ',xx%filled,' nonzeroes'
	  if (neq <= 16) then
            print*,''
	    print*,'left hand side'
	    call printm(xx)
	    print  '( '' right hand side:'',10(/10f8.2))',xy
	  endif
	endif

end subroutine

 subroutine setup_eq_rnorm
   integer :: i,j,k,l,t,addr_i,addr_j
   real (rh) ::  val=0.0

! Setup xx and xy for reaction norm model

   do t=1,2 
  ! t=1 for unknown covariate part and t=2 for theta part
     do irec=1,nrec
   ! FC keep data in memory 
       if (t == 1 .and. round == start) then 
         read(io_d,*,iostat=io) indata(irec,:)
         if (io.ne.0) then
           print*,'total number of record in data set: ', irec-1, &
              ' is not equal to the value you have entered: ', nrec
	   stop
         endif
       endif
       call decode_record ! necessary to predict missing values
       call find_addresses
      !FC if necessary changes r for ram or structural variances
       call decode_r
       call predict_missing_y
       call decode_record_rnorm(t) ! calculates y_h(t=1) or y_theta(t=2)
       do i=1,neff
	 if (t == 1 .and. effecttype(i) /= effuncov) then     !C(h) part
           cycle
	 elseif (t == 2 .and. effecttype(i) == effuncov) then !C(theta) part
           cycle
         endif 
         do j=1,neff
	   if (t == 1 .and. effecttype(j) /= effuncov) then     !C(h) part
             cycle
	   elseif (t == 2 .and. effecttype(j) == effuncov) then !C(theta) part
             cycle
           else
             do k=1,ntrait
 	       do l=1,ntrait
		 val=weight_cov(i,k)*weight_cov(j,l)*sqrt(weight_y(k))*sqrt(weight_y(l))*rinv(k,l)
		 call addm(val,address(i,k),address(j,l),xx)
               enddo
             enddo
           endif 
	 enddo
! XY part
	 do k=1,ntrait
	   do l=1,ntrait
	     xy(address(i,k))=xy(address(i,k))+sqrt(weight_y(k))*sqrt(weight_y(l))*rinv(k,l)*y(l)*weight_cov(i,k)
 	   enddo 
	 enddo     
       enddo
     enddo

     if (t == 1 .and. round == start) then
       call cpu_time(c_now)
       print*,'read ',nrec,' records in ',c_now-c_start,' s,  ',xx%filled,' nonzeroes'
     endif

! Random effects' contributions
     do i=1,neff
       if (t == 1 .and. effecttype(i) /= effuncov) then     !C(h) part
         cycle
       elseif (t == 2 .and. effecttype(i) == effuncov) then !C(theta) part
         cycle
       endif 
       select case (randomtype(i))
         case (g_fixed)
 	         continue                         ! fixed effect, do nothing
         case (g_diag)
 	         call add_g_diag(i)
         case (g_A, g_As, g_A_UPG,g_A_UPG_INB)
           rewind i+io_off
   	       call add_g_add(randomtype(i),i)
         case (g_A_MS)	                  	! FC add average A to parent equations
  	       call add_g_aav(i,maxnsires)
         case (g_A_MB, g_diag_MB)           ! FC add mbreed Ginv to xx
 	         if (max_msires==0) then
	           if (nbreed(i)==1) then
               rewind i+io_off
               call add_g_add(g_A,i)
	           else
	             call add_g_mb(i,maxnsires,randomnumb(i)*ntrait)
	           endif
 	          else
	           call add_g_mbs(i,maxnsires,randomnumb(i)*ntrait)
 	          endif
	        case (g_usr)
            rewind i+io_off
            call add_g_usr(i,'')
        	case (g_usr_inv)
            rewind i+io_off
            call add_g_usr(i,'inv')
          case(g_normal)             ! VJ add normal prior contributions
            rewind i+io_off
            call add_g_norm_prior1(i)
          case default
            print*,'unimplemented random type',randomtype(i)
       endselect
     enddo

     if (round == start) then
       call cpu_time(c_now)
       print*,'finished peds in ',c_now-c_start,' s,  ',xx%filled,' nonzeroes'
       if (neq <= 20) then
         print*,'left hand side'
         call printm(xx)
         print  '( '' right hand side:'',10(/10f8.2))',xy
       endif
     endif

     call gibbssol_rnorm(t) ! samples h(t=1) or theta(t=2)
   enddo 

 end subroutine


 subroutine decode_record
  integer :: i
  ! decodes data record into y's and weight

   do i=1,ntrait
      y(i)=indata(irec,pos_y(i))
  ! FC trait specific weights
      if (pos_weight(i) /= 0) then
        weight_y(i) = indata(irec,pos_weight(i))
        if (weight_y(i) <= 0) then
          print*,'Non-positive weight for trait ',i, ' record ', irec
          stop
        endif
      else
        weight_y(i)=1
      endif
   enddo

!if (round == start .and. irec==1) print*,'enter decode_record', irec , ' y ', y,'  weight_y ', weight_y, &
!                                         ' pos_y ', pos_y,'  pos_weight ', pos_weight

  end subroutine

  subroutine decode_record_rnorm(t)
    integer :: i,j,t
  ! adjusts data record for reaction norm model
    do i=1,ntrait
      do j=1,neff
	if (effecttype(j) == effrnorm) then ! don't change record just increase addr
	  cycle
	elseif (effecttype(j) == effuncov .and. t == 2) then ! adjust y(theta)
	  y(i) = y(i)-sol(address(j,i))  !weight_cov not used because includes sol(rnorm)
	elseif (effecttype(j) /= effuncov .and. t == 1) then ! adjust y(h)
	  y(i) = y(i)-weight_cov(j,i)*sol(address(j,i))
        endif
      enddo   
    enddo
  end subroutine

  subroutine decode_r
    integer :: i,j,k,addr
	!FC add mendelian variance to residuals
	!Single trait solution

	rinv=r; weight_r=1.
  ! FC structural variances
	if (neffr > 0) then
      do k=1,ntrait
		do i=1,neffr
		  if (nestedcov(pos_effr(i),k) > 0) then !allows reacnorm & rr heteroskedasticity base on covariate value
			  addr=address1(pos_effr(i),1,k)
			else
			  addr=address(pos_effr(i),k)
		  endif
		  weight_r=weight_r*(r_sol(addr)**weight_cov(pos_effr(i),k))
		enddo
      enddo
	endif
  ! FC robustness
    if (residualtype == r_slash .or. residualtype == r_struct_slash .or. &
        residualtype == r_studentt .or. residualtype == r_struct_studentt) then
      weight_r=weight_r/weight_robust(irec)
    endif
  ! FC reduced animal model
    if (any(randomtype==g_A_MS) .and. any(effecttype==effram)) then	
	  k=minloc(randomtype,dim=1,mask=randomtype==g_A_MS)
	  if (int(indata(irec,k)) > nlev(k)) then
	    rinv=rinv+g(k,1:ntrait,1:ntrait)/dinv(int(indata(irec,pos_eff(k,1))))
	  endif
	endif
	rinv=finverse_s(rinv*weight_r)
  end subroutine


  subroutine find_addresses 
  integer :: i,j,k,l

  weight_cov=0.0
  do j=1,ntrait
     naddr=0
     do i=1,neff
        if (pos_eff(i,j) /= miss) then
            dat_eff=indata(irec,pos_eff(i,j))
          else
            dat_eff=miss
        endif
        if (dat_eff == miss) then            !missing effect
              naddr=naddr+1
              address(naddr,j)=1    !dummy address
              cycle
        endif
        select case (effecttype(i))
           case (effcross)
              naddr=naddr+1
              address(naddr,j)=address1(i,int(dat_eff),j)
              weight_cov(naddr,j)=1.0
           case  (effcov)
              naddr=naddr+1
              weight_cov(naddr,j)=indata(irec,pos_eff(i,j))
              if (nestedcov(i,j) == 0) then
                  address(naddr,j)=address1(i,1,j)
                elseif (nestedcov(i,j) > 0) then
                  address(naddr,j)=address1(i,int(indata(irec,nestedcov(i,j))),j)
                else
                  print*,'wrong description of nested covariable'
                  stop
              endif
 	 ! FC unknown covariate for reaction norm model
           case (effuncov)
              naddr=naddr+1
              address(naddr,j)=address1(i,int(dat_eff),j)
              weight_cov(naddr,j)=weight_cov(naddr,j)+1.0 
	 ! FC reaction norm model with unknown covariate
           case (effrnorm)
              naddr=naddr+1
              address(naddr,j)=address1(i,int(indata(irec,nestedcov(i,j))),j)
	      l=0         
	      do k=1,neff ! seeks the unknown covariate address
		            l=l+1     ! keeps track of naddr 
                if (effecttype(k) == effuncov) then
		               if (pos_eff(k,j) == pos_eff(i,j)) then ! k is the unknown cov position
 	            ! uses solution of unknown covariate as weight for reaction norm
                      weight_cov(naddr,j)=sol(address1(k,int(dat_eff),j)) 
	            ! adds the solutions of reaction norm to xx of unknown covariate  
                      weight_cov(l,j)=weight_cov(l,j)+sol(address1(i,int(indata(irec,nestedcov(i,j))),j))
                      l=0 ! uses l as flag
		                  exit
                    endif 
                endif
	      enddo
	      if (l /= 0) then	 
                  print*,'unknown covariable address not found'
                  stop
        endif

 	 ! FC for the ram w/msires
           case (effram)
              if (int(dat_eff) <= nlev(i)) then  !parent
                naddr=naddr+1
                address(naddr,j)=address1(i,int(dat_eff),j)
                weight_cov(naddr,j)=1.0
              elseif (int(dat_eff) > nlev(i)) then !non-parent
		            if (ped(i,int(dat_eff),2) > 0) then
		              naddr=naddr+1
                  address(naddr,j)=address1(i,ped(i,int(dat_eff),2),j)
                  weight_cov(naddr,j)=0.5
		            endif
		            if (ped(i,int(dat_eff),3) > 0) then
	            	  naddr=naddr+1
                  address(naddr,j)=address1(i,ped(i,int(dat_eff),3),j)
                  weight_cov(naddr,j)=0.5
	            	elseif (ped(i,int(dat_eff),3) < 0) then !msired animal
		              do k=1,2*ms(p_ms(int(dat_eff))),2
		                naddr=naddr+1
                    address(naddr,j)=address1(i,int(ms(k+p_ms(int(dat_eff)))),j)
                    weight_cov(naddr,j)=0.5*ms(k+p_ms(int(dat_eff))+1)
		              enddo
	            endif
              else
                print*,'wrong description of parent in ram for animal: ', int(dat_eff)
                stop
	            endif
           case default
              print*,'unimplemented effect ',i
              stop
        end select
     enddo 
  enddo
  end subroutine
 
  function address1(e,l,t)
! returns address for level l of effect e and trait t
  integer :: e,l,t, address1,i
  logical,save::first=.true.
  integer,allocatable,save::offset(:)
  !
  if (first) then
     first=.false.
     allocate(offset(neff))
     do i=1,neff
        offset(i)=sum(nlev(1:i-1))*ntrait
     enddo
  endif
  !
  address1= offset(e)+(l-1)*ntrait+t
  !address1= sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t
  end function

  
  subroutine setup_g
! inverts g matrices
  interface
    subroutine ginv1(a,n,nmax,tol,rank)
     use kinds
     integer::n,nmax,rank 
     real (rh)::a(nmax,nmax),tol
    end subroutine
  end interface

  integer :: rank
  ginv=g  !FC keep g for mendelian sampling subroutine
  do i=1,neff
     if (randomnumb(i).ne.0) then
!         call printmat(ginv(i,1:randomnumb(i)*ntrait,1:randomnumb(i)*ntrait),'G')
       if (randomtype(i) /= g_A_MB .and. randomtype(i) /= g_diag_MB) & 
         call ginv1(ginv(i,:,:),randomnumb(i)*ntrait,maxcorr,real(1d-18,rh),rank)
     endif
  enddo
  end subroutine

  subroutine add_g_diag(eff)
! adds diagonal (IID) contributions to MME
  integer :: eff, i,j,k,l,m,t1,t2
  real (rh) :: one=1.0
  if (round == start) call zerom(ainv(eff),nlev(eff))

  do i=1,nlev(eff)
     if (round ==start) call addm(one,i,i,ainv(eff))
     do j=0,randomnumb(eff)-1
        do k=0,randomnumb(eff)-1
           do t1=1,ntrait
              do t2=1,ntrait
                 m=address1(eff+j,i,t1); l=address1(eff+k,i,t2)
!                 xx(m,l)=xx(m,l)+ginv(eff,t1+j*ntrait,t2+k*ntrait)
                 call addm(ginv(eff,t1+j*ntrait,t2+k*ntrait),m,l,xx)
              enddo
           enddo
        enddo
     enddo
  enddo
  if (round ==start) df_random(eff)=nlev(eff)

  end subroutine


  subroutine add_g_add(type,eff)
! generates contributions for additive sire or animal effects
  integer :: type,eff,i,j,t1,t2,k,l,m,n,io,p(4),ped_len,iped
  real (rh) :: w(3),res(4),mendel
!
 
  select case (type)
     case (g_A)
        w=(/1., -.5, -.5/)
        res=(/2., 4/3., 1., 0./)
        ped_len=3
     case (g_A_UPG)
        w=(/1., -.5, -.5/)
        res=(/2., 4/3., 1., 0./)
        !ped_len=4
        ped_len=3                ! VJ - pedigree is allocated as a three-dimensional matrix = ped(eff,length,3)
     case (g_A_UPG_INB)
        w=(/1., -.5, -.5/)
        ped_len=4
     case (g_As)
        w=(/1., -.5, -.25/)
        res=(/16/11., 4/3., 16/15., 1./)
        ped_len=3
  end select

  if (round == start) call zerom(ainv(eff),nlev(eff))

  do iped=1, lenped(eff)
     p(1:ped_len)=ped(eff,iped,1:ped_len)
     p(4)=ms_stat(iped)                     ! VJ - intergen only read the first three columns and we need a fourth with the number
                                            !      of known parents
     select case (type)
        case (g_A)
            i=1; if(p(2)==0) i=i+1; if (p(3) == 0) i=i+1
            mendel=res(i)
        case (g_A_UPG)
           mendel=res(p(4))
        case (g_A_UPG_INB)
           mendel=p(4)/1000.
        case (g_As)
            i=1; if(p(2)==0) i=i+2; if (p(3) == 0) i=i+1
            mendel=res(i)
     end select

     do i=0,randomnumb(eff) - 1
        do j=0,randomnumb(eff) - 1
           do t1=1,ntrait
              do t2=1,ntrait
                 if (pos_eff(eff+i,t1)/=0 .and. pos_eff(eff+j,t2)/=0) then
                    do k=1,3
                       do l=1,3	
                          if (p(k) /=0 .and. p(l) /=0) then
                             m=address1(eff+i,p(k),t1)
                             n=address1(eff+j,p(l),t2)
                             val=ginv(eff,t1+i*ntrait,t2+j*ntrait)&
                                                   *w(k)*w(l)*mendel
                             call addm(val,m,n,xx)
                          endif
                       enddo
                    enddo
                 endif
              enddo
            enddo
        enddo
     enddo

     ! generate A inverse
     if (round ==start) then
      do k=1,3
         do l=1,3	
            if (p(k) /=0 .and. p(l) /=0) then
               call addm(w(k)*w(l)*mendel,p(k),p(l),ainv(eff))
            endif
         enddo
      enddo
     endif

  enddo

    if (round ==start)   print*,' read ',lenped(eff),' additive pedigrees'
    if (lenped(eff) == 0) then
       print*,'Additive pedigree file for effect ',eff,' empty'
       stop
    endif
  
    if (round ==start) df_random(eff)=lenped(eff)
end subroutine

subroutine add_g_norm_prior(eff)
! VJ subroutine to add contributions on LHS.
! It is assumed that the effects are non-correlated

  integer   :: eff,i,j,k,l,m,n,row,col,t1,t2,type=1,maxrow=0,maxcol=0,maxdiag=0
  integer   :: upper=1 !,lower=2
  real (rh) :: val,x,ident(nlev(eff)),sol(nlev(eff)),xvb,b(nlev(eff)*ntrait)
  type (sparse_hashm) :: a_usr
  type (sparse_ija)   :: a_usr_ija
  
  if (round == start) then
     !call init(a_usr)
     call zerom(a_usr,nlev(eff)*ntrait)
     n=0

     do
        read(io_off+eff,*,iostat=io)t1,col,x,xvb
        if (io /= 0) exit
        n=n+1

        ! Save mean a priori values for normal prior effect
        b(n)=x
        row=col ! Assuming non-correlated normal prior effects !!!

        ! Evaluate if the file is not sorted correctelly
        if (t1 .gt. type) upper=1
        if ((type .gt. t1) .or. (col .ne. upper))then
          print*,''
          print*,'Mixed diagonals in norm_prior file:',t1,col,x
          print*,'Sort levels within traits in ascending order.'
          stop
        endif
        upper=upper+1
        type=t1
        ! Save varinace values in diagonals
        call addm(xvb,n,n,a_usr) 
     enddo

     print*,'norm_prior: read ',n,' elements'
     if (n == 0) then
       print*,'User defined file for effect ',eff,' empty'
       stop
     endif

     call zerom(ainv(eff),nlev(eff)*ntrait)

     ! invert a diagonal matrix
     do i=1,nlev(eff)*ntrait
      call addm(1/getm(i,i,a_usr),i,i,ainv(eff))
     enddo

     if (size(b) <= 12) then
       print*,'Inverted normal prior variance'
       call printm(ainv(eff))
       print*,''
       print*,'Mean prior effects'
       print*, b
     endif
       call reset(a_usr)     
  endif

  ! Calculate vinb*b for normal prior effect and add
  ! vinb contributions in X'RinvX
  j=1
  do t1=1,ntrait
    if (pos_eff(eff,t1)/=0) then
      do row=1,pos_eff(eff,t1)
        x=getm(j,j,ainv(eff))
        m=address1(eff,row,t1)
        call addm(x,m,m,xx)
        vb(m)=x*b(j)
        j=j+1
      enddo
    endif
  enddo
   ! Add contributions in xy due to normal prior effect
   !print*,xy
   xy=xy+vb
   !print*,''
   !print*,xy
end subroutine


subroutine add_g_usr(eff,text)
! generates contributions user defined matrix; if text is 'inv',
! this matrix is inverted 
! Originally code from remlf90 by M.P.E.
!
  integer :: eff,i,j,k,m,n,row,col,t1,t2,type=0,maxrow=0,maxcol=0,maxdiag=0
  integer,parameter::upper=1,lower=2
  real (rh) ::val,x,ident(nlev(eff)),sol(nlev(eff))
  type (sparse_hashm)::a_usr
  type (sparse_ija)::a_usr_ija
  character (*)::text
!
  if (round == start) then
     !call init(a_usr)
     call zerom(a_usr,nlev(eff))
     n=0

     do
        read(io_off+eff,*,iostat=io)row,col,x
        if (io /= 0) exit
        n=n+1

       ! Only upper or lower half stored matrices are accepted
        if (row /=col) then
           if (type == 0) then
               if (row < col) type = upper
               if (row > col) type = lower
           endif

           if ((type == upper .and. row > col) .or. &
               (type == lower .and. row < col)) then
              print*,'Mixed upper and lower diagonals in g_usr_inv:',row,col,x
              stop
           endif
        endif
        
        call addm(x,row,col,a_usr) 
        if (row /=col) call addm(x,col,row,a_usr) 
        maxrow=max(maxrow,row)
        maxcol=max(maxcol,col)
        if (row==col) maxdiag=max(maxdiag,row)
     enddo

     print*,'g_usr_inv: read ',n,' elements'
     if (n == 0) then
       print*,'User defined file for effect ',eff,' empty'
       stop
     endif

     df_random(eff)=nlev(eff)
     print*,'largest row, column,diagonal:',maxrow,maxcol,maxdiag,df_random(eff)
!  call printm(a_usr)
   !call init(ainv(eff))

   call zerom(ainv(eff),nlev(eff))

  if (text(1:3)=='inv') then ! invert a_usr
       a_usr_ija=a_usr
       call reset(a_usr) 
       call init(a_usr)
       call zerom(a_usr,nlev(eff))
       print*,'a_usr_ija'
       !call printm(a_usr_ija)
       do i=1,nlev(eff)
          ident=0; ident(i)=1
          call fspak90('solve',a_usr_ija,ident,sol)
          do j=1,nlev(eff)
	     if (sol(j) /= 0) call addm(sol(j),i,j,ainv(eff))
          enddo
       enddo
!       call fspak90('stat',a_usr_ija)
!       call reset(a_usr_ija)
    else
       ainv(eff)=a_usr 
  endif	   

  if (nlev(eff) <= 20) then
     print*,'user defined matrix'
     call printm(ainv(eff))
  endif

  call reset(a_usr)    	

  endif
  
  do i=0,randomnumb(eff) - 1 !loops for correlated effects
     do j=0,randomnumb(eff) - 1
        do t1=1,ntrait
           do t2=1,ntrait
              if (pos_eff(eff+i,t1)/=0 .and. pos_eff(eff+j,t2)/=0) then	 
	              do k=1,ainv(eff)%nel
                 row=ainv(eff)%x(1,k)
	               if (row /=0) then
                    col=ainv(eff)%x(2,k)
	                   x=ainv(eff)%x(3,k)
                     m=address1(eff+i,row,t1)
                     n=address1(eff+j,col,t2)
                     val=ginv(eff,t1+i*ntrait,t2+j*ntrait)*x 
                     call addm(val,m,n,xx)
     !		           if (m /=n) call addm(val,n,m,xx)
		                if (row /= col) call addm(val,n,m,xx)
		              endif
                 enddo
              endif
            enddo
        enddo
     enddo
   enddo  
  
end subroutine

subroutine add_g_aav(eff,maxnsires)
! FC generates contributions of additive average relationship for ram (App)
! update sires in the App part for gibbs cycle and updates App
  integer :: eff,i,j,t1,t2,k,l,m,n,nsires,iped,iani,sireold
  integer, intent(in) :: maxnsires
  real (rh) :: w(maxnsires+2)
!
  w=-.5; w(1)=1.;    

  if (round <= start+1) call zerom(ainv(eff),nlev(eff)) !use the average in round 1
  if (round == start) df_random(eff)=0

  do iped=1, nlev(eff)-ngrp(eff)
	probs=1.
    if (round == start) df_random(eff)=df_random(eff)+1
    ! Update sires if round > 1. Use average in round 1 so that all possible 
    ! positions in xx are initialized in round 1 -> avoid link reinitialization
    if (ms_stat(iped)>0 .and. round > start) then ! 1 for up & > 1 for twins w/ up 
	  sireold=ped(eff,iped,3)            !keeps the old sire
	  if (ms_alpha/=0) call sample_sire_prob(eff,iped)   !sample sire probabilities from dirichlet distr
	  call update_sire_app(eff,iped)   !updates sire for this cycle (eff=1)
	  ped(eff,iped:iped+ms_stat(iped)-1,3)=sirenew
    end if
    if (ped(eff,iped,3) < 0) then	  ! may use average for some assignments
	  nsires=2+ms(p_ms(iped))
    else
	  nsires=3
    endif
	! generate A inverse
    if (round <= start+1) then
      do k=1,nsires
        if (ped(eff,iped,3) < 0 .and. k >= 3) then	!Uncertain paternity
          m=int(ms((k-2)*2+p_ms(iped)-1))
	      else 
			m=ped(eff,iped,k)
		endif
        do l=1,nsires	
          if (ped(eff,iped,3) < 0 .and. l >= 3) then	!Uncertain paternity
	        probs(l)=ms((l-2)*2+p_ms(iped))
                n=int(ms((l-2)*2+p_ms(iped)-1))
	  else 
		n=ped(eff,iped,l)
          endif
          if (m /= 0 .and. n /= 0) then
            val=w(k)*w(l)*dinv(iped)*probs(k)*probs(l)
            call addm(val,m,n,ainv(eff))
          endif
        enddo
      enddo
	else if (ms_stat(iped)>0 .and. sireold /= sirenew) then !update App-1
	 ! subtract old sire
	 ped(eff,iped:iped+ms_stat(iped)-1,3)=sireold
     do iani=iped, iped+ms_stat(iped)-1  !for monozygous twins - sample sire jointly
      do k=1,3
        do l=1,3	
          if (ped(eff,iani,k) /= 0 .and. ped(eff,iani,l) /= 0) &
            call addm(-1*w(k)*w(l)*dinv(iani),ped(eff,iani,k),ped(eff,iani,l),&
			          ainv(eff))
        enddo
      enddo
	  ! add new sire
      ped(eff,iped:iped+ms_stat(iped)-1,3)=sirenew
      do k=1,3
        do l=1,3	
          if (ped(eff,iani,k) /= 0 .and. ped(eff,iani,l) /= 0) &
            call addm(w(k)*w(l)*dinv(iani),ped(eff,iani,k),ped(eff,iani,l),&
			          ainv(eff))
        enddo
      enddo
	 enddo
    endif

  enddo

 
  !create/update ainv_ija to be added to xx
  if (round <= start+1 .or. max_msires > 1) ainv_ija(eff)=ainv(eff)

  ! add contributions ginv*ainv to xx
  do i=0,randomnumb(eff) - 1
    do j=0,randomnumb(eff) - 1
      do t1=1,ntrait
        do t2=1,ntrait
          if (pos_eff(eff+i,t1)/=0 .and. pos_eff(eff+j,t2)/=0) then
            do iped=1, ainv_ija(eff)%n
  		      m=address1(eff+i,iped,t1)
              do k=ainv_ija(eff)%ia(iped),ainv_ija(eff)%ia(iped+1)-1
				n=address1(eff+j,ainv_ija(eff)%ja(k),t2)
                val=ginv(eff,t1+i*ntrait,t2+j*ntrait)*ainv_ija(eff)%a(k)
                call addm(val,m,n,xx)
                if (iped/=ainv_ija(eff)%ja(k)) call addm(val,n,m,xx)
              enddo
            enddo
          endif
        enddo
      enddo
    enddo
  enddo


  if (round == start)   print*,' read ',df_random(eff),' additive pedigrees'
  if (lenped(eff) == 0) then
     print*,'Additive pedigree file for effect ',eff,' empty'
     stop
  endif

  end subroutine


  subroutine add_g_mb(eff,maxnsires,dimw)
! FC generates contributions of multiple-breed additive or diagonal (co)variance
  integer :: eff,i,j,b1,b2,t1,t2,k,l,m,n,iped,dimw,rank
  integer, intent(in) :: maxnsires
  real (rh) :: w(maxnsires+2),winv(dimw,dimw)
!
  if (randomtype(eff) == g_diag_MB) then 
	w=0.		   !No parental contribution
	else  
      w=-.5
  endif	 
  w(1)=1.; probs=1.;   

  if (round == start) df_random(eff)=0

  do iped=1, lenped(eff)-ngrp(eff)	!id of upg should be > animal id's
    if (round == start) df_random(eff)=df_random(eff)+1
	winv=0
    do b1=1, nbreed(eff)	!Breed variances	  
	  !Individual contribution
      winv=winv+breed(eff,iped,b1)*g(eff,(b1-1)*dimw+1:b1*dimw,(b1-1)*dimw+1:b1*dimw)
	  if (randomtype(eff) /= g_diag_MB) then 
    	do k=2, 3	  	      ! Parents contribution
		  m=ped(eff,iped,k)
  	      if (m > 0 .and. m < ped(eff,iped,1)) & !Second condition for upg
		    winv=winv-.25*breed(eff,m,b1) &
		                 *g(eff,(b1-1)*dimw+1:b1*dimw,(b1-1)*dimw+1:b1*dimw)
 	    enddo
	  endif	
	  if (b1 < nbreed(eff)) then
	    do b2=b1+1, nbreed(eff)  ! Segregation variance between breeds
    	  do k=2, 3	  	         ! Parents contribution
			m=ped(eff,iped,k)
  	        if (m > 0) then
			  winv=winv+2*breed(eff,m,b1)*breed(eff,m,b2) &
 		                 *g(eff,(b1-1)*dimw+1:b1*dimw,(b2-1)*dimw+1:b2*dimw)
	          if (randomtype(eff) /= g_diag_MB .and. m < ped(eff,iped,1)) then 
				do l=2, 3	  	      ! Grandparents contribution
				  n=ped(eff,m,l)
				  if (n > 0) winv=winv &
				          -.5*breed(eff,n,b1)*breed(eff,n,b2) &
 		                     *g(eff,(b1-1)*dimw+1:b1*dimw,(b2-1)*dimw+1:b2*dimw)
				enddo
			  endif
			endif
 	      enddo
		enddo
	  endif
	enddo
    call pos_def(winv,'winv add_g_mb',1d-4)
    call ginv1(winv,dimw,dimw,real(1d-18,rh),rank)

  ! add contributions of ginv to xx
     do i=0,randomnumb(eff) - 1
        do j=0,randomnumb(eff) - 1
           do t1=1,ntrait
              do t2=1,ntrait
                 if (pos_eff(eff+i,t1)/=0 .and. pos_eff(eff+j,t2)/=0) then
                    do k=1,3
                       do l=1,3	
                          if (ped(eff,iped,k) /=0 .and. ped(eff,iped,l) /=0) then
                             m=address1(eff+i,ped(eff,iped,k),t1)
                             n=address1(eff+j,ped(eff,iped,l),t2)
                             val=winv(t1+i*ntrait,t2+j*ntrait)*w(k)*w(l)
                             call addm(val,m,n,xx)
                          endif
                       enddo
                    enddo
                 endif
              enddo
            enddo
        enddo
     enddo
  enddo

  if (round ==start)   print*,' read ',df_random(eff),' additive pedigrees'
  if (lenped(eff) == 0) then
     print*,'Additive pedigree file for effect ',eff,' empty'
     stop
  endif

  end subroutine


  subroutine add_g_mbs(eff,maxnsires,dimw)
! FC generates contributions of multiple-breed additive or diagonal (co)variance
  integer :: eff,i,j,b1,b2,t1,t2,k,l,m,n,o,p,iped,dimw,rank,nsires,ngsires
  integer, intent(in) :: maxnsires
  real (rh) :: w(maxnsires+2),probgs(maxnsires+2),winv(dimw,dimw)
!
  if (randomtype(eff) == g_diag_MB) then 
	w=0.		   !No parental contribution
	else  
      w=-.5
  endif	 
  w(1)=1.;    

  if (round == start) df_random(eff)=0

  do iped=1, lenped(eff)-ngrp(eff)	!id of upg should be > animal id's

	probs=1.
    if (round == start) df_random(eff)=df_random(eff)+1
	winv=0.

    if (ped(eff,iped,3) < 0) then	  ! may use average for some assignments
	  nsires=2+ms(p_ms(iped))
    else 
	  nsires=3
	endif

    do b1=1, nbreed(eff)	!Breed variances	  
	  !Individual contribution
      winv=winv+breed(eff,iped,b1)*g(eff,(b1-1)*dimw+1:b1*dimw,(b1-1)*dimw+1:b1*dimw)
	  if (randomtype(eff) /= g_diag_MB) then 
        do k=2,nsires		  ! Parents contribution
          if (ped(eff,iped,3) < 0 .and. k >= 3) then	!Uncertain paternity using average approach
            m=int(ms((k-2)*2+p_ms(iped)-1))
	        probs(k)=ms((k-2)*2+p_ms(iped))
	        else 
			  m=ped(eff,iped,k)
		  endif
		  !print*,'p contrib.',iped,m,probs(k)
  	      if (m > 0 .and. m < ped(eff,iped,1)) & !Second condition for upg
		    winv=winv-.25*probs(k)*breed(eff,m,b1) &
		                          *g(eff,(b1-1)*dimw+1:b1*dimw,(b1-1)*dimw+1:b1*dimw)
 	    enddo
	  endif	
	  if (b1 < nbreed(eff)) then
	    do b2=b1+1, nbreed(eff)  !Segregation variance between breeds
          do k=2,nsires		  ! Parents contribution
            if (ped(eff,iped,3) < 0 .and. k >= 3) then	!Uncertain paternity using average approach
              m=int(ms((k-2)*2+p_ms(iped)-1))
	          else 
			    m=ped(eff,iped,k)
		    endif
  	        if (m > 0) then
			  winv=winv+2*probs(k)*breed(eff,m,b1)*breed(eff,m,b2) &
 		                     *g(eff,(b1-1)*dimw+1:b1*dimw,(b2-1)*dimw+1:b2*dimw)
	          if (randomtype(eff) /= g_diag_MB .and. m < ped(eff,iped,1)) then
				if (ped(eff,m,3) < 0) then	  ! may use average for some assignments
				  ngsires=2+ms(p_ms(m))
				else 
				  ngsires=3
				endif
			  	probgs=1.
				do l=2, ngsires  	      ! Grandparents contribution
				  if (ped(eff,m,3) < 0 .and. l >= 3) then	!average approach
					n=int(ms((l-2)*2+p_ms(m)-1))
					probgs(l)=ms((l-2)*2+p_ms(m))
					else 
					  n=ped(eff,m,l)
				  endif
				  !print*,'gs contrib.',iped,n,probgs(l),probgs(l)*probs(k)
				  if (n > 0) &
				   winv=winv-.5*probs(k)*probgs(l)*breed(eff,n,b1)*breed(eff,n,b2) &
 		                     *g(eff,(b1-1)*dimw+1:b1*dimw,(b2-1)*dimw+1:b2*dimw)
				enddo
			  endif
			endif
 	      enddo
		enddo
	  endif
	enddo
    call ginv1(winv,dimw,dimw,real(1d-18,rh),rank)


  ! add contributions of ginv to xx
     do i=0,randomnumb(eff) - 1
        do j=0,randomnumb(eff) - 1
           do t1=1,ntrait
              do t2=1,ntrait
                 if (pos_eff(eff+i,t1)/=0 .and. pos_eff(eff+j,t2)/=0) then
					do k=1,nsires
					  if (ped(eff,iped,3) < 0 .and. k >= 3) then	!Uncertain paternity
					  	m=int(ms((k-2)*2+p_ms(iped)-1))
						else 
						  m=ped(eff,iped,k)
					  endif
					  do l=1,nsires	
					    if (ped(eff,iped,3) < 0 .and. l >= 3) then	!Uncertain paternity
						  n=int(ms((l-2)*2+p_ms(iped)-1))
						  else 
						    n=ped(eff,iped,l)
					    endif
					    if (m /= 0 .and. n /= 0) then
                          o=address1(eff+i,m,t1)
                          p=address1(eff+j,n,t2)
                          val=winv(t1+i*ntrait,t2+j*ntrait)*w(k)*w(l)*probs(k)*probs(l)
						  !print*,'xx contrib.',m,n,o,p,val
                          call addm(val,o,p,xx)
                          endif
                      enddo
                    enddo
                 endif
              enddo
            enddo
        enddo
     enddo

  enddo



  if (round ==start)   print*,' read ',df_random(eff),' additive pedigrees'
  if (lenped(eff) == 0) then
     print*,'Additive pedigree file for effect ',eff,' empty'
     stop
  endif

  end subroutine


  subroutine store_solutions
  ! store solutions in file 'solutions'
  integer e,i,j,l,t,nlevr
  open(io_s,file='solutions', status='replace')
  write(io_s,*)'Solution after: ', round,' rounds and burn-in of ',nburn 
  write(io_s,'(''trait effect level     mean            sd            last '')')

  j=0
  do e=1,neff
    do l=1,nlev(e)
       do t=1,ntrait
	     if (round > nburn) then
           if (pos_eff(e,t)/=0) write(io_s,'(2i3,i8,3f15.6)')t,e,l,&
              sumsol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)/((round-nburn)/skip),& 
               dsqrt((sssol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)- &
				sumsol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)**2/ &
				((round-nburn)/skip))/((round-nburn)/skip)),&
				sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)
         else
           if (pos_eff(e,t)/=0) write(io_s,'(2i3,i8,3f15.6)')t,e,l,&
              sumsol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t),& 
               sssol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t),&
				sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)
		 endif
       enddo
    enddo
  
  ! FC store solutions for nonparents in ram
    if (effecttype(e)==effram) then
      j=j+1
      do l=nlev(e)+1,lenped(minloc(effecttype,dim=1,mask=effecttype==effram))
        do t=1,ntrait
	     if (round > nburn) then
           if (pos_eff(e,t)/=0) write(io_s,'(2i3,i8,3f15.6)')t,e,l,&
              sumsol_n(l-nlev(e),t+(j-1)*ntrait)/((round-nburn)/skip),& 
               dsqrt((sssol_n(l-nlev(e),t+(j-1)*ntrait)- &
				sumsol_n(l-nlev(e),t+(j-1)*ntrait)**2/ &
				((round-nburn)/skip))/((round-nburn)/skip)),&
		         sol_n(l-nlev(e),t+(j-1)*ntrait)
		 else
           if (pos_eff(e,t)/=0) write(io_s,'(2i3,i8,3f15.6)')t,e,l,&
              sumsol_n(l-nlev(e),t+(j-1)*ntrait),& 
               sssol_n(l-nlev(e),t+(j-1)*ntrait),&
		        sol_n(l-nlev(e),t+(j-1)*ntrait)
		 endif
        enddo
      enddo
    endif

  enddo
  close(io_s)
  print*,'solutions stored in file: "solutions"'

				  
  !FC Multiple sires posterior probabilities
    if (max_msires>1) then
	  open(io_off+52,file='msiresppr', status='replace')
	  write(io_off+52,*)'Post. probs after: ', round,' rounds and burn-in of ',nburn 
	  write(io_off+52,'(''  animal    sire      post.prob.        sd           last       sample prop'')')
	  do e=1, lenped(1)
		if (ms_stat(e)>0) then
  		  do l=1, int(ms(p_ms(e)))
	       if (round > nburn) then
			write(io_off+52,'(2i8,4f14.6)') e,int(ms(l*2+p_ms(e)-1)),&
									summs(l*2+p_ms(e))/((round-nburn)/skip),&
                          dsqrt((ssms(l*2+p_ms(e))-summs(l*2+p_ms(e))**2/ &
 			               ((round-nburn)/skip))/((round-nburn)/skip)),ms(l*2+p_ms(e)),&
						   summs(l*2+p_ms(e)-1)/((round-nburn)/skip)
		   else
			write(io_off+52,'(2i8,4f14.6)') e,int(ms(l*2+p_ms(e)-1)),&
							summs(l*2+p_ms(e)),ssms(l*2+p_ms(e)),ms(l*2+p_ms(e)),summs(l*2+p_ms(e)-1)
		   endif
		  enddo
		endif
	  enddo
	  close(io_off+52)
	  print*,'msires ppr stored in file: "msiresppr"'
    endif

  
! FC Variance components
  open(io_off+53,file='varcomp', status='replace')
  write(io_off+53,*)'Variance components after: ', round,' rounds and burn-in of ',nburn 

! random effects (G)
  do i=1,neff
    if (randomnumb(i) > 0) then
	  if (randomtype(i)==g_A_MB .or. randomtype(i)==g_diag_MB) then
	    j=randomnumb(i)*ntrait*nbreed(i)
	    else
		  j=randomnumb(i)*ntrait
	  endif
      if (j > 0) then
        write(io_off+53,*)'Average variance components (G) for effect: ', i
		do e=1,j
	     if (round > nburn) then
          write(io_off+53,'(10f12.4)') sumg(i,e,1:j)/((round-nburn)/skip)
		 else
          write(io_off+53,'(10f12.4)') sumg(i,e,1:j)
		 endif
		enddo
        write(io_off+53,*)'Std deviation of variance components (G) for effect: ', i
		do e=1,j
	     if (round > nburn) then
          write(io_off+53,'(10f12.4)') dsqrt((ssg(i,e,1:j)-sumg(i,e,1:j)**2/((round-nburn)/skip))/ &
	                                ((round-nburn)/skip))
		 else
          write(io_off+53,'(10f12.4)') ssg(i,e,1:j)
		 endif
		enddo
        write(io_off+53,*)'Last variance components sample (G) for effect: ', i
		do e=1,j
          write(io_off+53,'(10f12.4)') g(i,e,1:j)
		enddo
     endif   
    endif 
  enddo


! residuals (R)
    
  write(io_off+53,*)'Average residuals (R)'
  do e=1,ntrait
   if (round > nburn) then
    write(io_off+53,'(10f12.4)') sumr(e,:)/((round-nburn)/skip)
   else
    write(io_off+53,'(10f12.4)') sumr(e,:)
   endif
  enddo
  write(io_off+53,*)'Std deviation of residuals (R)'
  do e=1,ntrait
   if (round > nburn) then
    write(io_off+53,'(10f12.4)') dsqrt((ssr(e,:)-sumr(e,:)**2/((round-nburn)/skip))/((round-nburn)/skip))
   else
    write(io_off+53,'(10f12.4)') ssr(e,:)
   endif
  enddo
  write(io_off+53,*)'Last residuals sample (R)'
  do e=1,ntrait
    write(io_off+53,'(10f12.4)') r(e,:)
  enddo
  close(io_off+53)
  print*,'variance components stored in file: "varcomp"'



! FC average likelihood and deviance at average solutions for each observation
  if (round > nburn) then
    open(io_off+54,file='loglike_obs', status='replace')
    write(io_off+54,*)'CPO,Loglike,residuals after: ', round,' rounds and burn-in of ',nburn 
    write(io_off+54,*)'obs       cpo             average logl    logl-thetabar   residuals'
    do i=1,nrec
      write(io_off+54,'(i10,10f18.12)') i, 1/(logl_obs(i,1)/((round-nburn)/skip)),logl_obs(i,2)/((round-nburn)/skip),logl_obs(i,3:3+ntrait)
    enddo
    close(io_off+54)
    print*,'loglike_obs stored in file: "loglike_obs"'
  endif


! FC stores sample of varcomps
  close(io_off+55)
  if (round < nround) open(io_off+55,file='varcompsam', status='old',position='append')

! FC average loglikelihood for each round
  close(io_off+56)
  if (round < nround) open(io_off+56,file='loglike_rnd', status='old',position='append')

! FC samples for two sire assignments
  if (ms_alpha/=0.) then 
    close(io_off+57)
    if (round < nround) open(io_off+57,file='msirespprsam',status='old',position='append')
  endif 

! FC proposal density degrees of belief tuning 
  if (any(randomtype == g_A_MB) .or. any(randomtype == g_diag_MB) .or. residualtype/=r_homo) then
    if(round <= .5*nburn) close(io_off+58)
    if(round < .5*nburn) open(io_off+58,file='mh_dbeliefchg',status='old',position='append')
  endif

  if (any(effectsave == 1)) then 
    close(io_off+61)
    if (round < nround) open(io_off+61,file='solutionsam',status='old',position='append')
  endif

! FC acceptance rates and degrees of belief for proposal density
  if (any(randomtype == g_A_MB) .or. any(randomtype == g_diag_MB) .or. residualtype/=r_homo) then
      open(io_off+59,file='mh_acceptance', status='replace')
	  do i=1,neff
		if (randomtype(i)==g_A_MB .or. randomtype(i)==g_diag_MB) then
		  do j=1,nbreed(i)
			do k=j,nbreed(i)
			  write(io_off+59,*) 'Acceptance rate & proposal dbelief for effect: ', i, 'Breed(s): ',j,k
			  write(io_off+59,'(f12.8,i10)') summb_accep(i,j,k)/((round-nburn)/skip), mb_dbelief(i,j,k)
			enddo
		  enddo
		endif   
	  enddo
      if (neffr>0) then
	   do i=1,neffr
		 e=pos_effr(i)
		 if (effecttype(e)==effcov .or. effecttype(e)==effrnorm) then
			write(io_off+59,*) 'Acceptance rate & proposal dbelief for structural R effect: ', e
			write(io_off+59,'(f12.8,i10)') sumr_accep(i)/((round-nburn)/skip), r_dbeliefprop(i)
		 elseif (effecttype(e)==effcross .and. randomnumb(e) > 0) then 
			write(io_off+59,*) 'Acceptance rate & proposal variance for alpha of structural R effect: ', e
			write(io_off+59,'(2f12.8)') sumr_alphaaccep(i)/((round-nburn)/skip), r_alpha_varprop(i)
         endif
	   enddo
      endif
	  if (residualtype == r_studentt .or. residualtype == r_struct_studentt) then
		write(io_off+59,*) 'Acceptance rate & proposal variance for alpha of robust R '
		write(io_off+59,'(2f12.8)') sumw_alphaaccep/((round-nburn)/skip), w_alpha_varprop
      endif
      close(io_off+59)
      print*,'acceptance rates and degrees of belief stored in file: "mh_acceptance"'
  endif
  !FC
  
! FC structural residual effects in file 'structural_r'
  if (neffr > 0) then
	  open(io_off+60,file='structural_r', status='replace')
	  write(io_off+60,*)'Structural R solution after: ', round,' rounds and burn-in of ',nburn 
	  write(io_off+60,'(''trait effect level     mean            sd            last '')')

	  do j=1,neffr
		e=pos_effr(j)
		if (effecttype(e)==effcov .or. effecttype(e)==effrnorm) then ! set nlev=1 even when nested covariate is present
			nlevr=1
          else
		    nlevr=nlev(e)
        endif
		if (effecttype(e)==effcross .and. randomnumb(e) > 0) then 
		  if (round > nburn) then
		    write(io_off+60,'(2i3,a8,3f15.6)')1,e,'  alpha',&
				  sumr_alpha(j)/((round-nburn)/skip), dsqrt((ssr_alpha(j)-sumr_alpha(j)**2/ &
					((round-nburn)/skip))/((round-nburn)/skip)),r_alpha(j)
		  else
		    write(io_off+60,'(2i3,a8,3f15.6)')1,e,'  alpha',sumr_alpha(j),ssr_alpha(j),r_alpha(j)
		  endif
		endif
		do l=1,nlevr
		   do t=1,ntrait
			 if (round > nburn) then
			   if (pos_eff(e,t)/=0) write(io_off+60,'(2i3,i8,3f15.6)')t,e,l,&
				  sumr_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)/((round-nburn)/skip),& 
				   dsqrt((ssr_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)- &
					sumr_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)**2/ &
					((round-nburn)/skip))/((round-nburn)/skip)),&
					r_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)
			 else
			   if (pos_eff(e,t)/=0) write(io_off+60,'(2i3,i8,3f15.6)')t,e,l,&
				  sumr_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t),& 
				   ssr_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t),&
					r_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)
			 endif
		   enddo
		enddo
	  enddo
	  close(io_off+60)
	  print*,'Structural R  stored in file: "Structural_r"'
  endif


! FC Robustness weights in file 'robustness_w'
  if (residualtype == r_slash .or. residualtype == r_struct_slash .or. &
      residualtype == r_studentt .or. residualtype == r_struct_studentt) then
	  open(io_off+62,file='robustness_w', status='replace')
	  write(io_off+62,*)'Robustness weights after: ', round,' rounds and burn-in of ',nburn 
	  write(io_off+62,'(''trait    record       mean            sd            last '')')

  	  do t=1,ntrait
		if (round > nburn) then
		  write(io_off+62,'(1i3,a10,3f15.6)')t,'    alpha',sumw_alpha/((round-nburn)/skip), &
		  dsqrt((ssw_alpha-sumw_alpha**2/((round-nburn)/skip))/((round-nburn)/skip)),w_alpha
		else
		  write(io_off+62,'(1i3,a10,3f15.6)')t,'    alpha',sumw_alpha,ssw_alpha,w_alpha
		endif
	    do j=1,nrec
		  if (round > nburn) then
		    write(io_off+62,'(1i3,i10,3f15.6)')t,j,sumweight_robust(j)/((round-nburn)/skip),& 
				    dsqrt((ssweight_robust(j)-sumweight_robust(j)**2/ &
					((round-nburn)/skip))/((round-nburn)/skip)),&
					weight_robust(j)
		  else
			write(io_off+62,'(1i3,i10,3f15.6)')t,j,sumweight_robust(j),& 
				    ssweight_robust(j),weight_robust(j)
		  endif
		enddo
	  enddo
	  close(io_off+62)
	  print*,'Robustness weights stored in file: "robustness_w"'
  endif
  end subroutine

  subroutine store_solutions_blup
  ! store solutions in file 'solutions' for blup
  integer e,l,t,j
  open(io_s,file='solutions', status='replace')

if(sol_se == 'se') then
    write(io_s,'(''trait/effect level  solution          s.e.'')')
    if (solv_method /= 'FSPAK') xx_ija=xx
    call fspak90('inverse',xx_ija)
    do e=1,neff
       do l=1,nlev(e)
          do t=1,ntrait
             j=address1(e,l,t)
	     write(io_s,fmt=fmt_sol_se)t,e,l,sol(j),sqrt(getm(j,j,xx_ija))
          enddo
       enddo
    enddo
    print*,'solutions and s.e. stored in file: "solutions"'

 else
    write(io_s,'(''trait/effect level  solution'')')
    do e=1,neff
       do l=1,nlev(e)
          do t=1,ntrait
             write(io_s,fmt=fmt_sol_se)t,e,l,sol(address1(e,l,t))
          enddo
       enddo
    enddo
    print*,'solutions stored in file: "solutions"'
 endif
 end subroutine



  subroutine predict_missing_y
  ! predicts missing traits 
  ! y1 ~ MVN[m1+V12 inv(V22)(y2-m2), MVN(0,V11-V12 inv(V22)V21]
  !                m=Xb+Zu
  ! 
  integer::i,j,k,nmiss,npres,ip(ntrait),im(ntrait),addr
  real (r8)::e(ntrait) 
  real (r8),allocatable,dimension(:,:)::v11,v12,v21,v22,v22inv
  !
  nmiss=0; npres=0
  do i=1,ntrait
     if (y(i) == 0) then
          nmiss=nmiss+1
	  im(nmiss)=i
	else
	  npres=npres+1  
	  ip(npres)=i
     endif
  enddo
   
   if (nmiss == ntrait) then
      print*,'All traits are missing; eliminate such records!'
      stop
   endif
         
   if (nmiss >0) then
     do i=1,ntrait
        e(i)=y(i)
	addr=0
        do j=1,neff
	   addr=addr+1
	   if (effecttype(j) == effuncov) then ! adjust y(theta)
	      e(i)=e(i)-sol(address(addr,i))  !weight_cov not used because x(h) includes sol(rnorm)
	   else ! x(theta) part
	      e(i)=e(i)-weight_cov(addr,i)*sol(address(addr,i))
	   endif
	enddo   
     enddo   
     allocate(v11(nmiss,nmiss),v22(npres,npres),v12(nmiss,npres),&
               v21(npres,nmiss),v22inv(npres,npres))

      v11=r(im(1:nmiss),im(1:nmiss))
      v22=r(ip(1:npres),ip(1:npres))
      v12=r(im(1:nmiss),ip(1:npres))
      v21=r(ip(1:npres),im(1:nmiss))
      ! FC adjust vij for trait specific weights
      if (sum(pos_weight) > 0) then 
        do i=1,nmiss
           do j=1,nmiss
              v11(i,j)=v11(i,j)/sqrt(weight_y(im(i)))/sqrt(weight_y(im(j)))
           enddo
           do j=1,npres
              v12(i,j)=v12(i,j)/sqrt(weight_y(im(i)))/sqrt(weight_y(ip(j)))
              v21(j,i)=v21(j,i)/sqrt(weight_y(ip(j)))/sqrt(weight_y(im(i)))
           enddo
        enddo  
        do i=1,npres
           do j=1,npres
              v22(i,j)=v22(i,j)/sqrt(weight_y(ip(i)))/sqrt(weight_y(ip(j)))
           enddo
        enddo  
      endif

      v22inv=finverse_s(v22)
      
      y(im(1:nmiss))=-e(im(1:nmiss))+matmul(matmul(v12,v22inv),e(ip(1:npres)))						
      y(im(1:nmiss))=gen_normal(real(y(im(1:nmiss)),r8),&
                                 v11-matmul(matmul(v12,v22inv),v21))	 
   endif				
 end subroutine


 subroutine gibbssol
! calculates MCMC estimates for location parameters
  integer :: i,j,k
  real (r8),dimension(ntrait,ntrait):: mydiag
  integer::firsteq,lasteq

  
  denseop_tol=1d-9 !reduce operational zero for dense Cholesky and inverse

  ! Use links for fast conversion from hash to ija formats
  call link_hash_ija(xx,xx_ija)
  ! Have to test if this link really saves cp time for multiple-sire problem
  !
  if (round == 1) then
	!FC use blup as starting values
    call default_iter(conv=1e-10,maxround=400,relax=1.4) !each parameter optional
    ! call solve_iterm(xx,xy,sol)
    call solve_pcg(xx,xy,sol)
    if (neq <= 20) print  '( '' start solutions: '',10(/10f8.2))',sol
    open(io_off+55,file='varcompsam', status='replace') !stores sample of varcomps
    open(io_off+56,file='loglike_rnd', status='replace')	 ! FC average loglikelihood for each round
    if (any(randomtype == g_A_MB) .or. any(randomtype == g_diag_MB) .or. residualtype/=r_homo) &
	    open(io_off+58,file='mh_dbeliefchg', status='replace')
    if (any(effectsave == 1)) open(io_off+61,file='solutionsam', status='replace')
  endif

  ! generate random samples for effects
   do i=1,neq,ntrait
      firsteq=i; lasteq=i+ntrait-1
      call solve_iterm_block(xx_ija,xy,sol,firsteq,lasteq,mydiag,'solve')
      sol(firsteq:lasteq)= &
			gen_normal(sol(firsteq:lasteq),finverse_s(mydiag))
      call solve_iterm_block(xx_ija,xy,sol,firsteq,lasteq,mydiag,'update')		
   enddo
  end subroutine

 subroutine gibbssol_rnorm(t)
! calculates MCMC estimates for location parameters in the reaction norms model
  integer :: i,j,k,l,t
  real (r8),dimension(ntrait,ntrait):: mydiag
  integer::firsteq,lasteq

  denseop_tol=1d-9 !reduce operational zero for dense Cholesky and inverse

  ! Use links for fast conversion from hash to ija formats
  call link_hash_ija(xx,xx_ija)
  !
  if (round == 1) then
    !FC use blup as starting values
    call default_iter(conv=1e-10,maxround=400,relax=1.4) !each parameter optional
    ! call solve_iterm(xx,xy,sol)
    call solve_pcg(xx,xy,sol)
    if (neq <= 20) print  '( '' start solutions: '',10(/10f8.2))',sol
    if (t == 1) then
      open(io_off+55,file='varcompsam', status='replace')  ! stores sample of varcomps
      open(io_off+56,file='loglike_rnd', status='replace') ! average loglikelihood for each round
      if (any(randomtype == g_A_MB) .or. any(randomtype == g_diag_MB) &
        .or. residualtype/=r_homo) open(io_off+58,file='mh_dbeliefchg', status='replace')
      if (any(effectsave == 1)) open(io_off+61,file='solutionsam', status='replace')
    endif
  endif

  ! generate random samples for effects
  do i=1,neff
    if (t == 1 .and. effecttype(i) /= effuncov) then     !C(h) part
      cycle
    elseif (t == 2 .and. effecttype(i) == effuncov) then !C(theta) part
      cycle
    elseif (t == 1 .and. effecttype(i) == effuncov .and. neffr > 0) then !C(h) part for struct model 
      !check if structural variance model is heteroskedastic on h
      l=0 ! uses l as flag
      do j=1,neffr
	k=pos_effr(j)
        if (effecttype(k) == effrnorm .and. pos_eff(k,1) == pos_eff(i,1)) then
          l=1
          exit
        endif 
      enddo
      if (l == 1) then ! use metropolis to sample h
        call sample_h_rnorm(i,k)
        cycle
      endif
    endif 
    do j=1,nlev(i)
      firsteq=sum(nlev(1:i-1))*ntrait+(j-1)*ntrait+1;
      lasteq=firsteq+ntrait-1
      call solve_iterm_block(xx_ija,xy,sol,firsteq,lasteq,mydiag,'solve')
      sol(firsteq:lasteq)=gen_normal(sol(firsteq:lasteq),finverse_s(mydiag))
      call solve_iterm_block(xx_ija,xy,sol,firsteq,lasteq,mydiag,'update')		
    enddo
  enddo
 end subroutine

 subroutine gibbsvar
! calculates MCMC estimates for dispersion parameters
  integer :: g1,g2,i,j,k,l,m,m1,m2,n,n1,n2,t1,t2,addr
  real (r8),dimension(neff,maxcorr,maxcorr)::  q_form,newg
  real (r8),dimension(ntrait,ntrait):: ee,mydiag,newr
  real (r8),dimension(nrec,ntrait)::  e
  integer::firsteq,lasteq,df

! generate random samples for variance components

  !  e'e for residuals (FC parents only in ram)
  ee=0
  do irec=1,nrec
    if (any(randomtype==g_A_MS) .and. any(effecttype==effram)) then	 ! for ram do parents only
      k=minloc(randomtype,dim=1,mask=randomtype==g_A_MS)
      if (int(indata(irec,pos_eff(k,1))) > nlev(k)) cycle
    endif
    call decode_record
    call find_addresses
    call decode_r
    call predict_missing_y

    do k=1,ntrait
      e(irec,k)=y(k)
      do i=1,neff
	    if (effecttype(i) == effuncov) then !weight_cov not used because C(h) includes sol(rnorm)
	      e(irec,k)=e(irec,k)-sol(address(i,k))  
	    else ! C(theta) part
	      e(irec,k)=e(irec,k)-weight_cov(i,k)*sol(address(i,k))
        endif
      enddo
      e(irec,k)=e(irec,k)*sqrt(weight_y(k))
    enddo

   !FC model choice criteria
    if (mod(round,skip)==0) call model_choice(irec,e(irec,:))

    do t1=1,ntrait
      e(irec,t1)=e(irec,t1)/dsqrt(weight_r) 
      do t2=t1,ntrait
        if (orig_r(t1,t2) /= 0) then
          ee(t1,t2)=ee(t1,t2)+e(irec,t1)*e(irec,t2)
          ee(t2,t1)=ee(t1,t2)
        endif   
      enddo                       
    enddo
    
  enddo

! generate random samples for variance components (except residuals)
   ! quadratic forms for random effects except residuals
   q_form=0
   do i=1,neff
     if (randomnumb(i) > 0) then
       if (randomtype(i)==g_A_MB .or. randomtype(i)==g_diag_MB) then
         if (max_msires==0) then
	       call sample_mb_g(i,randomnumb(i)*ntrait)
	     else
	       call sample_mbs_g(i,maxnsires,randomnumb(i)*ntrait)
	     endif
       else
         do j=0,randomnumb(i)-1
           do k=0,randomnumb(i)-1
             do t1=1,ntrait
               do t2=1,ntrait
                 m1=address1(i+j,1,t1); m2=address1(i+j,nlev(i),t1)
                 n1=address1(i+k,1,t2); n2=address1(i+k,nlev(i),t2)
                 g1= t1+j*ntrait; g2= t2+k*ntrait
                 if (orig_g(i,g1,g2) /=0 .and. g1 <= g2) then
                   q_form(i,g1,g2)=quadrf(sol(m1:m2:ntrait),&
                                       ainv(i),sol(n1:n2:ntrait))
                   q_form(i,g2,g1)=q_form(i,g1,g2)
                 endif
               enddo
             enddo
           enddo
         enddo
	! FC generate non parent solutions, ee and qform contributions for ram
         if (randomtype(i)==g_A_MS .and. effecttype(i)==effram) &
	   call mendelian(i,q_form,ee)
  	 j=randomnumb(i)*ntrait
 	 if (j > 0) then
	 ! FC add proper priors
	   df=df_random(i)+dbelief_g(i)
	   q_form(i,:j,:j)=q_form(i,:j,:j)+orig_g(i,:j,:j)*dbelief_g(i)
	   call pos_def(q_form(i,:j,:j),'qform not pos def, fixed',1d-4)
           newg(i,1:j,1:j)=gen_invwishart(finverse_s(q_form(i,1:j,1:j)),df) 
	 endif   
         where (orig_g(i,:,:) == 0)	! don't estimate newg(i,j) if orig_g(i,j)=0
           newg(i,:,:)=0
	     end where
         ! update variance components
         g(i,:,:)=newg(i,:,:)
       endif	     
     endif 
   enddo

!  sampling for residuals
    
! FC add proper priors
  df=nrec+dbelief_r
  ee=ee+orig_r*dbelief_r   

!if(round==start .or. mod(round,skip*10)==0) print*,'round, df, ee ',round , df, ee

!
  call pos_def(ee,'ee not pos def, fixed', 1d-4)
  newr= gen_invwishart(finverse_s(ee),df)  ! New version from IG but need to check first    
  where (orig_r == 0)	! don't estimate newr(i,j) if orig_r(i,j)=0
    newr=0
  end where      

! update variance components
  r=newr

! FC save variance components samples from the beginning of the chain each skip cycles
  if (round == nround .or. mod(round,skip)==0) then
    write(io_off+55,*) round
! random effects (G)
    do i=1,neff
      if (randomnumb(i) > 0) then
	    if (randomtype(i)==g_A_MB .or. randomtype(i)==g_diag_MB) then
	      j=randomnumb(i)*ntrait*nbreed(i)
 	    else
	      j=randomnumb(i)*ntrait
	    endif
	    if (j > 0) then
	      do k=1,j
	        do l=k,j
	        ! do not write samples if g(i,k,l) set to 0
	          if (orig_g(i,k,l)/=0) write(io_off+55,*) g(i,k,l) 
	        enddo
	      enddo
	    endif   
      endif 
    enddo

! residuals (R)
    do i=1,ntrait
      do j=i,ntrait
	if (orig_r(i,j)/=0) write(io_off+55,*) r(i,j)
      enddo
    enddo

! FC average loglikelihood for each round
   write(io_off+56,*) round,' ', logl_rnd

! FC save solutions samples for selected effects
    if (any(effectsave == 1)) then
      write(io_off+61,*) round
      do i=1,neff
	if (effectsave(i) == 1) then
	  do j=1, nlev(i)
	    do k=1, ntrait
	      if (pos_eff(i,k) /= 0) write(io_off+61,*) sol(address1(i,j,k))
	    enddo
	  enddo
	endif   
      enddo
    endif   
  endif   

! accumulate mean and variance of parameters after burning
  if (round > nburn .and. mod(round,skip)==0) then
! residual variance
  	sumr=sumr+newr
	ssr=ssr+newr**2
! genetic variances	
	sumg=sumg+g
	ssg=ssg+g**2
! solutions each cycle and write  mean, total and last
    sumsol=sumsol+sol
    sssol=sssol+sol**2
    if (randomtype(1)==g_A_MS) then
      sumsol_n=sumsol_n+sol_n
      sssol_n=sssol_n+sol_n**2
    endif
  endif	
! FC
 
! FC structural errors
  if (neffr > 0) call structural_r(e) 

! FC robustness
  if (residualtype == r_slash .or. residualtype == r_struct_slash) then
    call slash_weights(e) 
  elseif (residualtype == r_studentt .or. residualtype == r_struct_studentt) then
    call studentt_weights(e) 
  endif

  if(round < nround .and. mod(round,skip*100)==0) call store_solutions
 
  if(mod(round,skip*10)==0) then 
     call cpu_time(c_now)
	 print*,'round ',round , 'elapsed time ',c_now-c_start,' s ',&
	        '(',(skip*10)*60/(c_now-c_last), 'rounds/min )' 
     call cpu_time(c_last)
  endif

  end subroutine
     
  subroutine read_files
  ! FC read msires file, pedigree, dinv for ram
  integer :: i,j,k

  ! ms file; this can be generalize by adding one dim to files and looping over all eff
  if (max_msires > 0) then
  ! Allocate and read msfile
	allocate (ms(n_ms))
	open(io_ms,file=msfile)         ! msires file
	do i=1, n_ms
	  read(io_ms,*,iostat=io) ms(i)
	  if (io.ne.0) then
  		print*, 'Error reading the msires file at record: ', i
		stop
	  end if
	  if (i<10) then
		write (*,100) ms(i)
		100 format (1X, F10.4)
	  endif
	enddo
	print*,'read msfile',n_ms,' records '

	if (max_msires>1) then
	! allocate and initializes vectors to accumulate samples of sires
	  allocate ( summs(n_ms),ssms(n_ms),ms_alphas(n_ms) )
	  summs=0.; ssms=0.
	  if (ms_alpha > 0) then 
		ms_alphas=ms_alpha
	  else 
	        ms_alphas=ms
	  endif
	  if (rstart=='n' .and. ms_alpha/=0.) open(io_off+57,file='msirespprsam', status='replace')
	endif
  endif

  ! pedigree file, allow for more than one pedigree
  do i=1, neff
  
   select case (randomtype(i))

     case (g_fixed,g_diag)
        continue                    ! fixed,diag effect, do nothing
     case (g_A,g_As)		    ! additive 
		do j=1, lenped(i)
	   	    read(i+io_off,*,iostat=io) ped(i,j,:)
			if (io.ne.0) then
			  print*, 'Error reading the pedigree for effect: ',i,'file at animal: ',j
			  stop
			end if
			if (j<11) write (*,103) ped(i,j,:)
			103 format (1X, 3I10)
		end do
		print*,'read pedigree',i,'with ',lenped(i),' animals '

     case (g_A_UPG)		    ! additive `with unknown parent groups
		do j=1, lenped(i)
	   	    read(i+io_off,*,iostat=io) ped(i,j,:), ms_stat(j)
			if (io.ne.0) then
			  print*, 'Error reading the pedigree for effect: ',i,'file at animal: ',j
			  stop
			end if
			if (j<11) write (*,104) ped(i,j,:), ms_stat(j)
			104 format (1X, 4I10)
		end do
		print*,'read pedigree',i,'with ',lenped(i),' animals '

     case (g_A_MS)		            ! FC additive with msires
	allocate (dinv(lenped(i)))          ! Allows only one ms random_group
        if (effecttype(i)==effram) then	    ! allocate extra storage for ram effects
           allocate(sumSOL_N(lenped(i)-nlev(i),randomnumb(i)*ntrait),&
                    ssSOL_N(lenped(i)-nlev(i),randomnumb(i)*ntrait), &
                    sol_n(lenped(i)-nlev(i),randomnumb(i)*ntrait))
           sumSOL_N=0; ssSOL_N=0; sol_n=0
        endif	
	    n_ms=1 ! Reuse var to point positions in ms 
		do j=1, lenped(i)
			read(i+io_off,*,iostat=io) ped(i,j,:) !, dinv(j)        ! VJ - now the program calculates Dinv values
			if (io.ne.0) then
			  print*, 'Error reading the pedigree for effect: ',i,'file at animal: ',j
			  stop
			end if
			if (ped(i,j,3) < 0) then  
			  if (ms(n_ms) <= max_msires) then  
			                                  ! allows to use average relationship when there are too many sires
				ms_stat(j)=1              ! Animal has a multiple sires and use full model
			  endif
			  p_ms(j)=n_ms
			  n_ms=n_ms+2*ms(n_ms)+1
			end if
			if (j<11) then
			  !write (*,101) ped(i,j,:), dinv(j), ms_stat(j), p_ms(j)
                           write (*,101) ped(i,j,:)                               ! VJ - dinv's are compute on the next step
			  101 format (1X, 3I10, F9.4, 2I10)
			endif
		end do 
	 	print*,'read pedigree',i,'with ',lenped(i),' animals '

		! VJ - Calls the subroutine to calculate dinv using Intersires.f90 code
                call cpu_time(c_last)

	        call set_dinv(i)

 		call cpu_time(c_now)
		print*,'Elapsed time to compute dinv:',c_now-c_last,' seconds '
    		print*,""

     case (g_A_MB, g_diag_MB)   ! FC additive with mbreed 
        if (.not. allocated(breed)) then	
		  allocate (breed(neff,maxval(lenped),maxval(nbreed)),&
			        summb_accep(neff,maxval(nbreed),maxval(nbreed)),&
			        mb_dbelief(neff,maxval(nbreed),maxval(nbreed)))
		  summb_accep=0
	    endif
		mb_dbelief(i,:,:)=nlev(i)
	    n_ms=1 ! Reuse var to point positions in ms 
		do j=1, lenped(i)
		    read(i+io_off,*,iostat=io) ped(i,j,:), breed(i,j,:)
			if (io.ne.0) then
			  print*, 'Error reading the pedigree for effect: ',i,'file at animal: ',j
			  stop
			end if
			if (ped(i,j,3) < 0) then  
			  if (ms(n_ms) <= max_msires) then  
			  ! allows to use average relationship when there are too many sires
				ms_stat(j)=1  ! Animal has a multiple sires and use full model
			  endif
			  p_ms(j)=n_ms
			  n_ms=n_ms+2*ms(n_ms)+1
			end if
			if (j<11) write (*,102) ped(i,j,:), ms_stat(j), p_ms(j), breed(i,j,:)
			102 format (1X, 5I10, 10F9.4)
		end do
		!
		print*,'read pedigree',i,'with ',lenped(i),' animals ' 

     case default
       continue
   endselect

  enddo
  
  if (max_msires>0) print*,'nonparent starting position at ms ',n_ms
  
 end subroutine 


  subroutine mendelian(eff,q_form,ee)
! FC generates mendelian sampling contribution and non-parents solutions 
! for ram - the animal must be the fisrt effect
  integer :: eff,i,j,t1,t2,g1,g2,k,m,nsires,iped,i_ms
  real (r8),intent(inout):: q_form(neff,maxcorr,maxcorr),ee(ntrait,ntrait)
!
  real (rh) :: e(ntrait), varm
!
  probs=1.; sol_n=0; irec=0

  do iped=nlev(eff)+1, lenped(eff)

    probs=1;
  ! the record residual
	do 
	   irec=irec+1
	   if (irec > nrec) exit
	   if (iped==int(indata(irec,pos_eff(eff,1)))) then
    	      call decode_record
	      call find_addresses
   	      call decode_r
	      call predict_missing_y
	      e=y
	      do i=1,naddr
		 do j=1,ntrait
	            if (y(j) /= miss) &
		       e(j)=e(j)-weight_cov(i,j)*sol(address(i,j))
	            if (i == naddr) e(j)=e(j)*sqrt(weight_y(j))
		 enddo     
	      enddo
	      exit
	   elseif (iped<int(indata(irec,pos_eff(eff,1)))) then
	     irec=irec-1
	     exit
	   else
	     cycle
	   endif
	enddo


  ! calculate the mendelian sampling
    do i=0,randomnumb(eff) - 1
      do t1=1,ntrait
        if (e(t1)/=0 .and. pos_eff(eff+i,t1)/=0) then
		  if (mod(i,2)==0) then !Additive effect
		    varm=(1/((1/r(t1,t1))+dinv(iped)/g(eff,t1+i*ntrait,t1+i*ntrait)))
		    sol_n(iped-nlev(eff),t1+i*ntrait)= &
			  gen_normal(varm*(1/r(t1,t1))*e(t1),varm)
			e(t1)=e(t1)-sol_n(iped-nlev(eff),t1+i*ntrait)
		  else	! Maternal effect
			sol_n(iped-nlev(eff),t1+i*ntrait) = &
			  gen_normal(-(ginv(eff,t1+(i-1)*ntrait,t1+i*ntrait)/&
			  ginv(eff,t1+i*ntrait,t1+i*ntrait))*&
			  sol_n(iped-nlev(eff),t1+(i-1)*ntrait),&
			  (1/(ginv(eff,t1+i*ntrait,t1+i*ntrait)*dinv(iped))))
		  endif 
        endif
      enddo
    enddo

  ! use mendelian sampling to update ee q_form

  ! e'e for residuals of non-parents in ram
    
    do t1=1,ntrait
      do t2=t1,ntrait
!        if (orig_r(t1,t2) /= 0) then
          ee(t1,t2)=ee(t1,t2)+e(t1)*e(t2)
          ee(t2,t1)=ee(t1,t2)
!	    endif   
      enddo                       
    enddo
    
  !FC model choice criteria
    if (round > nburn .and. mod(round,skip)==0) call model_choice(irec,e)
  
  ! quadratic forms for random effects nonparents

    if (randomnumb(eff) > 0) then
      do i=0,randomnumb(eff)-1
        do j=i,randomnumb(eff)-1
          do t1=1,ntrait
            do t2=t1,ntrait
              g1= t1+i*ntrait; g2= t2+j*ntrait
!             if (orig_g(i,g1,g2) /=0 .and. g1 <= g2) then
              if (g1 <= g2) then
		        q_form(eff,g1,g2)=q_form(eff,g1,g2)+ &
				                  sol_n(iped-nlev(eff),g1)* &
				                  sol_n(iped-nlev(eff),g2)*dinv(iped)
		        q_form(eff,g2,g1)=q_form(eff,g1,g2)	
		      endif
            enddo                       
          enddo
        enddo
      enddo
	endif
    if (round == start) df_random(eff)=df_random(eff)+1



  ! add the parent average
    if (ped(eff,iped,3) < 0) then
	  nsires=2+ms(p_ms(iped))
    else 
      nsires=3
 	endif
    do i=0,randomnumb(eff) - 1
      do t1=1,ntrait
        if (pos_eff(eff+i,t1)/=0) then
          do k=2,nsires
            if (ped(eff,iped,3) < 0 .and. k >= 3) then	!Uncertain paternity
	          probs(k)=ms((k-2)*2+p_ms(iped))
              m=address1(eff+i,int(ms((k-2)*2+p_ms(iped)-1)),t1)
	          else 
				m=address1(eff+i,ped(eff,iped,k),t1)
			endif
            if (m /= address1(eff+i,0,t1)) then
			  sol_n(iped-nlev(eff),t1+i*ntrait) = &
			    sol_n(iped-nlev(eff),t1+i*ntrait) + .5*probs(k)*sol(m)
            endif
          enddo
        endif
      enddo
    enddo
	! update sire for next round
	if (ms_stat(iped)==1) then
	! add the previous sire solution and additive mendelian sampl to the residual
	  do t1=1,ntrait !not general for multiple trait yet
		e(t1)=e(t1)+sol_n(iped-nlev(eff),t1) !add non-parent solution 
		if (ped(eff,iped,2)>0) &  !subtract the dam back
		  e(t1)=e(t1)-.5*sol(address1(eff,ped(eff,iped,2),t1))
	  enddo
	  if (ms_alpha/=0) call sample_sire_prob(eff,iped)   !sample sire probabilities from dirichlet distr
      call update_sire_coeff(eff,iped,e) !updates sire for this cycle (eff=1)
	  ped(eff,iped,3)=sirenew
	end if

  enddo

  if (round ==start) &
    print*,' calculated ',(lenped(eff)-nlev(eff)),' mendelian samples and ', &
         'read ',irec,' records'
  
  end subroutine

  subroutine update_sire_app(eff,iped)

! Calculate probs for each possible sire of parent individuals (App)
! in the current cycle and updates the sire based on these probabilities
! For single trait ram with maternal effect
  integer :: eff,iped,i,j,k,t1,t2
  real (rh) :: mendeli,mendelj,scale,qform
!
  scale=0.;probs=0.

  do k=1, int(ms(p_ms(iped)))
    qform=0.
  ! calculate the mendelian sampling for each sire j
	do i=0,randomnumb(eff) - 1
	  do j=0,randomnumb(eff) - 1
		do t1=1,ntrait
		  do t2=1,ntrait
			mendeli=sol(address1(eff+i,iped,t1)) &
			   -.5*(sol(address1(eff+i,ped(eff,iped,2),t1))+ &
					sol(address1(eff+i,int(ms(k*2+p_ms(iped)-1)),t1)))
			mendelj=sol(address1(eff+j,iped,t2)) &
			   -.5*(sol(address1(eff+j,ped(eff,iped,2),t2))+ &
					sol(address1(eff+j,int(ms(k*2+p_ms(iped)-1)),t2)))
 			qform=qform+mendeli*mendelj*ginv(eff,t1+i*ntrait,t2+j*ntrait)
		  enddo
		enddo
	  enddo
	enddo
	probs(2+k)=probs(2+k)+ms(k*2+p_ms(iped))*dexp(-.5*dinv(iped)*qform)
  enddo
  
  scale=sum(probs) 

  do k=1, int(ms(p_ms(iped)))
    probs(2+k)=probs(2+k)/scale
  enddo

  call sample_sire(iped) ! Sampling a starting sire from the probs

  probs=1.   ! Reset probs vector
  
  end subroutine



  subroutine update_sire_coeff(eff,iped,e)
! Calculate probs of each possible msire of non-parent individuals (Z*)
! in the current cycle and updates the sire based on these probabilities
! For single trait ram with maternal effect and all non parents with record
  integer :: eff,iped,k,t1
  real (rh) :: e(ntrait),scale
!
  scale=0.;probs=0.

  do k=1, int(ms(p_ms(iped)))
  ! calculate the deviation for each sire k
    do t1=1,ntrait
      probs(2+k)=probs(2+k)+(rinv(t1,t1)*(e(t1) &
	             -.5*sol(address1(eff,int(ms(k*2+p_ms(iped)-1)),t1)))**2)
    enddo
    probs(2+k)=ms(k*2+p_ms(iped))*dexp(-.5*probs(2+k))
    scale=scale+probs(2+k)
  enddo

  do k=1, int(ms(p_ms(iped)))
    probs(2+k)=probs(2+k)/scale
  enddo

  call sample_sire(iped) ! Sampling a starting sire from the probs

  probs=1.   ! Reset probs vector

  end subroutine


  subroutine update_sire_mb(eff,iped,winv)
! Calculate probs for each possible sire of individuals 
! in the current cycle and updates the sire based on these probabilities
  integer :: eff,iped,i,j,k,t1,t2
  real (rh), intent(in) :: winv(randomnumb(eff)*ntrait,randomnumb(eff)*ntrait)
  real (rh) :: mendeli,mendelj,scale
!
  scale=0.;probs=0.


  do k=1, int(ms(p_ms(iped)))
  ! calculate the mendelian sampling for each sire j
    do i=0,randomnumb(eff) - 1
      do j=0,randomnumb(eff) - 1
        do t1=1,ntrait
		  do t2=1,ntrait
            mendeli=sol(address1(eff+i,iped,t1)) &
		       -.5*(sol(address1(eff+i,ped(eff,iped,2),t1))+ &
			        sol(address1(eff+i,int(ms(k*2+p_ms(iped)-1)),t1)))
            mendelj=sol(address1(eff+j,iped,t2)) &
		       -.5*(sol(address1(eff+j,ped(eff,iped,2),t2))+ &
			        sol(address1(eff+j,int(ms(k*2+p_ms(iped)-1)),t2)))
 	        probs(2+k)=probs(2+k)+mendeli*mendelj* &
		  			   winv(t1+i*ntrait,t2+j*ntrait)
		  enddo
		enddo
      enddo
    enddo
    probs(2+k)=ms(k*2+p_ms(iped))*dexp(-.5*probs(2+k))
    scale=scale+probs(2+k)
  enddo


  do k=1, int(ms(p_ms(iped)))
    probs(2+k)=probs(2+k)/scale
  enddo

  print*,iped,probs

  call sample_sire(iped) ! Sampling a starting sire from the probs

  probs=1.   ! Reset probs vector
  
  end subroutine


  subroutine sample_sire(iped)
! Sample the parent for this round using prob vector
  integer :: i,iped
  real (rh) :: cumprob,test ! cummulative probabilities of msires
!
  cumprob=0.; sirenew=0 ! sire assignement for this cycle
  test=gen_uniform() ! random uniform(0,1) number
  i=2	!position in the vector prob
  do while (sirenew==0)
    i=i+1
    cumprob=cumprob+probs(i)
	if (test<=cumprob) then 
	  sirenew=int(ms((i-2)*2+p_ms(iped)-1))
	  if (round > nburn .and. mod(round,skip)==0) then 
	    summs((i-2)*2+p_ms(iped)-1)=summs((i-2)*2+p_ms(iped)-1)+1 !accumulate sampled sire
      endif
	endif
  end do

  end subroutine


  subroutine sample_sire_prob(eff,iped)
! Sample sire assignment probabilities from a Dirichlet distr.
! uses alphas values as hyperparameters; update new values in ms;
  use ranlib
  integer :: eff,iped,k
  real (rh) :: scale
  real :: alphanew(int(ms(p_ms(iped))))
!
  scale=0.;alphanew=0.


  do k=1, int(ms(p_ms(iped)))
    alphanew(k)=ms_alphas(k*2+p_ms(iped))
    if (int(ms(k*2+p_ms(iped)-1))==ped(eff,iped,3) .or. round<=start+1) alphanew(k)=alphanew(k)+1
	! Add 1 to the alpha of last sampled sire & 
	! to all in rounds 1&2 (avoid floating error when alpha is small)
    alphanew(k)=sgamma(alphanew(k))
    scale=scale+alphanew(k)
  enddo

  do k=1, int(ms(p_ms(iped)))
    ms(k*2+p_ms(iped))=alphanew(k)/scale
	if (round > nburn .and. mod(round,skip)==0)  then
	  summs(k*2+p_ms(iped))=summs(k*2+p_ms(iped))+ms(k*2+p_ms(iped)) !accumulate sampled sire
	  ssms(k*2+p_ms(iped))=ssms(k*2+p_ms(iped))+ms(k*2+p_ms(iped))**2 !accumulate sampled sire
      if (p_ms(iped)==1) &
          write(io_off+57,*) round,iped,int(ms(k*2+p_ms(iped)-1)),ms(k*2+p_ms(iped))
	endif
  enddo


  end subroutine


  subroutine model_choice(irec,e)
! FC accumulates likelihood harmonic mean & loglikelihood for each observation &  
! the loglikelihood for each round
  integer :: irec,t1,t2
  real (r8) :: e(ntrait),logl
!
  logl=0
  do t1=1,ntrait
	do t2=1,ntrait
!	  if (orig_r(t1,t2) /= 0) then
  		logl=logl+e(t1)*e(t2)*rinv(t1,t2)
!	  endif   
	enddo                       
  enddo
  logl=-.5*ntrait*dlog(2*pi)+.5*dlog(fdet_s(rinv))-.5*logl
  logl_obs(irec,1)=logl_obs(irec,1)+1/dexp(logl)
  logl_obs(irec,2)=logl_obs(irec,2)+logl
  logl_rnd=logl_rnd+logl
 
  end subroutine

  subroutine deviance
! FC calculates log likelihood for each observation at posterior mean of parameters 
  integer :: i,t1,t2, addr
  real (r8) :: e(ntrait),logl
!
!FC points to nonparent starting position for ram(average)

  do irec=1,nrec
    logl=0;	rinv=sumr/((nround-nburn)/skip); weight_r=1.

    call decode_record
    call find_addresses

  ! FC structural variances
	if (neffr > 0) then
      do k=1,ntrait
		do i=1,neffr
		  if (nestedcov(pos_effr(i),k) > 0) then !allows reacnorm & rr heteroskedasticity base on covariate value
			  addr=address1(pos_effr(i),1,k)
			else
			  addr=address(pos_effr(i),k)
		  endif
		  weight_r=weight_r*((sumr_sol(addr)/((nround-nburn)/skip))&
		                       **weight_cov(pos_effr(i),k))
		enddo
      enddo
	endif
  ! FC robustness
    if (residualtype == r_slash .or. residualtype == r_struct_slash .or. &
        residualtype == r_studentt .or. residualtype == r_struct_studentt) then
      weight_r=weight_r/(sumweight_robust(irec)/((nround-nburn)/skip))
    endif
	rinv=finverse_s(rinv*weight_r)

    call predict_missing_y

    do t1=1,ntrait
      e(t1)=y(t1)
      do i=1,neff
	  ! CHECK THAT - SEEMS TO BE SUBTRACTING THE PARENT AVERAGE TO
        if (randomtype(i)==g_A_MS .and. int(indata(irec,pos_eff(i,1))) > nlev(i)) then	
		!subtract the animal effect instead of parents average for non-parents
		  weight_cov(i:(i+naddr-neff),t1)=0 !THIS COULD BE GENERALIZED
 		  e(t1)=e(t1)-sumsol_n(int(indata(irec,pos_eff(i,1)))-nlev(i),t1+(i-1)*ntrait)/&
																	((nround-nburn)/skip)
	endif
      enddo    
      addr=0
      do i=1,neff
	 addr=addr+1
	 if (effecttype(i) == effuncov) then !weight_cov not used because x(h) includes sol(rnorm)
	    e(t1)=e(t1)-sol(address(addr,t1))  
	 else ! x(theta) part
	    e(t1)=e(t1)-weight_cov(addr,t1)*(sumsol(address(addr,t1))/((nround-nburn)/skip))
	 endif
      enddo
    enddo

    do t1=1,ntrait
       logl_obs(irec,3+t1)=e(t1)*dsqrt(rinv(t1,t1))*sqrt(weight_y(t1))
       do t2=1,ntrait
!	   if (orig_r(t1,t2) /= 0) then
  	     logl=logl+e(t1)*e(t2)*rinv(t1,t2)*sqrt(weight_y(t1))*sqrt(weight_y(t2))
!	   endif   
       enddo                       
    enddo

    logl_obs(irec,3)=-.5*ntrait*dlog(2*pi)+.5*dlog(fdet_s(rinv))-.5*logl
             
  enddo
 
  end subroutine


  subroutine sample_mb_g(eff,dimw)
! FC sample multiple-breed and segregation variance via Metropolis Hasting alghorithm
  integer :: eff,i,j,b1,b2,t1,t2,k,l,m,n,iped,dimw,rank,t
!  integer, intent(in) :: maxnsires
  real (rh),dimension(neff,maxcorr,maxcorr):: newg
  real (rh) :: winv(dimw,dimw),winv_o(dimw,dimw),mendel(dimw), &
               logdet,logdet_o,qform,alpha,temp

  do t=1, mb_nround

   newg(eff,:,:)=g(eff,:,:)
	
   do i=1, nbreed(eff)	!Breed and segregation variances	  
    do j=i, nbreed(eff)
	  ! generate a candidate value with E(newg)=g
	  newg(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw) &
	      =gen_invwishart(finverse_s(g(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw) &
		                  *(mb_dbelief(eff,i,j)-dimw-1)),mb_dbelief(eff,i,j))
! 	  where (orig_g(eff,:,:) == 0)	! don't estimate newg(i,j) if orig_g(i,j)=0
!	    newg(eff,:,:)=0
!	  end where
	  ! calculate acceptance rates
      logdet=0.; logdet_o=0.; qform=0. ! Reset values 
	  do iped=1, lenped(eff)-ngrp(eff)
		winv=0.; winv_o=0.; mendel=0. ! Reset values 
        ! The mendelian sampling 
		do k=0,randomnumb(eff) - 1
		  do t1=1,ntrait
			mendel(t1+k*ntrait)=sol(address1(eff+k,iped,t1))
		    if (randomtype(eff) /= g_diag_MB) then 
    		  do l=2, 3	  	      ! Parents contribution
			    m=ped(eff,iped,l)
  			    if (m > 0) mendel(t1+k*ntrait)=mendel(t1+k*ntrait) &
				                       -.5*sol(address1(eff+k,m,t1))
			  enddo
			endif
		  enddo
		enddo
		do b1=1, nbreed(eff)	!Breed variances	  
		  !Individual contribution
		  winv=winv &
		      +breed(eff,iped,b1)*newg(eff,(b1-1)*dimw+1:b1*dimw,(b1-1)*dimw+1:b1*dimw)
		  winv_o=winv_o &
		        +breed(eff,iped,b1)*g(eff,(b1-1)*dimw+1:b1*dimw,(b1-1)*dimw+1:b1*dimw)
		  if (randomtype(eff) /= g_diag_MB) then 
    		do k=2, 3	  	      ! Parents contribution
			  m=ped(eff,iped,k)
  			  if (m > 0 .and. m < ped(eff,iped,1)) then
			    winv=winv-.25*breed(eff,m,b1) &
							 *newg(eff,(b1-1)*dimw+1:b1*dimw,(b1-1)*dimw+1:b1*dimw)
			    winv_o=winv_o-.25*breed(eff,m,b1) &
							   *g(eff,(b1-1)*dimw+1:b1*dimw,(b1-1)*dimw+1:b1*dimw)
			  endif
 			enddo
		  endif	
		  if (b1 < nbreed(eff)) then
			do b2=b1+1, nbreed(eff)  !Segregation variance between breeds
    		  do k=2, 3	  	      ! Parents contribution
			    m=ped(eff,iped,k)
  				if (m > 0) then
				  winv=winv+2*breed(eff,m,b1)*breed(eff,m,b2) &
 							 *newg(eff,(b1-1)*dimw+1:b1*dimw,(b2-1)*dimw+1:b2*dimw)
				  winv_o=winv_o+2*breed(eff,m,b1)*breed(eff,m,b2) &
 								 *g(eff,(b1-1)*dimw+1:b1*dimw,(b2-1)*dimw+1:b2*dimw)
				  if (randomtype(eff) /= g_diag_MB .and. m < ped(eff,iped,1)) then 
					do l=2, 3	  	      ! Grandparents contribution
					  n=ped(eff,m,l)
					  if (n > 0) then
					    winv=winv-.5*breed(eff,n,b1)*breed(eff,n,b2) &
 								 *newg(eff,(b1-1)*dimw+1:b1*dimw,(b2-1)*dimw+1:b2*dimw)
					    winv_o=winv_o-.5*breed(eff,n,b1)*breed(eff,n,b2) &
 								 *g(eff,(b1-1)*dimw+1:b1*dimw,(b2-1)*dimw+1:b2*dimw)
					  endif
					enddo
				  endif
				endif
 			  enddo
			enddo
		  endif
		enddo

        call ginv1(winv,dimw,dimw,real(1d-18,rh),rank)
        call ginv1(winv_o,dimw,dimw,real(1d-18,rh),rank)
		logdet=logdet+dlog(fdet_s(winv)); logdet_o=logdet_o+dlog(fdet_s(winv_o))

        ! quadratic form 
		do k=0,randomnumb(eff) - 1
		  do l=0,randomnumb(eff) - 1
	 		do t1=1,ntrait
	 		  do t2=1,ntrait
			    qform=qform+mendel(t1+k*ntrait)*mendel(t2+l*ntrait) &
				       *(winv_o(t1+k*ntrait,t2+l*ntrait)-winv(t1+k*ntrait,t2+l*ntrait))
			  enddo
			enddo
		  enddo
		enddo

	  enddo


	! acceptance probability
	  alpha = .5*logdet-.5*logdet_o+.5*qform+.5*(dbelief_g(eff)-2*mb_dbelief(eff,i,j)) &
	         *(dlog(fdet_s(g(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw))) &
		    -dlog(fdet_s(newg(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw)))) &
           +.5*sum(diag( &
	       matmul( &
	       (finverse_s(g(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw)) &
	       -finverse_s(newg(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw))) &
	       ,orig_g(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw) &
		   ) &
		   *(dbelief_g(eff)-dimw-1) &
		   )) & 
           +.5*(mb_dbelief(eff,i,j)-dimw-1)*sum(diag( &
	       matmul( &
		   finverse_s(newg(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw)) &
		   ,g(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw) &
		   ) &
	       -matmul( &
           finverse_s(g(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw)) &
		   ,newg(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw) &
		   ) &
		   )) 

      if (alpha<0) then
        alpha = dexp(alpha)
		else
		alpha=1
	  endif
	  temp=gen_uniform()
  	  if(temp < alpha) then
	    ! keep symmetry
	    if (i /= j) newg(eff,(j-1)*dimw+1:j*dimw,(i-1)*dimw+1:i*dimw) &
		           =newg(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw)
	    g(eff,:,:)=newg(eff,:,:)
		else
	      newg(eff,:,:)=g(eff,:,:)
  	  endif

	  ! Acceptance rates
      if(t == mb_nround) then
       if(round <= .5*nburn) then
	    summb_accep(eff,i,j)=summb_accep(eff,i,j)+alpha
		!tune MH alghorithm half-way through burnin
	    if(mod(round,mb_skip)==0) then
	      summb_accep(eff,i,j)=summb_accep(eff,i,j)/mb_skip
          if (summb_accep(eff,i,j) > .6) then 
	    	mb_dbelief(eff,i,j)=mb_dbelief(eff,i,j)*.9   
            elseif (summb_accep(eff,i,j) < .2) then
 	    	  mb_dbelief(eff,i,j)=mb_dbelief(eff,i,j)*1.05+1
		  endif
		  if ((mb_dbelief(eff,i,j)-dimw-1) < 4) mb_dbelief(eff,i,j) = dimw + 4 ! lower bound  	   
	      summb_accep(eff,i,j)=0.
		  write(io_off+58,*) mb_dbelief(eff,i,j)
        endif
	   elseif (round > nburn .and. mod(round,skip)==0) then
	    summb_accep(eff,i,j)=summb_accep(eff,i,j)+alpha
	   endif
	  endif

	enddo
   enddo

  enddo

  end subroutine


  subroutine sample_mbs_g(eff,maxnsires,dimw)
! FC sample multiple-breed and segregation variance via Metropolis Hasting alghorithm
! in the case of uncertain paternity
  integer :: eff,i,j,b1,b2,t1,t2,k,l,m,n,t,iped, dimw, rank,nsires,ngsires,sireold
  integer, intent(in) :: maxnsires
  real (rh),dimension(neff,maxcorr,maxcorr):: newg
  real (rh) :: winv(dimw,dimw),winv_o(dimw,dimw),mendel(dimw),logdet,logdet_o,qform &
			  ,probgs(maxnsires+2),alpha, temp

  do t=1, mb_nround

   newg(eff,:,:)=g(eff,:,:)
	
   do i=1, nbreed(eff)	!Breed and segregation variances	  
    do j=i, nbreed(eff)
	  ! generate a candidate value with E(newg)=g
	  newg(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw) &
	      =gen_invwishart(finverse_s(g(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw) &
		                  *(mb_dbelief(eff,i,j)-dimw-1)),mb_dbelief(eff,i,j))
	  ! calculate acceptance rates
      logdet=0.; logdet_o=0.; qform=0. ! Reset values 
	  do iped=1, lenped(eff)-ngrp(eff)
		winv=0.; winv_o=0.; mendel=0. ! Reset values 

		probs=1.
		if (ped(eff,iped,3) < 0) then	  ! may use average for some assignments
		  nsires=2+ms(p_ms(iped))
		else 
		  nsires=3
		endif

        ! The mendelian sampling 
		do k=0,randomnumb(eff) - 1
		  do t1=1,ntrait
			mendel(t1+k*ntrait)=sol(address1(eff+k,iped,t1))
		    if (randomtype(eff) /= g_diag_MB) then 
    		  do l=2, nsires	  	      ! Parents contribution
			    if (ped(eff,iped,3) < 0 .and. l >= 3) then	!Uncertain paternity using average approach
				  m=int(ms((l-2)*2+p_ms(iped)-1))
				  probs(l)=ms((l-2)*2+p_ms(iped))
				  else 
				    m=ped(eff,iped,l)
			    endif
			    !print*,'p contrib.',iped,m,probs(l)
  			    if (m > 0) mendel(t1+k*ntrait)=mendel(t1+k*ntrait) &
				                -.5*probs(l)*sol(address1(eff+k,m,t1))
			  enddo
			endif
		  enddo
		enddo


		do b1=1, nbreed(eff)	!Breed variances	  
		  !Individual contribution
		  winv=winv+breed(eff,iped,b1)*newg(eff,(b1-1)*dimw+1:b1*dimw,(b1-1)*dimw+1:b1*dimw)
		  winv_o=winv_o+breed(eff,iped,b1)*g(eff,(b1-1)*dimw+1:b1*dimw,(b1-1)*dimw+1:b1*dimw)
		  if (randomtype(eff) /= g_diag_MB) then 
            do k=2,nsires		  ! Parents contribution
              if (ped(eff,iped,3) < 0 .and. k >= 3) then	!Uncertain paternity using average approach
                m=int(ms((k-2)*2+p_ms(iped)-1))
	            else 
			      m=ped(eff,iped,k)
		      endif
  	          if (m > 0 .and. m < ped(eff,iped,1)) then
			    winv=winv-.25*probs(k)*breed(eff,m,b1) &
							 *newg(eff,(b1-1)*dimw+1:b1*dimw,(b1-1)*dimw+1:b1*dimw)
			    winv_o=winv_o-.25*probs(k)*breed(eff,m,b1) &
							   *g(eff,(b1-1)*dimw+1:b1*dimw,(b1-1)*dimw+1:b1*dimw)
			  endif
 			enddo
		  endif	
		  if (b1 < nbreed(eff)) then
			do b2=b1+1, nbreed(eff)  !Segregation variance between breeds
              do k=2,nsires		  ! Parents contribution
                if (ped(eff,iped,3) < 0 .and. k >= 3) then	!Uncertain paternity using average approach
                  m=int(ms((k-2)*2+p_ms(iped)-1))
	              else 
			        m=ped(eff,iped,k)
		        endif
  	            if (m > 0) then
				  winv=winv+2*probs(k)*breed(eff,m,b1)*breed(eff,m,b2) &
 							 *newg(eff,(b1-1)*dimw+1:b1*dimw,(b2-1)*dimw+1:b2*dimw)
				  winv_o=winv_o+2*probs(k)*breed(eff,m,b1)*breed(eff,m,b2) &
 								 *g(eff,(b1-1)*dimw+1:b1*dimw,(b2-1)*dimw+1:b2*dimw)
				  if (randomtype(eff) /= g_diag_MB .and. m < ped(eff,iped,1)) then
					if (ped(eff,m,3) < 0) then	  ! may use average for some assignments
					  ngsires=2+ms(p_ms(m))
					else 
					  ngsires=3
					endif
			  		probgs=1.
					do l=2, ngsires  	      ! Grandparents contribution
					  if (ped(eff,m,3) < 0 .and. l >= 3) then	!average approach
						n=int(ms((l-2)*2+p_ms(m)-1))
						probgs(l)=ms((l-2)*2+p_ms(m))
						else 
						  n=ped(eff,m,l)
					  endif
					  !print*,'gs contrib.',iped,n,probgs(l)
					  if (n > 0) then
					    winv=winv-.5*probs(k)*probgs(l)*breed(eff,n,b1)*breed(eff,n,b2) &
 								 *newg(eff,(b1-1)*dimw+1:b1*dimw,(b2-1)*dimw+1:b2*dimw)
					    winv_o=winv_o-.5*probs(k)*probgs(l)*breed(eff,n,b1)*breed(eff,n,b2) &
 								 *g(eff,(b1-1)*dimw+1:b1*dimw,(b2-1)*dimw+1:b2*dimw)
					  endif
					enddo
				  endif
				endif
 			  enddo
			enddo
		  endif
		enddo

        call ginv1(winv,dimw,dimw,real(1d-18,rh),rank)
        call ginv1(winv_o,dimw,dimw,real(1d-18,rh),rank)
		logdet=logdet+dlog(fdet_s(winv)); logdet_o=logdet_o+dlog(fdet_s(winv_o))

        ! quadratic form 
		do k=0,randomnumb(eff) - 1
		  do l=0,randomnumb(eff) - 1
	 		do t1=1,ntrait
	 		  do t2=1,ntrait
			    qform=qform+mendel(t1+k*ntrait)*mendel(t2+l*ntrait) &
				       *(winv_o(t1+k*ntrait,t2+l*ntrait)-winv(t1+k*ntrait,t2+l*ntrait))
			  enddo
			enddo
		  enddo
		enddo

		! sample sire probabilities and assignments
		if (ms_stat(iped)==1 .and. t==1 .and. i==1 .and. j==1) then 
		  sireold=ped(eff,iped,3)            !keeps the old sire
		  if (ms_alpha/=0) call sample_sire_prob(eff,iped)   !sample sire probabilities from dirichlet distr
		  call update_sire_mb(eff,iped,winv)   !updates sire for this cycle (eff=1)
		  ped(eff,iped,3)=sirenew
		end if

	  enddo


	! acceptance probability
	  alpha = .5*logdet-.5*logdet_o+.5*qform+.5*(dbelief_g(eff)-2*mb_dbelief(eff,i,j)) &
	         *(dlog(fdet_s(g(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw))) &
		    -dlog(fdet_s(newg(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw)))) &
           +.5*sum(diag( &
	       matmul( &
	       (finverse_s(g(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw)) &
	       -finverse_s(newg(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw))) &
	       ,orig_g(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw) &
		   ) &
		   *(dbelief_g(eff)-dimw-1) &
		   )) & 
           +.5*(mb_dbelief(eff,i,j)-dimw-1)*sum(diag( &
	       matmul( &
		   finverse_s(newg(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw)) &
		   ,g(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw) &
		   ) &
	       -matmul( &
           finverse_s(g(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw)) &
		   ,newg(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw) &
		   ) &
		   )) 

      if (alpha<0) then
        alpha = dexp(alpha)
		else
		alpha=1
	  endif
	  temp=gen_uniform()
  	  if(temp < alpha) then
	    ! keep symmetry
	    if (i /= j) newg(eff,(j-1)*dimw+1:j*dimw,(i-1)*dimw+1:i*dimw) &
		           =newg(eff,(i-1)*dimw+1:i*dimw,(j-1)*dimw+1:j*dimw)
	    g(eff,:,:)=newg(eff,:,:)
		else
	      newg(eff,:,:)=g(eff,:,:)
  	  endif

	  ! Acceptance rates
      if(t == mb_nround) then
       if(round <= .5*nburn) then
	    summb_accep(eff,i,j)=summb_accep(eff,i,j)+alpha
		!tune MH alghorithm half-way through burnin
	    if(mod(round,mb_skip)==0) then
	      summb_accep(eff,i,j)=summb_accep(eff,i,j)/mb_skip
          if (summb_accep(eff,i,j) > .6) then 
	    	mb_dbelief(eff,i,j)=mb_dbelief(eff,i,j)*.9   
            elseif (summb_accep(eff,i,j) < .2) then
 	    	  mb_dbelief(eff,i,j)=mb_dbelief(eff,i,j)*1.05+1
		  endif
		  if ((mb_dbelief(eff,i,j)-dimw-1) < 4) mb_dbelief(eff,i,j) = dimw + 4 ! lower bound  	   
	      summb_accep(eff,i,j)=0.
		  write(io_off+58,*) mb_dbelief(eff,i,j)
        endif
	   elseif (round > nburn .and. mod(round,skip)==0)  then
	    summb_accep(eff,i,j)=summb_accep(eff,i,j)+alpha
	   endif
	  endif

	enddo
   enddo

  enddo

  end subroutine

  subroutine restart(start)
! FC use previous values to continue chain
  integer, intent(in) :: start
  integer :: e,i,j,l,t,io,tempi,nlevr
  real :: tempr
  character (20) :: tempc           !FC dummy character variable


  open(io_s,file='solutions', status='old')
  read(io_s,*,iostat=io) tempc
  read(io_s,*,iostat=io) tempc

  j=0
  do e=1,neff
    do l=1,nlev(e)
       do t=1,ntrait
	     if (start-1 > nburn) then
           if (pos_eff(e,t)/=0) read(io_s,*,iostat=io) tempi,tempi,tempi, &
              sumsol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t), & 
              sssol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t), &
 			  sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)
           sumsol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)= &
		     sumsol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)*((start-nburn-1)/skip) 
           sssol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)= &
		     sssol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)**2*((start-nburn-1)/skip) &
			 +sumsol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)**2/((start-nburn-1)/skip)
		 else
           if (pos_eff(e,t)/=0) read(io_s,*,iostat=io) tempi,tempi,tempi,tempr,tempr, &
 				sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)
           sumsol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)= 0
           sssol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)= 0
		 endif
		 if (io.ne.0) then
		   print*,'error reading solutions at effect: ', e, ' level :', l,&
			      ' trait: ', t
		   stop
		 endif
       enddo
    enddo


  ! FC store solutions for nonparents in ram
    if (effecttype(e)==effram) then
      j=j+1
      do l=nlev(e)+1,lenped(minloc(effecttype,dim=1,mask=effecttype==effram))
        do t=1,ntrait
	     if (start-1 > nburn) then
           if (pos_eff(e,t)/=0) read(io_s,*,iostat=io) tempi,tempi,tempi, &
              sumsol_n(l-nlev(e),t+(j-1)*ntrait),& 
              sssol_n(l-nlev(e),t+(j-1)*ntrait), &
		      sol_n(l-nlev(e),t+(j-1)*ntrait)
              sumsol_n(l-nlev(e),t+(j-1)*ntrait)= &
			    sumsol_n(l-nlev(e),t+(j-1)*ntrait)*((start-nburn-1)/skip) 
              sssol_n(l-nlev(e),t+(j-1)*ntrait)= &
			    sssol_n(l-nlev(e),t+(j-1)*ntrait)**2*((start-nburn-1)/skip) &
				+sumsol_n(l-nlev(e),t+(j-1)*ntrait)**2/((start-nburn-1)/skip)
		 else
           if (pos_eff(e,t)/=0) read(io_s,*,iostat=io) tempi,tempi,tempi,tempr,tempr, &
		      sol_n(l-nlev(e),t+(j-1)*ntrait)
              sumsol_n(l-nlev(e),t+(j-1)*ntrait)= 0
              sssol_n(l-nlev(e),t+(j-1)*ntrait)= 0
		 endif
 		  if (io.ne.0) then
		   print*,'error reading solutions at effect: ', e, ' level :', l,&
			      ' trait: ', t
		   stop
		  endif
         enddo
      enddo
    endif

  enddo

  close(io_s)
  print*,'read solutions from file: "solutions"'

				  
  !FC Multiple sires posterior probabilities
    if (max_msires>1) then
	  open(io_off+52,file='msiresppr', status='old')
      read(io_off+52,*,iostat=io) tempc
      read(io_off+52,*,iostat=io) tempc
	  do e=1, lenped(1)
		if (ms_stat(e)>0) then
  		  do l=1, int(ms(p_ms(e)))
	       if (start-1 > nburn) then
            read(io_off+52,*,iostat=io) tempi, tempi,&
			         summs(l*2+p_ms(e)),ssms(l*2+p_ms(e)),ms(l*2+p_ms(e)),summs(l*2+p_ms(e)-1)
			         summs(l*2+p_ms(e))=summs(l*2+p_ms(e))*((start-nburn-1)/skip)
					 ssms(l*2+p_ms(e))=ssms(l*2+p_ms(e))**2*((start-nburn-1)/skip) &
					   +summs(l*2+p_ms(e))**2/((start-nburn-1)/skip)
					 summs(l*2+p_ms(e)-1)=summs(l*2+p_ms(e)-1)*((start-nburn-1)/skip)
		   else
            read(io_off+52,*,iostat=io) tempi, tempi,&
			         summs(l*2+p_ms(e)),ssms(l*2+p_ms(e)),ms(l*2+p_ms(e)),summs(l*2+p_ms(e)-1)
			         summs(l*2+p_ms(e))=0
					 ssms(l*2+p_ms(e))=0
					 summs(l*2+p_ms(e)-1)=0
		   endif
		  enddo
		endif
 		if (io.ne.0) then
		  print*,'error reading msires ppr at effect/animal: ', e
		  stop
		endif
	  enddo
	  close(io_off+52)
	  print*,'read msires ppr from file: "msiresppr"'
    endif

  
  
! FC Variance components
  open(io_off+53,file='varcomp', status='old')
  read(io_off+53,*,iostat=io) tempc

! random effects (G)
  do i=1,neff
    if (randomnumb(i) > 0) then
	  if (randomtype(i)==g_A_MB .or. randomtype(i)==g_diag_MB) then
	    j=randomnumb(i)*ntrait*nbreed(i)
	    else
		  j=randomnumb(i)*ntrait
	  endif
      if (j > 0) then
        read(io_off+53,*,iostat=io) tempc
		do e=1,j
          read(io_off+53,*,iostat=io) sumg(i,e,1:j)
		  if (io.ne.0) then
		    print*,'error reading G at effect: ', i
		    stop
		  endif
	      if (start-1 > nburn) then
 		    sumg(i,e,1:j)=sumg(i,e,1:j)*((start-nburn-1)/skip)
		  else
		    sumg(i,e,1:j)=0
		  endif
		enddo

        read(io_off+53,*,iostat=io) tempc
		do e=1,j
          read(io_off+53,*,iostat=io) ssg(i,e,1:j)
		  if (io.ne.0) then
		    print*,'error reading G at effect: ', i
		    stop
		  endif
	      if (start-1 > nburn) then
 		    ssg(i,e,1:j)=ssg(i,e,1:j)**2*((start-nburn-1)/skip)+sumg(i,e,1:j)**2/((start-nburn-1)/skip)
		  else
		    ssg(i,e,1:j)=0
		  endif
		enddo

        read(io_off+53,*,iostat=io) tempc
		do e=1,j
          read(io_off+53,*,iostat=io) g(i,e,1:j)
		  if (io.ne.0) then
		    print*,'error reading G at effect: ', i
		    stop
		  endif
		enddo

     endif   
    endif 
  enddo


! residuals (R)
    
  read(io_off+53,*,iostat=io) tempc
  do e=1,ntrait
	  read(io_off+53,*,iostat=io) sumr(e,:)
	  if (io.ne.0) then
		print*,'error reading R at trait: ', e
		stop
	  endif
	  if (start-1 > nburn) then
 		sumr(e,:)=sumr(e,:)*((start-nburn-1)/skip)
	  else
		sumr(e,:)=0
	  endif
  enddo

  read(io_off+53,*,iostat=io) tempc
  do e=1,ntrait
	  read(io_off+53,*,iostat=io) ssr(e,:)
	  if (io.ne.0) then
		print*,'error reading R at trait: ', e
		stop
	  endif
	  if (start-1 > nburn) then
 		ssr(e,:)=ssr(e,:)**2*((start-nburn-1)/skip)+sumr(e,:)**2/((start-nburn-1)/skip)
	  else
		ssr(e,:)=0
	  endif
  enddo

  read(io_off+53,*,iostat=io) tempc
  do e=1,ntrait
	  read(io_off+53,*,iostat=io) r(e,:)
	  if (io.ne.0) then
		print*,'error reading R at trait: ', e
		stop
	  endif
  enddo
  close(io_off+53)
  print*,'read variance components from file: "varcomp"'


! FC average likelihood and deviance at average solutions for each observation
  if (start-1 > nburn) then
	open(io_off+54,file='loglike_obs', status='old')
	read(io_off+54,*,iostat=io) tempc
	read(io_off+54,*,iostat=io) tempc
	do i=1,nrec
!	  read(io_off+54,*,iostat=io) tempi, logl_obs(i,:)
	  read(io_off+54,*,iostat=io) tempi, logl_obs(i,1:2)
	  if (io.ne.0) then
		print*,'error reading loglilke_obs at record: ', i
		stop
	  endif
	  logl_obs(i,1)=((start-nburn-1)/skip)/logl_obs(i,1)
	  logl_obs(i,2)=((start-nburn-1)/skip)*logl_obs(i,2)
	enddo
	close(io_off+54)
	print*,'read loglike_obs from file: "loglike_obs"'
  else
	logl_obs=0
  endif

! FC stores sample of varcomps
  open(io_off+55,file='varcompsam', status='old',position='append')

! FC average loglikelihood for each round
  open(io_off+56,file='loglike_rnd', status='old',position='append')

! FC samples for two sire assignments
  if (ms_alpha/=0.) open(io_off+57,file='msirespprsam',status='old',position='append')

! FC proposal density degrees of belief tuning 
  if (any(randomtype == g_A_MB) .or. any(randomtype == g_diag_MB) .or. residualtype/=r_homo) then
    if(start <= .5*nburn) open(io_off+58,file='mh_dbeliefchg',status='old',position='append')
  endif
  if (any(effectsave == 1)) open(io_off+61,file='solutionsam',status='old',position='append')


! FC acceptance rates and degrees of belief for proposal density
  if (any(randomtype == g_A_MB) .or. any(randomtype == g_diag_MB) .or. residualtype/=r_homo) then
	  open(io_off+59,file='mh_acceptance', status='old')
	  do i=1,neff
		if (randomtype(i)==g_A_MB .or. randomtype(i)==g_diag_MB) then
		  do j=1,nbreed(i)
			do k=j,nbreed(i)
			  read(io_off+59,*,iostat=io) tempc
			  read(io_off+59,*,iostat=io) summb_accep(i,j,k), mb_dbelief(i,j,k)
			  if (io.ne.0) then
				print*,'error reading acceptance rate at effect: ', i
				stop
			  endif
 			  summb_accep(i,j,k)=summb_accep(i,j,k)*((start-nburn-1)/skip)
			enddo
		  enddo
          if (start-1 < nburn) summb_accep(i,:,:)=0
		endif   
	  enddo
      if (neffr>0) then
	do i=1,neffr
	  e=pos_effr(i)
	  if (effecttype(e)==effcov .or. effecttype(e)==effrnorm) then
	    read(io_off+59,*,iostat=io) tempc
	    read(io_off+59,*,iostat=io) sumr_accep(i), r_dbeliefprop(i)
	    if (io.ne.0) then
	      print*,'error reading acceptance rate at structural effect: ', pos_effr(i)
	      stop
	    endif
 	    sumr_accep(i)=sumr_accep(i)*((start-nburn-1)/skip)
	  elseif (effecttype(e)==effcross .and. randomnumb(e) > 0) then 
	    read(io_off+59,*,iostat=io) tempc
	    read(io_off+59,*,iostat=io) sumr_alphaaccep(i), r_alpha_varprop(i)
	    if (io.ne.0) then
	      print*,'error reading acceptance rate for alpha at structural effect: ', pos_effr(i)
	      stop
	    endif
 	    sumr_alphaaccep(i)=sumr_alphaaccep(i)*((start-nburn-1)/skip)
          endif
	enddo
        if (start-1 < nburn) sumr_accep=0
      endif
	  if (residualtype == r_studentt .or. residualtype == r_struct_studentt) then
		read(io_off+59,*,iostat=io) tempc
		read(io_off+59,*,iostat=io) sumw_alphaaccep, w_alpha_varprop
		if (io.ne.0) then
			print*,'error reading acceptance rate for alpha robust'
			stop
		endif
 		sumw_alphaaccep=sumw_alphaaccep*((start-nburn-1)/skip)
      endif
	  close(io_off+59)
	  print*,'Read acceptance rates and proposal dbelief from file: "mh_acceptance"'
  endif
! FC

! FC acceptance rates and degrees of belief for proposal density


! FC structural residual effects in file 'solutions'
  if (neffr > 0) then
	  open(io_off+60,file='structural_r', status='old')
	  read(io_off+60,*,iostat=io) tempc
	  read(io_off+60,*,iostat=io) tempc

	  do j=1,neffr
		e=pos_effr(j)
		if (effecttype(e)==effcov .or. effecttype(e)==effrnorm) then ! set nlev=1 even when nested covariate is present
			nlevr=1
          else
		    nlevr=nlev(e)
        endif
		if (effecttype(e)==effcross .and. randomnumb(e) > 0) then 
			if (start-1 > nburn) then 
			      read(io_off+60,*,iostat=io) tempi,tempi,tempc,&
				  sumr_alpha(j),ssr_alpha(j),r_alpha(j)
				  sumr_alpha(j)=sumr_alpha(j)*((start-nburn-1)/skip)
				  ssr_alpha(j)=ssr_alpha(j)**2*((start-nburn-1)/skip) &
				              +sumr_alpha(j)**2/((start-nburn-1)/skip)
			 else
			      read(io_off+60,*,iostat=io) tempi,tempi,tempc,tempr,tempr,r_alpha(j)
				  sumr_alpha(j)=0
				  ssr_alpha(j)=0
			 endif
			 if (io.ne.0) then
			   print*,'error reading structural R alpha at effect: ', e
			   stop
			 endif
		endif
		do l=1,nlevr
		   do t=1,ntrait
			 if (start-1 > nburn) then
			   if (pos_eff(e,t)/=0) read(io_off+60,*,iostat=io) tempi,tempi,tempi, &
				  sumr_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t), & 
				  ssr_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t), &
 				  r_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)
			      sumr_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)= &
				    sumr_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)*((start-nburn-1)/skip) 
			      ssr_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)= &
				    ssr_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)**2*((start-nburn-1)/skip) &
				    +sumr_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)**2/((start-nburn-1)/skip)
			 else
			   if (pos_eff(e,t)/=0) read(io_off+60,*,iostat=io) tempi,tempi,tempi,tempr,tempr, &
 				  r_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)
				  sumr_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)=0 
				   ssr_sol(sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t)=0
			 endif
			 if (io.ne.0) then
			   print*,'error reading structural R at effect: ', e, ' level :', l,&
					  ' trait: ', t
			   stop
			 endif
		   enddo
		enddo
	  enddo
	  close(io_off+60)
	  print*,'read structural R from file: "Structural_r"'
  endif

! FC Robustness weights in file 'robustness_w'
  if (residualtype == r_slash .or. residualtype == r_struct_slash .or. &
      residualtype == r_studentt .or. residualtype == r_struct_studentt) then
	  open(io_off+62,file='robustness_w', status='old')
	  read(io_off+62,*,iostat=io) tempc
	  read(io_off+62,*,iostat=io) tempc

	  do t=1,ntrait
		if (start-1 > nburn) then 
		  read(io_off+62,*,iostat=io) tempi,tempc,sumw_alpha,ssw_alpha,w_alpha
		  sumw_alpha=sumw_alpha*((start-nburn-1)/skip)
		  ssw_alpha=ssw_alpha**2*((start-nburn-1)/skip) &
			       +sumw_alpha**2/((start-nburn-1)/skip)
		else
		  read(io_off+62,*,iostat=io) tempi,tempc,tempr,tempr,w_alpha
		  sumw_alpha=0
		  ssw_alpha=0
		endif
		if (io.ne.0) then
		  print*,'error reading robustness w_alpha'
		  stop
		endif
	    do j=1,nrec
		  if (start-1 > nburn) then
		    read(io_off+62,*,iostat=io) tempi,tempi, sumweight_robust(j),& 
				    ssweight_robust(j), weight_robust(j)
				    sumweight_robust(j)=sumweight_robust(j)*((start-nburn-1)/skip) 
				    ssweight_robust(j)=ssweight_robust(j)**2*((start-nburn-1)/skip) &
				    +sumweight_robust(j)**2/((start-nburn-1)/skip)
		  else
		    read(io_off+62,*,iostat=io) tempi,tempi,tempr,tempr,weight_robust(j)
				    sumweight_robust(j)=0 
				    ssweight_robust(j)=0
		  endif
		  if (io.ne.0) then
		    print*,'error reading robustness w at record: ', j, ' trait: ', t
		    stop
		  endif
		enddo
	  enddo
	  close(io_off+62)
	  print*,'Read robustness weights from file: "robustness_w"'
  endif


  end subroutine

 subroutine structural_r(e)
! sample structural effects on r
! for single trait... so far
  integer :: i,j,k,n,o,t,t1,t2
  real (r8):: ee(neq,ntrait,ntrait),r_solold(neq),e(nrec,ntrait)
  real (r8) :: prop,ratio,alpha,temp,suminv,sumlog
  real (r4) :: a,b

  rinv=finverse_s(r)
  r_solold=1.

  do t=1, r_nround
   ee=0
   do i=1,neffr
     n=pos_effr(i) !address in pos_eff
     if (effecttype(n)==effcov .or. effecttype(n)==effrnorm) then  ! MH sampling required
       o=address1(n,1,1) 
!       r_solold(o)=r_sol(o)
     ! generate a candidate value with E(r_solnew)=r_sol
!       r_sol(o)=gen_invwishart(1/(r_sol(o)*(r_dbeliefprop(i)-2)),r_dbeliefprop(i))
!       ratio=r_sol(o)/r_solold(o)
       prop=gen_invwishart(1/(r_sol(o)*(r_dbeliefprop(i)-2)),r_dbeliefprop(i))
       ratio=prop/r_sol(o)
       if (t == 1 .and. effecttype(n)==effrnorm) r_nrec(o)=0.0  ! reset sum of weight_cov
!      print*, 'eff:', n, 'proposal:',r_sol(o), 'current:', r_solold(o), 'address:', o
     endif
     do irec=1,nrec
       call find_addresses
       ! adds the sum of weight_cov for each scaling factor
       if (round == start .and. t == 1 .and. effecttype(n)/=effrnorm) then
	 o=address(n,1)
         if (effecttype(n)==effuncov) then
           r_nrec(o)=r_nrec(o)+1 ! weight_cov for effuncov includes rnorm solutions          
         else
           r_nrec(o)=r_nrec(o)+weight_cov(n,1) ! 
         endif   
       ! allows reacnorm environmental heteroskedasticity base on covariate value
       elseif (t == 1 .and. effecttype(n)==effrnorm) then 
         ! update sum of weight_cov 
	 o=address1(n,1,1) 
	 r_nrec(o)=r_nrec(o)+weight_cov(n,1)  
       endif
      ! updates e according to r_sol changes
       do t1=1,ntrait
 	 ! the last effect from previous MH cycle
	 if (i == 1 .and. t > 1 .and. weight_cov(pos_effr(neffr),1)/=0) then
           ! allows reacnorm & rr heteroskedasticity base on covariate value
           if (effecttype(pos_effr(neffr))==effrnorm) then 
             o=address1(pos_effr(neffr),1,1) 
	   else
	     o=address(pos_effr(neffr),1)
	   endif
    	   if (r_sol(o)/=r_solold(o)) e(irec,t1)=e(irec,t1)*dsqrt((r_solold(o)/ &
				    r_sol(o))**weight_cov(pos_effr(neffr),1))
         ! the previous effect
	 elseif (i > 1 .and. weight_cov(pos_effr(i-1),1)/=0) then 
	   ! allows reacnorm & rr heteroskedasticity base on covariate value
	   if (effecttype(pos_effr(i-1))==effrnorm) then 
	     o=address1(pos_effr(i-1),1,1) 
	   else
	     o=address(pos_effr(i-1),1)
	   endif
	   if (r_sol(o)/=r_solold(o)) e(irec,t1)=e(irec,t1)*dsqrt((r_solold(o)/ &
				    r_sol(o))**weight_cov(pos_effr(i-1),1))
	 endif
	 do t2=t1,ntrait
!		if (orig_r(t1,t2) /= 0 .and. weight_cov(n,1) /=0) then
	   if (weight_cov(n,1) /=0) then
	   ! qform difference for MH 
	     if (effecttype(n)==effcov .or. effecttype(n)==effrnorm) then 
 	       o=address1(n,1,1) 
	       ee(o,t1,t2)=ee(o,t1,t2)+(e(irec,t1)*e(irec,t2)- &
	         (e(irec,t1)*e(irec,t2)/(ratio**weight_cov(n,1))))*rinv(t1,t2)/2
             else ! effcross qform for inverted scaled chi-square
	       o=address(n,1)
	       ee(o,t1,t2)=ee(o,t1,t2)+(e(irec,t1)*e(irec,t2)* &
		          (r_sol(address(n,1))**weight_cov(n,1)))*rinv(t1,t2)
! weight_cov remove because should be always 1 and is not for effuncov
!		          (r_sol(address(n,1))))*rinv(t1,t2)
             endif
	     if (t2 > t1) ee(address(n,1),t2,t1)=ee(address(n,1),t1,t2)
	   endif   
	 enddo                       
       enddo
     enddo

!if (round == start .and. t == 1) print*, 'effect', n

     if (effecttype(n)==effcov .or. effecttype(n)==effrnorm) then  ! MH sampling required
       o=address1(n,1,1) 
     ! acceptance probability
       r_solold(o)=r_sol(o)
       alpha = -.5*(r_nrec(o)-2*r_dbeliefprop(i)-2)*dlog(ratio)&
	          +.5*(r_dbeliefprop(i)-2)*(1/ratio-ratio)+ee(o,1,1)
       if (alpha<0) then
	 alpha = dexp(alpha)
       else
         alpha=1
       endif
       temp=gen_uniform()
       if(alpha >= temp) r_sol(o)=prop

! if(round==start .or. mod(round,skip*10)==0) print*,'round, effect, addr, r_nrec, ee', &
!                  n, o, r_nrec(o), ee(o,1,1), alpha, temp, r_solold(o), r_sol(o), prop 

     ! Acceptance rates
       if(t == r_nround) then
	 if(round <= .5*nburn) then
	   sumr_accep(i)=sumr_accep(i)+alpha
 	   !tune MH alghorithm half-way through burnin
	   if(mod(round,r_skip)==0) then
	     sumr_accep(i)=sumr_accep(i)/r_skip
	     if (sumr_accep(i) > .6) then 
	       r_dbeliefprop(i)=r_dbeliefprop(i)*.9   
	     elseif (sumr_accep(i) < .2) then
 	       r_dbeliefprop(i)=r_dbeliefprop(i)*1.05+1
	     endif
	     if ((r_dbeliefprop(i)-2) < 4) r_dbeliefprop(i) = 6 ! lower bound  	   
	     sumr_accep(i)=0.
	     write(io_off+58,*) r_dbeliefprop(i)
	   endif
	 elseif (round > nburn .and. mod(round,skip)==0) then
	   sumr_accep(i)=sumr_accep(i)+alpha
	 endif
       endif
       if (t == r_nround .and. mod(round,skip)==0 .and. effrsave(i)==1) &
           write(io_off+55,*) r_sol(o)
     elseif (effecttype(n)==effcross .or. effecttype(n)==effuncov) then
       if (randomnumb(n) > 0) then
	 suminv=0.; sumlog=0.
         do j=1, nlev(n)
	   o=address1(n,j,1) 
!print*,'address in solutions and r_nrec of i', n, j, o, r_sol(o), r_nrec(o)
          !call pos_def(ee,'ee not pos def, fixed', 1d-4)
           r_solold(o)=r_sol(o)
	   a=r_nrec(o)/2+r_alpha(i)
	   b=ee(o,1,1)/2+r_alpha(i)-1
           r_sol(o)= b/sgamma(a) 
	   if (r_sol(o)<=0) print*, 'eff:', n, 'level:',j, 'r_nrec:', &
	     r_nrec(o), 'ee: ', ee(o,1,1), 'a: ', a, 'b: ', b, 'r_alpha', &
             r_alpha(i), 'sol: ', r_sol(o) 
           suminv=suminv+1/r_sol(o)
	   sumlog=sumlog+dlog(r_sol(o))
           if (t == r_nround .and. mod(round,skip)==0 .and. effrsave(i)==1) & 
             write(io_off+55,*) r_sol(o)
         enddo
         call sample_r_alpha_ggam(i,suminv,sumlog)
       elseif (randomnumb(n)==0) then
         do j=1, nlev(n)-1
	   o=address1(n,j,1) 
           !call pos_def(ee,'ee not pos def, fixed', 1d-4)
           r_solold(o)=r_sol(o)
	   a=r_nrec(o)/2-1
	   b=ee(o,1,1)/2
           r_sol(o)= b/sgamma(a) 
           if (t == r_nround .and. mod(round,skip)==0 .and. effrsave(i)==1) &
             write(io_off+55,*) r_sol(o)
         enddo
       endif
     else
       print*, 'wrong effecttype: ', n
       stop
     endif
   enddo
  enddo

! accumulate mean and variance of dispersion parameters after burning
  if (round > nburn .and. mod(round,skip)==0) then
    sumr_sol=sumr_sol+r_sol
    ssr_sol=ssr_sol+r_sol**2
  endif	

  if (residualtype == r_struct_slash .or. residualtype == r_struct_studentt) then
  ! Update e with the last sample effect on r
    do irec=1,nrec
      call find_addresses
      do t1=1,ntrait
        !update e according to r_sol changes
        if (weight_cov(pos_effr(neffr),1)/=0) then
        ! the last effect from previous MH cycle
        ! allows reacnorm & rr heteroskedasticity base on covariate value
	  if (nestedcov(pos_effr(neffr),1) > 0) then 
	    o=address1(pos_effr(neffr),1,1) 
	  else
	    o=address(pos_effr(neffr),1)
	  endif
	  if (r_sol(o)/=r_solold(o)) e(irec,t1)=e(irec,t1)*dsqrt((r_solold(o)/ &
				    r_sol(o))**weight_cov(pos_effr(neffr),1))
        endif
      enddo
    enddo
  endif	

 end subroutine

  subroutine sample_r_alpha(i,suminv,sumlog)
! sample heterogeneity parameter for ith random structural effects on r
  integer, parameter:: r_alpha_nround=5, r_alpha_skip=10
  integer :: i,j,k,n,t
  real (r8) :: ratio,alpha,temp,suminv,sumlog,r_alphanew,gamma_log,priormed=6
 
  n=pos_effr(i)
  
  do t=1, r_alpha_nround
	 ! generate a candidate value with E(r_alphanew)=r_alpha
      r_alphanew=dexp(gen_normal(dlog(r_alpha(i)),r_alpha_varprop(i)))

	! acceptance probability
	  alpha  = nlev(n)*( &
	         (r_alphanew/2)*dlog((r_alphanew-2)/2) &
	         -(r_alpha(i)/2)*dlog((r_alpha(i)-2)/2) &
             -gamma_log(r_alphanew/2)+gamma_log(r_alpha(i)/2)) & 
			 -2*(dlog(priormed+r_alphanew)-dlog(priormed+r_alpha(i))) &
			 -.5*(r_alphanew-r_alpha(i))*(suminv+sumlog) &
			 +dlog(r_alphanew)-dlog(r_alpha(i))
      if (alpha<0) then
        alpha = dexp(alpha)
		temp=gen_uniform()
	    if(temp < alpha) r_alpha(i)=r_alphanew
	  else
		alpha=1
		r_alpha(i)=r_alphanew
	  endif
  enddo

! Acceptance rates
  if(round <= .5*nburn) then
	sumr_alphaaccep(i)=sumr_alphaaccep(i)+alpha
	!tune MH alghorithm half-way through burnin
	if(mod(round,r_alpha_skip)==0) then
	  sumr_alphaaccep(i)=sumr_alphaaccep(i)/r_alpha_skip
      if (sumr_alphaaccep(i) > .6) then 
	    r_alpha_varprop(i)=r_alpha_varprop(i)*1.05   
      elseif (sumr_alphaaccep(i) < .2) then
 	    r_alpha_varprop(i)=r_alpha_varprop(i)*.9
	  endif
	  sumr_alphaaccep(i)=0.
	  write(io_off+58,*) r_alpha_varprop(i)
    endif
  elseif (round > nburn .and. mod(round,skip)==0) then
	sumr_alphaaccep(i)=sumr_alphaaccep(i)+alpha
  endif
! Save samples
  if (mod(round,skip)==0) write(io_off+55,*) r_alpha(i)

! accumulate mean and variance of r_alpha after burning
  if (round > nburn .and. mod(round,skip)==0) then
  	sumr_alpha(i)=sumr_alpha(i)+r_alpha(i)
  	ssr_alpha(i)=ssr_alpha(i)+r_alpha(i)**2
  endif	

  end subroutine


  subroutine sample_r_alpha_gam(i,suminv,sumlog)
! sample heterogeneity parameter for ith random structural effects on r
  integer, parameter:: r_alpha_nround=5, r_alpha_skip=10
  integer :: i,j,k,n,t
  real (r8) :: ratio,alpha,temp,suminv,sumlog,r_alphanew,gamma_log,r_alpha_r=.06,r_alpha_s=.01
 
  n=pos_effr(i)
  
  do t=1, r_alpha_nround
	 ! generate a candidate value with E(r_solnew)=r_sol
      r_alphanew=dexp(gen_normal(dlog(r_alpha(i)),r_alpha_varprop(i)))

	! acceptance probability
	  alpha = nlev(n)*( &
	         (r_alphanew/2)*dlog((r_alphanew-2)/2) &
	         -(r_alpha(i)/2)*dlog((r_alpha(i)-2)/2) &
             -gamma_log(r_alphanew/2)+ gamma_log(r_alpha(i)/2)) &
			 -r_alpha_s*(r_alphanew-r_alpha(i)) &
			 -.5*(r_alphanew-r_alpha(i))*(suminv+sumlog) &
			 +r_alpha_r*(dlog(r_alphanew)-dlog(r_alpha(i)))

      if (alpha<0) then
        alpha = dexp(alpha)
		temp=gen_uniform()
	    if(temp < alpha) r_alpha(i)=r_alphanew
	  else
		alpha=1
	    r_alpha(i)=r_alphanew
	  endif
  enddo

! Acceptance rates
  if(round <= .5*nburn) then
	sumr_alphaaccep(i)=sumr_alphaaccep(i)+alpha
	!tune MH alghorithm half-way through burnin
	if(mod(round,r_alpha_skip)==0) then
	  sumr_alphaaccep(i)=sumr_alphaaccep(i)/r_alpha_skip
      if (sumr_alphaaccep(i) > .6) then 
	    r_alpha_varprop(i)=r_alpha_varprop(i)*1.05   
      elseif (sumr_alphaaccep(i) < .2) then
 	    r_alpha_varprop(i)=r_alpha_varprop(i)*.9
	  endif
	  sumr_alphaaccep(i)=0.
	  write(io_off+58,*) r_alpha_varprop(i)
    endif
  elseif (round > nburn .and. mod(round,skip)==0) then
	sumr_alphaaccep(i)=sumr_alphaaccep(i)+alpha
  endif
! save samples
  if (mod(round,skip)==0) write(io_off+55,*) r_alpha(i)

! accumulate mean and variance of r_alpha after burning
  if (round > nburn .and. mod(round,skip)==0) then
  	sumr_alpha(i)=sumr_alpha(i)+r_alpha(i)
  	ssr_alpha(i)=ssr_alpha(i)+r_alpha(i)**2
  endif	

  end subroutine


  subroutine sample_r_alpha_ggam(i,suminv,sumlog)
! sample heterogeneity parameter for ith random structural effects on r
  integer, parameter:: r_alpha_nround=5, r_alpha_skip=10
  integer :: i,j,k,n,t
  real (r8) :: ratio,alpha,temp,suminv,sumlog,r_alphanew,gamma_log,r_alpha_r=.03,r_alpha_s=.01
 
  n=pos_effr(i)
  
  do t=1, r_alpha_nround
	 ! generate a candidate value with E(r_solnew)=r_sol >1;
      r_alphanew=dexp(gen_normal(dlog(r_alpha(i)),r_alpha_varprop(i)))
     if (r_alphanew > 1) then

	! acceptance probability
	  alpha = nlev(n)*( &
	         r_alphanew*dlog(r_alphanew-1) &
	         -r_alpha(i)*dlog(r_alpha(i)-1) &
             -gamma_log(r_alphanew)+gamma_log(r_alpha(i)) &
			 ) &
			 -r_alpha_s*(r_alphanew-r_alpha(i)) &
			 -(r_alphanew-r_alpha(i))*(suminv+sumlog) &
			 +r_alpha_r*(dlog(r_alphanew)-dlog(r_alpha(i)))

      if (alpha<0) then
        alpha = dexp(alpha)
		temp=gen_uniform()
	    if(temp < alpha) r_alpha(i)=r_alphanew
	  else
		alpha=1
	    r_alpha(i)=r_alphanew
	  endif

	 else
	  alpha=0  
     endif
  enddo

! Acceptance rates
  if(round <= .5*nburn) then
	sumr_alphaaccep(i)=sumr_alphaaccep(i)+alpha
	!tune MH alghorithm half-way through burnin
	if(mod(round,r_alpha_skip)==0) then
	  sumr_alphaaccep(i)=sumr_alphaaccep(i)/r_alpha_skip
      if (sumr_alphaaccep(i) > .6) then 
	    r_alpha_varprop(i)=r_alpha_varprop(i)*1.05   
      elseif (sumr_alphaaccep(i) < .2) then
 	    r_alpha_varprop(i)=r_alpha_varprop(i)*.9
	  endif
	  sumr_alphaaccep(i)=0.
	  write(io_off+58,*) r_alpha_varprop(i)
    endif
  elseif (round > nburn .and. mod(round,skip)==0) then
	sumr_alphaaccep(i)=sumr_alphaaccep(i)+alpha
  endif
! save samples
  if (mod(round,skip)==0) write(io_off+55,*) r_alpha(i)

! accumulate mean and variance of r_alpha after burning
  if (round > nburn .and. mod(round,skip)==0) then
  	sumr_alpha(i)=sumr_alpha(i)+r_alpha(i)
  	ssr_alpha(i)=ssr_alpha(i)+r_alpha(i)**2
  endif	

  end subroutine


  subroutine slash_weights(e)
! sample robustness weights with slash specification on y's
! for single trait... so far
  integer :: i,j,k,n,o,t,t1,t2,status
  real (r8):: e(nrec,ntrait),ee(ntrait,ntrait)
  real (r8):: a,b,p,uplim,bound,sumlog,w_alpha_r=.015,w_alpha_s=.01

 
  rinv=finverse_s(r)
  a=w_alpha+.5
  uplim=1.
  sumlog=0.

  do irec=1,nrec

	do t1=1,ntrait
	  do t2=t1,ntrait
!		if (orig_r(t1,t2) /= 0) then
		  ee(t1,t2)=e(irec,t1)*e(irec,t2)*rinv(t1,t2)
		  if (t2 > t1) ee(t2,t1)=ee(t1,t2)
!		endif   
	  enddo                       
	enddo
    b=ee(1,1)/(2*weight_robust(irec))
	p=0.
	do while (p==0)
      call cdfgam ( 1, p, 1-p, uplim, a, b, status, bound )
	  if (p==0) b=b*10 ! approx trunc gamma when a large b very small with a large b
    enddo
	p=p*gen_uniform() 
    call cdfgam ( 2, p, 1-p, weight_robust(irec), a, b, status, bound )
	if (weight_robust(irec)>1 .or. weight_robust(irec)<=0) then
	  print *,'round',round,'rec',irec,'status',status,'a',a,'b',b,'bound',bound,'w',weight_robust(irec)
	  stop
    endif
    sumlog=sumlog+dlog(weight_robust(irec))
    if (mod(round,skip)==0 .and. nrec==22717) then ! to dg_pwg
      if (irec==10581 .or. irec==10656 .or. irec==10765 .or. irec==19483 .or. irec==19754 &
	      .or. irec==19286) write(io_off+55,*) weight_robust(irec)
    endif
  enddo

  ! generate a sample of the robustness parameter w_alpha
  a=nrec+w_alpha_r
  b=w_alpha_s-sumlog
  p=gen_uniform() 
  call cdfgam ( 2, p, 1-p, w_alpha, a, b, status, bound )

  if (mod(round,skip)==0) write(io_off+55,*) w_alpha

! accumulate mean and variance of dispersion parameters after burning
  if (round > nburn .and. mod(round,skip)==0) then
  	sumweight_robust=sumweight_robust+weight_robust
	ssweight_robust=ssweight_robust+weight_robust**2
	sumw_alpha=sumw_alpha+w_alpha
	ssw_alpha=ssw_alpha+w_alpha**2
  endif	

  end subroutine

  subroutine studentt_weights(e)
! sample robustness weights with slash specification on y's
! for single trait... so far
  integer :: i,j,k,n,o,t,t1,t2,status
  real (r8):: e(nrec,ntrait),ee(ntrait,ntrait)
  real (r8):: a,b,p,uplim,bound,sumlog,sum,w_alpha_r=.04,w_alpha_s=.01

  rinv=finverse_s(r)
  a=(w_alpha+1)/2
  sumlog=0.; sum=0.

  do irec=1,nrec

	do t1=1,ntrait
	  do t2=t1,ntrait
!		if (orig_r(t1,t2) /= 0) then
		  ee(t1,t2)=e(irec,t1)*e(irec,t2)*rinv(t1,t2)
		  if (t2 > t1) ee(t2,t1)=ee(t1,t2)
!		endif   
	  enddo                       
	enddo
    b=(w_alpha+ee(1,1)/weight_robust(irec))/2
	p=gen_uniform() 
    call cdfgam ( 2, p, 1-p, weight_robust(irec), a, b, status, bound )
	if (weight_robust(irec)<=0 .or. status /= 0) then
	  print *,'round', round,'rec',irec,'status',status,'a',a,'b',b,'bound',bound, 'w',weight_robust(irec)
	  stop 
    endif
    sumlog=sumlog+dlog(weight_robust(irec))
    sum=sum+weight_robust(irec)
    if (mod(round,skip)==0 .and. nrec==22717) then ! to dg_pwg
      if (irec==10581 .or. irec==10656 .or. irec==10765 .or. irec==19483 .or. irec==19754 &
	      .or. irec==19286) write(io_off+55,*) weight_robust(irec)
    endif
  enddo

! accumulate mean and variance after burning
  if (round > nburn .and. mod(round,skip)==0) then
  	sumweight_robust=sumweight_robust+weight_robust
	ssweight_robust=ssweight_robust+weight_robust**2
  endif	
  
  call sample_w_alpha_gam(sum,sumlog)

  end subroutine

  subroutine sample_w_alpha_gam(sum,sumlog)
! sample robustness parameter for student t with gamma(r,s) prior
  integer, parameter:: w_alpha_nround=5, w_alpha_skip=10
  integer :: i,j,k,n,t
  real (r8) :: ratio,alpha,temp,sum,sumlog,w_alphanew,gamma_log,w_alpha_r=.06,w_alpha_s=.01
  real (r8) :: two
 
  two=2.

  do t=1, w_alpha_nround
	 ! generate a candidate value with E(r_solnew)=r_sol
      w_alphanew=dexp(gen_normal(dlog(w_alpha),w_alpha_varprop))

	! acceptance probability
	  alpha = nrec*( &
	         (w_alphanew/2)*dlog(w_alphanew) &
	         -(w_alpha/2)*dlog(w_alpha) &
             -gamma_log(w_alphanew/2)+gamma_log(w_alpha/2)&
			 -((w_alphanew-w_alpha)/2)*dlog(two)) &
			 +.5*(w_alphanew-w_alpha)*(sumlog-sum-w_alpha_s) &
			 +w_alpha_r*(dlog(w_alphanew)-dlog(w_alpha))

      if (alpha<0) then
        alpha = dexp(alpha)
		temp=gen_uniform()
	    if(temp < alpha) w_alpha=w_alphanew
	  else
		alpha=1
	    w_alpha=w_alphanew
	  endif
  enddo

! Acceptance rates
  if(round <= .5*nburn) then
	sumw_alphaaccep=sumw_alphaaccep+alpha
	!tune MH alghorithm half-way through burnin
	if(mod(round,w_alpha_skip)==0) then
	  sumw_alphaaccep=sumw_alphaaccep/w_alpha_skip
      if (sumw_alphaaccep > .6) then 
	    w_alpha_varprop=w_alpha_varprop*1.05   
      elseif (sumw_alphaaccep < .2) then
 	    w_alpha_varprop=w_alpha_varprop*.9
	  endif
	  sumw_alphaaccep=0.
	  write(io_off+58,*) w_alpha_varprop
    endif
  elseif (round > nburn .and. mod(round,skip)==0) then
	sumw_alphaaccep=sumw_alphaaccep+alpha
  endif
! save samples
  if (mod(round,skip)==0) write(io_off+55,*) w_alpha

! accumulate mean and variance of w_alpha after burning
  if (round > nburn .and. mod(round,skip)==0) then
  	sumw_alpha=sumw_alpha+w_alpha
  	ssw_alpha=ssw_alpha+w_alpha**2
  endif	

  end subroutine

 subroutine sample_h_rnorm(i,k)
! sample unknown covariate (h) for heteroskedastic reaction norms model
! i = unknown covariate effect and k = reaction norm effect
! for single trait... 
  integer :: i,j,k,l,m,o,t,firsteq,lasteq
  real (r8):: ee_diff(nlev(i)),solold(neq) 
  real (r8) :: accept,alpha,temp

  o=address1(k,1,1) ! address of rnorm heteroskedasticity parameter on r_sol
  firsteq=sum(nlev(1:i-1))*ntrait+1;
  lasteq=sum(nlev(1:i))*ntrait
  do t=1, r_nround
    do j=1,nlev(i)
      m=address1(i,j,1) 
      solold(m)=sol(m)
      sol(m)= gen_normal(sol(m),h_varprop(i))  ! random walk
    enddo
    ee_diff=0
    do irec=1,nrec
      call decode_record
      call find_addresses
      call decode_r
      call predict_missing_y
      call decode_record_rnorm(1) ! calculates y_h(t=1) 
      j=address(i,1) ! address in solutions and r_nrec of i
      if (round == start .and. t == 1) then
	r_nrec(j)=r_nrec(j)+1.
      endif
      l=indata(irec,pos_eff(i,1))  ! level of i for current record
      ee_diff(l)=ee_diff(l)+((((y(1)-weight_cov(i,1)*solold(j))**2) &
                 /(r_sol(o)**solold(j))) &
                 -(((y(1)-weight_cov(i,1)*sol(j))**2) &
                 /(r_sol(o)**sol(j))))*weight_y(1)
    enddo
    ! acceptance probability
    do j=1,nlev(i)
      m=address1(i,j,1) 
      alpha = 0.5*(r_nrec(m)*dlog(r_sol(o)))*(solold(m)-sol(m))  &
             +0.5*(1/g(i,1,1))*(solold(m)**2-sol(m)**2)  &
	     +0.5*(1/orig_r(1,1))*ee_diff(j)
      if (alpha<0) then
	alpha = dexp(alpha)
	temp=gen_uniform()
	if(temp > alpha) sol(m)=solold(m) !retain current sample
      else
        alpha=1
      endif
     ! acceptance rates
      sumh_accep(m)=sumh_accep(m)+alpha
    enddo
  enddo  

! Acceptance rates
  if(round <= .5*nburn) then
    !tune MH alghorithm half-way through burnin
    if(mod(round,r_skip)==0) then
      accept=sum(sumh_accep(firsteq:lasteq))/nlev(i)/r_skip/r_nround
      if (accept > .6) then 
	h_varprop(i)=h_varprop(i)*1.05   
      elseif (accept < .2) then
 	h_varprop(i)=h_varprop(i)*.9
      endif
      sumh_accep(firsteq:lasteq)=0.
      write(io_off+58,*) i, h_varprop(i), accept
    endif
  endif

  !print acceptance rates
  if(round == nround) then
    sumh_accep(firsteq:lasteq)=sumh_accep(firsteq:lasteq)/(nround-.5*nburn)/r_nround
    accept=sum(sumh_accep(firsteq:lasteq))/nlev(i)
    print*,' MH global acceptance rate for unknown covariate: ',i,' = ', accept
    print*,' level        acceptance '
    do j=firsteq,lasteq  
      print*,j-firsteq+1,sumh_accep(j)
    enddo
  endif

 end subroutine

 ! VJ - Code from Intersires.f90
 subroutine set_dinv(eff) ! Famula (1992) JAS algorithm
 integer :: i, k, eff, iped
 ! u = contem os elementos da diagonal de A 
 ! v = contem os elementos da diagonal da matriz de Henderson L
 ! w = contem soma dos produtos das probabilidades dos touros e vacas para todos os animais

 real :: v(lenped(eff))
 real :: w(lenped(eff))  ! pode ser declarado com w pa não existe essa variavel global
 real :: u(lenped(eff)) 
 real:: ts, td           ! se é a posição, não seria inteiro? , 
 	                ! mesmo tentando com integer deu o mesmo erro.

 w=0; v=0; u=0
 ts=0  
 td=0                    ! aqui esta indicando o numero de touros candidatos 

 do iped=1, nlev(eff)-ngrp(eff)    !iped=1, lenped(eff)
     if (ped(eff,iped,2)== 0 .and. ped(eff,iped,3)== 0)  then ! Animal com pae e mae desconhecido
        v(iped)=1   
     else 
        v(iped)=SQRT( 1-u(iped)+.5*w(iped) )
     endif
     ! Terceiro passo do algoritmo, inicia-se o loop
     do i=iped+1,nlev(eff)-ngrp(eff)
        if (ped(eff,i,2) >= 0 .and. ped(eff,i,3) < 0) then    ! somente o pai é incerto
           if (ped(eff,i,2) >= iped) ts=v(ped(eff,i,2))
           do l=1,2*ms(p_ms(i))-1,2                           ! (,2)= colocar para iniciar na 2 coluna!
	      if ( int(ms(p_ms(i)+l)) >= iped ) td=td+(v(int(ms(p_ms(i)+l)))*ms(p_ms(i)+l+1))
           enddo
          else 
          if (ped(eff,i,2) >= iped) ts=v(ped(eff,i,2)) 
          if (ped(eff,i,3) >= iped) td=v(ped(eff,i,3)) 
       endif
       w(i)=w(i)+(ts*td)
       v(i)=.5*(ts+td)
       ts=0     
       td=0
     enddo
      do i=iped,nlev(eff)-ngrp(eff) 
        u(i)=u(i)+(v(i)**2)
        dinv(i)=(1/(v(i)**2))
     enddo
 enddo
! VJ - modified 
 print*, " "
 print*, "Dinv elements calculated"
 ! Write file with elements for all animals in pedigree
 open(io_dinv,file='dinv_elements', status='replace')
 write(io_dinv,"(2X,A,10X,A,9X,A,10X,A,9X,A)"), "u","v","w","dinv","p_ms"
 do i=1,nlev(eff)-ngrp(eff)
    ! VJ - print on screen
    !if(i .lt. 11) write (*,"(3i3,4f10.4,i3)") ped(eff,i,:), u(i), v(i), w(i), dinv(i), p_ms(i)
    ! VJ - Save on file
    write (io_dinv,103) u(i), v(i), w(i), dinv(i), p_ms(i)
    103 format (4f10.4,i3)
 enddo

 close(io_dinv)
 end subroutine

! VJ Subroutine adapted from blupf90.f90 to compute yhat, adjusted y (y*) and residuals
  subroutine write_adjusted_data
  integer::i,j,k,l
  character (40):: form
  real :: y2(ntrait),&
          y3(ntrait)!,&       ! variable to compute y - Xb = Zu + e
!          yhat

  write(*,*) 'Writing yhat values. Please wait...'
  form='(100f14.6)'
  rewind io_d
  open(io_s+1,file='yhat_residual')
  !write(io_s+1,'(''yhat adjusted_y residual '')')
  do irec=1,nrec
    read(io_d,*,iostat=io) indata(irec,:)
    call decode_record
    call find_addresses

    y2=y; y3=y
    do k=1,ntrait
       if (y(k) /= miss) then
          do i=1,neff
            y(k) = y(k) - weight_cov(i,k)*sol(address(i,k))
	        ! Use only fixed effects solutions to adjust phenotype
            if ( randomtype(i) .eq. g_fixed ) y3(k) = y3(k)-weight_cov(i,k)*sol(address(i,k))
          enddo
       endif
    enddo
    write(io_s+1, form) y2-y, y3, y
  enddo

  print*,'wrote yhat, adjusted_y and residual in file "yhat_residual"'
  print*,''
  end subroutine write_adjusted_data


subroutine add_g_norm_prior1(eff)
! VJ subroutine to add contributions on LHS and RHS.
! It is assumed that the effects are non-correlated
! 05-2016

  integer   :: eff,i,j,k,l,m,n,row,col,t1,t2,type=1,maxrow=0,maxcol=0,maxdiag=0
  integer   :: upper=1 !,lower=2
  real (rh) :: val,x,ident(nlev(eff)),sol(nlev(eff)),xvb,b(nlev(eff)*ntrait)
  type (sparse_hashm) :: a_usr
  type (sparse_ija)   :: a_usr_ija
  
  if (round == start) then
     !call init(a_usr)
     call zerom(a_usr,nlev(eff)*ntrait)
     n=0

     print*
     print*, 'Adding contribution due to normal prior effect'
     do
        read(io_off+eff,*,iostat=io)t1,col,x,xvb
        if (io /= 0) exit
        n=n+1

        ! Evaluate if the file is not sorted correctelly
        if (t1 .gt. type) upper=1
        if ((type .gt. t1) .or. (col .ne. upper))then
          print*,''
          print*,'Mixed diagonals in norm_prior file:',t1,col,x
          print*,'Sort levels within traits in ascending order.'
          stop
        endif
        upper=upper+1
        
        ! Add prior values
        m=address1(eff,col,t1)
        call addm(1/xvb,m,m,xx)
        xy(m)=xy(m)+x*(1/xvb)
     enddo

     print*,'norm_prior: read ',n,' elements'
     if (n == 0) then
       print*,'User defined file for effect ',eff,' empty'
       stop
     endif
  endif
end subroutine


! VJ Subroutine adapted from CBLUPTHR

subroutine threshold
  ! calculates pseudo R and Y, as in Quaas, BIF Guidelines
 
  integer,parameter::mxc=10 !Maximum number of categorories

  real(8),parameter::xpmin=1e-15
  real(8)::xp(mxc+1),eta(mxc+1),xint(nthr+2),xdis(mxc+1),p(mxc+1),l(mxc+1) 
  real(8)::t(mxc+1,mxc+1)  
  real (8) :: xl1,xl2,y2      
  real(8)::xmi,bea,sige,sig2,w,v,xrhs2
  real(8)::zik,zik1
  integer:: cat,nco
  real (8), allocatable :: stc(:),alpha(:)
  logical,save::first=.true.
  real (r8),dimension(ntrait)::  e

  t=0
  
  if (first) then
      sol(neq+1:neq+nthr)=(/((i-1.)/nthr,i=1,nthr)/) ! set priors for thresholds
      if (r(1,1) /= 1.0) then
         print*,'The residual variance for the categorical trait must be 1!'
         stop
      endif  
      first=.false.
  endif

  xmi=0
  xmi2=0

  do i=1,neff
    do j=1,ntrait
      xmi2(j)=xmi2(j)+(weight_cov(i,j)*sol(address(i,j)))
    enddo
  enddo 
   
 if (y(1) /= miss) then     
   do i=1,neff
     xmi=xmi+sol(address(i,1))*weight_cov(i,1)
   enddo   
   nco=0
   do i=2,ntrait
      if (y(i) /= miss) then
        nco=nco+1
        place(i)=nco
      else
        place(i)=0          
      endif
   enddo

   allocate (rc(nco,nco),stc(nco),alpha(nco))
   if (.not. allocated(rci)) allocate(rci(nco,nco))

   rci=0
   alpha=0
   stc=0
   rc=0

   do i=2,ntrait
     if (place(i)==0) cycle
     stc(place(i))=r(1,i)
     do j=2,ntrait
       if (place(j)==0) cycle
       rc(place(i),place(j))=r(i,j)
     enddo             
   enddo

   alpha=fsolve_s(rc,stc)
   rci=finverse_s(rc)

   do i=2,ntrait
     if (y(i)/=miss) then
       b(i)=alpha(place(i))
     else
       b(i)=0 
     endif
   enddo

   if (fthr == 1) then
     b=0   ! treat traits as uncorrelated at round 1 for better convergence
   endif

   do i=2,ntrait
     if (y(i) /= miss) then
       bea=y(i)-sum(sol(address(:,i))*weight_cov(:,i))
       bea=bea*b(i)
       xmi=xmi+bea  
     else
       bea=0
       b(i)=0      
     endif
   enddo
 endif

 sig2=0
 do i=2,ntrait
   do j=2,ntrait
     sig2=sig2+b(i)*b(j)*r(i,j)
   enddo
 enddo
  
 sige=sqrt(1-sig2)
       
 if (y(1) /= miss) then
   xint(1)=0
   xint(nthr+2)=1 
   xdis(1)=0
   xdis(nthr+2)=0 
   do j=2,nthr+1   
     eta(j)=(sol(neq+j-1)-xmi)/sige 
     xint(j)=normalcdf(eta(j))
     xdis(j)=normal(eta(j))
   enddo
    
   do j=2,nthr+2   
     p(j-1)=0   
     xp(j)=xint(j)-xint(j-1)       
     if (xp(j).lt.xpmin) xp(j)=xpmin
   enddo
 
   cat=y(1) 
   if (cat.gt.(nthr+1) .or. cat.lt.0) then
   	 print*,'category from record', irec, 'too small or large:', cat
   	 stop
   endif

   v=-(-xdis(cat+1)+xdis(cat))/(xp(cat+1)*sige)
   if (cat==1) then             
     p(1)=xdis(2)/(xp(2)*sige)    
   elseif (cat==nthr+1) then          
	 p(nthr)=-xdis(nthr+1)/(xp(nthr+2)*sige)    
   else                                       
	 p(cat-1)=-xdis(cat)/(xp(cat+1)*sige)      
	 p(cat)=xdis(cat+1)/(xp(cat+1)*sige)      
   endif                                
	  
	  w=0   
	  do j=2,nthr+2
	    w=w+(xdis(j-1)-xdis(j))**2/xp(j)/sige**2    
	  enddo
	  if (w < xpmin) w=xpmin

      do j=2,nthr+1                       
        xl1=(xdis(j)-xdis(j-1))/xp(j)       
        xl2=(xdis(j+1)-xdis(j))/xp(j+1)   
        l(j)=-xdis(j)*(xl1-xl2)/sige**2   
      enddo                           
  
      do j=2,nthr+1                                                  
        t(j,j)=(xp(j)+xp(j+1))/(xp(j)*xp(j+1)*sige**2)*xdis(j)**2    
        k=j+1                  
      if (k>nthr+1) cycle              
      t(j,k)=-xdis(j)*xdis(k)/(xp(k)*sige**2)                      
      t(k,j)=t(j,k)                                                
      enddo                  
     
      do i=1,nthr   
        weight_cov(neff+i,1)=l(i+1)/w 
      enddo         
   
      y2=0
      lt=0
      do i=1,nthr
        lt=lt+weight_cov(neff+i,1)*sol(neq+i) 
      enddo 
  
      y2=xmi2(1)-v/w    
      y(1)=y2
    
      do i=1,nthr                 
         xrhs(i)=xrhs(i)+p(i)  
         do j=1,nthr                       
           xrhs(i)=xrhs(i)+t(i+1,j+1)*sol(neq+j)
           xlhs(i,j)=xlhs(i,j)+t(i+1,j+1)-weight_cov(neff+i,1) &  
                            *weight_cov(neff+j,1)*w 
         enddo                
      enddo                 
  endif

  xrhs2=0

  do i=1,nthr
    xrhs2=xrhs2+weight_cov(neff+i,1)*w*sol(neq+i)
  enddo

  bwb=0

  rinv(1,1)=w
  do i=2,ntrait
    rinv(1,i)=-b(i)*w
    rinv(i,1)=rinv(1,i)
    do j=i,ntrait
      if (rci(place(i),place(j))==0) cycle
      bwb(i,j)=b(i)*b(j)*w
      rinv(i,j)=bwb(i,j)+rci(place(i),place(j))
      rinv(j,i)=rinv(i,j)  
    enddo
  enddo
 
  do i=1,ntrait
    if (y(i) == miss) then
      rinv(i,:)=0
      rinv(:,i)=0
    endif
  enddo  

  deallocate(rc,stc,alpha)

end subroutine threshold


! VJ Subroutine adapted from CBLUPTHR

subroutine setup_eq_threshold

!set addresses of thresholds and weigths
address=0; weight_cov=0
address(neff+1:neff+nthr,1)=(/(i,i=neq+1,neq+nthr)/)
address(neff+1:neff+nthr,2:ntrait)=1


do fthr=1,rthr

 rewind io_d

 call setup_g                    ! invert G matrices
 
 call zerom(xx,neq+NTHR); xy=0

 do irec=1,nrec
   read(io_d,*,iostat=io) indata(irec,:)
   if (io.ne.0) then
       print*,'total number of record in data set: ', irec-1, ' is not equal &
       &  to the value you have entered: ', nrec
       stop
   endif

   call decode_record
   call find_addresses
   !call find_rinv

   !  Threshold-model computations of pseudo RINV and Y are done here
   !call threshold
   call threshold_iod

   do i=1,neff+NTHR
       do j=1,neff+NTHR
          do k=1,ntrait
             do l=1,ntrait
                 val=weight_cov(i,k)*weight_cov(j,l)*sqrt(weight_y(k))*sqrt(weight_y(l))*rinv(k,l)
                 call addm(val,address(i,k),address(j,l),xx)
             enddo
          enddo     
       enddo
       do k=1,ntrait
          if (k==1) then
             if (i>neff) then
               do l=1,ntrait
                 xy(address(i,k))=xy(address(i,k))+rinv(k,l)*xmi2(l) & 
                                *weight_cov(i,k)*sqrt(weight_y(k))*sqrt(weight_y(l))
               enddo
             else 
               do l=1,ntrait
                 if (l==1) then          
                   xy(address(i,k))=xy(address(i,k))+rinv(k,l)*(y(l)+lt)&
                                *weight_cov(i,k)*sqrt(weight_y(k))*sqrt(weight_y(l))
                 else
                     xy(address(i,k))=xy(address(i,k))+rinv(k,l)*xmi2(l) &
                          *weight_cov(i,k)*sqrt(weight_y(k))*sqrt(weight_y(l))
                 endif
               enddo
             endif
          else  
            do l=1,ntrait
              if (l==1) then
                xy(address(i,k))=xy(address(i,k))+rinv(k,l)*(y(l)+lt) &
                        *weight_cov(i,k)*sqrt(weight_y(k))*sqrt(weight_y(l))
              else 
                xy(address(i,k))=xy(address(i,k))+(rci(place(k),place(l))*y(l) &
                 +bwb(k,l)*xmi2(l))*weight_cov(i,k)*sqrt(weight_y(k))*sqrt(weight_y(l))
              endif 
            enddo 
          endif
       enddo     
   enddo
  enddo

  ! Random effects' contributions
  do i=1,neff
     select case (randomtype(i))
       case (g_fixed)
          continue                ! fixed effect, do nothing
       case (g_diag)
          call add_g_diag(i)
       case (g_A, g_As, g_A_UPG,g_A_UPG_INB)
          rewind i+io_off
          call add_g_add(randomtype(i),i)
       case (g_A_MS)    ! FC add average A to parent equations
          call add_g_aav(i,maxnsires)
       case (g_A_MB, g_diag_MB)  ! FC add mbreed Ginv to xx
          if(max_msires==0) then
            if(nbreed(i)==1) then
                  rewind i+io_off
                  call add_g_add(g_A,i)
            else
              call add_g_mb(i,maxnsires,randomnumb(i)*ntrait)
            endif
          else
            call add_g_mbs(i,maxnsires,randomnumb(i)*ntrait)
          endif
       case (g_usr)
          rewind i+io_off
          call add_g_usr(i,'')
       case (g_usr_inv)
          rewind i+io_off
          call add_g_usr(i,'inv')
       case(g_normal)                ! VJ add contributions due to normal prior effects
            rewind i+io_off
            call add_g_norm_prior1(i)
       case default
         print*,'unimplemented random type',randomtype(i)
     endselect
  enddo
 
  ! add modifications to threshold parts: XLHS AND XRHS
  do i=1,nthr
     xy(i+neq)=xy(i+neq)+xrhs(i)
     xrhs(i)=0
     do j=1,nthr
         call addm(xlhs(i,j),i+neq,j+neq,xx)
         xlhs(i,j)=0
     enddo
  enddo

   if(solv_method == 'SOR') then
      call default_iter(conv=conv_crit,maxround=maxround_thr,relax=1.3,zerosol=.false.)
      call solve_iterm(xx,xy,sol)
   elseif(solv_method == 'PCG') then
      call default_iter(conv=conv_crit,maxround=maxround_thr,relax=1.3,zerosol=.false.)
      call solve_pcg(xx,xy,sol,blksize)
   else
      print*,'Wrong solve_method ',solv_method
      stop
   endif

   print*,'round ',fthr 
   print*,'thresholds:'
   do i=1,nthr      
     print*,sol(neq+i)    
     write(io_thr,fmt='(2i5,f15.7)') fthr,i,sol(neq+i)    
   enddo        

   call setup_g

end do

if (neq <= 8) then
  print*,''
  print*,'left hand side'
  call printm(xx)
  print  '( '' right hand side:'',10(/10f8.2))',xy
endif

end subroutine setup_eq_threshold



! Subroutine from CBLUPIOD

subroutine threshold_iod
  ! calculates pseudo R and Y, as in Quaas, BIF Guidelines
 
  integer,parameter::mxc=10 !Maximum number of categorories

  real(8),parameter::xpmin=1e-15
  real(8) :: xp(mxc+1),eta(mxc+1),xint(nthr+2),xdis(mxc+1),p(mxc+1),l(mxc+1) 
  real(8) :: t(mxc+1,mxc+1)  
  real(8) :: xl1,xl2,y2			
  real(8) :: xmi,bea,sige,sig2,w,v,xrhs2
  integer :: cat,nco
  real (8), allocatable :: stc(:),alpha(:),rc(:,:)
  logical,save :: first=.true.

  t=0;w=0.
  
!if(nrec==15525) call printmat(rci,'RCI')

  if (first) then
  	  sol(neq+1:neq+nthr)=(/((i-1.)/nthr,i=1,nthr)/) ! set priors for thresholds
      if (r(1,1) /= 1.0) then
         print*,'The residual variance for the categorical trait must be 1!'
	 stop
      endif	 
      first=.false.
  endif

  xmi=0
  xmi2=0

  do i=1,neff
    do j=1,ntrait
      xmi2(j)=xmi2(j)+(weight_cov(i,j)*sol(address(i,j)))
    enddo
  enddo 
   
 if (y(1) /= miss) then
     
     do i=1,neff
       xmi=xmi+sol(address(i,1))*weight_cov(i,1)
     enddo
   
   nco=0;place=0
   do i=2,ntrait
      if (y(i) /= miss) then
        nco=nco+1
	      place(i)=nco
      else
        place(i)=0          
      endif
   enddo

if(nco > 0) then
!   allocate (rc(nco,nco),stc(nco),alpha(nco),rci(nco,nco))
  
   if (.not. allocated(rc)) allocate(rc(nco,nco))
   if (.not. allocated(stc)) allocate(stc(nco))
   if (.not. allocated(alpha)) allocate(alpha(nco))
   if (.not. allocated(rci)) allocate(rci(0:nco,0:nco))

   rci=0
   alpha=0
   stc=0
   rc=0
 
   do i=2,ntrait
      if (place(i)/=0) then
        stc(place(i))=r(1,i)
        do j=2,ntrait
          if (place(j)/=0) then
            rc(place(i),place(j))=r(i,j)
          endif
        enddo             
      endif
   enddo

   alpha=fsolve_s(rc,stc)

   rci(1:nco,1:nco)=finverse_s(rc(1:nco,1:nco))

   do i=2,ntrait
     if (y(i)/=miss) then
       b(i)=alpha(place(i))
     else
       b(i)=0 
     endif
   enddo

   !if (round < 4 .and. ncont < 1) then
   if (fthr.eq.1) then
      b=0   ! treat traits as uncorrelated at round 1 for better convergence
   endif

   do i=2,ntrait
      if (y(i) /= miss) then
	    bea=y(i)-sum(sol(address(:,i))*weight_cov(:,i))
	    bea=bea*b(i)
	    xmi=xmi+bea  
	  else
            bea=0
            b(i)=0      
      endif
   enddo

   deallocate (rc,stc,alpha)
endif
   
endif

!if(nrec==15525) call printmat(r,'R')
 sig2=0
     do i=2,ntrait
       do j=2,ntrait
        sig2=sig2+b(i)*b(j)*r(i,j)
       enddo
    enddo
  
  sige=sqrt(1-sig2)
       
  if (y(1) /= miss) then
      xint(1)=0
      xint(nthr+2)=1 
      xdis(1)=0
      xdis(nthr+2)=0 
      do j=2,nthr+1   
        eta(j)=(sol(neq+j-1)-xmi)/sige 
	    xint(j)=normalcdf(eta(j))
        xdis(j)=normal(eta(j))
      enddo
      do j=2,nthr+2   
        p(j-1)=0   
        xp(j)=xint(j)-xint(j-1)       
	if (xp(j).lt.xpmin) xp(j)=xpmin
      enddo
 
      cat=y(1) 
      if (cat.gt.(nthr+1) .or. cat.lt.0)  print*,'category too small or large:', nthr,y(1)
     
      v=-(-xdis(cat+1)+xdis(cat))/(xp(cat+1)*sige)
      if (cat==1) then             
          p(1)=xdis(2)/(xp(2)*sige)    
        elseif (cat==nthr+1) then          
          p(nthr)=-xdis(nthr+1)/(xp(nthr+2)*sige)    
        else                                       
	  p(cat-1)=-xdis(cat)/(xp(cat+1)*sige)      
          p(cat)=xdis(cat+1)/(xp(cat+1)*sige)      
      endif                                
      
    w=0   
    do j=2,nthr+2
      w=w+(xdis(j-1)-xdis(j))**2/xp(j)/sige**2    
    enddo
      if (w < xpmin) w=xpmin

      do j=2,nthr+1                       
        xl1=(xdis(j)-xdis(j-1))/xp(j)       
	xl2=(xdis(j+1)-xdis(j))/xp(j+1)   
	l(j)=-xdis(j)*(xl1-xl2)/sige**2   
      enddo	                          
  
      do j=2,nthr+1                                                  
        t(j,j)=(xp(j)+xp(j+1))/(xp(j)*xp(j+1)*sige**2)*xdis(j)**2    
        k=j+1							     
	if (k>nthr+1) cycle					     
	t(j,k)=-xdis(j)*xdis(k)/(xp(k)*sige**2)                      
	t(k,j)=t(j,k)                                                
      enddo							     
     
      do i=1,nthr		
 	weight_cov(neff+i,1)=l(i+1)/w	
      enddo					
   
    y2=0
    lt=0
    do i=1,nthr
     lt=lt+weight_cov(neff+i,1)*sol(neq+i) 
    enddo 
  
     y2=xmi2(1)-v/w		
     y(1)=y2
    
     do i=1,nthr							    
         xrhs(i)=xrhs(i)+p(i)  
         do j=1,nthr							         
	 xrhs(i)=xrhs(i)+t(i+1,j+1)*sol(neq+j)
       	 xlhs(i,j)=xlhs(i,j)+t(i+1,j+1)-weight_cov(neff+i,1) &	
	  	                      *weight_cov(neff+j,1)*w	
       enddo								
     enddo  								
  endif

xrhs2=0

 do i=1,nthr
  xrhs2=xrhs2+weight_cov(neff+i,1)*w*sol(neq+i)
 enddo

bwb=0
     rinv(1,1)=w
     do i=2,ntrait
        rinv(1,i)=-b(i)*w
        rinv(i,1)=rinv(1,i)
        do j=i,ntrait
	      if (rci(place(i),place(j))/=0) then
            !if (place(i)/=0.and.place(j)/=0) then
            !if(nrec==15525) print*,i,j,place(i),place(j),rci(place(i),place(j))
	        bwb(i,j)=b(i)*b(j)*w
	        rinv(i,j)=bwb(i,j)+rci(place(i),place(j))
            rinv(j,i)=rinv(i,j)	 
          endif
        enddo
     enddo

   do i=1,ntrait
     if (y(i) == miss) then
       rinv(i,:)=0
       rinv(:,i)=0
     endif
   enddo  

end subroutine threshold_iod


subroutine copyright
implicit none
!integer :: un=99999999
!open(unit=un,file='copyright')
write(*,'(a)') ' '
write(*,'(a)') '================================================================================='
write(*,'(a)') 'Copyright (C) 2008, 2012, 2016 Fernando Cardoso, Vinicius Junqueira, Rodrigo Mota'
write(*,'(a)') ' '
write(*,'(a)') '    This program is free software: you can redistribute it and/or modify'
write(*,'(a)') '    it under the terms of the GNU General Public License as published by'
write(*,'(a)') '    the Free Software Foundation, either version 3 of the License, or'
write(*,'(a)') '    (at your option) any later version.'
write(*,'(a)') ' '
write(*,'(a)') '    This program is distributed in the hope that it will be useful,'
write(*,'(a)') '    but WITHOUT ANY WARRANTY; without even the implied warranty of'
write(*,'(a)') '    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the'
write(*,'(a)') '    GNU General Public License for more details.'
write(*,'(a)') ' '
write(*,'(a)') '    You should have received a copy of the GNU General Public License'
write(*,'(a)') '    along with this program.  If not, see <http://www.gnu.org/licenses/>.'
write(*,'(a)') '================================================================================='
write(*,'(a)') ' '
!close(un)
end subroutine


end program intergen

