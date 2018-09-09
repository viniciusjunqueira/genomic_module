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

! blups1.f90
! Modules and subroutines to read the parameter file and specify the model
! FC, 1/1/2001 - 9/1/2007 - multisires, breeds and reaction norms modifications
! IM, 2/23/2000 - 2/24/2000


module textop
   use model
! contains subroutines for
!       - reading and writing parameters
!       - for decomposing one line into alphanumeric and numeric 
!         components separated by spaces
!       - reading optional parameters from parameter files
! FC, 1/1/2001 - 9/1/2007 - multisires, breeds and reaction norms modifications
! IM, 2/23/2000 - 12/15/2002
   integer :: max_string_readline = 800 

  CONTAINS

  subroutine read_parameters()
  ! reads and decodes parameter file
   implicit none
   integer n, x(40),ef,stat,curr_eff,i,j,k
   character xc(40)*40
   real xf(40)                            ! VJ - variable to receive convergence criteria value
        

   write(*,'(a)',advance='no')'  name of parameter file? '
   read '(a)',parfile
   parfile=adjustl(parfile)               ! ignore leading spaces
   !write(*,'(/a)') trim(parfile)
   open(io_p,file=parfile)

   call chkfmtnewline(io_p)
 
   call readline(io_p,x,xc,n)
   if (n == -1) then
       print*,'INPUT FILE EMPTY'
       stop
   endif
! FC - HGC - VJ - MCMC and BLUP information
   if (xc(1) /= 'MCMC_CHAIN:' .and. xc(1) /= 'BLUP:') then
       print*,'MCMC_CHAIN or BLUB: not found'
       stop
   endif
   ! FC - HGC - parameters for MCMC only
   if (xc(1) .eq. 'MCMC_CHAIN:') then
      call readline(io_p,x,xc,n)
      nround=x(1)
      nburn=x(2)
      skip=x(3)

   ! FC SEED for random number generation 
      call readline(io_p,x,xc,n)
      if (xc(1) /= 'SEED') then
         print*,'SEED not found'
         stop
      endif
      call readline(io_p,x,xc,n) 
      seed=x(1)

   ! FC allows restart the chain
      call readline(io_p,x,xc,n)
      if (xc(1) /= 'RESTART:') then
         print*,'RESTART: not found'
         stop
      endif
      call readline(io_p,x,xc,n) 
      rstart=xc(1)
      if (rstart == 'y') then
         start=x(2)
      else 
         start=1
      endif
   else
      call readlinef(io_p,xf,xc,n) ! VJ - function that return x as floating
      nround=0
      ! VJ - changes default values if present on parameter file
      if (n > 0) then
         solv_method=xc(1)
         if(xf(2).gt.0) conv_crit=xf(2)
         sol_se=xc(3)
      endif      
   endif
  ! FC
   call readline(io_p,x,xc,n)
   if (xc(1) /= 'DATAFILE') then
       print*,'DATAFILE not found'
       stop
   endif
   call readline(io_p,x,xc,n) 
   datafile=xc(1)
   nrec=x(2)          !FC read the number of records

   call readline(io_p,x,xc,n)
   if (xc(1) /= 'NUMBER_OF_TRAITS') then
       print*,'NUMBER_OF_TRAITS not found'
       stop
   endif
   call readline(io_p,x,xc,n) 
   ntrait=x(1)

   call readline(io_p,x,xc,n)
   if (xc(1) /= 'NUMBER_OF_EFFECTS') then
       print*,'NUMBER_OF_EFFECTS not found'
       stop
   endif
   call readline(io_p,x,xc,n) 
   neff=x(1)

   allocate(pos_weight(ntrait),pos_y(ntrait),pos_eff(neff,ntrait),&
            nlev(neff),effecttype(neff),effectsave(neff),nestedcov(neff,ntrait),&
            randomtype(neff),randomnumb(neff),randomfile(neff),dbelief_g(neff),&
            r(ntrait,ntrait),rinv(ntrait,ntrait),g(neff,maxcorr,maxcorr),&
            lenped(neff),ngrp(neff),nbreed(neff))
   
   randomtype=g_fixed
   randomnumb=0
   nestedcov=0
   effectsave=0
   lenped=0; ngrp=0; nbreed=0

   call readline(io_p,x,xc,n)
   if (xc(1) /= 'OBSERVATION(S)') then
       print*,' OBSERVATION(S) not found'
       stop
   endif
   call readline(io_p,x,xc,n) 
   if (n == ntrait) then
        pos_y(1:ntrait)=x(1:ntrait)
       else
         print*,ntrait,' numbers expected after OBSERVATION(S)'
	 stop
   endif

   call readline(io_p,x,xc,n)
   if (xc(1) /= 'WEIGHT(S)') then
       print*,' WEIGHT(S) not found'
       stop
   endif
   call readline(io_p,x,xc,n) 
! FC trait specific weights
   if (n == ntrait) then
        pos_weight(1:ntrait)=x(1:ntrait)
       elseif (n == 0) then
        pos_weight(1:ntrait)=0
       else
         print*,'None or ', ntrait,' numbers expected after WEIGHT(S)'
	 stop
   endif

   call readline(io_p,x,xc,n)
   if (xc(1) /= 'EFFECTS:') then
       print*,' EFFECTS: not found'
       stop
   endif
   do ef=1,neff
      call readline(io_p,x,xc,n)
      if (n == ntrait+3 .or. n == 2*ntrait+3 .or. n == ntrait+4 .or. n == 2*ntrait+4) then
          pos_eff(ef,:)=x(1:ntrait)
          nlev(ef)=x(ntrait+1)
          select case (xc(ntrait+2)(1:3))
               case ('cro')         !crossclassified
                  effecttype(ef)=effcross
               case ('cov')         !covariate
                  effecttype(ef)=effcov
                  if (n==2*ntrait+2) then
                       nestedcov(ef,:)=x(ntrait+3:n)
                  elseif (n==2*ntrait+3) then
                       nestedcov(ef,:)=x(ntrait+4:n)
                  endif
               case ('ram')         ! FC reduced animal model
                  effecttype(ef)=effram
               case ('unk')         ! FC unknown covariable
                  effecttype(ef)=effuncov
               case ('rno')         ! FC reaction norm
                  effecttype(ef)=effrnorm
                  if (n==2*ntrait+2) then
                       nestedcov(ef,:)=x(ntrait+3:n)
                  elseif (n==2*ntrait+3) then
                       nestedcov(ef,:)=x(ntrait+4:n)
                  endif
               case default
                  print*,'Unknown type of effect: ',xc(ntrait+2)
          end select
        if (xc(ntrait+3)(1:1) == 'y') effectsave(ef)=1
        else
          print*,'too few or too many numbers for effect ',ef,'. Was ',n
          stop
      endif 
   enddo

! FC allow for heterogeneous residuals
   call readline(io_p,x,xc,n)
   if (xc(1) /= 'RANDOM_RESIDUAL:') then
       print*,' RANDOM_RESIDUAL: not found'
       stop
   endif

   call readline(io_p,x,xc,n) 
   select case (xc(1)(1:16))
     case ('homogeneous')
        residualtype=r_homo
     case ('structural')
        residualtype=r_structural ! structural heterogeneous variance
     case ('slash')
        residualtype=r_slash 	  ! robust slash errors
     case ('struct_slash')
        residualtype=r_struct_slash ! structural heterogeneous variance and slash errors
     case ('student_t')
        residualtype=r_studentt 	 ! robust student t errors
     case ('struct_student_t')
        residualtype=r_struct_studentt ! structural heterogeneous variance and student t errors
     case default
        print*,'unknown RANDOM_RESIDUAL-TYPE: ', xc(1)
        print*,' Should be homogeneous, structural, slash, struct_slash, student_t, struct_student_t'
        stop
   end select  
   dbelief_r=x(2)

! VJ - Stops if is BLUP with heterogeneous residual
   if (nround .eq. 0) then
      if (residualtype .ne. r_homo ) then
          print*,'ERROR: only homoskedastic residuals are possible for BLUP'
          stop
      endif
   endif

  if (nround .eq. 0 .and. any(effecttype == effrnorm)) then
        print*,'ERROR: BLUP not implemented for reaction norms with unknown covariates'
        stop
  endif

 ! FC for structural residual effects 
     
	if (residualtype==r_structural .or. residualtype==r_struct_slash &
	    .or. residualtype==r_struct_studentt) then

	   call readline(io_p,x,xc,n)
	   if (xc(1) /= 'METROPOLIS_STEP_OF_STRUCTURAL_EFFECTS:') then
		   print*,'METROPOLIS_STEP_OF_STRUCTURAL_EFFECTS: not found'
		   stop
	   endif
	   call readline(io_p,x,xc,n) 
	   r_nround=x(1)
       r_skip=x(2)

	   call readline(io_p,x,xc,n)
	   if (xc(1) /= 'NUMBER_OF_STRUCTURAL_EFFECTS') then
		   print*,'NUMBER_OF_STRUCTURAL_EFFECTS not found'
		   stop
	   endif
	   call readline(io_p,x,xc,n) 
	   neffr=x(1)

	   allocate(pos_effr(neffr),effrsave(neffr))

	   call readline(io_p,x,xc,n)
	   if (xc(1) /= 'STRUCTURAL_EFFECTS:') then
		   print*,' STRUCTURAL_EFFECTS: not found'
		   stop
	   endif
	   do ef=1,neffr
		  call readline(io_p,x,xc,n)
		  if (n == 1 .or. n == 2) then
		      pos_effr(ef)=x(1)
			  if (xc(2)(1:1) == 'y') effrsave(ef)=1
			else
			  print*,'too few or too many numbers for effect ',ef,'. Was ',n
			  stop
		  endif 
	   enddo
     endif
! FC

   call readline(io_p,x,xc,n)
   if (xc(1) /= 'RESIDUAL_PRIOR_(CO)VARIANCES') then
      print*,' RESIDUAL_PRIOR_(CO)VARIANCES not found'
      stop
   endif

   read(io_p,*,iostat=stat)((r(i,j),i=1,ntrait),j=1,ntrait)
   if (stat.ne.0) then
      print*,'error reading RANDOM_RESIDUAL variances'
      stop
   endif

   do
       call readline(io_p,x,xc,n)
          if (n == -1) exit
          if (xc(1) /='RANDOM_GROUP') exit

          call readline(io_p,x,xc,n)
          if (n<1 .or. n>neff) then
             print*,'line after RANDOM_GROUP has too few or too many numbers'
             stop
          endif
          if (x(1)<1 .or. x(n)> neff) then
              print*,'number ', x(1:n),' after RANDOM_GROUP out of range'
              stop
          endif
          do i=1,n
             if (i+x(1)-1 /= x(i)) then
                print*,' correlated effects:',x(1:n),' should be consecutive'
                stop
             endif
          enddo
          curr_eff=x(1)
          randomnumb(curr_eff)=n
          call readline(io_p,x,xc,n)
          if (xc(1) /= 'RANDOM_TYPE') then
              print*,' RANDOM_TYPE not found'
              stop
          endif

          call readline(io_p,x,xc,n)
          select case (xc(1)(1:11))
             case('diagonal')
                randomtype(curr_eff)=g_diag
             case('add_animal')
                randomtype(curr_eff)=g_A
             case('add_an_upg')
                randomtype(curr_eff)=g_A_UPG
             case('add_an_upgi')
                randomtype(curr_eff)=g_A_UPG_INB
             case('add_sire')
                randomtype(curr_eff)=g_As
             case('par_domin')
                randomtype(curr_eff)=g_PD
             case('add_an_ms')	                   ! FC am/ram w/msires
                randomtype(curr_eff)=g_A_MS
             case('add_an_mb')	                   ! FC mbreed additive
                randomtype(curr_eff)=g_A_MB
             case('diag_mb')	                     ! FC ram diagonal
                randomtype(curr_eff)=g_diag_MB
             case('user_file')
                randomtype(curr_eff)=g_usr
             case('user_file_i')
                randomtype(curr_eff)=g_usr_inv
             case('norm_prior')
                randomtype(curr_eff)=g_normal      ! VJ normal prior
             case default
                print*,'unknown RANDOM-TYPE: ', xc(1)
                print*,' Should be diagonal, add_animal, add_sire',&
                       ' add_an_upg, add_an_upginb, or par_domin',&
                       ' add_an_ms, add_an_mb, diag_mb','norm_prior',&
		                   ' user_file, user_file_inv'
                stop
          end select  
          dbelief_g(curr_eff)=x(2)

          call readline(io_p,x,xc,n)
          if (xc(1) /= 'PEDIGREEFILE:') then
              print*,' PEDIGREEFILE for random effect ',curr_eff,' not found'
              stop
          endif
          call readline(io_p,x,xc,n)
          randomfile(curr_eff)=xc(1)
          if (randomtype(curr_eff) /= g_diag) open(curr_eff+io_off,file=xc(1))

! FC read extra info for ms/ram/mb/mtbs
          select case (randomtype(curr_eff))
             case (g_A,g_As,g_A_UPG)            ! FC - extra effect added
                lenped(curr_eff)=x(2)
             case (g_A_MS)
				        lenped(curr_eff)=x(2)
				        ngrp(curr_eff)=x(3)
                call readline(io_p,x,xc,n)
				        if (xc(1) /= 'MULTIPLE_SIRES:') then
				           print*,' MULTIPLE_SIRES: for random effect ',curr_eff,' not found'
				           stop
				        endif
				        call readline(io_p,x,xc,n)
				        max_msires=x(1)	                  ! maximum # of sires to infer upon paternity (1 -> average/0 -> all known)
				        msfile=xc(2)
				        n_ms=x(3)
				        if (max_msires > 1) ms_alpha=x(4)	! dirichlet hyperpar: 1 -> flat; 
													                        ! 0 -> no dirichlet;  < 0 -> user defined from ms vector 
             case (g_A_MB,g_diag_MB)
        				lenped(curr_eff)=x(2)
        				ngrp(curr_eff)=x(3)
        				nbreed(curr_eff)=x(4)
        				if (randomtype(curr_eff) == g_A_MB) then 
        					call readline(io_p,x,xc,n)
        					if (xc(1) /= 'MULTIPLE_SIRES:') then
        					  print*,' MULTIPLE_SIRES: for random effect ',curr_eff,' not found'
        					  stop
        					endif
        					call readline(io_p,x,xc,n)
        					max_msires=x(1)	                  ! maximum # of sires to infer upon paternity (1 -> average/0 -> all known)
        					msfile=xc(2)
        					n_ms=x(3)
        					if (max_msires > 1) ms_alpha=x(4)	! dirichlet hyperpar; 
          			endif
                ! FC metropolis step of multibreed model
      			    call readline(io_p,x,xc,n)
      			    if (xc(1)(1:20) /= 'METROPOLIS_STEP_OF_M') then
      				   print*,'METROPOLIS_STEP_OF_MULTIBREED_(CO)VARIANCES: not found'
      				   stop
      			    endif
        			  call readline(io_p,x,xc,n) 
      			    mb_nround=x(1)
      			    mb_skip=x(2)
! VJ normal prior
             case(g_normal)
                lenped(curr_eff)=x(2) ! dimension = number of traits * number of levels
             case default
                continue
          end select  
! VJ read (co)variances matrix only if random type is different of normal_prior effect
        if(randomtype(curr_eff) /= g_normal) then
! FC
          call readline(io_p,x,xc,n)
          if (xc(1) /= '(CO)VARIANCES') then
              print*,' (CO)VARIANCES not found'
              stop
          endif
		      if (randomtype(curr_eff)==g_A_MB .or. & 
		        randomtype(curr_eff)==g_diag_MB) then
			      k=randomnumb(curr_eff)*ntrait*nbreed(curr_eff)
			    else
			      k=randomnumb(curr_eff)*ntrait
		      endif
          if (k <= maxcorr) then
             read(io_p,*,iostat=stat)((g(curr_eff,i,j),i=1,k),j=1,k)
             if (stat /= 0) then
                print*,'error reading variances for effect ',curr_eff
                stop
             endif
          else
              print*,'maxcorr should be increased to at least ',k
          endif
        endif
   enddo

   ! Categorical parameters
   call getoption('cat',n,xf)
   if(n.ne.-1) then
     if(nround.gt.0) then   ! stops if MCMC with threshold model
      stop 'Error: Threshold models can only be perfomed on BLUP analysis.'      
     end if
     if(n.lt.2) then
       stop 'Parameters of threshold model: OPTION cat [number of categories] [number of rounds].'
     endif
     cat_value=n;nthr=xf(1)-1;rthr=xf(2)
     if(solv_method=='FSPAK') then
       stop 'Error: Threshold models solve method must be set as PCG or SOR.'       
     endif
     open(io_thr,file='thresholds', status='replace')
   else 
     if (sol_se.eq.'se') solv_method='FSPAK' ! VJ FSPAK if "se" is on to quantitative traits
   endif

end subroutine


subroutine readline(unit,x,xc,n,nohash)
! Reads a line from unit and decomposes it into numeric fields in x and
! alphanumeric fields in xc; the number of fields is stored in n.
! If optional variable nohash is absent, 
!  all characters after character # are ignored.

  integer unit,x(:),n,stat,i
  integer,optional::nohash
  character xc(:)*(*)
  character(len=max_string_readline) :: a   
 
  do
     read(unit,'(a)',iostat=stat)a
     if (stat /= 0) then
         n=-1
         return
     endif
     a=adjustl(a)               ! ignore leading spaces
     
     if (present(nohash)) then 
        exit
       else	!ignore characters after #
        i=scan(a,'#')
        if (i == 0) exit
        if (i>1) then
            a=a(1:i-1)
            exit
        endif
     endif  	
  enddo
  call nums(a,x,xc,size(x),n)
  return
end subroutine


subroutine readline1(unit,x,xc,n)
! as readline but does not ignore entries after character #

  integer unit,x(:),n,stat,i,nohash
  character xc(:)*(*)

  call readline(unit,x,xc,n,nohash)
end subroutine  



subroutine readlinef(unit,x,xc,n)
! As readline but decodes x to floating point

  integer unit,n,stat,i
  character (*)::xc(:)
  character(len=max_string_readline) :: a
  real::x(:)

  read(unit,'(a)',iostat=stat)a
  if (stat == 0) then
      a=adjustl(a)               ! ignore leading spaces
      call nums2(a,n,x,xc)
    else
      n=-1
  endif
end subroutine



subroutine print_parameters(title)
  ! Prints parameters of the model
   use denseop
   implicit none
   integer :: i,j,k,l,b1,b2
   character (20):: name
   character (40):: fmt
   character (*):: title
   logical::not_pd

   print*,'    '
   print*,'    ',title
   print '(/'' Parameter file:'',t30,a)',trim(parfile)
   ! VJ - Print only if MCMC will be perfomed
   if(nround .gt. 0) then
     ! FC MCMC information
     print '('' Method: '',t30, ''MCMC'')'
     print '('' Total number of MCMC cycles:'',t30,i8)',nround
     print '('' Burn-in period:'',t30,i8)',nburn
     print '('' Thinning interval:'',t30,i8)',skip
     print '('' Seed:'',t30,i8)',seed
     print '('' Is this a restart?:'',t30,a)',rstart
     if (rstart == 'y') print '(''Restarting cycle'',t30,i8)',start
   else
     ! VJ
     print '('' Method: '',t30, ''BLUP'')'                              ! VJ blup
     if(sol_se.eq.'se') then 
       print '('' Standard error: '',t30, ''Yes'')'
     else 
       print '('' Standard error: '',t30, ''No'')'
     end if
     if(cat_value.ne.-1) then
       print '('' Model: '',t30, ''Threshold'')'                        ! VJ threshold parameters
       print '('' Number of thresholds: '',t29,i2)', nthr
       print '('' Number of iteration rounds: '',t29,i3)', rthr         
     else
       print '('' Model: '',t30, ''Quantitative'')'
     endif
     print '('' Solve Method:'',t30,a)', solv_method           
     print '('' Convergence Criteria:'',t29,ES9.2)', conv_crit
   endif

! FC
   print '('' Data file:'',t30,a)',trim(datafile)
   print '('' Number of Traits'',t29,i2)',ntrait                          ! VJ Space modifications
   print '('' Number of Effects'',t29,i2)',neff                           ! VJ Space modifications
   print '('' Position of Observations'',t28,20i3)',(pos_y(i),i=1,ntrait) ! VJ Space modifications
   print '('' Position of Weights'',t28,20i3)',(pos_weight(i),i=1,ntrait) ! FC trait specific weights. VJ Space modifications
   print '('' Value of Missing Trait/Observation'',t40,i8)',miss

   print '(/,'' EFFECTS'')'
   print '(a,t25,a,t42,a)',' #  type','position (2)','levels  save  [positions for nested]'
   do i=1,neff
      select case (effecttype(i))
           case (effcross)
              name='cross-classified'
           case (effcov)
              name='covariable'
		! FC animal effect for reduced animal model
           case (effram)
              name='ram'
		! FC unknown covariable for reaction norm model
           case (effuncov)
              name='unknown-covariable'
		! FC reaction norm
           case (effrnorm)
              name='reaction-norm'
           case default
              name='???'
      end select
 
      ! writing with number of fields dependent on ntrait and presence of cov
      write(fmt,'(''(i2,2x,a19,'',i2,''i3,t37,i10,t51,i2,t55,10i3)'')')ntrait
      if (nestedcov(i,1) /= 0) then
           print fmt,i,name,pos_eff(i,:),nlev(i),effectsave(i),nestedcov(i,:)
          else
           print fmt,i,name,pos_eff(i,:),nlev(i),effectsave(i)
      endif
   enddo	

   if (any(effecttype == effram) .and. any(effecttype == effrnorm)) then
     print*,'ERROR: Joint implementation of RAM and reaction-norm effects not possible'
     stop
   endif

! FC Residuals
   select case (residualtype)
       case (r_homo)
            print '('' Type of Residuals:'',t30,a)','homogeneous'
       case (r_structural)
            print '('' Type of Residuals:'',t30,a)','structural'
            print '('' Metropolis steps within round:'',t40,i8)',r_nround
            print '('' Cycles to tune MH during burn-in:'',t40,i8)',r_skip
            print '('' Number structural effects:'',t40,i3)',neffr
            print '('' Structural effects:'',t40,50i3)',pos_effr(1:neffr)
            print '('' Structural effects save:'',t40,50i3)',effrsave(1:neffr)
       case (r_slash)
            print '('' Type of Residuals:'',t30,a)','slash'
       case (r_struct_slash)
            print '('' Type of Residuals:'',t30,a)','struct_slash'
            print '('' Metropolis steps within round:'',t30,i8)',r_nround
            print '('' Cycles to adjust metropolis variance during burn-in:'',t30,i8)',r_skip
            print '('' Number structural effects:'',t30,i3)',neffr
            print '('' Structural effects:'',t30,50i3)',pos_effr(1:neffr)
            print '('' Structural effects save:'',t30,50i3)',effrsave(1:neffr)
       case (r_studentt)
            print '('' Type of Residuals:'',t30,a)','student_t'
       case (r_struct_studentt)
            print '('' Type of Residuals:'',t30,a)','struct_student_t'
            print '('' Metropolis steps within round:'',t30,i8)',r_nround
            print '('' Cycles to adjust metropolis variance during burn-in:'',t30,i8)',r_skip
            print '('' Number structural effects:'',t30,i3)',neffr
            print '('' Structural effects:'',t30,50i3)',pos_effr(1:neffr)
            print '('' Structural effects save:'',t30,50i3)',effrsave(1:neffr)
       case default
            print*,'unknown RESIDUAL-TYPE'
            stop
   end select  
   if (nround .gt. 0) print '('' R prior degrees of belief:'',t30,i12)',dbelief_r
   if (ntrait > 1 .and. residualtype /= r_homo) then
     print*,'ERROR: only homoskedastic residuals are possible for multitrait models'
     stop
   endif
   if (.not. checksym(r) ) then
      print*,'R Matrix not symmetric'
      stop 
   endif 
   print '(/'' Residual (co)variance Matrix'')'
   call pos_def(r,'    Matrix not positive definite: corrected',1d-4)
   call printmat(r,fmt='(20g12.5)')

!   do i=1,ntrait
!      print 	       ,r(i,:)
!   enddo

   do i=1,neff
      if (randomnumb(i) == 0) then
          cycle
      else if (randomnumb(i) == 1) then
          print '(/,'' Random Effect(s)'',t20,10i2)',(/(j,j=i,i+randomnumb(i)-1)/)
      else if(randomnumb(i) > 1) then
           print '(/,'' correlated random effects'',t30,10i3)',&
                     (/(j,j=i,i-1+randomnumb(i)) /)
      endif
      select case (randomtype(i))
         case (g_diag)
            print '('' Type of Random Effect:'',t30,a)','diagonal'
            ! VJ - Print only in BLUP
            if(nround .gt. 0) print '('' Prior degrees of belief:'',t30,i12)',dbelief_g(i)
         case(g_A)
            print '('' Type of Random Effect:'',t30,a)','additive animal'
            print '('' Pedigree File:'',t30,a)',randomfile(i)
            ! VJ - Print only in BLUP
            if(nround .gt. 0) print '('' Prior degrees of belief:'',t30,i12)',dbelief_g(i)
         case(g_As)
            print '('' Type of Random Effect:'',t30,a)','additive sire'
            print '('' Pedigree File: '',t30,a)',randomfile(i)
            ! VJ - Print only in BLUP
            if(nround .gt. 0) print '('' Prior degrees of belief:'',t30,i12)',dbelief_g(i)
         case(g_A_UPG)
            print '('' Type of Random Effect:'',t30,a,a)','additive animal',&
                     ' with unknown parent groups'
            print '('' Pedigree File: '',t30,a)',randomfile(i)
            ! VJ - Print only in BLUP
            if(nround .gt. 0) print '('' Prior degrees of belief:'',t30,i12)',dbelief_g(i)
         case(g_A_UPG_INB)
            print '('' Type of Random Effect:'',t30,a,a)','additive animal',&
                      ' with unknown parent groups and inbreeding'
            print '('' Pedigree File: '',t30,a)',randomfile(i)
            ! VJ - Print only in BLUP
            if(nround .gt. 0) print '('' Prior degrees of belief:'',t30,i12)',dbelief_g(i)
         case(g_PD)
            print '('' Type of Random Effect:'',t30,a,a)','parental dominance'
            print '('' Pedigree File: '',t30,a)',randomfile(i)
	! FC ram w/ msires
         case(g_A_MS)
            print '('' Type of Random Effect:'',t30,a,a)','(reduced) animal model',&
                      ' with msires and parent groups'
                      ! VJ - Print only in BLUP
            if(nround .gt. 0) print '('' Prior degrees of belief:'',t30,i12)',dbelief_g(i)
            print '('' Pedigree File: '',t30,a)',randomfile(i)
            print '('' Number of Animals: '',t30,i12)', lenped(i)
            print '('' Number of groups: '',t30,i12)', ngrp(i)
            if (max_msires > 0) then
              print '('' Msires File: '',t30,a)', msfile 
              print '('' Dimension of Msfile: '',t30,i12)', n_ms 
              print '('' Max_N of Sires to infer: '',t30,i12)', max_msires 
              if (max_msires > 1) print '('' Dirichlel hyperparameters: '',t30,f12.4)', ms_alpha
			      endif
         case(g_A_MB)
            print '('' Type of Random Effect:'',t30,a,a)','multiple breed animal model',&
                      ' with additive effects'
            print '('' Metropolis steps within round:'',t30,i8)',mb_nround
            print '('' Cycles to adjust metropolis variance during burn-in:'',t30,i8)',mb_skip
            ! VJ - Print only in BLUP
            if(nround .gt. 0) print '('' Prior degrees of belief:'',t30,i12)',dbelief_g(i)
            print '('' Pedigree File: '',t30,a)',randomfile(i)
            print '('' Number of Animals: '',t30,i12)', lenped(i)
            print '('' Number of breeds: '',t30,i12)', nbreed(i)
            print '('' Number of groups: '',t30,i12)', ngrp(i)
            if (max_msires > 0) then
			        print '('' Msires File: '',t30,a)', msfile 
              print '('' Dimension of Msfile: '',t30,i12)', n_ms 
              print '('' Max_N of Sires to infer: '',t30,i12)', max_msires 
              if (max_msires > 1) print '('' Dirichlel hyperparameters: '',t30,f12.4)', ms_alpha
			      endif
         case(g_diag_MB)
            print '('' Type of Random Effect:'',t30,a,a)','multiple breed uncorrelated effect'
            print '('' Metropolis steps within round:'',t30,i8)',mb_nround
            print '('' Cycles to adjust metropolis variance during burn-in:'',t30,i8)',mb_skip
            ! VJ - Print only in BLUP
            if(nround .gt. 0) print '('' Prior degrees of belief:'',t30,i12)',dbelief_g(i)
            print '('' Pedigree File: '',t30,a)',randomfile(i)
            print '('' Number of Animals: '',t30,i12)', lenped(i)
            print '('' Number of breeds: '',t30,i12)', nbreed(i)
            print '('' Number of groups: '',t30,i12)', ngrp(i)
	! FC
	       case(g_usr)
            print '('' Type of Random Effect:'',t30,a,a)',&
	            'user defined from file'
            print '('' User File: '',t30,a)',randomfile(i)
	       case(g_usr_inv)
            print '('' Type of Random Effect:'',t30,a,a)',&
	           'user defined from file and inverted'
            print '('' User File: '',t30,a)',randomfile(i)
  ! VJ Normal prior effect
         case(g_normal)
            print '('' Type of Random Effect:'',t30,a,a)',&
             'normal prior from file'
            print '('' User File: '',t30,a)',randomfile(i)
         case default
            print*,'unknown RANDOM-TYPE'
            stop
      end select  

      ! VJ Evaluate only if the randomtype is different of g_normal
      if(randomtype(i) .ne. g_normal) then
        j=ntrait*randomnumb(i)
        if (.not. checksym(g(i,:j,:j) )) then 
          print*,'G Matrix not symmetric'
          stop
        endif       
        call pos_def(g(i,:j,:j),'    Matrix not positive definite: corrected',1d-4)
      endif

	  if (randomtype(i)==g_A_MB .or. & 
			randomtype(i)==g_diag_MB) then
		  print '('' breeds  trait   effect    (CO)VARIANCES'')'
		  do b1=1,nbreed(i)
			 do b2=b1,nbreed(i)
			  do k=0,randomnumb(i)-1
				 do j=1,ntrait
					 l=(b1-1)*randomnumb(i)*ntrait+k*ntrait+j
					 print '(x,2i3,2x,i4,i8,2x,20g12.5)',b1,b2,j,k+i, &
							   g(i,l,(b2-1)*randomnumb(i)*ntrait+1:b2*randomnumb(i)*ntrait)
				 enddo
			  enddo
			 enddo
		  enddo
	  else if(randomtype(i) .ne. g_normal) then
		  print '('' trait   effect    (CO)VARIANCES'')'
		  do k=0,randomnumb(i)-1
			 do j=1,ntrait
				 l=j+k*ntrait
				 print '(i3,i8,2x,20g12.5)',j,k+i,g(i,l,1:randomnumb(i)*ntrait)
			 enddo
		  enddo
	  endif
   enddo

   print '(/'' REMARKS'')'
   print*,' (1) Weight position 0 means no weights utilized'
   print*,' (2) Effect positions of 0 for some effects and traits means that such'
   print*,'     effects are missing for specified traits'
   print*
end subroutine



subroutine nums2(a,n,x,xc)
! separates array a into items delimited by blanks. character elements are
! put into optional character vector xc, decoded numeric values 
! into optional real vector x, and n contains the number of items. The 
! dimension of x and xc can be lower than n.
! A modification of nums() from f77 to f90
! Now accepts real numbers
! 2/23/2000

character (*)::a
character (*),optional::xc(:)
real,optional::x(:)
integer::n,curr,first,last,lena,stat,i

curr=1;lena=len(a);n=0
if (present(xc)) xc=' '
if (present(x)) x=0
do 
  ! search for first nonspace
  first=0
  do i=curr,lena
     if (a(i:i) /= ' ') then
        first=i
	exit
     endif
  enddo
  if (first == 0) exit
  
  
  ! search for first space
  curr=first+1
  last=0
  do i=curr,lena
      if (a(i:i) == ' ') then
        last=i
        exit
      endif
  enddo
  
  if (last == 0) last=lena
  
  n=n+1
  if (present(xc)) then
     if (size(xc) >= n) then
           xc(n)=a(first:last)

         else
              print*, "Error in nums2 splitting string starting: "
              print*,a(1:80)
              print*,'into m:',size(xc),'items'
           stop
     endif
  endif
  if (present(x)) then
     if (size(x) >= n) then
            read(a(first:last),'(f12.0)',iostat=stat)x(n)
            if (stat /=0) x(n)=0

         else
              print*, "Error in nums2 splitting string starting: "
              print*,a(1:80)
              print*,'into m:',size(xc),'items'            
              stop
     endif
  endif
  curr=last+1
enddo
end subroutine

subroutine getoption(name,n,x,xc)
! In unit io_p=40, which reads parameter line in BLUPF90, locates line:
! OPTION name str1 str2
! where str1, str2 are strings separated by spaces.
! Then, it assigns: xc(1)=str1,  xc(2)=str2,...
! and attempts to decode strings into real values: x1=value(str1),....
!
! n contains the number of strings.  x and xc are optional and their 
! dimensions may be smaller than n in which case some strings/values are
! not stored.
!
! Upon exit, unit io_p=40 points to line next to the one located.
!
! If the line cannot be located, n=-1
! 

character (*)::name
integer::n
real,optional::x(:)
character (*),optional::xc(:)

real::x1(1000)
integer::stat,m
character (50)::xc1(1000)
character (400)::a

rewind io_p
n=-1

do
   read (io_p,'(a)',iostat=stat)a
   if (stat /= 0) exit
      call nums2(a,m,x1,xc1)       
      if ( m>1 .and. xc1(1) == 'OPTION' .and. xc1(2) == name) then
         n=m-2
	 if (present(xc)) xc=xc1(3:size(xc)+2)
	 if (present(x)) x=x1(3:size(x)+2)
	 exit
      endif
enddo
end subroutine

subroutine chkgetoption(xa)
! check all arguments to be in the list of valid ones. 
character (*) :: xa(:)
real::x1(1000)
integer::stat,m,i
character (50)::xc1(1000)
character (400)::a
logical :: found
rewind io_p
do
   read (io_p,'(a)',iostat=stat)a
   if (stat /= 0) exit
      call nums2(a,m,x1,xc1)       
      if ( m>1 .and. xc1(1) == 'OPTION') then
         found=.false.
         do i=1,size(xa)
            if (trim(xa(i))==trim(xc1(2)) ) then
               found=.true. 
               exit
            endif
         enddo 
         if (.not. found) then 
            print '(/a)'  ,' ***** OPTION NOT valid for this program !!!! '
            print '(a,a,a)' ,'       "',trim(adjustl(xc1(2))),'"'
            print '(a,a)' ,'       Check for typing errors (e.g. upper/lower cases) in the key'
            stop
         endif
      endif

enddo
end subroutine

subroutine getoption_r8(name,n,x,xc)
! As getoption() but with x in r8 precision
! 
character (*)::name
integer::n
real (r8)::x(:)
real::xr4(size(x))
character (*),optional::xc(:)

!
call getoption(name,n,xr4,xc)
x=xr4
end subroutine

subroutine chkfmtnewline(un)
character(50000) :: a,i_name ! big number to allow check genotypes files
integer :: i,un
!
inquire(un,name=i_name)
read(un,'(a)') a 
do i=1,len(a) 
   if (ichar(a(i:i)) == 13) then ! Mac files
      print '(/a,a,a)', ' Mac old newline format detected in file: "',trim(i_name(index(i_name,'/',back=.true.)+1:)),'"'
      print '(a)'  , ' Convert it to Unix newline format, e.g. "flip -u <file>"'
      stop
   endif
enddo

if  (scan(achar(9),a)/=0) then 
   print '(/a,a,a)',' TAB was found, Use spaces to separeate columns in file: "',trim(i_name(index(i_name,'/',back=.true.)+1:)),'"'
   print '(a)'     ,' Convert to TAB to spaces, e.g. "expand file > newfile"'
   stop
endif 
rewind(un)
endsubroutine
end module

