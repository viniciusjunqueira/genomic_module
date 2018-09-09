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

module model
use kinds
implicit none

! FC MCMC chain controls and inputs
integer:: nround,nburn,round,start,&  ! number of rounds, length of burn-in, starting round
          seed,skip,&                 ! save each "skip" cycles
          mb_nround,mb_skip,&  ! # of MH rounds within cycle, and number of round between tuning process of MH
          r_nround,r_skip   ! # of MH rounds within cycle, and number of round between tuning process of MH

character (1) :: rstart     ! choice of new chain or restart

! HGC BLUP controls and inputs
real :: r_factor=1.4,conv_crit=1e-12
character (20):: solv_method='PCG',sol_se= ' '
integer :: maxrounds=5000,&
           blksize=1
logical :: lzerosol=.true.
character(50) ::      fmt_sol_se=        '(2i4,i10,2f20.8)'


!        Types of effects
integer,parameter::effcross=0,&   !effects can be cross-classified 
                   effcov=1, &    !or covariables
                   effram=2, &    !FC reduced animal model effect
                   effuncov=3, &  !FC unknown covariables
                   effrnorm=4     !FC reaction norm regressed on unknown covariables

!        Types of random effects
integer,  parameter ::  g_fixed=1,&       ! fixed effect
                        g_diag=2, &       ! diagonal
                        g_A=3, &          ! additive animal
                        g_A_UPG=4, &      ! additive animal with unknown parent groups (upg)
                        g_A_UPG_INB=5, &  ! additive animal with upg and inbreeding
                        g_As=6,&          ! additive sire
                        g_PD =7, &        ! parental dominance
                        g_A_MS=8, &       ! FC am/ram w/msires and genetic groups
                        g_A_MB=9, &       ! FC multiple breed additive am 
                        g_diag_MB=10,&    ! FC multiple breed uncorrelated
	                      g_usr=11,&	      ! user defined
		                    g_usr_inv=12,&	  ! user defined and inverted 
                        g_last=13,&       ! last type
                        g_normal=14       ! VJ normal prior

! FC     Types of random residual effects
integer,  parameter ::  r_homo=1,&            ! homogeneous
                        r_structural=2,&      ! structural heterogeneous variance
                        r_slash=3,&		        ! robust slash errors
                        r_struct_slash=4,&    ! structural heterogeneous variance and slash errors
                        r_studentt=5,&	      ! robust student t errors
                        r_struct_studentt=6,& ! structural heterogeneous variance and student t errors
                        r_last=5              ! last type

character (250)     ::   parfile, &     ! name of parameter file
                         datafile,&     ! name of data set
                         msfile         ! FC name of file w/msires info

integer :: ntrait,&                    ! number of traits
           neff,&                      ! number of effects
           miss=0,&                    ! value of missing trait/effect
           nrec=0,&                    ! FC number of records
           max_msires=0, &             ! FC max # msires to use full model
           n_ms, &				             ! FC dimension of msires vector
           neffr=0  				           ! FC number of effects (structutal) on residual var

real :: ms_alpha=0			! Dirichlet alpha values equal for all sire assignments 
							          ! If alpha=1 => flat prior; if alpha<0 get values from ms    

integer,allocatable :: pos_y(:)        ! positions of observations
integer,allocatable :: pos_weight(:)   ! position of weight of records; zero if none ! FC trait specific weights


integer,allocatable :: pos_eff(:,:),&   ! positions of effects for each trait
                       nlev(:),&        ! number of levels
                       lenped(:),&      ! length of pedigree files
                       effecttype(:),&  ! type of effects
                       nestedcov(:,:),& ! position of nesting effect for each trait
                                        ! if the effect is nested covariable
                     & randomtype(:),&  ! status of each effect, as above
                       randomnumb(:),&  ! number of consecutive correlated effects
                       nbreed(:),&      ! FC number of breeds in population
		                   dbelief_g(:),&   ! FC degrees of belief on prior g (co)variances
                       ngrp(:),&        ! FC # of unknown parent groups
                       effectsave(:),&  ! FC if = 1 save samples of effect
		                   pos_effr(:),&    ! FC position of struc eff on r in the pos_eff vector
                       effrsave(:)      ! FC if = 1 save samples of effect on r

integer ::             residualtype,&   ! FC type of residual
		                   dbelief_r        ! FC degrees of belief on prior r (co)variances


character (250),allocatable::  randomfile(:)   ! name of file associated with given effect

real (rh), allocatable :: r(:,:),&      ! residual (co)variance matrix
                          rinv(:,:),&   ! and its inverse
                          g(:,:,:),&    ! The random (co)variance matrix for each trait
                          ginv(:,:,:)   ! FC and its inverse

integer,parameter::   maxcorr=50,&     ! maximum number of correlated effects; used for g 
                      io_off=50,&      ! unit number = i+io_off for effect i
                      io_d=50,&        ! unit number for data file
                      io_s=30,&        ! unit number for solution file
		                  io_p=40,&        ! unit number for parameter file
                      io_ms=49,&       ! unit number for msires parameter file
                      io_dinv=69,&     ! unit number for dinv file
                      io_thr=79,&      ! unit number for thresholds solution file
                      maxnsires=50     ! FC maximum number of msires; must be at least 1 for ram 

! Option variables
integer  ::  cat_value=-1,&    ! variable to inform categorical analysis
             nthr=0,&          ! number of thresholds
             rthr=0,&          ! number of thresholds iterations
             fthr=0,&          ! first iteration
             maxround_thr=40
!real     ::  cat_nthr(1)=0.0
                      
end module model
