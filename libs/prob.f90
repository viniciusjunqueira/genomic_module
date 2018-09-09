
module prob
! Collection of statistical subroutines helpful to implement threshold
!    model and Gibbs sampling
!   
!  Version using the RANLIB library with help and several
!  subroutines form by Deukhwan Lee
!
!  Ignacy Misztal, 04/15/98-12/7/2001

use denseop; use ranlib
implicit none


!=====================================================================
!Interfaces
!=====================================================================

  interface gen_uniform
      module procedure gen_uniform_real,&
                       gen_uniform_integer
  end interface


  interface gen_normal
      module procedure gen_normal_scalar,&
                       gen_normal_multi
  end interface

  interface gen_trunc_normal
      module procedure gen_trunc_normal_scalar,gen_trunc_normal_mult,&
                       gen_trunc_normal_mult_mt
  end interface
  
  
  interface gen_wishart
      module procedure gen_chisq,&
                       gen_Wishart_multi
  end interface


  interface gen_invwishart
      module procedure gen_invchisq,&
                       gen_invWishart_multi
  end interface


!----------------------------
! Module variables
!----------------------------



   CONTAINS
   
!==================================
! Module subroutines and functions
!==================================   

!----------------------------------------------------------------
! Seeting seed 
!----------------------------------------------------------------

   subroutine set_seed(n)
   ! sets seed for all random number generators in this module 
   
   integer::n
   real::x
   call inrgcm()
   x=ranf()	!need to generate at least one sample before setting the seed
   call setsd(n,n)		
   end subroutine
   
   
!----------------------------------------------------------------
! Uniform random number generators
!----------------------------------------------------------------   

   function gen_uniform_base() result (x)
   ! generates uniform(0,1)
   
   real (r8)::x
   real ,save::y(100)
   real::z
   integer::i,n
   
   x=ranf()
   end function   
	 
	 
	 
   function gen_uniform_real(a,b) result (x)	 
   ! generates unform real random number in interval <a,b)
   ! if a is missing, it is replaced by 0
   ! if b is missing, it is replaced by 1
   
   real (r8)::x,a1,b1
   real (r8),optional::a,b
   
   a1=0; b1=1
   if (present(a)) a1=a
   if (present(b)) b1=b
   
   x=gen_uniform_base()*(b1-a1)+a1
   end function
   
   
   
   function gen_uniform_integer(a,b) result (x)	 
   ! generates unform integer random number between a and b 
   
   integer::x,a,b
   
   x=floor(gen_uniform_base()*(b-a+1))+a
   end function

!----------------------------------------------------------------
!  generate normal univariate and multivariate samples with 
!  given mean and variance
!----------------------------------------------------------------

  function  gen_normal_scalar(mean,var) result (value)
  ! univariate
  real (r8),optional::mean,var
  real (r8)::value,p,mean1,var1
  
  mean1=0; var1=1
  if (present(mean)) mean1=mean
  if (present(var))  var1=var
 
  if (var1 <= 0) then
      value=mean1
    else
      value=mean1+sqrt(var1)*snorm()
  endif
  end function

   
  function  gen_normal_multi(mean,var) result (value)
  ! multivariate
  real (r8)::mean(:),var(:,:)
  real (r8),dimension(size(mean))::value
  real (r8),save, allocatable :: u(:),v(:,:)
  integer::i,j,m
  integer, save :: nmx=0
  
  m=size(mean)
  if (m>nmx) then 
     if (allocated(v)) then 
        deallocate(v,u)
     endif
     allocate(v(m,m),u(m))
     nmx=m
  endif      
  do i=1,m
     u(i)=snorm()
  enddo
  !value=mean+matmul(fchol(var),u)
  v(1:m,1:m)=fchol(var) 
  do i=1,m
     value(i)=mean(i)
     do j=1,m
        value(i)=value(i)+v(i,j)*u(j)
    enddo
  enddo
  end function

!----------------------------------------------------------------
!  generate (regular and inverted) chi square or wishart samples with 
!  given inverted quadratic forms and degrees of freedom
!----------------------------------------------------------------

function  gen_chisq(v,df) result (value)
  ! chi square

  real (r8)::v,value
  integer :: df
  
  value=2*sgamma(df/2.)*v
  
  end function


  function  gen_invchisq(v,df) result (value)
  ! chi square: i

  real (r8)::v,value
  integer :: df
  
  value=1/gen_chisq(v,df)  
  
  end function


 
  function  gen_wishart_multi(v,df) result (value)
  ! wishart distribution

  real (r8)::v(:,:)
  real (r8),dimension(size(v,1),size(v,2))::value,t,a,l
  integer :: df,n,i,j
  
  n=size(v,dim=1)
  do i=1,n
      t(i,i)=sqrt(2*sgamma((df-i+1)/2.))
       do  j=1,i-1
          t(i,j)=0.
          t(j,i)=snorm()
       enddo
  enddo
  
  a=matmul(t,transpose(t))
  l=fchol(v)
  value=matmul(matmul(l,a),transpose(l))
  end function


  function  gen_invwishart_multi(v,df) result (value)
  ! inverted wishart distribution

  real (r8)::v(:,:)
  real (r8),dimension(size(v,1),size(v,2))::value
  integer :: df
  
  value=finverse_s(gen_wishart(v,df))
  
  end function

 

!================================================================
! This is to generate IW distr. with restriction of categorical
! effects to 1 as referred by Korsgaard et.al.
! (1999, Genet. Sel. Evol.31:177-181)
! By Deukhwan Lee
! Modified to account for nonzero covariances among
! restricted effects
!================================================================
subroutine  gen_invwishart_restrict(LC,CC,V,df,value,SCAL)
!================================================================
 integer:: LC(:),CC(:),df
 real(r8)::V(:,:),VALUE(size(V,1),size(V,2)),SCAL(size(CC))
 real(r8)::R11(size(LC),size(LC)),R12(size(LC),size(CC)), &
           R22(size(CC),size(CC)),T11(size(LC),size(LC)), &
           T12(size(LC),size(CC)),T22(size(CC),size(CC)), &
           T2var(size(LC)*size(CC),size(LC)*size(CC)),    &
           T2mean(size(LC)*size(CC))
 integer::n,i1,i2,j1,j2,ii,jj,nlc,ncc

  n=size(V,1); nlc=size(LC); ncc=size(CC)
  R11=V(LC(:),LC(:)); R12=V(LC(:),CC(:)); R22=V(CC(:),CC(:))
  
  T11=gen_invwishart(R11,df)
  R11=finverse_s(R11)
  R22=R22-matmul(matmul(transpose(R12),R11),R12)
  R12=matmul(R11,R12)

  !Sampling covariances between un_restr. & restr. from normal
  do i1=1,ncc
     do j1=1,nlc
        ii=(i1-1)*nlc+j1
        T2mean(ii)=R12(j1,i1)
     enddo
  enddo

  ! Kronecker product
  do i1=1,ncc
     do j1=1,nlc
        ii=(i1-1)*nlc+j1
        do i2=1,ncc
           do j2=1,nlc
              jj=(i2-1)*nlc+j2
              T2var(ii,jj)=T11(j1,j2)*R22(i1,i2)
           enddo
        enddo
     enddo
  enddo

  T2mean=gen_normal(T2mean,T2var)

  do i1=1,ncc
     do j1=1,nlc
        ii=(i1-1)*nlc+j1
        T12(j1,i1)=T2mean(ii)
     enddo
  enddo
  T22=gen_invwishart(R22,df-nlc)
  do i1=1,NCC
     SCAL(i1)=sqrt(T22(i1,i1))
     T22(i1,:)=T22(i1,:)/SCAL(i1)
     T22(:,i1)=T22(:,i1)/SCAL(i1)
  enddo

  T11=T11+matmul(matmul(T12,T22),transpose(T12))
  T12= -1.d0*matmul(T12,T22)
  VALUE(LC(:),LC(:))=T11
  VALUE(LC(:),CC(:))=T12
  VALUE(CC(:),LC(:))=transpose(T12)
  VALUE(CC(:),CC(:))=T22

  end subroutine

!================================================================
 subroutine  gen_invwishart_restrict1(LC,CC,V,df,value,SCAL)
!================================================================
! a simpler version of the above procedure where a sample from a regular
! inverted Wishart distribution is scaled to 1 by linear transformation and 
! the scaling values put into variable scal
!
 integer:: LC(:),CC(:),df,i,j,n
 real(r8)::V(:,:),VALUE(size(V,1),size(V,2)),SCAL(size(v,1))

 value=gen_invwishart(v,df)
 n=size(v,1); scal=1


 do i=1,size(cc)
   scal(cc(i))=sqrt(value(cc(i),cc(i)))
 enddo

 do i=1,n
   do j=1,n
      value(i,j)=value(i,j)/scal(i)/scal(j)
   enddo
 enddo

 end subroutine         

!--------------------------------------------------------------
! normal probability, cdf, inverse cdf
!--------------------------------------------------------------

function normal(x)
!
! Returns probability density of normal density function
!
real (r8):: normal,x
real (r8), parameter::pi=3.14159265358979
!
normal=1/sqrt(2.*pi)*dexp(-x*x/2)
end function


function normalcdf(x,mean,sd) result(p)
! cumulative normal density function

real(r8)::x,p,q,mean1=0,sd1=1,bound
real (r8),optional::mean,sd
integer:: status
!
  if (present(mean)) mean1=mean
  if (present(sd)) sd1=sd
  
  call cdfnor(1,p,q,x,mean1,sd1,status,bound)
  
  if (status /= 0) then
     print*,'Parameters for NORMALCDF incorrect'
     stop
  endif   
end function



function normal_invcdf(p,mean,sd) result(x)
! return inverse of CDF, i.e., such x: p=cdf(x)
!   for stability with low/high p, always returns |invcdf(x)|<6.36.

real(r8)::x,p,q,mean1=0,sd1=1,bound
real (r8),optional::mean,sd
integer:: status
real(r8),parameter::epsCDF=1E-10	!,epsX=6.36134090240406

!
  if (present(mean)) mean1=mean
  if (present(sd)) sd1=sd
  
  if(p>(1.d0-epsCDF)) then
     p=1.d0-epsCDF
  elseif(p<epsCDF) then
     p=epsCDF
  endif 

  q=1-p
  
  call cdfnor(2,p,q,x,mean1,sd1,status,bound)
  
  if (status /= 0) then
     print*,'Parameters for NORMAL_INVCDF incorrect'
     stop
  endif   

end function


!-------------------------------------------------------------------
! generates truncated normal distribution
!-------------------------------------------------------------------

function gen_trunc_normal_scalar(a,b,mean,var) result (x)
! Generates normal distribution N(mean,var) truncated to <a,b>;
!   missing mean defaults to 0, missing var defaults to 1.

real(r8)::a,b,cdfa,cdfb,x
real (r8),optional::mean,var
real (r8)::mean1,var1

if (b < a) then
   Print*,'GEN_TRUNC_NORMAL: bound b =',b,' < bound a =',a
   stop
endif   


mean1=0; var1=1
if (present(mean)) mean1=mean
if (present(var))  var1=var
   
if (var1 <= 0) then
   Print*,'GEN_TRUNC_NORMAL: variance =',var1,' <= 0'
   stop
endif
     
      
cdfa=normalcdf(a,mean1,sqrt(var1))
cdfb=normalcdf(b,mean1,sqrt(var1)) 
x=sqrt(var1)*normal_invcdf(cdfa+(cdfb-cdfa)*gen_uniform())+mean1
if(x>b)x=b
if(x<a)x=a
end function



function gen_trunc_normal_mult(a,b,mean,Var) result (x)
! Generates multivariate normal distribution N(0,V) with v(1,1)=1 and
! the first variable truncated to <a,b>. f

real (r8)::a,b,var(:,:),mean(:)
real (r8)::x(size(mean))
real (r8)::l(size(var,dim=1),size(var,dim=1)),u(size(mean))
integer::i

l=fchol(var)
u(1)=gen_trunc_normal((a-mean(1))/sqrt(var(1,1)), (b-mean(1))/sqrt(var(1,1)))
do i=2,size(mean)
  u(i)=gen_normal(0d0,1d0)
enddo

x=matmul(l,u)+mean
end function


function gen_trunc_normal_mult_mt(lbound,ubound,mean,var) result(Ui)
! generates MVN(mean,var) samples where variables corresponding to dimensions 
! where ubound>lbound = .true. are trunctated by lower bounds given by lbound
! and upper bounds given by ubound
! Written by Deukhwan Lee and modified by IM
!
  real(r8):: lbound(:),ubound(:),mean(:),var(:,:)
  real(r8):: wm(1),wv(1,1),vinv(size(var,1),size(var,2)),Ui(size(mean))
  integer:: i,j,n
  logical::is_trunc(size(lbound))
  
  is_trunc = ubound > lbound
  
  n=size(mean)
  if(is_trunc(1)) then
     Ui(1)=gen_trunc_normal(lbound(1),ubound(1),mean(1),var(1,1))
   else  
     Ui(1)=gen_normal(mean(1),var(1,1))
  endif
 
  if(n>1) then
    do i=2,n
       j=i-1
       vinv(:j,:j)=finverse_s(var(:j,:j))
       wm(1:1)=mean(i:i)+matmul(matmul(var(i:i,:j), &
               vinv(:j,:j)),(Ui(:j)-mean(:j)))
       wv(1:1,1:1)=var(i:i,i:i)-matmul(matmul(var(i:i,:j), &
               vinv(:j,:j)),var(:j,i:i))
       if(is_trunc(i)) then
           Ui(i)=gen_trunc_normal(lbound(i),ubound(i),wm(1),wv(1,1))
	  else 
           Ui(i)=gen_normal(wm(1),wv(1,1))
       endif
    enddo
  endif
  end function



function gen_trunc_normal_mult_mt1(lbound,ubound,mean,var) result(Ui)
! Function as above but modified to generate all the initial unrestricted
! dimensions jointly rather than sequentially

  real(r8):: lbound(:),ubound(:),mean(:),var(:,:)
  real(r8):: wm(1),wv(1,1),vinv(size(var,1),size(var,2)),Ui(size(mean))
  integer:: i,j,k,n
  logical::is_trunc(size(lbound))
!  
  is_trunc=ubound>lbound
  n=size(mean)
  
  if(is_trunc(1)) then
     Ui(1)=gen_trunc_normal(lbound(1),ubound(1),mean(1),var(1,1))
     k=2
   else 
     ! generate all initial untruncated dimensions
     do i=2,n
        if (is_trunc(i)) then
	   k=i-1
	   exit
	endif
     enddo
     Ui(1:k)=gen_normal(mean(1:k),var(1:k,1:k))
     k=k+1
  endif
 
  if(n>1) then
    do i=k,n
       j=i-1
       vinv(:j,:j)=finverse_s(var(:j,:j))
       wm(1:1)=mean(i:i)+matmul(matmul(var(i:i,:j), &
               vinv(:j,:j)),(Ui(:j)-mean(:j)))
       wv(1:1,1:1)=var(i:i,i:i)-matmul(matmul(var(i:i,:j), &
               vinv(:j,:j)),var(:j,i:i))
       if(is_trunc(i)) then
           Ui(i)=gen_trunc_normal(lbound(i),ubound(i),wm(1),wv(1,1))
	  else 
           Ui(i)=gen_normal(wm(1),wv(1,1))
       endif
    enddo
  endif
  end function




function gen_trunc_normal_mult_mt2(lbound,ubound,mean,var) result(Ui)
! 
! Function as above but modified to generate unrestricted
! dimensions jointly 
!
!
  real(r8):: lbound(:),ubound(:),mean(:),var(:,:),Ui(size(mean))
  integer:: i,ind(size(lbound)),last,first
!    
! find permutation in ind that would align all the unrestricted dimensions 
!   first  

  first=1; last=size(mean)
   
  do i=1,size(mean)
     if (ubound(i)>lbound(i)) then
          ind(last)=i
	  last=last-1
	else
	  ind(first)=i
	  first=-first+1
     endif
  enddo
  
  ui(ind)=gen_trunc_normal_mult_mt1(lbound(ind),ubound(ind),&
                      mean(ind),var(ind,ind))    	  
  end function

end module





