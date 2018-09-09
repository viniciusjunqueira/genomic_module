! 20060418 flag initialization bugs fixed
! 20060328 flag incorporated into chols4 & chols8
!     Last change:  ST   15 Aug 2001    4:38 pm
MODULE DENSEOP
! 5/Oct/1998-16/Apr/03 by Tomasz Strabel & Ignacy Misztal, University of Georgia.
! Modification 8/15/2001 by Shogo Tsuruta
! Last modified 12/2/2003
USE KINDS; USE LAPACK90
implicit none

INTERFACE SOLVE_S
 MODULE PROCEDURE SOLVS4F, SOLVS8F, SOLVS4T, SOLVS8T
END INTERFACE

INTERFACE CHOL
 MODULE PROCEDURE CHOLS4, CHOLS8, CHOLS4T, CHOLS8T
END INTERFACE

INTERFACE FDET_S
 MODULE PROCEDURE FDETS4, FDETS8, FDETS4P, FDETS8P
END INTERFACE 

INTERFACE INVERSE_S
 MODULE PROCEDURE INVCS4, INVCS8, INVCS4T, INVCS8T
END INTERFACE

INTERFACE EIGEN
 MODULE PROCEDURE SYEV4, SYEV8
END INTERFACE

INTERFACE SOLVE
 MODULE PROCEDURE GESV8, GESV4
END INTERFACE

INTERFACE INVERSE
 MODULE PROCEDURE GEINV8, GEINV4
END INTERFACE

INTERFACE EYE
 MODULE PROCEDURE EYE4, EYE8
END INTERFACE

INTERFACE FSOLVE
 MODULE PROCEDURE FGESV8, FGESV4
END INTERFACE

INTERFACE FCHOL
 MODULE PROCEDURE FCHOLS4, FCHOLS8, FCHOLS4T, FCHOLS8T
END INTERFACE 

INTERFACE FINVERSE_S
 MODULE PROCEDURE FINVCS4, FINVCS8, FINVCS4T, FINVCS8T
END INTERFACE 

INTERFACE FSOLVE_S
 MODULE PROCEDURE FSOLVS4, FSOLVS8, FSOLVS4T, FSOLVS8T
END INTERFACE 

INTERFACE FINVERSE
 MODULE PROCEDURE FGEINV4, FGEINV8
END INTERFACE 

INTERFACE PRINTMAT
 MODULE PROCEDURE PRINTM4, PRINTM8, PRINTM4P, PRINTM8P
END INTERFACE

INTERFACE PACKIT
 MODULE PROCEDURE PACKIT4, PACKIT8
END INTERFACE

INTERFACE UNPACKIT
 MODULE PROCEDURE UNPACKIT4, UNPACKIT8
END INTERFACE

interface diag
  module procedure diag_vec_mat8, diag_mat_vec8, diag_vec_mat4, diag_mat_vec4
end interface  

interface kron
  module procedure kron4,kron8
end interface

interface mat_mul_vec_symm
  module procedure mat_mul_vec_symm4, mat_mul_vec_symm8
end interface

interface checksym
  module procedure checksym4,checksym8
end interface 


real(r8)::denseop_tol=1.d-10


CONTAINS

!subroutines with 4/8 in their names work with r4/r8 arguments


SUBROUTINE INVCS4T (x,rank)
!SUBROUTINE INVCS4 TO WORK WITH MATRICES STORED IN PACKED FORM
integer ::i,j,n
integer, intent(out), optional::rank
real(r4)::x(:)
real(r4), allocatable:: x_o(:,:)
n=(-1+sqrt(1.0+8*size(x,1)))/2
allocate(x_o(n,n))

call unpackit(x,x_o)
call inverse_s(x_o,rank)
call packit(x_o,x)
END SUBROUTINE INVCS4T

FUNCTION FINVCS4T(X) RESULT(Y)
 real(r4)::X(:),y(size(x,1))
 y=x
 call inverse_s(y)
END FUNCTION FINVCS4T

FUNCTION FINVCS8T(X) RESULT(Y)
 real(r8)::X(:),y(size(x,1))
 y=x
 call inverse_s(y)
END FUNCTION FINVCS8T

SUBROUTINE INVCS8T (x,rank)
!SUBROUTINE INVCS8 TO WORK WITH MATRICES STORED IN PACKED FORM
integer ::i,j,n
integer, intent(out), optional::rank
real(r8)::x(:)
real(r8), allocatable:: x_o(:,:)
n=(-1+sqrt(1.0+8*size(x,1)))/2
allocate(x_o(n,n))

call unpackit(x,x_o)
call inverse_s(x_o,rank)
call packit(x_o,x)
deallocate(x_o)
END SUBROUTINE INVCS8T

SUBROUTINE SOLVS4T (x,v,sol,rank,flag)
!SUBROUTINE INVCS4 TO WORK WITH MATRICES STORED IN PACKED FORM
integer ::i,j,n
integer, intent(out), optional::rank,flag
real(r4)::x(:),v(:),sol(:)
real(r4), allocatable:: x_o(:,:)
n=(-1+sqrt(1.0+8*size(x,1)))/2
allocate(x_o(n,n))
sol=v
call unpackit(x,x_o)
call solvs4(x_o,sol,rank,flag)
deallocate(x_o)
END SUBROUTINE SOLVS4T

FUNCTION FSOLVS4(X,RHS) RESULT(SOL)
 real(r4)::X(:,:),RHS(:),sol(size(x,1))
 sol=0
 call solve_s(x,rhs,sol)
END FUNCTION FSOLVS4

FUNCTION FSOLVS4T(X,RHS) RESULT(SOL)
 real(r4)::X(:),RHS(:),sol(size(rhs,1))
 sol=0
 call solve_s(x,rhs,sol)
END FUNCTION FSOLVS4T

FUNCTION FSOLVS8T(X,RHS) RESULT(SOL)
 real(r8)::X(:),RHS(:),sol(size(rhs,1))
 sol=0
 call solve_s(x,rhs,sol)
END FUNCTION FSOLVS8T

FUNCTION FSOLVS8(X,RHS) RESULT(SOL)
 real(r8)::X(:,:),RHS(:),sol(size(x,1))
 sol=0
 call solve_s(x,rhs,sol)
END FUNCTION FSOLVS8

SUBROUTINE SOLVS8T (x,v,sol,rank,flag)
!SUBROUTINE CHOLS8 TO WORK WITH MATRICES STORED IN PACKED FORM
integer ::i,j,n
integer, intent(out), optional::rank,flag
real(r8)::x(:),v(:),sol(:)
real(r8), allocatable:: x_o(:,:)
n=(-1+sqrt(1.0+8*size(x,1)))/2
allocate(x_o(n,n))
sol=v
call unpackit(x,x_o)
call solvs8(x_o,sol,rank,flag)
deallocate(x_o)
END SUBROUTINE SOLVS8T


SUBROUTINE CHOLS4T (x,rank,flag)
!SUBROUTINE CHOLS4 TO WORK WITH MATRICES STORED IN PACKED FORM
integer ::i,j,n
integer, intent(out), optional::rank,flag
real(r4)::x(:)
real(r4), allocatable:: x_o(:,:)
n=(-1+sqrt(1.0+8*size(x,1)))/2
allocate(x_o(n,n))
call unpackit(x,x_o)
call chols4(x_o,rank,flag)
call packit(x_o,x)
deallocate(x_o)
END SUBROUTINE CHOLS4T

FUNCTION FCHOLS4T(X) RESULT(Y)
 real(r4)::X(:),y(size(x,1))
 y=x
 call chol(y)
END FUNCTION FCHOLS4T

FUNCTION FCHOLS8T(X) RESULT(Y)
 real(r8)::X(:),y(size(x,1))
 y=x
 call chol(y)
END FUNCTION FCHOLS8T

SUBROUTINE CHOLS8T (x,rank,flag)
!SUBROUTINE CHOLS4 TO WORK WITH MATRICES STORED IN PACKED FORM
integer ::i,j,n
integer, intent(out), optional::rank,flag
real(r8)::x(:)
real(r8), allocatable:: x_o(:,:)
n=(-1+sqrt(1.0+8*size(x,1)))/2
allocate(x_o(n,n))

call unpackit(x,x_o)
call chol(x_o,rank,flag)
call packit(x_o,x)

END SUBROUTINE CHOLS8T

SUBROUTINE SOLVS4F(l,rhs,sol,rank,flag)
integer, intent(out), optional::rank,flag
real(r4)::l(:,:),rhs(:),sol(:)
real(r4),allocatable::l_o(:,:)
allocate(l_o(size(l,1),size(l,1)))
l_o=l
sol=rhs
call solvs4(l_o,sol,rank,flag)
deallocate(l_o)
END SUBROUTINE SOLVS4F

SUBROUTINE SOLVS8F(l,rhs,sol,rank,flag)
integer, intent(out), optional::rank,flag
real(r8)::l(:,:),rhs(:),sol(:)
real(r8),allocatable::l_o(:,:)
allocate(l_o(size(l,1),size(l,1)))
l_o=l
sol=rhs
call solvs8(l_o,sol,rank,flag)
deallocate(l_o)
END SUBROUTINE SOLVS8F


SUBROUTINE SOLVS4(l,rhs,rank,flag)
! solve: cholesky & LL'x=rhs
integer::n,i,j,rank_o,flag_o
integer, intent(out), optional::rank,flag
real(r4)::l(:,:),rhs(:)
real(r4),allocatable::z(:)
allocate(z(size(l,1)))
z=0

flag_o=0
rank_o=0
n=size(l,1)

call chols4(l,rank_o,flag_o) ! cholesky factorization
if(flag_o/=0) then
   print*, 'Solution not possible'
   if(present(flag)) flag=-1
   RETURN
endif

z=0			   ! Lz=y  z=?
do i=1,n
  if (l(i,i)>0)then
    z(i)=(rhs(i)-dot_product(l(i,1:i-1),z(1:i-1)))/l(i,i)
  end if 
end do

rhs=0			! L'x=z x=?
do i=n,1,-1
  if (l(i,i)>0)then
    rhs(i)=(z(i)-dot_product(l(i+1:n,i),rhs(i+1:n)))/l(i,i)
  end if 
end do

if (present(rank)) then
 rank=rank_o
end if

if (present(flag)) flag=flag_o

deallocate(z)

END SUBROUTINE SOLVS4

FUNCTION FCHOLS4(X) RESULT(Y)
 real(r4)::X(:,:),y(size(x,1),size(x,1))
 y=x
 call chol(y)
END FUNCTION FCHOLS4

FUNCTION FCHOLS8(X)  RESULT(Y)
 real(r8)::X(:,:),y(size(x,1),size(x,1))
 y=x
 call chol(y)
END FUNCTION FCHOLS8


SUBROUTINE CHOLS4(x,rank,flag)
! cholesky decomposition
integer :: n,i,j,k,rank_o, flag_o
integer, intent(out), optional::rank, flag
real(r4)::diagsq
real(r4)::x(:,:)
rank_o=0
flag_o=0
n=size(x,1)
do i=1,n
   diagsq=x(i,i)-dot_product(x(i,1:i-1),x(i,1:i-1))
   if (abs(diagsq).lt.denseop_tol) then
        x(i,:)=0;x(:,i)=0       !zero row and column
      elseif (diagsq.lt.0) then
        print*,' Matrix not semipositive-definite, row ',i
        flag_o=-1
        EXIT
        !stop
!   endif
   else
    rank_o=rank_o+1
   x(i,i)=sqrt(diagsq)
   do j=i+1,n     
      x(j,i)=(x(j,i)-dot_product(x(j,1:i-1),x(i,1:i-1)))/x(i,i)
      x(i,j)=x(j,i)
   enddo
   end if
enddo

! zero upper-diagonals
do i=1,n
   x(i,i+1:n)=0
enddo   

if (present(rank))then
 rank=rank_o
end if

if (present(flag)) then
 flag=flag_o
end if
END SUBROUTINE CHOLS4


SUBROUTINE INVCS4(l,rank,flag)
! calculates inverse of LL', where L is cholesky decomposition 
!USE LA_PRECISION, ONLY: WP => DP
integer  ::n,i,j,rank_o,flag_o
integer, intent(out), optional::rank,flag
real(r4) :: l(:,:)
real(r4) ::w(size(l,1))
n=size(l,1)

flag_o=0
rank_o=0

call chols4(l,rank_o,flag_o) 	! cholesky factorization
if(flag_o/=0) then
   print*, 'Inverse not possible'
   if(present(flag)) flag=-1
   RETURN
endif

do i=1,n
   w(i:n)=0
   if (abs(l(i,i)).gt.denseop_tol) w(i)=1/l(i,i)
   ! forward substitution
   do j=i+1,n
    if (abs(l(j,j)).gt.denseop_tol)  w(j)=-dot_product(l(j,i:j-1),w(i:j-1))/l(j,j)
   enddo
   !backward substitution
   do j=n,i,-1
    if (abs(l(j,j)).gt.denseop_tol)  w(j)=(w(j)-dot_product(l(j+1:n,j),w(j+1:n)))/l(j,j)
   enddo
   l(i:n,i)=w(i:n)
   l(i,i:n)=w(i:n)
enddo   

if (present(rank)) then
 rank=rank_o
end if

END SUBROUTINE INVCS4

FUNCTION FINVCS4(X) RESULT(Y)
 real(r4)::X(:,:),y(size(x,1),size(x,1))
 y=x
 call inverse_s(y)
END FUNCTION FINVCS4

FUNCTION FINVCS8(X) RESULT(Y)
 real(r8)::X(:,:),y(size(x,1),size(x,1))
 y=x
 call inverse_s(y)
END FUNCTION FINVCS8


SUBROUTINE CHOLS8(x,rank,flag)
! cholesky decomposition
integer :: n,i,j,k,rank_o,flag_o
integer, intent(out), optional::rank,flag
real(r8)::diagsq,d_p
real(r8)::x(:,:)
rank_o=0
flag_o=0
n=size(x,1)
do i=1,n
   d_p=0
   do k=1,i-1
      d_p=d_p+x(i,k)*x(i,k)
   enddo
!   diagsq=x(i,i)-dot_product(x(i,1:i-1),x(i,1:i-1))
   diagsq=x(i,i)-d_p
  if (abs(diagsq).lt.denseop_tol) then
        x(i,:)=0;x(:,i)=0       !zero row and column
      elseif (diagsq.lt.0) then
        print*,' Matrix not semipositive-definite, row ',i
        flag_o=-1
        EXIT
        !stop
!   endif
   else
    rank_o=rank_o+1
   x(i,i)=sqrt(diagsq)
   do j=i+1,n     
      d_p=0
      do k=1,i-1
         d_p=d_p+x(j,k)*x(i,k)
      end do
!      x(j,i)=(x(j,i)-dot_product(x(j,1:i-1),x(i,1:i-1)))/x(i,i)
      x(j,i)=(x(j,i)-d_p)/x(i,i)
      x(i,j)=x(j,i)
   enddo
   end if
enddo
! zero upper-diagonals

do i=1,n
   x(i,i+1:n)=0
enddo   

if (present(rank))then
 rank=rank_o
end if

if (present(flag)) then
 flag=flag_o
end if
END SUBROUTINE CHOLS8


SUBROUTINE SOLVS8(l,rhs,rank,flag)
! solve: cholesky & LL'x=rhs
integer::n,i,j,rank_o,flag_o
integer, intent(out), optional::rank,flag
real(r8)::l(:,:),rhs(:)
real(r8)::z(size(l,1))
z=0

flag_o=0
rank_o=0
n=size(l,1)

call chols8(l,rank_o,flag_o) 	! cholesky factorization
if(flag_o/=0) then
   print*, 'Solving not possible'
   if(present(flag)) flag=-1
   RETURN
endif

z=0			! Lz=y  z=?
do i=1,n
  if (l(i,i)>0)then
    z(i)=(rhs(i)-dot_product(l(i,1:i-1),z(1:i-1)))/l(i,i)
  end if 
end do

rhs=0			! L'x=z x=?
do i=n,1,-1
  if (l(i,i)>0)then
    rhs(i)=(z(i)-dot_product(l(i+1:n,i),rhs(i+1:n)))/l(i,i)
  end if 
end do

if (present(rank)) then
 rank=rank_o
end if

if (present(flag)) then
 flag=flag_o
end if

END SUBROUTINE SOLVS8

FUNCTION FDETS8(matrix) result(det)
integer::i
real(r8)::matrix(:,:),matrix_o(size(matrix,1),size(matrix,1)),det
matrix_o=matrix
call chol(matrix_o)
det=1
do i=1,size(matrix,1)
   if (matrix_o(i,i)/=0.0) det=det*matrix_o(i,i)
end do
det=det**2
END FUNCTION FDETS8

FUNCTION FDETS4(matrix) result(det)
real(r4)::matrix(:,:),det
real(r8)::matrix_o(size(matrix,1),size(matrix,1))
matrix_o=matrix
det=fdets8(matrix_o)
END FUNCTION FDETS4

FUNCTION FDETS4P(matrix) result(det)
integer::n
real(r4)::matrix(:),det 
real(r4),allocatable::matrix_o(:,:)
n=(-1+sqrt(1.0+8*size(matrix,1)))/2
allocate(matrix_o(n,n))
call unpackit(matrix,matrix_o)
det=fdets4(matrix_o)
END FUNCTION FDETS4P

FUNCTION FDETS8P(matrix) result(det)
integer::n
real(r8)::matrix(:),det 
real(r8),allocatable::matrix_o(:,:)
n=(-1+sqrt(1.0+8*size(matrix,1)))/2
allocate(matrix_o(n,n))
call unpackit(matrix,matrix_o)
det=fdets8(matrix_o)
deallocate(matrix_o)
END FUNCTION FDETS8P

SUBROUTINE INVCS8(l,rank,flag)
! calculates inverse of LL', where L is cholesky decomposition 
!USE LA_PRECISION, ONLY: WP => DP
integer  ::n,i,j,k,rank_o,flag_o
integer, intent(out), optional::rank,flag
real(r8) :: l(:,:),d_p
real(r8) ::w(size(l,1))
n=size(l,1)

flag_o=0
rank_o=0

call chols8(l,rank_o,flag_o) 	! cholesky factorization
if(flag_o/=0) then
   print*, 'Inverse not possible'
   if(present(flag)) flag=-1
   RETURN
endif

do i=1,n
   w(i:n)=0
   if (abs(l(i,i)).gt.denseop_tol) w(i)=1/l(i,i)
   ! forward substitution
   do j=i+1,n
      d_p=0
      do k=i,j-1
         d_p=d_p+l(j,k)*w(k)
      end do
!      if (abs(l(j,j)).gt.denseop_tol) w(j)=-dot_product(l(j,i:j-1),w(i:j-1))/l(j,j)
      if (abs(l(j,j)).gt.denseop_tol) w(j)=-d_p/l(j,j)
   enddo
   !backward substitution
   do j=n,i,-1
      d_p=0
      do k=j+1,n
         d_p=d_p+l(k,j)*w(k)
      end do
!      if (abs(l(j,j)).gt.denseop_tol) w(j)=(w(j)-dot_product(l(j+1:n,j),w(j+1:n)))/l(j,j)
      if (abs(l(j,j)).gt.denseop_tol) w(j)=(w(j)-d_p)/l(j,j)
   enddo
   l(i:n,i)=w(i:n)
   l(i,i:n)=w(i:n)
enddo   

if (present(rank)) then
 rank=rank_o
end if

END SUBROUTINE INVCS8


FUNCTION POS_SYMM(i,j,n) result (address)
 !finds position of (i,j) element of a n*n matrix upper-stored
 integer :: i,j,n,address
 if (j >= i) then
      address=(i-1)*n-(i*(i-3))/2+j-i
   else
      address=(j-1)*n-(j*(j-3))/2+i-j
 endif
END FUNCTION POS_SYMM

function eye4(m) 
integer::i,n
real(r4)::m(:,:)
real(r4)::eye4(size(m,1),size(m,1))
n=size(m,1)
eye4=0
do i=1,n
 eye4(i,i)=1
end do
end function eye4

function eye8(m) 
integer::i,n
real(r8)::m(:,:)
real(r8)::eye8(size(m,1),size(m,1))
n=size(m,1)
eye8=0
do i=1,n
 eye8(i,i)=1
end do
end function eye8

!-------------------------------------------------------!
!interfaces to work with subroutines of lapack90        !
!-------------------------------------------------------!

FUNCTION FGESV4(X,RHS) RESULT(SOL)
 real(r4)::X(:,:),RHS(:),sol(size(x,1))
 sol=0
 call solve(x,rhs,sol)
END FUNCTION FGESV4

FUNCTION FGESV8(X,RHS) RESULT(SOL)
 real(r8)::X(:,:),RHS(:),sol(size(x,1))
 sol=0
 call solve(x,rhs,sol)
END FUNCTION FGESV8

FUNCTION FGEINV4(X) RESULT(Y)
 real(r4)::X(:,:),y(size(x,1),size(x,1))
 y=x
 call inverse(y)
END FUNCTION FGEINV4

FUNCTION FGEINV8(X) RESULT(Y)
 real(r8)::X(:,:),y(size(x,1),size(x,1))
 y=x
 call inverse(y)
END FUNCTION FGEINV8

SUBROUTINE GEINV8 (x)
USE LAPACK90
 real(r8)::x(:,:)
 real(r8),allocatable::x_o(:,:)
 allocate(x_o(size(x,1),size(x,1)))
 x_o=eye(x_o)
 CALL dgesv_f90(x,x_o)
 x=x_o
END SUBROUTINE GEINV8

SUBROUTINE GEINV4(x)
USE LAPACK90
real(r8), allocatable:: x_o(:,:),x_o1(:,:)
real(r4)::x(:,:)
allocate (x_o(size(x,1),size(x,1)))
allocate (x_o1(size(x,1),size(x,1)))
x_o=eye(x_o1)
x_o1=x
CALL dgesv_f90(x_o1,x_o)
x=x_o
END SUBROUTINE GEINV4

SUBROUTINE GESV8(x,v,sol)
USE LAPACK90
integer::i
real(r8)::x(:,:),v(:),sol(:)
real(r8),allocatable::v_o(:,:),x_o(:,:)

!convert vector v(:) -> v_o(:,1)
allocate(v_o(size(v,1),1))
allocate(x_o(size(x,1),size(x,1)))
x_o=x
 do i=1,size(v,1)
   v_o(i,1)=v(i)
 end do

call dgesv_f90(x_o,v_o)

!rewrite vector v_o(:,1) -> sol(:)
 do i=1,size(v,1)
   sol(i)=v_o(i,1)
 end do

END SUBROUTINE GESV8

SUBROUTINE GESV4(x,v,sol)
USE LAPACK90
integer::i
real(r4)::x(:,:),v(:),sol(:)
real(r8),allocatable::v_o(:,:)
real(r8),allocatable::x_o(:,:)
allocate(x_o(size(x,1),size(x,1)))
x_o=x

!convert vector v(:) -> v_o(:,1)
allocate(v_o(size(v,1),1))
 do i=1,size(v,1)
   v_o(i,1)=v(i)
 end do

call dgesv_f90(x_o,v_o)

!rewrite vector v_o(:,1) -> v(:)
 do i=1,size(v,1)
   sol(i)=v_o(i,1)
 end do
END SUBROUTINE GESV4

SUBROUTINE SYEV8(x,v,d)
USE LAPACK90
real(r8)::x(:,:),v(:),d(:,:)
d=x
call dsyev_f90(d,v,'v')
END SUBROUTINE SYEV8

SUBROUTINE SYEV4(x,v,d)
USE LAPACK90
real(r8), allocatable:: x_o(:,:), v_o(:)
real(r4)::x(:,:),v(:),d(:,:)
allocate (x_o(size(x,1),size(x,1)))
allocate (v_o(size(v,1)))
x_o=x
v_o=v
call dsyev_f90(x_o,v_o,'v')
d=x_o
v=v_o
END SUBROUTINE SYEV4

SUBROUTINE PRINTM4 (matrix,text,fmt,un)
! prints matrix and text(optionaly) using form (optionally) to standard
! output or a unit(optionally)
character(*), optional::text*(*),fmt
character(40)::fmt1
integer :: n,i
integer,optional::un
real (r4):: matrix(:,:)
if (present(text)) then
   if (present(un)) then
      write(un,'(a)')text
     else
      print*,text
    end if
endif
if (present(fmt)) then; fmt1=fmt;else
fmt1= '(10g12.4/)' ;end if
n=(size(matrix,1))
do i=1,n
   if (present(un)) then
         write(un,fmt1)matrix(i,1:n)
      else
         print fmt1,matrix(i,1:n)
   endif 
enddo
END SUBROUTINE PRINTM4

SUBROUTINE PRINTM8 (matrix,text,fmt,un)
! prints matrix and text(optionaly) using form (optionally) to standard
! output or a unit(optionally)
character(*), optional::text*(*),fmt
character(40)::fmt1
integer :: n,i
integer,optional::un
real (r8):: matrix(:,:)
if (present(text)) then
   if (present(un)) then
         write(un,'(a)')text
      else
         print*,text; 
   endif
end if
if (present(fmt)) then; fmt1=fmt;else
fmt1= '(10g12.4/)';end if
n=(size(matrix,1))
do i=1,n
   if (present(un)) then
         write(un,fmt1)matrix(i,1:n)
      else
         print fmt1,matrix(i,1:n)
   endif
enddo
END SUBROUTINE PRINTM8

SUBROUTINE PRINTM4P (matrix,text,fmt,un)
! prints matrix and text(optionally) using form (optionally) to standard
! output or a unit(optionally)
character(*), optional::text*(*),fmt
integer,optional::un
integer :: n,i
real (r4):: matrix(:)
real (r4), allocatable::matrix_f(:,:)
n=(-1+sqrt(1.0+8*size(matrix,1)))/2
allocate(matrix_f(n,n))
call unpackit(matrix,matrix_f)
call printmat(matrix_f,text,fmt,un)
END SUBROUTINE PRINTM4P

SUBROUTINE PRINTM8P (matrix,text,fmt,un)
! prints matrix and text(optionally) using form (optionally) to standard
! output or a unit(optionally)
character(*), optional::text*(*),fmt
integer,optional::un
integer :: n,i
real (r8):: matrix(:)
real (r8), allocatable::matrix_f(:,:)
n=(-1+sqrt(1.0+8*size(matrix,1)))/2
allocate(matrix_f(n,n))
call unpackit(matrix,matrix_f)
call printmat(matrix_f,text,fmt,un)
END SUBROUTINE PRINTM8P


SUBROUTINE PACKIT4(fullmx,packedmx)
!conversion n*n matrix to triangular storage
integer::i,j,n
real(r4)::fullmx(:,:),packedmx(:)
n=size(fullmx,1)
 do i=1,n
  do j=1,i
    packedmx(pos_symm(i,j,n))=fullmx(i,j)
  end do
 end do
END SUBROUTINE PACKIT4

SUBROUTINE PACKIT8(fullmx,packedmx)
!conversion n*n matrix to triangular storage
integer::i,j,n
real(r8)::fullmx(:,:),packedmx(:)
n=size(fullmx,1)
 do i=1,n
  do j=1,i
    packedmx(pos_symm(i,j,n))=fullmx(i,j)
  end do
 end do
END SUBROUTINE PACKIT8


SUBROUTINE UNPACKIT4(packedmx,fullmx)
!conversion from packed to quadratic form
integer::i,j,n
real(r4)::fullmx(:,:),packedmx(:)
fullmx=0
n=size(fullmx,1)
 do i=1,n
  do j=1,i
     fullmx(i,j)=packedmx(pos_symm(i,j,n))
     fullmx(j,i)=packedmx(pos_symm(i,j,n))
  end do
 end do
END SUBROUTINE UNPACKIT4

SUBROUTINE UNPACKIT8(packedmx,fullmx)
!conversion from packed to quadratic form
integer::i,j,n
real(r8)::fullmx(:,:),packedmx(:)
fullmx=0
n=size(fullmx,1)
 do i=1,n
  do j=1,i
     fullmx(i,j)=packedmx(pos_symm(i,j,n))
     fullmx(j,i)=packedmx(pos_symm(i,j,n))
  end do
 end do
END SUBROUTINE UNPACKIT8

function diag_mat_vec8(x) result (y)
! y=diag(x) if x is a matrix and y a vector
real (r8)::x(:,:)
real (r8)::y(size(x,1))
integer::i
!
do i=1,size(y)
   y(i)=x(i,i)
enddo
end function   

function diag_mat_vec4(x) result (y)
! single-precision version
real (r4)::x(:,:)
real (r4)::y(size(x,1))
integer::i
!
do i=1,size(y)
   y(i)=x(i,i)
enddo
end function   



function diag_vec_mat8(x) result (y)
! y=diag(x) if x is a vector and y a matrix
real (r8)::x(:)
real (r8)::y(size(x),size(x))
integer::i
!
y=0
do i=1,size(x)
   y(i,i)=x(i)
enddo
end function   


function diag_vec_mat4(x) result (y)
! single-precision version
real (r4)::x(:)
real (r4)::y(size(x),size(x))
integer::i
!
y=0
do i=1,size(x)
   y(i,i)=x(i)
enddo
end function   

  

RECURSIVE SUBROUTINE POS_DEF(x,text,min_eig,stat,showeig)
! Ensures that symmetric X is positive definite with condition number
! >= 1/min_eig; ignores columns/rows with 0 values.
!
! Decomposes X = VDV', where D are eigenvalues and V eigenvectors while 
! ignoring rows/columns with 0 only. 
!
! If min(d) < min_eig * max(d), then .
!	a) if optional variable text is present, it is printed,

!	c) if optional min_eig is present, then:
!		if |d(i)| < min_eig * max(d) then d(i)=min_eig * max(d)
!	   if optional min_eig is missing , then:
!		if |d(i)| < min_eig_const * max(d) then &
!                                        d(i)=min_eig_const * max(d)
!	d) X = VDV'					
!	e) if optional stat is present, it is set to .true.
!       f) if optional showeig is present print eigenvalues (IA)
!

real (r8)::x(:,:)
character (*),optional::text
real (r8), optional::min_eig
real (r8)::v(size(x,1),size(x,1)),d(size(x,1))
real (r8),allocatable::y(:,:)
integer::nzero(size(x,1)),nnz,n,i
logical,optional::stat,showeig
real (r8), parameter::min_eig_const=1e-5
real (r8) :: min_eig_c
!
if (present(stat)) stat=.false.

! if some rows/columns zero, skip them
 n=size(x,1); nnz=0
 do i=1,n
    if (sum(x(i,:)**2)+sum(x(:,i)**2) == 0) cycle
     nnz=nnz+1
     nzero(nnz)=i
 enddo
 
 if (nnz .ne. n) then
    allocate(y(nnz,nnz))
    y=x(nzero(1:nnz),nzero(1:nnz))
    call pos_def(y,text,min_eig,stat,showeig)    
    x(nzero(1:nnz),nzero(1:nnz))=y
    deallocate(y)
    return
 endif 

call eigen(x,d,v)

! Set minimum "correct" eigenvalus to a small negative number rather than
! 0; subroutine eigen may return small negative eigenvalue(s) for
! semi-positive matrices. 
! min_eig = relative eigenvalues = eigenvalues / maximum eigenvalues
    min_eig_c=min_eig_const
    if (present(min_eig)) min_eig_c=min_eig
  if (minval(d) < min_eig_c) then
     if (present(showeig))  print '("  eigenvalues"/100(9g13.5))',d	
     if (present(stat)) stat=.true.
     if (present(min_eig)) then
         where (d/maxval(d) < min_eig) d= min_eig * maxval(d)
       else
         where (d/maxval(d) < min_eig_const) d= min_eig_const * maxval(d)
     endif 
     if (present(text)) print*,text
     if (present(showeig))  print '("  eigenvalues after"/100(9g13.5))',d
     x=0
     x=matmul(matmul(v,diag(d)),transpose(v))
  else
    if (present(showeig))  print '("  eigenvalues"/100(9g13.5))',d
  endif
  
  ! check for compiler bug in Absoft f90 version 6
  !if (maxval(matmul(v,transpose(v))-eye(v)) > 1d-6) then
  !   print*,'eigenvalues incorrect; recompile lapack90r without optimization'
  !   stop
  !endif   
end subroutine  


  function kron8(x,y) result (z)
! z = x Kronecker_product y for real r8

  real (r8)::x(:,:),y(:,:)
  real (r8)::z(size(x,dim=1)*size(y,dim=1), size(x,dim=2)*size(y,dim=2))
  integer::i,j,nx1,nx2,ny1,ny2
 
  nx1=size(x,dim=1); nx2= size(x,dim=2)
  ny1=size(y,dim=1); ny2= size(y,dim=2)
 
  do i=1,ny1
     do j=1,ny2
     !print*,'k',i,j,1+(i-1)*nx1,i*nx1,1+(j-1)*nx2,j*nx2,x,y(i,j)
        z(1+(i-1)*nx1:i*nx1,1+(j-1)*nx2:j*nx2)=x*y(i,j)
     enddo
  enddo
  end function   	


 function kron4(x,y) result (z)
 !z = x Kronecker_product y for real r4

 real (r4)::x(:,:),y(:,:)					       
 real (r4)::z(size(x,dim=1)*size(y,dim=1), size(x,dim=2)*size(y,dim=2))
 integer::i,j,nx1,nx2,ny1,ny2

 nx1=size(x,dim=1); nx2= size(x,dim=2)
 ny1=size(y,dim=1); ny2= size(y,dim=2)

 do i=1,ny1
    do j=1,ny2
    !print*,'k',i,j,1+(i-1)*nx1,i*nx1,1+(j-1)*nx2,j*nx2,x,y(i,j)
       z(1+(i-1)*nx1:i*nx1,1+(j-1)*nx2:j*nx2)=x*y(i,j)
    enddo
 enddo
 end function


function mat_mul_vec_symm8(x,y) result (z)
! z = X*y where X symm in packed storage for r8
real(r8) :: x(:),y(:),z(size(y))
integer i,j,k,n
n=size(y)
z=0
k=0
i=1
do
   z(i) = z(i) + x(i+k)*y(i)
   do j=i+1,n
      z(i) = z(i) + x(j+k)*y(j)
      z(j) = z(j) + x(j+k)*y(i)
   enddo
i=i+1
k=k+n-i+1
if (i>n) exit
enddo
end function


function mat_mul_vec_symm4(x,y) result (z)
! z = X*y where X symm in packed storage for r4
real(r4) :: x(:),y(:),z(size(y))
integer i,j,k,n
n=size(y)
z=0
k=0
i=1
do
   z(i) = z(i) + x(i+k)*y(i)
   do j=i+1,n
      z(i) = z(i) + x(j+k)*y(j)
      z(j) = z(j) + x(j+k)*y(i)
   enddo
i=i+1
k=k+n-i+1
if (i>n) exit
enddo
end function

function checksym4(x) result (a) 
! check that matrix x is symmetric
real(r4) :: x(:,:) 
integer :: i,j
logical :: a 
!
a=.true. 
do i=1,size(x,1)
   do j=i,size(x,1)
      if (x(i,j) /= x(j,i) ) then 
         a = .false.
         exit 
      endif 
   enddo 
enddo 
end function

function checksym8(x) result (a) 
! check that matrix x is symmetric
real(r8) :: x(:,:) 
integer :: i,j
logical :: a 
!
a=.true. 
do i=1,size(x,1)
   do j=i,size(x,1)
      if (x(i,j) /= x(j,i) ) then 
         a = .false.
         exit 
      endif
   enddo
enddo
end function
END MODULE DENSEOP

