!============================================================================================
module genomic
! Collection of subroutines and functions to implement  
!  single-step evaluation using full pedigree and genomic information
!
!  Methodology applied:
!  Misztal et al 2009, Legarra et al 2009 and Aguilar et al 2010
!
!  Ignacio Aguilar and Ignacy Misztal, University of Georgia
!   01/29/09 - 01/31/12
!============================================================================================
!
!  

use kinds; use model; use textop; use sparsem
use sparseop; use denseop; use prob; use pcg
implicit none
! 
character :: version*5= "1.108"
real(r8) ::  wGimA22i =  1.d0,&     
             wGi      =  1.d0,&    
             wA22i    = -1.d0      
logical :: sameweights
integer :: whichH=0             
!
! storage matrices
!
type(densem),save :: A22i,Gi,GimA22i


!=============================================================
! Overloading interfaces
!=============================================================
 
interface add_g_usr_gs
   module procedure add_g_usr_gs_hash_densem,&
                    add_g_usr_gs_hash_dense_symm,&
                    add_g_usr_gs_densem_densem,&
                    add_g_usr_gs_densem_dense_symm                   
end interface 

interface AddGimA22i_ainv
   module procedure AddGimA22i_ainv_densem,&
                    AddGimA22i_ainv_dense_symm
end interface

interface createGmA22
 module procedure createGmA22_densem ,&
                  createGmA22_dense_symm ,&
                  createGmA22_sparse_hashm
end interface
contains
!==================================================================
! Subroutines & functions
!==================================================================
  
subroutine setup_genomic(eff, unit_eff)
   ! 
   ! Create/Store/Read G matrix, A22 matrix and GmA22 matrix and its inverses
   !  
   implicit none
   integer :: eff,unit_eff
   !
   print*,''
   print '(a)'      ,' *--------------------------------------------------------------*'
   print '(3a)'     ,' *                 Genomic Library: Version ',version,'               *'
   print '(a)'      ,' *                                                              *' 
   print '(a,i3,a)' ,' *  Modified relationship matrix (H) created for effect: ',eff,'    *'
   print '(a)'      ,' *--------------------------------------------------------------*'     
   print*,''
   print*,          '  * Not available in this distribution'
  stop   
end subroutine

subroutine defaults_and_checks_GS()
   ! check optional parameters for Genomic Selection
   integer::n
   character (100):: xc(50)
   real :: x(50)
   print*,''
   do
      call getoption('SNP_file',n,x,xc)
      if (n > 0) then ! if SNP_file is present 
            call setup_genomic(0,0) 
            stop
      endif ! if SNP_file present
      
      call getoption('end',n,x,xc)
      if(n <= 0) exit
   enddo
end subroutine


!=========================================================
! Overloading Subroutines
!=========================================================
!
subroutine  AddGimA22i_ainv_densem(x,w,a)
   ! add a densem matrix to hash matrix with indexes from permutation vector genID
   type(densem) :: x
   type(sparse_hashm) :: a 
   real(r8) :: w ,val
   integer :: i,j
   !
   print*,'AddGimA22i_ainv_densem: not implemented'
   stop
end subroutine

subroutine  AddGimA22i_ainv_dense_symm(x,w,a)
   ! add a dense_symm matrix to hash matrix with indexes from permutation vector genID
   type(dense_symm) :: x
   type(sparse_hashm) :: a 
   real(r8) :: w 
   !
   print*,'AddGimA22i_ainv_dense_symm: not implemented'
   stop
end subroutine

!-------------------------------------------------------------
! add_g_usr_gs(xx,i,nlevsum,randomnumb(i),ntrait,Gi,1.d0)
!-------------------------------------------------------------

subroutine add_g_usr_gs_hash_densem(xx,eff,nlevsum,erandomnumb,entrait,A,weight)
   ! add user matrix in dense format to a matrix in hash format by correlated effects and traits
   type(densem):: A
   type(sparse_hashm):: xx
   integer :: nlevsum(:),eff,erandomnumb,entrait
   real(r8) :: weight,val
   integer :: i,j,k,l,m,n,t1,t2
   !
   print*,'add_g_usr_gs_hash_densem: not implemented'
   stop
end subroutine

subroutine add_g_usr_gs_hash_dense_symm(xx,eff,nlevsum,erandomnumb,entrait,A,weight)
   ! add user matrix in dense_symm format to matrix in hash format
   type(dense_symm):: A
   type(sparse_hashm):: xx
   integer :: nlevsum(:),eff,erandomnumb,entrait
   real(r8) :: weight
   !
   print*,'add_g_usr_gs_densem_densem: not implemented'
   stop
end subroutine

subroutine add_g_usr_gs_densem_densem(xx,eff,nlevsum,erandomnumb,entrait,A,weight)
   ! add user matrix in dense format to matrix in dense format
   type(densem):: A
   type(densem):: xx
   integer :: nlevsum(:),eff,erandomnumb,entrait
   real(r8) :: weight,val
   integer :: i,j,k,l,m,n,t1,t2
   !
   print*,'add_g_usr_gs_densem_densem: not implemented'
   stop
end subroutine

subroutine add_g_usr_gs_densem_dense_symm(xx,eff,nlevsum,erandomnumb,entrait,A,weight)
   ! add user matrix in densem format to matrix in densem_sym format
   type(dense_symm):: A
   type(densem):: xx
   integer :: nlevsum(:),eff,erandomnumb,entrait
   real(r8) :: weight
   !
   print*,'add_g_usr_gs_densem_dense_symm: not implemented'
   stop
end subroutine

!--------------------------------------------------------
! Create G minus A22  regular or inverse matrices
!--------------------------------------------------------
subroutine createGmA22_dense_symm(X,Y,XmY,w2,wx2,wy2)
   ! creates Matrix Adelta (wx*X + wy*Y)*w in dense_symm format
   !
   type(densem) :: X,Y
   type(dense_symm) :: XmY
   integer  :: i,j,n
   real(r8) :: val,w,wx,wy
   real(r8), optional :: w2,wx2,wy2
   real :: t1
   !
endsubroutine
subroutine createGmA22_densem(X,Y,XmY,w2,wx2,wy2)
   ! creates Matrix Adelta (wx*X + wy*Y)*w in densem format
   !
   type(densem) :: X,Y,XmY
   integer  :: i,j,n
   real(r8) :: val,w,wx,wy
   real(r8), optional :: w2,wx2,wy2
   real :: t1
endsubroutine

subroutine createGmA22_sparse_hashm(X,Y,XmY,w2,wx2,wy2)
   ! creates Matrix Adelta (wx*X + wy*Y)*w in hashm format
   ! 
   type(densem) :: X,Y
   type(sparse_hashm) :: XmY
   integer  :: i,j,n
   real(r8) :: val,w,wx,wy
   real(r8), optional :: w2,wx2,wy2
   real :: t1
   !
endsubroutine

end module
