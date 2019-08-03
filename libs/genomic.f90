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
!use omp_lib
use kinds; use model; use textop; use sparsem
use sparseop; use denseop; use prob; use pcg

implicit none
! 
character :: version*5= "1.000"
real(r8) ::  wGimA22i =  1.d0,&     
             wGi      =  1.d0,&    ! tau
             wA22i    = -1.d0,&    ! omega
             wG       =  0.95,&    ! alpha
             wA22     =  0.05      ! beta
logical :: sameweights
integer :: whichH=0

!
! storage matrices
!
type(densem),save :: A22i !,Gi, GimA22i
type(sparse_ija), save :: Gi_ija, A22i_ija, GimA22i_ija
type(sparse_hashm), allocatable :: Hinv, ainv(:)

integer :: io_g=10, io_i, io, i, j, l, begSNP, lastSNP, len_gen
integer :: num_snp=80000, len_char=50
integer :: unit_xref = 123
integer, allocatable :: xref(:)

real (rh) :: temp(1),t1,t2,val

real (rh), allocatable :: gen(:,:), W(:,:), Z(:,:), k(:), Gi(:,:), GimA22i(:,:)

character (len=100) :: a,          &
                       gen_file,   & ! = "genotypes.gen",&
                       g_xref,     & !="genotypes.gen_XrefID",&
                       ped_file      ! = "pedigree.ped"






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
   !print*,          '  * Not available in this distribution'
   
   gen_file  = "genotypes.gen_clean"
   g_xref    =  trim(gen_file) // "_XrefID"  !"genotypes.gen_XrefID"
   !ped_file = "pedigree.ped"
   !print*, g_xref
   !eff=1
   
   open(unit_xref, file=g_xref)
   open(io_g, file=gen_file)

   call read_xref_file(len_gen)
   call read_first_and_last_snp(len_gen,begsnp,lastsnp)

   print '(a35,i10)','Number of genotyped individuals = ',len_gen
   print '(a35,i10)','Position of first SNP = ',begsnp
   print '(a35,i10)','Number of genotypes = ',lastsnp
   print*,''
   
   allocate( gen(len_gen,lastsnp) )
   call create_Gi_dense_symm(gen,len_gen,begsnp,lastsnp)
   call createGmA22_densem(X=Gi,Y=ainv(neff),XmY=GimA22i,wx2=wGi,wy2=wA22i)
   call create_Hinv(n=lenped(eff), eff=eff)

end subroutine


subroutine create_Hinv(n,eff)
! ! n = number of genotyped individuals
! ! X = Hinv / sparse
! ! Y = Ainv / sparse
! ! Z = GimA22 / dense

integer :: n, i, j, eff, lengen
real(rh) :: val

!allocate( Hinv(n )
call zerom(Hinv,n)

Hinv=ainv(eff)
lengen = size(GimA22i,1)

print*,'Creating Hinv matrix'
do i=1,lengen
   do j=1,lengen
      val = GimA22i(i,j)
      call addm(val, xref(i), xref(j), Hinv)
   enddo
enddo

print*,'Saving Hinv...'
open(101012,file="Hinv",status='replace')
do i=1,n
   do j=i,n
      write(101012,fmt="(2i,f11.3)") xref(i), xref(j), getm(i,j,Hinv)
   enddo
enddo
close(101012)

end subroutine



subroutine dense_inversion_mkl(A,n)
   integer :: n
   real :: t1,t2
   real(rh), intent(inout) :: A(:,:)
   real(rh) :: WORK(n)
   integer :: IPIV(n)
   integer :: info,error

   print*,''
   print '(a)', '------------------------------------------------------------'
   print '(a)', 'Inversion performed using DGETR subroutines from MKL Library'
   print '(a)', '------------------------------------------------------------'
   print*,''

   !print '(a)', 'Factoriz G matrix...'
   call cpu_time(t1)
   call dgetrf(n,n,A,n,IPIV,info)
   if(info .ne. 0) then
      write(*,*) "Factorization failed!"
      stop
   end if
   call cpu_time(t2)
   print '(a35,f10.2)', 'Time to factorize = ',t2-t1

   call cpu_time(t1)
   call dgetri(n,A,n,IPIV,WORK,n,info)
   if(info .ne. 0) then
      write(*,*) "Inversion failed!"
      stop
   end if
   call cpu_time(t2)
   print '(a35,f10.2)', 'Time for inversion = ',t2-t1
end subroutine



subroutine create_Gi_dense_symm(x,n,pos1,pos2)
   implicit none
   real(rh) :: x(:,:)
   integer :: n, pos1, pos2, irank
   integer  :: i,j,l
   real (rh) :: val, ALPHA, BETA
   integer :: io_save=30

   allocate( w(0:2,pos2), Z(n,pos2), k(pos2) )
   allocate( Gi(n,n), GimA22i(n,n) )
   !call zerom(Gi,n)
   !call zerom(GimA22i,n)

   call read_markers(x, n, pos1, pos2)
   call createK(x,n)
   call createW(x,n,pos2)

   !open(io_save,file="G",status='replace')
   print '(a)', ' > Creating Z matrix...'
   call cpu_time(t1)
   Z=0
   do l=1,pos2
      Z(:,l)=W(x(:,l),l) ! centering column by gene content
   end do
   
   print '(a)', ' > Creating G matrix...'
   Gi = matmul(Z, transpose(Z))
   !CALL DGEMM('N','N',n,n,pos2,ALPHA,Gi,n,transpose(Gi),pos2,BETA,Gi,n)
   
   !do i=1,n
   !   do j=i,n
   !      do l=1,pos2
   !         !val = getm(i,j,Gi) + Z(i,l)*Z(j,l)
   !         !call setm(val,i,j,Gi)
   !         val = Gi(i,j) + Z(i,l)*Z(j,l)
   !         Gi(i,j) = val
   !      end do
   !   end do
   !end do

   do i=1,n
      do j=i,n
         !val = getm(i,j,Gi)/sqrt(k(i)*k(j))
         !call setm(val,i,j,Gi)
         val = Gi(i,j)/sqrt(k(i)*k(j))
         Gi(i,j) = val
         if(i .ne. j) Gi(j,i) = val
         !if (i<=j) write(io_save,fmt="(2i,f11.3)") i,j, getm(i,j,Gi)
      end do
   end do
   !call cpu_time(t2)
   !close(io_save)
   !print '(a35,f10.2)', 'Time to build G = ',t2-t1

   !print*, 'Scaling G matrix'
   !print '(a)', ' > Scaling G matrix...'
   !call cpu_time(t1)
   do i=1,n
       !do j=1,len_gen
       !    val = getm(i,j,Gi) * wG
       !    call setm(val,i,j,Gi)
       !    if (i==j) then
       !     val = getm(i,j,Gi) + wA22
       !     call setm(val,i,j,Gi)
       !    endif
       !enddo
       !val = getm(i,i,Gi) + 0.001
       !call setm(val,i,i,Gi)
       val = Gi(i,i) + 0.001
       Gi(i,i) = val
   enddo
   call cpu_time(t2)
   !print '(a35,f10.2)', 'Time to scale G = ',t2-t1
   print '(a35,f10.2)', 'Time to build G = ',t2-t1

   ! Inverting genomic matrix (ija format)
   !call cpu_time(t1)
   !Gi_ija = Gi
   !call fspak90('invert', Gi_ija)
   !Gi = Gi_ija
   !call cpu_time(t2)
   !print '(a35,f10.2)', 'Time to invert G = ',t2-t1
   !call reset(Gi_ija)
   !if(n <= 10) call printm(Gi)
   
   !do i=1,n
   !   print'(<n>f7.2)', Gi(i,:)
   !enddo
   !stop

   call dense_inversion_mkl(Gi,n)
   
   !print*,''
   !do i=1,n
   !   print'(<n>f7.2)', Gi(i,:)
   !enddo
   !stop
   
   print*,'Saving Gi...'
   open(io_save,file="Gi",status='replace')
   do i=1,n
      do j=i,n
         write(io_save,fmt="(2i,f11.3)") xref(i), xref(j), Gi(i,j) !getm(i,j,Gi)
      enddo
   enddo
   

end subroutine


! Function that compute variance of a vector
function var(x) result(c)
implicit none
real (rh) :: x(:),c,mu
integer :: i

! Mean
mu=sum(x)/size(x)
! Sum of Squares
do i=1,size(x)
  c=c+(x(i)-mu)**2
end do 
c=c/size(x)
end function


real function freq(values,n) result(x)
integer :: i,j,n
real (rh) :: f,values(:),np,npq,nmiss=0
np = count(values==0.)
npq = count(values==1.)
nmiss = count(values==5.)
x = (np + npq*.5)/(n-nmiss)
end function


! Subroutine to create a vector of 2pq
subroutine createk(m,n)
implicit none
real (rh) :: m(:,:), p, q, pq, val
integer :: i,j,n,io_save=20
double precision :: t1, t2


!open(io_save,file="allele_freq")
print '(a)', 'Creating scaling (k) vector'
!t1 = omp_get_wtime()
call cpu_time(t1)
do i=1,n
    p=0;val=0
    do j=1,size(m,2)
        p=freq(m(:,j),size(m,1))
        q=1-p
        !write(io_save,fmt="(5f15.3)"), p,q,p*q,2*p*q
        val=val+2*p*q
    enddo
    !k(i) = val
    k=val
    exit
enddo
!t2 = omp_get_wtime()
call cpu_time(t2)
print '(a35,f10.2)', 'Time to compute k  = ',t2-t1
close(io_save)
end subroutine


subroutine createW(m,n,snp)
implicit none
real (rh) :: m(:,:), val
integer :: n,snp,i

do i=1,snp
    val = freq(m(:,i),n)
    W(0,i)=-2*val; W(1,i)=1-2*val; W(2,i)=2-2*val
enddo
end subroutine


subroutine defaults_and_checks_GS()
   ! check optional parameters for Genomic Selection
   integer :: n
   character (100):: xc(50)
   real :: x(50)
   print*,''
   !do
      !call getoption('SNP_file',n,x,xc)
      !if (n > 0) then ! if SNP_file is present 
            call setup_genomic(1,0) 
      !      stop
      !endif ! if SNP_file present
      !
      !call getoption('end',n,x,xc)
      !if(n <= 0) exit
   !enddo
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
   ! X   = Gi;  wx2 = wGi
   ! Y   = A22; wy2 = wA22i
   ! XmY = GimA22i

   !type(densem) :: X,Y,XmY
   !type(densem) :: X,XmY
   type(sparse_hashm) :: Y

   real(rh) :: X(:,:)
   real(rh), intent(inout) :: XmY(:,:)
   real(r8) :: val,w,wx,wy
   real(r8), optional :: w2,wx2,wy2
   real :: t1
   integer  :: i,j,n

   print*,'Creating GimA22i matrix'
   do i=1,len_gen
      do j=i,len_gen         
         val = X(i,j)*wx2 + getm(xref(i),xref(j),Y)*wy2
         XmY(i,j) = val
         if(i/=j) XmY(j,i)=val
      enddo
   enddo

   print*,'Saving GimA22i...'
   open(1541,file="GimA22i",status='replace')
   do i=1,len_gen
      do j=i,len_gen
         write(1541, fmt="(2i5, f15.3)"), xref(i), xref(j), XmY(i,j)
      enddo
   enddo

endsubroutine

subroutine createGmA22_sparse_hashm(X,Y,XmY,w2,wx2,wy2)
   ! creates Matrix Adelta (wx*X + wy*Y)*w in hashm format
   type(densem) :: X,Y
   type(sparse_hashm) :: XmY
   integer  :: i,j,n
   real(r8) :: val,w,wx,wy
   real(r8), optional :: w2,wx2,wy2
   real :: t1
   !

endsubroutine

! Subroutine that computes the number of markers and samples
subroutine read_first_and_last_snp(n,beg,bot)
implicit none
character(len=num_snp) :: line
integer :: n, beg, bot, curr_pos

do i=1, n
  read(io_g, fmt='(a)', iostat=io) line
  if (i == 1) then
    beg = extract_snp(line, index( line(1:1000), " ")+1)
    bot = len(trim(line(beg:)))
  else
    curr_pos = extract_snp(line, index( line(1:1000), " ")+1)
    if (curr_pos /= beg) then
        print*, 'Position of first SNP is different! Check genotypic file!'
        print '( 3(a,i5) )', 'Position level     1 = ', beg, ' / Position level ',i,' = ', curr_pos
        stop
    endif
  endif
end do
rewind io_g
end subroutine read_first_and_last_snp


integer function extract_snp(x,curPos) result(updatePOS)
character(len_char) :: x
integer :: curPos
do
    if( x(curPos:curPos) /= " ") exit
    curPos = curPos+1
enddo
updatePOS = curPos
end function

! Read cross reference file
subroutine read_xref_file(n)
integer :: i,io,n,tmp(2)

xref=0; n=0
do
  read(unit_xref,*,iostat=io) i
  if( io .ne. 0) exit
  n=n+1
enddo
rewind unit_xref
allocate(xref(n))
do i=1,n
   read(unit_xref,*,iostat=io) xref(i)
enddo

end subroutine

! Subroutine that reads genotype file as character and converts it into real
subroutine read_markers(x,n,first_snp,nsnp)
implicit none

real (rh) :: x(:,:)
integer :: n,i,j,pos,snp,first_snp,nsnp
character(len=num_snp) :: line
character(len=1) :: a

print '(a)', 'Reading marker data...'
call cpu_time(t1)
do i=1,n
    !if(mod(i,1000)==0) print*,i
    pos=first_snp
    read(io_g,fmt='(a)') line
    do j=1,nsnp    
        a=line(pos:pos)
        read(a,'(i)') snp
        gen(i,j)=snp
        pos=pos+1
    end do
enddo
call cpu_time(t2)
print '(a35,f10.2)', 'Time to read genotypic data = ',t2-t1
end subroutine



end module
