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
character :: out_fmt*11 = '(2i,f20.12)'
real(r8) ::  wGimA22i =  1.d0,&     
             wGi      =  1.d0,&    ! tau
             wA22i    = -1.d0,&    ! omega
             wG       =  1.00,&    ! alpha
             wA22     =  0.00,&    ! beta
             wGamma   =  0.05      ! gamma

! Options read from parameter file:
!  * SNP file: genotypes.gen_clean
!  * SNP Xref file: genotypes.gen_clean_XrefID
!  * Save all G Matrices (default=.false.)
!  * Save G Inverse matrix (default=.false.)
!  * Save GmA22 matrix (default=.false.)
!  * Save H Inverse matrix (default=.false.)
!  * No Quality Control Checks !!!!! (default .false.):  T
!  * AlphaBeta defaults alpha=0.95, beta=0.05) :  1.00  0.00
!  * GammaDelta (defaults gamma=0.0, delta=0.0) :  0.05  0.00
!  * Matrix in Ascii format(default=binary)

logical :: sameweights, saveG=.true., saveGi=.true., saveGimA22i=.true., saveHinv=.true., saveBinary=.true.
logical, allocatable :: is_ancestral(:)
integer :: whichH=0
integer :: io_save=30123456


! storage matrices
! type(densem),save :: A22i !,Gi, GimA22i
type(sparse_ija), save :: Gi_ija, Aija, A22i_ija, GimA22i_ija
type(sparse_hashm), allocatable :: ainv(:)
type(sparse_hashm) :: Hinv

integer :: io_g=10, io_i, io, i, j, l, begSNP, lastSNP, lengen
integer :: num_snp=80000, len_char=50
integer :: unit_xref = 123
integer, allocatable :: xref(:)

real (rh) :: temp(1),t1,t2,val

real (rh), allocatable :: gen(:,:), W(:,:), Z(:,:), k(:), Gi(:,:), GimA22i(:,:), A22(:,:), A22i(:,:)

character (len=100) :: a,          &
                       gen_file,   & ! = "genotypes.gen",&
                       g_xref,     & ! = "genotypes.gen_XrefID",&
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

interface createGimA22i
 module procedure createGimA22i_densem !,&
                  !createGmA22_dense_symm ,&
                  !createGmA22_sparse_hashm
end interface

interface print_basic_stats
 module procedure print_basic_stats_dense_symm, &
                  print_basic_stats_hashm,&
                  print_basic_stats_ija
end interface

interface check_symmetry
 module procedure check_symmetry_dense,&
                  check_symmetry_densem,&
                  check_symmetry_hashm,&
                  check_symmetry_ija
end interface

interface save_matrix
 module procedure save_dense_,&
                  save_hashm,&
                  save_ija
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
   real(rh) :: val
   
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

   open(unit_xref, file=g_xref)
   open(io_g, file=gen_file)

   call read_xref_file(lengen)
   call read_first_and_last_snp(lengen,begsnp,lastsnp)

   print '(a35,i10)','Number of genotyped individuals = ',lengen
   print '(a35,i10)','Position of first SNP = ',begsnp
   print '(a35,i10)','Number of genotypes = ',lastsnp
   print*,''

   allocate( gen(lengen,lastsnp), A22(lengen,lengen), A22i(lengen,lengen) )


   call get_A22
   call check_symmetry(A22, 'Statistics of A22 Matrix')
   call print_basic_stats(A22, 'A22')
   A22i = A22
   print*,''
   call dense_inversion_mkl(A22i, lengen)
   call check_symmetry(A22i, 'A22i')
   call print_basic_stats(A22i, 'Statistics of A22 Inverse Matrix')
   print*,''
   ! stop
   
   call create_Gi_dense_symm(gen, lengen, begsnp, lastsnp, ainv(eff))
   call createGimA22i(eff)
   call create_Hinv( maxval(lenped), eff )

   if (saveGi) call save_matrix(Gi, 'Gi', 'Gi.txt')
   if (saveGimA22i) call save_matrix(GimA22i, 'GimA22i', 'GimA22i.txt')

   if (saveHinv) then
      print '(a)', 'Saving Hinv'
      open(io_save, file="Hinv.txt", status='replace')
      do i=1,lenped(eff)
         do j=1,lenped(eff)
            val = getm(i, j, Hinv)
            if( val /= 0.0) write(io_save, fmt=out_fmt) i, j, val
         enddo
      enddo
      close(io_save)
   endif
end subroutine



subroutine set_A22_orig_ped(perm,iperm,xref)

real(rh) :: val
real(rh), allocatable :: tmp(:,:)
integer, allocatable :: genid(:)
integer :: perm(:), iperm(:), xref(:), i, j, n, io, jo

n = size(A22,1)
allocate( tmp(n,n), genid(n) )

genid = perm(xref)
do i=1,n
   io = findloc(xref, iperm(genid(i)), dim = 1)
   do j=i,n
      jo = findloc(xref, iperm(genid(j)), dim = 1)
      ! print*, genid(j), iperm(genid(j)), jo
      ! if ( j > 2) stop
      tmp(io,jo) = A22(i,j)
      if (i /= j) tmp(jo,io) = tmp(io,jo)
   enddo
enddo
A22 = tmp

deallocate(tmp)



end subroutine



subroutine get_A22

integer :: i, j, top, bot, id, sire, dam, miss=0, Beg, curr, nval, pos, n
integer, allocatable :: gen(:), tgen(:), perm(:), iperm(:), igen(:), indices(:)
logical, allocatable :: nindex(:)
integer :: ped(maxval(lenped) , 3)


!
! reads pedigree (again) -> Need to solve in a better way!
n = maxval(lenped)
allocate( gen(n),tgen(n),perm(n),iperm(n),igen(n),indices(n), nindex(n) )

gen=1;tgen=1;perm=1;iperm=1;indices=0
nindex=.false.

rewind(io_off)
do i=1, n
   read(io_off,*,iostat=io) ped(i,:)
enddo

! Step that computes 'generic' indices for each animal based on progeny. This is useful
! to sort individuals
print'(a)', 'Building A22 '
print '(a)', "   Calculating A22 Matrix by Colleau"
do
    top=sum(gen)
    do i=1,n

        id=ped(i,1)
        sire=ped(i,2)
        dam=ped(i,3)

        !sire
        if(sire>0) then
            if(gen(id) >= gen(sire)) gen(sire)=gen(sire)+1
        endif
        !dam
        if(dam>0) then
            if(gen(id) >= gen(dam)) gen(dam)=gen(dam)+1
        endif
    enddo
    bot=sum(gen)
    if(top==bot) exit
    miss=miss+1
    if(miss>50) then
        print*,"Circular pedigree!"
        stop
    endif
enddo
tgen=gen(ped(:,1))

!
! Organize indices to be in needed order
!
! print '(a)',"Sorting indices..."
Beg=maxval(tgen)
curr=1
do
    if(all(nindex)) exit
    nval=count(tgen==beg); !print*,nval
    do i=1,nval
        pos=findloc(tgen, VALUE=beg, DIM=1, MASK=nindex==.false.); !print*,pos
        nindex(pos)=.true.        
        indices(pos)=curr
        curr=curr+1
    enddo
    Beg=Beg-1
enddo

call sort_ped(ped,indices,n)

!
! Recode pedigree based on pre-defined order
do i=1,n
  id=ped(i,1)
  iperm(i)=id ! Original identification
  perm(id)=i  ! Recoded identification
enddo

call get_new_ped(ped,perm,n)
call compute_A22(ped,n,perm,iperm)

call set_A22_orig_ped(perm, iperm, xref)


end subroutine



! Algorithm from Ignacio JABS (2011)
subroutine A_times_v(w,v,ped,f)
integer :: i,s,d,n,ped(:,:)
real :: w(:),v(:),f(:),tmp,di
real, allocatable :: q(:)

n=size(w)
allocate(q(n))
q=0
do i=n, 1, -1
    q(i)=q(i)+v(i)
    s=ped(i,2); d=ped(i,3)
    if(s/=0) q(s)=q(s)+q(i)*0.5
    if(d/=0) q(d)=q(d)+q(i)*0.5
enddo

do i=1,n
    s=ped(i,2); d=ped(i,3)
    di=(count(ped(i,2:3)==0)+2)/4.d0 - 0.25*(f(s)+f(d))
    tmp=0
    if(s/=0)tmp=tmp+w(s)
    if(d/=0)tmp=tmp+w(d)
    w(i)=0.5*tmp
    w(i)=w(i)+di*q(i)
enddo
deallocate(q)
end subroutine




! Algorithm from Ignacio JABS (2011)
subroutine compute_A22(ped,n,perm, iperm)
real :: f(n),w(n),v(n)
integer :: i,irank,ped(:,:),ngen,n,perm(:),iperm(:)
integer, allocatable :: genid(:)

f=0;
A22=0

allocate( genid(lengen) )
ngen=size(A22,1)
genid = perm(xref)

do i=1,ngen
    v=0;w=0
    v(genid(i))=1
    call A_times_v(w,v,ped,f)
    A22(:,i) = w(genid)
enddo
end subroutine


subroutine get_new_ped(ped, perm, n)
integer, intent(inout) :: ped(:,:), perm(:)
integer :: i,n
do i=1,n
    ped(i,1)=perm(ped(i,1))
    if(ped(i,2)>0) ped(i,2)=perm(ped(i,2))
    if(ped(i,3)>0) ped(i,3)=perm(ped(i,3))
enddo
end subroutine


subroutine sort_ped(ped,index, n)
integer, intent(inout) :: ped(:,:), index(:)
integer :: pedt(n,3),i,n
pedt=0
do i=1,n
    pedt(index(i),:)=ped(i,:)
enddo
ped=pedt
end subroutine



subroutine check_symmetry_densem(x,fn)
implicit none
type(densem) :: x
integer :: i,j,n
real(rh) :: val
character(*) :: fn

n=x%n
do i=1,n
   do j=1,n
      if( getm(i,j,x) /= getm(j,i,x) ) then
         print '(2a)', fn,'is not symetric'
         print '(a,2i,2f10.4)', '', i, j, getm(i,j,x), getm(j,i,x)
      endif
   enddo
enddo
end subroutine 


subroutine check_symmetry_hashm(x,fn)
implicit none
type(sparse_hashm) :: x
integer :: i,j,n
real(rh) :: val
character(*) :: fn

n=x%n
do i=1,n
   do j=1,n
      if( getm(i,j,x) /= getm(j,i,x) ) then
         print '(2a)', fn,'is not symetric'
         print '(a,2i,2f10.4)', '', i, j, getm(i,j,x), getm(j,i,x)
      endif
   enddo
enddo
end subroutine 


subroutine check_symmetry_ija(x,fn)
implicit none
type(sparse_ija) :: x
integer :: i,j,n
real(rh) :: val
character(*) :: fn

n=x%n
do i=1,n
   do j=1,n
      if( getm(i,j,x) /= getm(j,i,x) ) then
         print '(2a)', fn,'is not symetric'
         print '(a,2i,2f10.4)', '', i, j, getm(i,j,x), getm(j,i,x)
      endif
   enddo
enddo
end subroutine 


subroutine check_symmetry_dense(x,fn)
implicit none
real(rh) :: x(:,:), val
integer :: i,j,n
character(*) :: fn

n=size(x,1)
do i=1,n
   do j=1,n
      if( x(i,j) /= x(j,i) ) then
         print '(2a)', fn,' is not symmetric'
         print '(a,2i,2f10.4)', '', i, j, x(i,j), x(j,i) 
      endif
   enddo
enddo
end subroutine


subroutine save_hashm(x,text,fn)
implicit none
type(sparse_hashm) :: x
character(*) :: text, fn
integer :: i,j,n,k

n = x%n
print '(2a)', 'Saving ', text
print*,''

k=index(fn,'.')
open(io_save, file=fn, status='replace')
do i=1,n
   do j=i,n
      val = getm(i,j,x)
      if( val .ne. 0.0) write(io_save,fmt="(2i,f25.15)") i, j , val  ! saves values in the same order as Xref file
   enddo
enddo
close(io_save)
end subroutine


subroutine save_ija(x,text,fn)
implicit none
type(sparse_ija) :: x
character(*) :: text, fn
integer :: i,j,n,k

print '(2a)', 'Saving ', text
print*,''

n = x%n
k=index(fn,'.')
open(io_save, file=fn, status='replace')
do i=1,n
   do j=i,n
      val = getm(i,j,x)
      if( val /= 0.0) write(io_save,fmt="(2i,f25.15)") i, j , val  ! saves values in the same order as Xref file
   enddo
enddo
close(io_save)
end subroutine


subroutine save_dense_(x,text,fn)
implicit none
real(rh) :: x(:,:)
character(*) :: text, fn
integer :: i,j,n,k


print '(2a)', 'Saving ', text
print*,''

n = size(x,1)
k=index(fn,'.')
open(io_save, file=fn, status='replace')
!print*,fn(1:k-1)
!if( saveBinary ) open(unit, file=fn(1:k-1), FORM='UNFORMATTED', access='stream')
do i=1,n
   do j=i,n
      val = x(i,j)
      if( val .ne. 0.0) write(io_save,fmt=out_fmt) i, j , val  ! saves values in the same order as Xref file
      !if( val .ne. 0.0) write(io_save,fmt="(2i,f25.20)") xref(i), xref(j), val
   enddo
enddo
close(io_save)
end subroutine


subroutine create_Hinv(n,eff)
integer :: n, i, j, eff, ngen
real(rh) :: val

call zerom(Hinv, n)
ngen = size(GimA22i,1)

print*, ''
print '(a)', 'Building H Inverse matrix'
print*, 'wGimA22i ', wGimA22i

! ainv elements
do i=1,n
   do j=i,n
      val = getm(i,j,ainv(eff))
      call addm(val,i,j,Hinv)
      !if (i/=j) call addm(val,j,i,Hinv)
   enddo
enddo

! ainv - gima22i*wGimA22i
do i=1,ngen
   do j=i,ngen
      val = GimA22i(i,j)*wGimA22i
      call addm(val, xref(i), xref(j), Hinv)
      if (i /= j) call addm(val, xref(j), xref(i), Hinv)
   enddo
enddo
call print_basic_stats(Hinv, 'Statistics of H Inverse Matrix')
end subroutine


subroutine dense_inversion_mkl(A,n)
   implicit none
   integer :: n
   real :: t1, t2
   real(rh), intent(inout) :: A(:,:)
   real(rh) :: WORK(n)
   integer :: IPIV(n)
   integer :: info,error,i,j
   
   ! print*,''
   ! print '(a)', 'Calculating G Inverse matrix'
   print '(a)', ' Inversion performed using dpotrf/i subroutines from MKL Library'

   call cpu_time(t1)
   ! call dgetrf(n,n,A,n,IPIV,info)
   call dpotrf ('L', n, A, n, info) ! Computes the Cholesky factorization of a symmetric positive-definite matrix
  
   if(info /= 0) then
      write(*,*) "Factorization failed!"
      print '(a,i)', 'info =',info
      stop
   end if
   call cpu_time(t2)
   print '(a35,f10.2)', 'Time to factorize = ',t2-t1

   call cpu_time(t1)
   ! call dgetri(n , A, n, IPIV, WORK, n, info)
   call dpotri( 'L', n, A, n, info) ! Computes the inverse of a symmetric positive-definite matrix
   if(info /= 0) then
      write(*,*) "Inversion failed!"
      print '(a,i)', 'info =',info
      stop
   end if
   call cpu_time(t2)
   print '(a35,f10.2)', 'Time for inversion = ',t2-t1

   ! Set Lower triangular values to upper triangular values. DPOTRI overwrittes lower triangle of the inverse of A
   do i=1,n
      do j=1,i
         if ( i /= j ) A(j,i) = A(i,j)
      enddo
   enddo
end subroutine



subroutine print_basic_stats_ija(x,name)
implicit none
type(sparse_ija) :: x
integer :: i,j,k,n,count
real(rh) :: val, max_val, min_val, std_val, sumSQR, sum,&
            max_val_off, min_val_off, std_val_off, sumSQR_off, sum_off
character(*) :: name

sum = 0
sumSQR = 0
max_val = 0
min_val = 1

n=x%n
! diagonal
do i=1, n
   val = getm(i,i,x)
   sum = sum + val
   max_val = max(max_val,val)
   min_val = min(min_val,val)
   sumSQR = sumSQR + val*val
enddo

sum_off = 0
sumSQR_off = 0
max_val_off = 0
min_val_off = 1
count=0
! off-diagonal
do i=1,n
   do j=1,n
      if (i==j) cycle
      val = getm(i,j,x)
      sum_off = sum_off + val
      max_val_off = max(max_val_off,val)
      min_val_off = min(min_val_off,val)
      sumSQR_off = sumSQR_off + val*val
      count=count+1
   enddo
enddo

print*,''
print '(a)', name
!print '(a)' ,  '==========================================='
print '(a15, 5a12)', '','N','Min','Mean','Max','Var'
print '(a15, i12, 4f12.3)',adjustl('Diagonal'), n, min_val, sum/n, max_val, (sumSQR - sum*sum/n)/(n-1)
print '(a15, i12, 4f12.3)',adjustl('Off-diagonal'), count, min_val_off, sum_off/count, max_val_off, (sumSQR_off - sum_off*sum_off/count)/(count-1)
!print '(a)' ,  '==========================================='
print*,''
end subroutine




subroutine print_basic_stats_hashm(x,name)

type(sparse_hashm) :: x
integer :: i,j,k,n,count
real(rh) :: val, max_val, min_val, std_val, sumSQR, sum,&
            max_val_off, min_val_off, std_val_off, sumSQR_off, sum_off
character(*) :: name

sum = 0
sumSQR = 0
max_val = 0
min_val = 1
n=x%n

! diagonal
do i=1, n
   val = getm(i,i,x)
   sum = sum + val
   max_val = max(max_val,val)
   min_val = min(min_val,val)
   sumSQR = sumSQR + val*val
enddo

sum_off = 0
sumSQR_off = 0
max_val_off = 0
min_val_off = 1
count=0
! off-diagonal
do i=1,n
   do j=1,n
      if (i==j) cycle
      val = getm(i,j,x)
      sum_off = sum_off + val
      max_val_off = max(max_val_off,val)
      min_val_off = min(min_val_off,val)
      sumSQR_off = sumSQR_off + val*val
      count=count+1
   enddo
enddo

print*,''
print '(a)', name
!print '(a)' ,  '==========================================='
print '(a15, 5a12)', '','N','Min','Mean','Max','Var'
print '(a15, i12, 4f12.3)',adjustl('Diagonal'), n, min_val, sum/n, max_val, (sumSQR - sum*sum/n)/(n-1)
print '(a15, i12, 4f12.3)',adjustl('Off-diagonal'), count, min_val_off, sum_off/count, max_val_off, (sumSQR_off - sum_off*sum_off/count)/(count-1)
!print '(a)' ,  '==========================================='
print*,''
end subroutine



subroutine print_basic_stats_dense_symm(x,name)

real(rh) :: x(:,:)
integer :: i,j,k,n,count
real(rh) :: val, max_val, min_val, std_val, sumSQR, sum,&
            max_val_off, min_val_off, std_val_off, sumSQR_off, sum_off
character(*) :: name

sum = 0
sumSQR = 0
max_val = 0
min_val = 1

n=size(x,1)
! diagonal
do i=1, n
   val = x(i,i)
   sum = sum + val
   max_val = max(max_val,val)
   min_val = min(min_val,val)
   sumSQR = sumSQR + val*val
enddo

sum_off = 0
sumSQR_off = 0
max_val_off = 0
min_val_off = 1
count=0
! off-diagonal
do i=1,n
   do j=1,n
      if (i==j) cycle
      val = x(i,j)
      sum_off = sum_off + val
      max_val_off = max(max_val_off,val)
      min_val_off = min(min_val_off,val)
      sumSQR_off = sumSQR_off + val*val
      count=count+1
   enddo
enddo


print*,''
print '(a)', name
!print '(a)' ,  '==========================================='
print '(a15, 5a12)', '','N','Min','Mean','Max','Var'
print '(a15, i12, 4f12.3)',adjustl('Diagonal'), n, min_val, sum/n, max_val, (sumSQR - sum*sum/n)/(n-1)
print '(a15, i12, 4f12.3)',adjustl('Off-diagonal'), count, min_val_off, sum_off/count, max_val_off, (sumSQR_off - sum_off*sum_off/count)/(count-1)
!print '(a)' ,  '==========================================='
!print*,''
end subroutine 



subroutine create_Gi_dense_symm(x,n,pos1,pos2,a)
   implicit none
   real(rh) :: x(:,:)
   !type(densem) :: a22i_densem
   type(sparse_hashm) :: a
   type(sparse_hashm) :: A22_to_check
   type(sparse_ija) :: A22_to_check_ija
   integer :: n, pos1, pos2, irank
   integer  :: i,j,l
   real (rh) :: val, ALPHA, BETA
   real(rh) :: sumDiag, sumOffDiag, nDiag, nOffDiag, MeanDiagG, MeanOffDiagG, MeanDiagA, MeanOffDiagA


   allocate( w(0:2,pos2), Z(n,pos2), k(pos2) )
   allocate( Gi(n,n), GimA22i(n,n) )
   
   ! Inverting ainv / inv( ainv )
   ! aija = ainv(1)
   ! call fspak90('invert', aija)
   
   call read_markers(x, n, pos1, pos2)
   call createK(x,n)
   call createW(x,n,pos2)

   print '(a)', 'Building genomic (G) matrix'
   print '(a)', ' > Creating Z matrix...'
   call cpu_time(t1)
   Z=0
   do l=1,pos2
      Z(:,l)=W(x(:,l),l) ! centering column by gene content
   end do
   
   print '(a)', " > Creating ZZ' matrix..."
   Gi = matmul(Z, transpose(Z))
   ! CALL DGEMM('N','N',n,n,pos2,ALPHA,Gi,n,Gi,pos2,BETA,Gi,n)

   print '(a)', ' > Scaling G matrix'
   do i=1,n
      do j=i,n
         val = Gi(i,j)/sqrt(k(i)*k(j))
         Gi(i,j) = val
         if (i /= j) Gi(j,i) = val
      end do
   end do   
   !call save_matrix(Gi, 'Gini', 'Gini.txt')
   call check_symmetry(Gi, 'Gini')


   if (.true.) then
      ! Calculating diag(A) and off-diag(A)
      sumDiag=0; sumOffDiag=0; nDiag=0; nOffDiag=0
      do i=1,lengen
         do j=i,lengen

            ! val = getm(xref(i),xref(j),aija)
            val = A22(i,j)
            if (val <= 0.) then
               val = 0.
            else
               val = nint(val * 1000.0) * 1E-3
            endif

            if (i == j) then
               sumDiag=sumDiag+val
               nDiag=nDiag+1
            else 
               sumOffDiag=sumOffDiag+val
               nOffDiag=nOffDiag+1
            endif
         enddo
      enddo
      MeanDiagA=sumDiag/nDiag
      MeanOffDiagA=sumOffDiag/nOffDiag
      
      ! Calculating diag(G) and off-diag(G)
      sumDiag=0; sumOffDiag=0; nDiag=0; nOffDiag=0
      do i=1,n
         do j=i,n
            val = Gi(i,j)
            if (i == j) then
               sumDiag=sumDiag+val
               nDiag=nDiag+1
            else 
               sumOffDiag=sumOffDiag+val
               nOffDiag=nOffDiag+1
            endif
         enddo
      enddo
      MeanDiagG=sumDiag/nDiag
      MeanOffDiagG=sumOffDiag/nOffDiag
      
      print '(a)', ' > Tuning G matrix'
      do i=1,n
         do j=i,n
            if (i == j) then
               val = MeanDiagG - MeanDiagA
               Gi(i,j) = Gi(i,j) - val
            else
               val = MeanOffDiagG - MeanOffDiagA
               Gi(i,j) = Gi(i,j) - val
               Gi(j,i) = Gi(i,j)
            endif
         end do
      end do

      print '(a,f8.4,a,f8.4)', '  - mean(diag(A22)) ', MeanDiagA, ' / diag(G) adjusted by', MeanDiagG - MeanDiagA
      print '(a,f8.4,a,f8.4)', '  - mean(off_diag(A22))', MeanOffDiagA, ' / off_diag(G) adjusted by', MeanOffDiagG - MeanOffDiagA      
   endif

   if (.true.) then
      ! G*alpha + A22*beta
      print*,''
      print '(a)', ' G*alpha + A22*beta'
      do i=1,n
         do j=i,n

            ! val = getm(xref(i), xref(j), Aija)
            val = A22(i,j)
            if (val <= 0.) then
               val = 0.
            else
               val = nint(val * 1000.0) * 1E-3
            endif

            val = Gi(i,j)*wG + val*wA22
            Gi(i,j) = val
            if (i /= j) Gi(j,i) = val
         enddo
      enddo
      call save_matrix(Gi,'G_Alpha_Beta','G_Alpha_Beta.txt')
      call check_symmetry(Gi, 'G_Alpha_Beta')
   endif

   if (.true.) then
      ! G + I*gamma
      print '(a)', ' G + I*gamma'
      do i=1,n
            val = Gi(i,i) + wGamma
            Gi(i,i) = val
      enddo
      call save_matrix(Gi,'G_Gamma','G_Gamma.txt')
      call check_symmetry(Gi, 'G_Gamma')
   endif

   call cpu_time(t2)
   print '(a35,f10.2)', 'Time to build G = ', t2-t1

   call print_basic_stats(Gi, 'Statistics of Genomic Matrix')
   
   print*,''
   if (saveG) call save_matrix(Gi, 'G', 'G.txt')
   call check_symmetry(Gi, 'G')

   !if (.false.) then
      ! checking max value in matrix / will be removed latter on
   !    print*,''
   !    print '(a)', 'Maxval pos tunning'
   !    val = maxval(Gi)
   !    print*,'looking for = ', val
   !    do i=1,n
   !       do j=1,n
   !          if( val == Gi(i,j) ) then
   !             print*,i,j,val
   !          endif
   !       enddo
   !    enddo
   ! endif

   ! if (.false.) then
   !    print '(a)', 'First 10 rows from Gi - pre mkl'
   !    do i=1, 10
   !       print '(10f7.3)', Gi(i,:10)
   !    enddo
   ! endif
   
   print '(a)', 'Calculating G Inverse matrix'
   call dense_inversion_mkl(Gi, n)
   
   ! if (.false.) then
   !    print '(a)', 'First 10 rows from Gi - pos mkl'
   !    do i=1, 10
   !       print '(10f7.3)', Gi(i,:10)
   !    enddo
   ! endif

   ! if (.false.) then
   !    ! checking max value in matrix / will be removed latter on
   !    print*,'Maxval pos mkl'
   !    val = maxval(Gi)
   !    print*,'looking for = ', val
   !    do i=1,n
   !       do j=1,n
   !          if( val == Gi(i,j) ) then
   !             print*,i,j,val
   !          endif
   !       enddo
   !    enddo
   ! endif

   call print_basic_stats(Gi, 'Statistics of Inv. Genomic Matrix')
   call check_symmetry(Gi, 'G Inverse')

end subroutine


! Function that compute variance of a vector
function var_dense_symm(x) result(c)
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
np = count(values==2.0)
npq = count(values==1.0)
nmiss = count(values==5.0)
x = (np + npq*.5)/(n-nmiss)
end function


! Subroutine to create a vector of 2pq
subroutine createk(m,n)
implicit none
real (rh) :: m(:,:), p, q, pq, val
integer :: i,j,n,io_save=20
double precision :: t1, t2


!open(io_save,file="sum2pq")
print '(a)', 'Creating scaling (k) vector'
!t1 = omp_get_wtime()
call cpu_time(t1)
do i=1,n
    p=0;val=0
    do j=1, size(m,2)
        p=freq(m(:,j),size(m,1))
        q=1-p
!        write(io_save,fmt="(5f15.3)"), p,q,p*q,2*p*q
        val = val + p*q
    enddo
    !k(i) = 2*val
    k=2*val
    exit
enddo
!t2 = omp_get_wtime()
call cpu_time(t2)
print '(a35, f10.2)', 'Time to compute k = ', t2-t1
print '(a35, f22.10)', 'Scale by Sum(2pq) = ', k(1)
!close(io_save)
end subroutine


subroutine createW(m,n,snp)
implicit none
real (rh) :: m(:,:), val
integer :: n, snp, i
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


subroutine createGimA22i_densem(eff)
   implicit none
   real(rh) :: val
   real :: t1, t2
   integer  :: i,j,n,eff

   n=lengen
   print*,''
   print'(a)','Building GimA22i matrix'
   do i=1,n
      do j=i,n
         ! val = getm(xref(i),xref(j),ainv(eff))*wA22i
         ! if (val <= 0.) val = 0.
         ! if (val > 0 ) val = nint(val * 1000.0) * 1E-3
         ! val = Gi(i,j)*wGi + getm(xref(i),xref(j),ainv(eff))*wA22i
         val = Gi(i,j)*wGi + A22i(i,j)*wA22i
         GimA22i(i,j) = val
         if (i/=j) GimA22i(j,i) = val
      enddo
   enddo
   call cpu_time(t2)
   print '(a35,f10.2)', 'Time to build GimA22i = ',t2-t1
   call print_basic_stats(GimA22i,'Statistics of Inv. Genomic - A22 Matrix')
endsubroutine


subroutine createGimA22i_sparse_hashm(X,Y,XmY,w2,wx2,wy2)
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
   !print*, xref(i)
enddo
!stop
end subroutine


! Subroutine that reads genotype file as character and converts it into real
subroutine read_markers(x,n,first_snp,nsnp)
implicit none

real (rh) :: x(:,:)
integer :: n,i,j,pos,snp,first_snp,nsnp
character(len=num_snp) :: line
character(len=1) :: a

print '(a)', 'Reading SNP file'
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
