program genomic
use mkl_service

implicit none


character :: version*5= "1.000"
character (len=100) :: a, gen_file = "../examples/ex1/deltag.gen_clean", gen_index="../examples/ex1/deltag_50K_imp.gen_clean_XrefID"
real     ::  wGimA22i =  1.d0,&     
             wGi      =  1.d0,&    
             wA22i    = -1.d0,&
             wOmega    =  0.95  
logical :: sameweights

integer :: whichH=0, io_g=10, io_i, io, i, j, l, BegSNP, LastSNP, nid
integer :: num_snp=80000, len_char=50
real, allocatable :: G(:,:), gen(:,:), W(:,:), Z(:,:), k(:), Gi(:,:)

! mkl
integer :: num_thread,nthrg

integer :: lwork,info,lda
real :: temp(1),t1,t2,val
character :: jobz='N',uplo='U'
real, allocatable :: ww(:),rand(:)
real, allocatable :: work(:,:)
logical :: wantz
INTEGER :: NTHREADS, TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM,&
           omp_get_max_threads
integer, allocatable :: ipiv(:,:)




open(io_g,file=gen_file)
call read_first_and_last_snp(nid,BegSNP,LastSNP)

print*,''
print*,'*******************************'
print*,'Genomic module - Version ',version
print*,'*******************************'

!$omp parallel
 !$ TID = OMP_GET_THREAD_NUM()
 !$ TID = omp_get_max_threads()
 !$ call omp_set_num_threads(TID)
!$omp end parallel
!$ print'(a35,i10,a)','OMP activated using = ',TID,' Threads'

print '(a35,i10)','Number of genotyped individuals = ',nid
print '(a35,i10)','Position of first SNP = ',BegSNP
print '(a35,i10)','Number of genotypes = ',LastSNP
print*,''


allocate( G(nid,nid), Gi(nid,nid), gen(nid,LastSNP), W(3,LastSNP), Z(nid,LastSNP), k(LastSNP) )
allocate( rand(nid) )
lda=max(1,nid)
allocate( ipiv(nid,1) )
allocate( work(nid,1) )


call read_markers(gen,nid,BegSNP,LastSNP)
call createk(gen,nid,k)
call createW(gen,W,nid,LastSNP)
call createG_dense_symm(gen,Z,G,W,k,nid,LastSNP)


! Store A in Ainv to prevent it from being overwritten by LAPACK
Gi=G

call cpu_time(t1)
! DGETRF computes an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges.
!call DGETRF(n, n, xinv, n, ipiv, info)
print '(a)', 'Decomposing genomic matrix...'
call DPOTRF( UPLO, nid, Gi, LDA, INFO )

if (info /= 0) then
 stop 'Matrix is numerically singular!'
end if

! DGETRI computes the inverse of a matrix using the LU factorization computed by DGETRF.
!call DGETRI(n, Gi, n, ipiv, work, n, info)
print '(a)', 'Inverting genomic matrix...'
call DPOTRI( UPLO, nid, Gi, LDA, INFO )
call cpu_time(t2)
print '(a35,f10.2)', 'Time to compute genomic inverse  = ',t2-t1

if (info /= 0) then
 stop 'Matrix inversion failed!'
end if

open(40, file="Gi", status='replace')
do i=1,nid
    do j=1,nid
        if(i<=j) write(40,fmt="(2i,f11.3)") i,j,Gi(i,j)
    enddo
end do



contains

subroutine createGi(G)
implicit none

real :: G(:,:)




end subroutine



real function freq(values,n) result(x)
integer :: i,j,n
real :: f,values(:),np,npq,nmiss=0
np = count(values==2.)
npq = count(values==1.)
nmiss = count(values==5.)
x = (np + npq*.5)/(n-nmiss)
end function


subroutine createk(m,n,k)
implicit none
real :: m(:,:), k(:), p,q,pq,val=0
integer :: i,j,n,io_save=20

open(io_save,file="allele_freq")
print '(a)', 'Creating scaling (k) vector'
call cpu_time(t1)
!do i=1,n
    p=0;val=0
    do j=1,size(m,2)
        p=freq(m(:,j),size(m,1)); q=1-p
        write(io_save,fmt="(5f15.3)"), p,q,p*q,2*p*q
        val=val+p*q
    enddo
    !k(i) = val
    k=2*val
    !write(io_save,fmt="(5f15.3)"), p,(1-p),p*(1-p),2*p*(1-p)!,k(i)
!enddo
call cpu_time(t2)
print '(a35,f10.2)', 'Time to compute k  = ',t2-t1
close(io_save)
end subroutine


subroutine createW(m,W,n,snp)
implicit none
real :: m(:,:), W(:,:), val
integer :: n,snp,i

do i=1,snp
    val = freq(m(:,i),n)
    W(1,i)=-2*val; W(2,i)=1-2*val; W(3,i)=2-2*val
enddo
end subroutine


subroutine createG_dense_symm(m,Z,G,w,k,n,snp)
implicit none
integer  :: i,j,n,l,snp
real :: val, W(:,:), m(:,:), G(:,:), Z(:,:), k(:)
integer :: io_save=30

open(io_save,file="G",status='replace')
print '(a)', 'Creating genomic matrix...'
call cpu_time(t1)
! ====================================================================
! Ignacio JABG paper (2011)
! ====================================================================
G=0
do l=1,snp
    Z(:,l)=W(M(:,l)+1,l) ! centering column by gene content
end do
do i=1,n
    do j=1,n
        do l=1,snp
            if (i==j) then
                G(i,j)=G(i,j)+Z(i,l)*Z(j,l)*wOmega+(1.-wOmega)*1.
            else
                G(i,j)=G(i,j)+Z(i,l)*Z(j,l)*wOmega
            endif
        end do
    end do
end do
do i=1,n
    do j=1,n
        G(i,j)= G(i,j)/sqrt(k(i)*k(j))
        if(i<=j) write(io_save,fmt="(2i,f11.3)") i,j,G(i,j)
    end do
    if(n < 11) print '(10f10.2)', G(i,:)
end do
! ====================================================================
call cpu_time(t2)
print '(a35,f10.2)', 'Time to construct G = ',t2-t1
close(io_save)
end subroutine


! Subroutine that reads genotype file as character and converts it into reale
subroutine read_markers(x,n,first_snp,nsnp)
implicit none

real :: x(:,:)
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


! Subroutine that computes the number of markers and samples
subroutine read_first_and_last_snp(n,beg,bot)
implicit none
character(len=num_snp) :: line
integer :: tmp1, tmp2, n, beg, bot, curr_pos

n=0
do
  read(io_g, fmt='(a)', iostat=io) line
  if(io /= 0) exit
  n=n+1
  if (n == 1) then
    beg = extract_snp(line, index( line(1:1000), " ")+1)
    bot = len(trim(line(beg:)))
  else
    curr_pos = extract_snp(line, index( line(1:1000), " ")+1)
    if (curr_pos /= beg) then
        print*, 'Position of first SNP is different! Check genotypic file!'
        print '( 3(a,i5) )', 'Position level     1 = ', beg, ' / Position level ',n,' = ', curr_pos
        stop
    endif
  endif
end do
rewind io_g
end subroutine read_first_and_last_snp


integer function extract_snp(x,curPos) result(UpdatePOS)
character(len_char) :: x
integer :: curPos
do
    if( x(curPos:curPos) /= " ") exit
    curPos = curPos+1
enddo
UpdatePOS = curPos
end function



end program