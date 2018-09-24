program t
!use mkl_service
implicit none

integer :: num_thread,i,j,nthrg
real(8), dimension(400,400) :: x,xinv

integer :: lwork,info,n,lda
real :: temp(1),t1,t2,val
character :: jobz='N',uplo='U'
real, allocatable :: w(:),rand(:)
real, allocatable :: work(:,:)
logical :: wantz
INTEGER :: NTHREADS, TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM,&
           omp_get_max_threads
integer, allocatable :: ipiv(:,:)

allocate( rand(size(x,1)) )
n = size(x,1)
lda=max(1,n)
allocate( ipiv(n,1) )
allocate( work(n,1) )

!$omp parallel
 !$ TID = OMP_GET_THREAD_NUM()
 !$ TID = omp_get_max_threads()
 !$ call omp_set_num_threads(TID)
!$omp end parallel
!$ print'(a,i,a)','OMP activated using = ',TID,' Threads'

x=0;val=0.
do i=1,size(x,1)
  do j=1,size(x,2)
    val=(i+j)/2.
    x(i,j)=val
    if(i==j) x(i,j)=x(i,j)+10-i*j**i
  end do
  !print*, x(i,:)
end do

! Store A in Ainv to prevent it from being overwritten by LAPACK
xinv=x

!print*,''
!do i=1,n
!  print '(100f10.5)', xinv(i,:)
!end do

call cpu_time(t1)

! DGETRF computes an LU factorization of a general M-by-N matrix A
! using partial pivoting with row interchanges.
!call DGETRF(n, n, xinv, n, ipiv, info)
call DPOTRF( UPLO, N, xinv, LDA, INFO )

if (info /= 0) then
 stop 'Matrix is numerically singular!'
end if

! DGETRI computes the inverse of a matrix using the LU factorization
! computed by DGETRF.
!call DGETRI(n, xinv, n, ipiv, work, n, info)
call DPOTRI( UPLO, N, xinv, LDA, INFO )
call cpu_time(t2)

if (info /= 0) then
 stop 'Matrix inversion failed!'
end if

print*,''
do i=1,n
  do j=i,n
    if(i<j) xinv(j,i)=xinv(i,j)
  end do
  print '(100f10.5)', xinv(i,:)
end do

print*,'Elapsed time ',t2-t1,' in seconds'

!print*,''
!do i=1,n
!  print*, xinv(i,:)
!end do

!n=size(x,1); lda=max(1,n)
!allocate( w(n) )
!call dsyev('N', 'U', n, x, lda, w, temp, -1, info)
!lwork = temp(1)
!allocate( work(lwork) )
!call dsyev('N', 'U', n, x, lda, w, work, lwork, info)
!deallocate(work)

!IF( INFO.GT.0 ) THEN
!  WRITE(*,*)'The algorithm failed to compute eigenvalues.'
!  STOP
!END IF

end program
