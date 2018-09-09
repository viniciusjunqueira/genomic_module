program ev_omp
!use mkl_service
use model; use sparsem; use sparseop; use denseop; use textop; use iounf
implicit none

type (sparse_hashm) :: ai

integer :: io_a22i=1,io,io_out=2
integer :: i,j,row,col,x,la,nthr,nthrg
real(rh) :: val,val1
real(rh), allocatable :: work(:), dd(:,:)
character (30) :: fname,ioname,fname_out


!$omp parallel
!$ print'(a)','OMP working!'
!$omp end parallel


print*,'Input file name:'
read*,fname
fname_out= trim(fname) // '2'

print*,'Output file name:'
read*,ioname
ioname=trim(ioname)

open(io_a22i,file=fname)
!open(io_a22i,file='test')
open(io_out,file=ioname)

print*,''
write(*,'(3a)') 'Computing eigenvalue of ',fname,' file'
write(*,'(a)') 'Counting number of elemets'
x=0;la=0;val=0.;val1=0.
do
  read(io_a22i,*,iostat=io) row,col,val
  if(io/=0) exit
  la=max(la,max(row,col))
  x=x+1
end do
rewind io_a22i

! Initialization
call zerom(ai,la)

write(*,'(a,i10,a)') 'Storaging',x,' elements'
do
  read(io_a22i,*,iostat=io) row,col,val
  if(io/=0) exit
  call addm(val,row,col,ai)
  !if(row/=col) call addm(val,col,row,ai)
end do
print*,''


allocate( work(ai%n), dd(ai%n,ai%n) )
work=0;dd=0

write(*,'(a,i10)')'Matrix dimension (square matrix):', size(dd,1)
write(*,'(a)')'Converting hash matrix into dense format'
do i=1,ai%n
  do j=i,ai%n
      val=getm(i,j,ai)
      if(val==0.) cycle    ! Next if is missing data
      dd(i,j)=val; dd(j,i)=val
  end do
end do

! Print if has less than 9 rows
if(size(dd,1)<=9) then
  print*,'Matrix'
  call printm(ai)
end if

write(*,'(a)')'Computing eigenvalues...'
call DSYEV_F90(dd,work)

print*,''
write(*,'(a)') 'Eigenvalues'
write(*,'(a,f15.6)') 'Minimum = ', minval(work)
write(*,'(a,f15.6)') 'Average = ', sum(work)/size(work)
write(*,'(a,f15.6)') 'Maximum = ', maxval(work)
write(*,'(a,f15.6)') 'SD      = ', sqrt(var(work))

! Save eigenvalues
write(*,'(a)') 'Saving eigenvalues...'
do i=1,ai%n
  write(io_out,'(f15.6)') work(i)
end do
print*,''



contains

! Function that compute variance of a vector
function var(x) result(c)
implicit none
real(rh) :: x(:),c,mu
integer :: i

! Mean
mu=sum(x)/size(x)
! Sum of Squares
do i=1,size(x)
  c=c+(x(i)-mu)**2
end do 
c=c/size(x)

end function



end program







