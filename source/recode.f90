program recode
implicit none

character (len=100) :: ped_file != "unsort.ped"
!integer, parameter :: n=6
integer :: n, io_i, io, i, j, nid, io_off=10, top, bot, id, sire, dam, miss=0, Beg, curr, nval, pos
integer, allocatable :: ped(:,:), gen(:),tgen(:),perm(:),iperm(:),igen(:),indices(:)!=0
logical,allocatable :: nindex(:) !=.false.


print '(a)', 'Pedigree file name: '
read*, ped_file

!
! pedigree file, allow for more than one pedigree
!
print '(a)',"Reading pedigree..."
open(io_off,file=ped_file)
n=0;
do
    read(io_off,*,iostat=io) ped(j,:)
    if( io /= 0) exit
    n=n+1
enddo
print '(a,i10)','Number of individuals = ',n
rewind(io_off)

!
! Allocating arrays
!
allocate( ped(n,3), gen(n),tgen(n),perm(n),iperm(n),igen(n),indices(n), nindex(n) )
gen=1;tgen=1;perm=1;iperm=1;indices=0
nindex=.false.

!
! Loading pedigree data
!
do j=1,n
    read(io_off,*) ped(j,:)    
    !print*,ped(j,:)
enddo

! Step that computes 'generic' indices for each animal based on progeny. This is useful
! to sort individuals
!
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
        print*,"circular pedigree"
        stop
    endif
enddo
tgen=gen(ped(:,1))

!
! Organize indices to be in needed order
!
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
enddo; !print*, indices
ped=ped(indices,:)
!do i=1,n
!    print*, ped(indices(i),:)
!enddo

!
! Recode pedigree based on pre-defined order
!
do i=1,n
  id=ped(i,1)
  iperm(i)=id ! Original identification
  perm(id)=i  ! Recoded identification
enddo



end program