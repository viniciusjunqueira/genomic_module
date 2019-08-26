program genomic_module
!use omp_lib
use mkl_service
use kinds; use model;use sparsem; use sparseop; use pcg; use textop; use genomic
!use omp_lib

implicit none


type(sparse_hashm) :: xx !, GimA22
type(sparse_ija) :: xx_ija!, aija !, GimA22_ija
type(sparse_ija), allocatable :: ainv_ija(:)
integer,allocatable:: ped(:,:,:)        ! pedigree data (could be multiple pedigrees)
integer,allocatable:: address(:,:),&    ! start and address of each effect
                      df_random(:)

logical :: is_genomic=.false.,i_genotyped=.false.,j_genotyped=.false.

! mkl
integer :: num_thread,nthrg
integer :: NTHREADS, TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM,&
           omp_get_max_threads


print*,''
!$omp parallel
 !$ TID = OMP_GET_THREAD_NUM()
 !$ TID = omp_get_max_threads()
 !$ call omp_set_num_threads(TID)
!$omp end parallel
!$ print'(a35,i10,a)','OMP activated using = ',TID,' Threads'
print*,''


allocate( lenped(1) )
neff=1
lenped(neff)=12755
allocate ( ainv(neff), ainv_ija(neff), df_random(neff), ped(neff,maxval(lenped),3) )


ped_file = 'pedigree.ped'


! pedigree file, allow for more than one pedigree
open(io_off,file=ped_file)
do i=1, neff
    do j=1, lenped(i)
        read(io_off,*,iostat=io) ped(i,j,:)
        ! print*, ped(i,j,:)
    enddo
enddo

! Random effects' contributions
do i=1,neff
    call add_g_add(g_A,i)
enddo

is_genomic=.true.



call defaults_and_checks_GS





contains


subroutine add_g_add(type,eff)
! generates contributions for additive sire or animal effects
integer :: type,eff,i,j,t1,t2,k,l,m,n,io,p(4),ped_len,iped
real (rh) :: w(3),res(4),mendel
!

select case (type)
    case (g_A)
    w=(/1., -.5, -.5/)
    res=(/2., 4/3., 1., 0./)
    ped_len=3
    case (g_A_UPG)
    w=(/1., -.5, -.5/)
    res=(/2., 4/3., 1., 0./)
    !ped_len=4
    ped_len=3                ! VJ - pedigree is allocated as a three-dimensional matrix = ped(eff,length,3)
    case (g_A_UPG_INB)
    w=(/1., -.5, -.5/)
    ped_len=4
    case (g_As)
    w=(/1., -.5, -.25/)
    res=(/16/11., 4/3., 16/15., 1./)
    ped_len=3
end select

if (round == start) call zerom(ainv(eff),maxval(lenped))

do iped=1, lenped(eff)
    p(1:ped_len)=ped(eff,iped,1:ped_len)
    select case (type)
    case (g_A)
        i=1; if(p(2)==0) i=i+1; if (p(3) == 0) i=i+1
        mendel=res(i)
    case (g_A_UPG)
        mendel=res(p(4))
    case (g_A_UPG_INB)
        mendel=p(4)/1000.
    case (g_As)
        i=1; if(p(2)==0) i=i+2; if (p(3) == 0) i=i+1
        mendel=res(i)
    end select

    ! generate A inverse
    if (round ==start) then
    do k=1,3
        do l=1,3
        if (p(k) /=0 .and. p(l) /=0) then
            call addm(w(k)*w(l)*mendel,p(k),p(l),ainv(eff))
        endif
        enddo
    enddo
    endif
enddo

if (round == start .and. fthr <= 1) print*,' read ',lenped(eff),' additive pedigrees'
if (lenped(eff) == 0) then
    print*,'Additive pedigree file for effect ',eff,' empty'
    stop
endif

if (round ==start) df_random(eff)=lenped(eff)
end subroutine



function address1(e,l,t)
! returns address for level l of effect e and trait t
integer :: e,l,t, address1,i
logical,save::first=.true.
integer,allocatable,save::offset(:)
!
if (first) then
    first=.false.
    allocate(offset(neff))
    do i=1,neff
    offset(i)=sum(nlev(1:i-1))*ntrait
    enddo
endif
!
address1= offset(e)+(l-1)*ntrait+t
!address1= sum(nlev(1:e-1))*ntrait+(l-1)*ntrait+t
end function


end program
