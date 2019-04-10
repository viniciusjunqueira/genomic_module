program genomic_module
use mkl_service
use kinds; use model;use sparsem; use sparseop; use pcg; use textop; use genomic
!use omp_lib

implicit none




type(sparse_hashm) :: xx,GimA22
type(sparse_ija) :: xx_ija,GimA22_ija
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
lenped(neff)=4
allocate ( ainv(neff), ainv_ija(neff), df_random(neff), ped(neff,maxval(lenped),3) )


ped_file = 'pedigree.ped'

! pedigree file, allow for more than one pedigree
open(io_off,file=ped_file)
do i=1, neff
    do j=1, lenped(i)
        read(io_off,*,iostat=io) ped(i,j,:)
    enddo
    print*, ped(i,:,:)
enddo

! Random effects' contributions
do i=1,neff
    call add_g_add(g_A,i)
enddo


call defaults_and_checks_GS
is_genomic=.true.


! Store A in Ainv to prevent it from being overwritten by LAPACK
!do i=1,lenped(neff)
!    do j=i,lenped(neff)
        !i_genotyped = xref(i)
        !j_genotyped = xref(j)
        !print*, i, j, getm(i,j,ainv(neff))
        !val = Gi(i,j) - getm(xref(i),xref(j),ainv(neff))
        !print*,xref(i),xref(j),val
        !call setm(val,xref(i),xref(j),GimA22)
        !print '(3f7.3)',Gi(i,j) , getm(xref(i),xref(j),ainv(neff)) , getm(xref(i),xref(j),GimA22)
!    enddo
!enddo

!print*,''

!do i=1,lenped(neff)
!    do j=i,lenped(neff)
!        i_genotyped = any(xref == i)
!        j_genotyped = any(xref == j)
!        if(i_genotyped .and. j_genotyped) then
!            val = getm(i,j,ainv(neff)) + getm(i,j,GimA22)
!        else 
!            val = getm(i,j,ainv(neff))
!        endif
!        !print '(2i,3f7.3)', i,j,getm(i,j,ainv(neff)),getm(i,j,GimA22),val 
!    enddo
!enddo









contains

subroutine create_GimA22(m)
implicit none

real(rh) :: m(:,:)
integer :: irank,i

!do i=1,size(m,1)
!    print*, m(i,:)
!enddo
!print*,''

call ginv1(m,size(m,1),size(m,1),real(1d-12,rh),irank)

!print*,'Inverting genomic matrix'
!do i=1,size(m,1)
!   print*, m(i,:)
!enddo

end subroutine



! Algorithm from Ignacio JABS (2011)
!subroutine compute_A22(ped)
!real(rh) :: f(nped),w(nped),v(nped),tmp(nsamples,nsamples)
!integer :: i,irank,ped(:,:)
!call zerom(A22i,nsamples)
!f=0;A22=0
!n=nped
!do i=1,nsamples
!    v=0;w=0
!    v(perm(xref(i,1)))=1
!    call A_times_v(w,v,ped,f)
!    A22(:,i)=w(perm(xref(:,1)))    
!enddo
!call printmat(A22)
!tmp=A22
!call ginv1(tmp,size(tmp,1),size(tmp,1),real(1d-12,rh),irank)
!do i=1,nsamples
!    do j=i,nsamples
!        val=tmp(i,j)
!        if(val/=0) then
!            call setm(val,i,j,A22i); call setm(val,j,i,A22i)
!        endif
!    enddo
!enddo
!call printm(A22i)
!end subroutine

! Algorithm from Ignacio JABS (2011)
subroutine A_times_v(w,v,ped,f)
integer :: i,s,d,n,ped(:,:)
real(rh) :: w(:),v(:),f(:),tmp,di
real(rh), allocatable :: q(:)

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

subroutine createG_dense_symm(m,Z,G,w,k,n,snp)
implicit none
integer  :: i,j,n,l,snp
real (rh) :: val, W(:,:), m(:,:), G(:,:), Z(:,:), k(:),diag,ofdiag
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
            G(i,j)=G(i,j)+Z(i,l)*Z(j,l)
        end do
    end do
end do
do i=1,n
    do j=1,n
        G(i,j)= G(i,j)/sqrt(k(i)*k(j))
        if (i<=j) write(io_save,fmt="(2i,f11.3)") i,j,G(i,j)
    end do
    if (n < 11) print '(10f10.2)', G(i,:)
end do
! ====================================================================
call cpu_time(t2)

! scaling
!print*, 'Scaling G matrix'
!do i=1,n
!    do j=1,n
!        G(i,j) = G(i,j) * .95
!        if (i==j) G(i,j) = G(i,j) + .05
!    enddo
!enddo
print '(a35,f10.2)', 'Time to construct G = ',t2-t1
close(io_save)
end subroutine







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
