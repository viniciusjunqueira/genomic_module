program genomic
use mkl_service
use model; use sparsem
!use gibbs;
!use prob; use pcg; use textop; use ranlib
!use omp_lib
implicit none


!character (len=100) :: a, gen_file = "../examples/ex1/deltag_50K_imp.gen_clean", g_xref="../examples/ex1/deltag_50K_imp.gen_clean_XrefID",&
!                          ped_file = "../examples/ex1/renadd01.ped"
character (len=100) :: a, gen_file = "genotypes.gen", g_xref="genotypes.gen_XrefID",&
                          ped_file = "pedigree.ped"


type(sparse_hashm) :: xx,GimA22
type(sparse_ija) :: xx_ija,GimA22_ija
type(sparse_hashm), allocatable :: ainv(:)
type(sparse_ija), allocatable :: ainv_ija(:)
integer,allocatable:: ped(:,:,:)        ! pedigree data (could be multiple pedigrees)
integer,allocatable:: address(:,:),&    ! start and address of each effect
                      df_random(:)


character(5) :: version= "1.000"
logical :: sameweights
logical :: is_genomic=.false.,i_genotyped=.false.,j_genotyped=.false.

real     ::  wG       =  0.95,&
             wA22     =  0.05,&
             wGimA22i =  1.d0,&     
             wGi      =  1.d0,&    
             wA22i    = -1.d0,&
             wOmega   =  1.d0
real (rh), allocatable :: gen(:,:), W(:,:), Z(:,:), k(:), Gi(:,:)
real (rh) :: temp(1),t1,t2,val

integer :: whichH=0, io_g=10, io_i, io, i, j, l, begSNP, lastSNP, nid
integer :: num_snp=80000, len_char=50
integer :: unit_xref = 123
integer, allocatable :: xref(:)

! mkl
integer :: num_thread,nthrg
integer :: NTHREADS, TID, OMP_GET_NUM_THREADS, OMP_GET_THREAD_NUM,&
           omp_get_max_threads


!integer :: lwork,info,lda
!character :: uplo='U'
!real (rh), allocatable :: ww(:), rand(:)
!real (rh), allocatable :: work(:,:)!, work1(:)
!integer, allocatable :: ipiv(:,:)



print*,''
print*,'*******************************'
print*,'Genomic module - Version ',version
print*,'*******************************'

is_genomic=.true.

!$omp parallel
 !$ TID = OMP_GET_THREAD_NUM()
 !$ TID = omp_get_max_threads()
 !$ call omp_set_num_threads(TID)
!$omp end parallel
!$ print'(a35,i10,a)','OMP activated using = ',TID,' Threads'


open(unit_xref, file=g_xref)
open(io_g, file=gen_file)


call read_first_and_last_snp(nid,begSNP,lastsnp)

print '(a35,i10)','Number of genotyped individuals = ',nid
print '(a35,i10)','Position of first SNP = ',BegSNP
print '(a35,i10)','Number of genotypes = ',LastSNP
print*,''

allocate( lenped(1) )
neff=1
lenped(neff)=4 !12668
allocate ( ainv(neff), ainv_ija(neff), df_random(neff), ped(neff,maxval(lenped),3) )
allocate( Gi(nid,nid), gen(nid,lastsnp), W(3,lastsnp), Z(nid,lastsnp), k(lastsnp) )
!allocate( gen(nid,lastsnp), W(3,lastsnp), Z(nid,lastsnp), k(lastsnp) )
!call zerom(Gi,nid)
allocate( xref(nid) )

call read_xref_file(nid)
!print*, xref


!allocate( rand(nid) )
!lda=max(1,nid); lwork=max(1,3*nid-1)
!allocate( ipiv(nid,1) )
!allocate( work(nid,1), work1(max(1,lwork)) )
!work1=0

! pedigree file, allow for more than one pedigree
open(io_off,file=ped_file)
do i=1, neff
    do j=1, lenped(i)
        read(io_off,*,iostat=io) ped(i,j,:)
    enddo
enddo

! Random effects' contributions
do i=1,neff
    call add_g_add(g_A,i)
enddo

call read_markers(gen, nid, BegSNP, LastSNP)
call createk(gen,nid,k)
call createW(gen,W,nid,LastSNP)
call createG_dense_symm(gen,Z,Gi,W,k,nid,LastSNP)

! Store A in Ainv to prevent it from being overwritten by LAPACK
!Gi=Genomic
!call create_Gi(n=nid, snp=LastSNP, X=Gi)
call create_GimA22(gi)


call zerom(GimA22,lenped(neff))
do i=1,nid
    do j=i,nid
        !i_genotyped = xref(i)
        !j_genotyped = xref(j)
        val = Gi(i,j) - getm(xref(i),xref(j),ainv(neff))
        !print*,xref(i),xref(j),val
        call setm(val,xref(i),xref(j),GimA22)
        print '(3f7.3)',Gi(i,j) , getm(xref(i),xref(j),ainv(neff)) , getm(xref(i),xref(j),GimA22)
    enddo
enddo

print*,''

do i=1,lenped(neff)
    do j=i,lenped(neff)
        i_genotyped = any(xref == i)
        j_genotyped = any(xref == j)
        if(i_genotyped .and. j_genotyped) then
            val = getm(i,j,ainv(neff)) + getm(i,j,GimA22)
        else 
            val = getm(i,j,ainv(neff))
        endif
        !print '(2i,3f7.3)', i,j,getm(i,j,ainv(neff)),getm(i,j,GimA22),val 
    enddo
enddo



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



! Read cross reference file
subroutine read_xref_file(n)
integer :: i,io,n,tmp(2)
xref=0
do i=1,n
    read(unit_xref,*,iostat=io) xref(i)
enddo
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



!subroutine create_Gi(n,snp,X)
!implicit none

!real (rh), intent(inout) :: X(:,:)
!integer :: n,snp,info
!real (rh) ::  work(n)

!work=0

!print*,'X row and column dim = ',size(X,1),size(X,2)
!print*,'work array dim = ',size(work)

!call cpu_time(t1)
!write(*,'(a)') 'Computing eigenvalues...'
!call DSYEV_F90(X, work)
!call cpu_time(t2)
!print '(a35,f10.2)', 'Time to compute eigenvalues  = ',t2-t1

!print*,''
!write(*,'(a)') 'Eigenvalues'
!write(*,'(a,f15.6)') 'Minimum = ', minval(work)
!write(*,'(a,f15.6)') 'Average = ', sum(work)/size(work)
!write(*,'(a,f15.6)') 'Maximum = ', maxval(work)
!write(*,'(a,f15.6)') 'SD      = ', sqrt(var(work))

! Save eigenvalues
!open(100,file="eigen")
!write(*,'(a)') 'Saving eigenvalues...'
!do i=1,n
!  write(100,'(f15.6)') work(i)
!end do
!print*,''

!call cpu_time(t1)
! DGETRF computes an LU factorization of a general M-by-N matrix A using partial pivoting with row interchanges.
!print '(a)', 'Decomposing genomic matrix...'
!call DGETRF(nid, nid, Gi, nid, ipiv, info)
!call DPOTRF( UPLO, nid, Gi, LDA, INFO )
!call cpu_time(t2)
!print '(a35,f10.2)', 'Time to factorize matrix  = ',t2-t1

!if (info /= 0) then
!    print*,' INFO variable =',info
!    stop ' Matrix is numerically singular!'
!end if
!print*,'LU G(1,1) = ',G(1,1),' Gi(1,1) = ',Gi(1,1)
!print*,'LU G(1,2) = ',G(1,2),' Gi(1,2) = ',Gi(1,2)

!call cpu_time(t1)
! DGETRI computes the inverse of a matrix using the LU factorization computed by DGETRF.
!print '(a)', 'Inverting genomic matrix...'
!call DGETRI(nid, Gi, nid, ipiv, work, nid, info)
!call DPOTRI( UPLO, nid, Gi, LDA, INFO )
!call cpu_time(t2)
!print '(a35,f10.2)', 'Time to compute inverse  = ',t2-t1

!if (info /= 0) then
! stop 'Matrix inversion failed!'
!end if
!print*,'Inv G(1,1) = ',G(1,1),' Gi(1,1) = ',Gi(1,1)
!print*,'Inv G(1,2) = ',G(1,2),' Gi(1,2) = ',Gi(1,2)

!open(40, file="Gi")
!do i=1,nid
!    do j=1,nid
!        if(i<=j) write(40,*) i,j,Gi(i,j)
!    enddo
!end do
!end subroutine



real function freq(values,n) result(x)
integer :: i,j,n
real (rh) :: f,values(:),np,npq,nmiss=0
np = count(values==0.)
npq = count(values==1.)
nmiss = count(values==5.)
x = (np + npq*.5)/(n-nmiss)
end function


! Subroutine to create a vector of 2pq
subroutine createk(m,n,k)
implicit none
real (rh) :: m(:,:), k(:), p, q, pq, val
integer :: i,j,n,io_save=20

!open(io_save,file="allele_freq")
print '(a)', 'Creating scaling (k) vector'
call cpu_time(t1)
do i=1,n
    p=0;val=0
    do j=1,size(m,2)
        p=freq(m(:,j),size(m,1))
        q=1-p
        !write(io_save,fmt="(5f15.3)"), p,q,p*q,2*p*q
        val=val+2*p*q
    enddo
    k(i) = val
enddo
call cpu_time(t2)
print '(a35,f10.2)', 'Time to compute k  = ',t2-t1
close(io_save)
end subroutine


subroutine createW(m,W,n,snp)
implicit none
real (rh) :: m(:,:), W(:,:), val
integer :: n,snp,i

do i=1,snp
    val = freq(m(:,i),n)
    W(1,i)=-2*val; W(2,i)=1-2*val; W(3,i)=2-2*val
enddo
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
print*, 'Scaling G matrix'
do i=1,n
    do j=1,n
        G(i,j) = G(i,j) * .95
        if (i==j) G(i,j) = G(i,j) + .05
    enddo
enddo

! matrix statistics
!diag=0; ofdiag=0
!do i=1,n
!    do j=i,n
!        if(i==j) then
!            diag=diag+G(i,j)
!        else
!            ofdiag=ofdiag+G(i,j)
!        endif
!    enddo
!enddo

!print '(a,f5.2)','Diagonal mean     =',diag/n
!print '(a,f5.2)','Off Diagonal mean =',ofdiag/((n**2-n)/2)

print '(a35,f10.2)', 'Time to construct G = ',t2-t1
close(io_save)
end subroutine


! Subroutine that reads genotype file as character and converts it into reale
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


integer function extract_snp(x,curPos) result(updatePOS)
character(len_char) :: x
integer :: curPos
do
    if( x(curPos:curPos) /= " ") exit
    curPos = curPos+1
enddo
updatePOS = curPos
end function



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
