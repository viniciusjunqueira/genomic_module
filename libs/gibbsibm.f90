! Collection of data manipulation subroutines helpful to implement 
!   Gibbs sampling
! IM  04/15/98-07/03/06

module gibbs
use sparsem; use denseop
implicit none



! vector for links between SPARSE_HASHM and SPARSE_IJA formats
   integer, allocatable,save::linkhi(:)

  contains

!=====================================================================
!Module subroutines
!=====================================================================

function mmatmul(x,y) result(z)
! z=X*y
! This is a matrix x vector multiplication to circumvent a bug in the
! xlf90 (IBM) compiler
!
real (r8)::x(:,:),y(:),z(size(y))
integer::i,j,k,n
!
n=size(y)
z=0
do i=1,n
   do j=1,n
      z(i)=z(i)+x(i,j)*y(j)
   enddo
enddo
end function



!----------------------------------------------------------------
!  Convert hash to ija fromats and create link for fast  
!  subsequent update
!----------------------------------------------------------------

  recursive subroutine link_hash_ija(xx,xx_ija)
  type (sparse_hashm)::xx
  type (sparse_ija)  ::xx_ija
  integer::i,j,k,l
  
  if (.not. allocated(linkhi)) then	
!     first entry, construct links
      xx_ija=xx
      xx=xx_ija	 !back conversion to restore xx
      allocate(linkhi(xx%filled))
      do i=1,xx_ija%n
         do j=xx_ija%ia(i),xx_ija%ia(i+1)-1
	    k=xx_ija%ja(j)
 	    l=hashvr((/real(i,rh),real(k,rh)/),2,xx%x,xx%nel,3,0,xx%filled)
            if (l <1) then 
                print*,'link_hash_ija: hash conversion problem'
		stop
              else
	        linkhi(j)=l
	    endif   
         enddo
      enddo
      
    elseif( size(linkhi) == xx%filled) then
!     subsequent entry, update values
      do i=1,size(linkhi)
         xx_ija%a(i)=xx%x(3,linkhi(i))
      enddo	   
      
    else
!     dimensions do not agree,  reinitialize
      print*,' link_hash_ija: reinitialized'
      deallocate(linkhi)
      call link_hash_ija(xx,xx_ija)
  endif
  
  end subroutine             
    

!----------------------------------------------------------------
! solve blocks of equations and update them    
!----------------------------------------------------------------

  subroutine solve_iterm_block(xx,xy,sol,first,last,diag,op)
  ! if op='solve' : diag=xx(first:last,first:last)
  !                  sol(first:last)=inverse(diag)*xy(first:last)
  !
  ! if op='update': adjust xy for current equations to account for xx 
  !                 triangual storage
  !xx is of type sparse_ija, vectors are of type (rh).
  !                 
  type (sparse_ija)::xx
  real (r8)::xy(:),sol(:)
  integer::first,last,i,j,k,l
  real(r8)::diag(:,:)
  character (*)::op
  !
  if (op == 'solve') then
  
    ! extract block and solve by finite method  
     diag=0
     do i=first,last
       do j=xx%ia(i),xx%ia(i+1)-1
         if (xx%ja(j) >=first .and. xx%ja(j) <=last) then
	     k=i-first+1;   l=xx%ja(j)-first+1
	     diag(k,l)=xx%a(j)
	     diag(l,k)= diag(k,l)
	   else
	     xy(i)=xy(i)-xx%a(j)*sol(xx%ja(j))
	 endif
       enddo
     enddo 	 
     sol(first:last)=fsolve_s(diag,xy(first:last))
    
  elseif (op == 'update') then
  
    ! adjust right hand sides because of triangular storage
    do i=first,last
       do j=xx%ia(i),xx%ia(i+1)-1
          if (xx%ja(j) > last) xy(xx%ja(j)) = xy(xx%ja(j))-xx%a(j)*sol(i)
       enddo
    end do
  endif
  
  end subroutine


!------------------------------------------------------------------
!  Convert matrix in hash to ija format but fully stored 
!------------------------------------------------------------------
 subroutine convert_hash_ija_full(y,x)
 ! converts symmetric hash matrix x into ija matrix y, fully stored
 
 type (sparse_ija), intent(out) :: y
 type (sparse_hashm), intent(in) :: x
 
   call convert_hash_ija_general(y,x,conv_upper_full)
 end subroutine
   

!-----------------------------------------------------------------
! adjust RHS for one block of linear equations; left hand side is 
! block stored as xx [x] v, where [x] is kronecker product.
!-----------------------------------------------------------------
   
  subroutine adjust_rhs_block(n,xx,v,arhs,sol,diagonal,offset1,offset2)
! arhs(block n)=arhs(block n)-
!          sum[xx(n+offset1,i+offset2)*v*sol(block i)]
! diagonal=diagonal+v* xx(n+offset1,n+offset2)
!
! block size is determined by the size of square matrix v
!
! xx must be of type sparse_ija and fully storred

  implicit none
  integer::n,offset1,offset2
  type(sparse_ija)::xx
  real (r8)::v(:,:),arhs(:),sol(:),diagonal(:,:)
  integer::block,i,j,eq1,eq2
  real (r8)::a

  if (offset1 >=n) return	!matrix is zero
  if (offset1+xx%n > size(sol) .or. offset2+xx%n > size(sol)) then
      print*,'ADJUST_RHS_BLOCK: matrix with offsets too large'; stop
  endif    

  block=size(v,dim=1)
  if (size(diagonal,dim=1) /= block) then
     print*,'ADJUST_RHS_BLOCK: v and diag have different dimensions'; stop
  endif

  eq1=(n-1)*block+1; eq2=n*block
  do i=xx%ia(n-offset1), xx%ia(n-offset1+1)-1
     j=xx%ja(i)+offset2
     a=xx%a(i)
     !arhs(eq1:eq2)=arhs(eq1:eq2)-a*matmul(v,sol((j-1)*block+1:j*block))
     ! In IBM, matmul here causes bad answers; mmatmul substituted
     arhs(eq1:eq2)=arhs(eq1:eq2)-a*mmatmul(v,sol((j-1)*block+1:j*block))
     if (n == j) then
        diagonal=diagonal+v*a
     endif
  enddo

end subroutine


!-----------------------------------------------------------------
! As above but no diagonal (nd)
!-----------------------------------------------------------------
   
  subroutine adjust_rhs_block_nd(n,xx,v,arhs,sol,offset1,offset2)
! arhs(block n)=arhs(block n)-
!          sum[xx(n+offset1,i+offset2)*v*sol(block i)]
!
! block size is determined by the size of square matrix v
!
! xx must be of type sparse_ija and fully storred

  implicit none
  integer::n,offset1,offset2
  type(sparse_ija)::xx
  real (r8)::v(:,:),arhs(:),sol(:)
  integer::block,i,j,eq1,eq2
  real (r8)::a

  if (offset1 >=n) return	!matrix is zero
  if (offset1+xx%n > size(sol) .or. offset2+xx%n > size(sol)) then
      print*,'ADJUST_RHS_BLOCK: matrix with offsets too large'; stop
  endif    

  block=size(v,dim=1)

  eq1=(n-1)*block+1; eq2=n*block
  do i=xx%ia(n-offset1), xx%ia(n-offset1+1)-1
     j=xx%ja(i)+offset2
     a=xx%a(i)
     arhs(eq1:eq2)=arhs(eq1:eq2)-a*matmul(v,sol((j-1)*block+1:j*block))
  enddo

end subroutine



end module
