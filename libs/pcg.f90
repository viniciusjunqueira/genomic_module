module pcg
! Solves the system of equations A*sol=rhs by preconditioned gradient algorithm. 
!  The default preconditioner is diagonal but can be blockdiagonal for
!  blksize>1. The block implementation uses the sparse_hashm structure and
!  fspak90, and is designed primarly for testing and not efficiency. More
!  efficient implementation would use an array of dense symmetric matrices. 
!  BiCGSTAB and CGSQ solvers for Unsymmetric A  1/2009 by I.Aguilar

   use denseop; use sparseop
   implicit none

    interface solve_pcg
       module procedure solve_densem_pcg,&
                        solve_sparse_hashm_pcg!,&
!                        solve_sparse_ija_pcg
    end interface

    interface solve_cgsq
       module procedure solve_densem_cgsq,&
                        solve_sparse_hashm_cgsq
    end interface


    interface solve_bicgstab
       module procedure solve_densem_bicgstab,&
                        solve_sparse_hashm_bicgstab
    end interface    
   contains
   
   
   
   subroutine solve_densem_pcg(A,rhs,sol,blksize,verbose)
   ! solve A sol = rhs by preconditioned conjugate gradient for densem A
   !   Preconditioner can be block-diagonal if blksize>1
   type (densem)::A
   real (rh)::rhs(:),sol(:)
   integer,optional::blksize
   integer,optional:: verbose
   integer :: verb
   real (r8),allocatable::m(:),r(:),p(:),z(:),w(:)
   real (r8)::alpha,beta,tau,oldtau,conv,val
   integer::i,j,k,block
   real :: t1,t2
   !variables for blockdiagonal if any
   type (sparse_hashm) :: m_hash
   type (sparse_ija) :: m_ija
   
   if (present(verbose)) then 
      verb=verbose
   else 
      verb=2
   endif

   if (present(blksize)) then
         block=blksize
      else
         block=1
    endif
    
    if (block == 1) then
        allocate(m(a%n))
      else
        call zerom(m_hash,a%n) !,int(a%n*block*(block+1)/2.0*1.3))
   endif       
   
   ! find the preconditioner 
   if (block == 1 ) then 
        !diagonal preconditioner
        m=0
        do i=1,a%n
!           if (a%x(i,i)/=0) m(i)=1/a%x(i,i)
            val=getm(i,i,a)
            if (val /= 0) m(i)=1/val
        enddo
     else   
        ! block preconditioner; not efficiently implemented but simple
        do i=1,a%n,block
           do j=0,block-1
              do k=0,block-1
                 call addm(getm(i+j,i+k,a),i+j,i+k,m_hash)
              enddo
           enddo
        enddo
        m_ija=m_hash
        call reset(m_hash)
        call fspak90('invert',m_ija)
        call fspak90('reset',m_ija)
     endif
   
   allocate (r(a%n),p(a%n),z(a%n),w(a%n))
   
   if (default_zerosol) then
      sol=0
      r=rhs
   else 
      r=rhs-mult(a,sol)
   endif
   
   call cpu_time(t1)  
   do k=1,default_rounds
      if (block==1) then
            z=m*r
         else
            z=mult(m_ija,r)
      endif
      tau=dot_product(z,r)
      if (k == 1) then
         beta=0; p=z
       else
         beta=tau/oldtau
         p=z+beta*p
      end if
      w=mult(a,p)
      alpha=tau/dot_product(p,w)
      sol=sol+alpha*p
      if (mod(k,100) /= 0) then
           r=r-alpha*w
         else
           r=rhs-mult(a,real(sol,r8))
      endif
      conv=dot_product(r,r)/dot_product(rhs,rhs)
      
      if (verb>=2) write(*,'(a8,i5,a8,g11.4)') 'round = ',k,'  convergence = ',conv
      if (conv < default_conv) exit
      oldtau=tau
   enddo
   call cpu_time(t2)   
   if (block == 1) then
         deallocate(m)
      else
         call reset(m_ija)
   endif
   deallocate (r,p,z,w)
   if (verb>=1) print '(i0,a,g10.3,a,g10.3)',k,' iterations,   convergence criterion=',conv,' sec per round=',(t2-t1)/k

  end subroutine
  
     
   
   subroutine solve_sparse_hashm_pcg(A,rhs,sol,blksize,verbose)
   ! solve A sol = rhs by preconditioned conjugate gradient for sparse_hashm A
   !   Preconditioner can be block-diagonal if blksize>1
   type (sparse_hashm)::A
   real (rh)::rhs(:),sol(:)
   integer,optional::blksize
   integer, optional:: verbose
   integer :: verb
   real (rh),allocatable::m(:),r(:),p(:),z(:),w(:)
   real (r8)::alpha,beta,tau,oldtau,conv,val
   integer::i,j,k,block
   
   !variables for blockdiagonal if any
   type (sparse_hashm) :: m_hash
   type (sparse_ija) :: m_ija


   if (A%storage_type=='f') then 
       print*,'SOLVE_SPARSE_HASHM_PCG: not work Non-Symetric full-stored matrix'
    endif   

   if (present(verbose)) then 
      verb=verbose
   else 
      verb=2
   endif

   if (present(blksize)) then
         block=blksize
      else
         block=1
    endif
    
    if (block == 1) then
        allocate(m(a%n))
      else
        call zerom(m_hash,a%n) !,int(a%n*block*(block+1)/2.0*1.3))
   endif       
   
   ! find the preconditioner 
   if (block == 1 ) then 
        !diagonal preconditioner
        m=0
        do i=1,a%n
!           if (a%x(i,i)/=0) m(i)=1/a%x(i,i)
            val=getm(i,i,a)
            if (val /= 0) m(i)=1/val
        enddo
     else   
        ! block preconditioner; not efficiently implemented but simple
        do i=1,a%n,block
           do j=0,block-1
              do k=0,block-1
                 call addm(getm(i+j,i+k,a),i+j,i+k,m_hash)
              enddo
           enddo
        enddo
        m_ija=m_hash
        call reset(m_hash)
        call fspak90('invert',m_ija)
        call fspak90('reset',m_ija)
     endif
   
   allocate (r(a%n),p(a%n),z(a%n),w(a%n))
   if (default_zerosol) then
      sol=0
      r=rhs
   else 
      r=rhs-mult(a,sol)
   endif
  
   do k=1,default_rounds
      if (block==1) then
            z=m*r
         else
            z=mult(m_ija,real(r,r8))
      endif
      tau=dot_product(z,r)
      if (k == 1) then
         beta=0; p=z
       else
         beta=tau/oldtau
         p=z+beta*p
      end if
      w=mult(a,p)
      alpha=tau/dot_product(p,w)
      sol=sol+alpha*p
      if (mod(k,100) /= 0) then
           r=r-alpha*w
         else
           r=rhs-mult(a,sol)
      endif
      conv=dot_product(r,r)/dot_product(rhs,rhs)
      if (verb>=2) write(*,'(a,i5,a,g11.4)') 'round = ',k,'  convergence = ',conv
      if (conv < default_conv) exit
      oldtau=tau
   enddo
   
   if (block == 1) then
         deallocate(m)
      else
         call reset(m_ija)
   endif
   deallocate (r,p,z,w)
   if (verb>=1) write(*,'(i5,a,g11.4)') k,' iterations,   convergence criterion=',conv

  end subroutine
  
   subroutine solve_densem_cgsq(A,rhs,sol,blksize)
   ! solve A sol = rhs by Conjugate Gradient Squared for densem A
   !   Preconditioner can be block-diagonal if blksize>1
   ! memory share phat qhat -> pqh
   !              vhat uhat -> vuh
   ! TODO reuse more variables
   ! 01/29/09 IA
   type (densem)::A
   real (rh)::rhs(:),sol(:)
   integer,optional::blksize
   real (r8),allocatable:: m(:),r(:),p(:),u(:),q(:),pqh(:),vuh(:)
   real (r8)::alpha,beta,tau,oldtau,conv,val
   integer::i,j,k,block
   real :: t1,t2
   !variables for blockdiagonal if any
   type (sparse_hashm) :: m_hash
   type (sparse_ija) :: m_ija
   print*,'Solving by SOLVE_DENSEM_CGS'
   if (present(blksize)) then
         block=blksize
      else
         block=1
    endif

    if (block == 1) then
        allocate(m(a%n))
      else
        call zerom(m_hash,a%n) !,int(a%n*block*(block+1)/2.0*1.3))
   endif
   ! find the preconditioner
   if (block == 1 ) then
        !diagonal preconditioner
        m=0
        do i=1,a%n    
            val=getm(i,i,a)
            if (val /= 0) m(i)=1/val
        enddo
     else
        ! block preconditioner; not efficiently implemented but simple
        do i=1,a%n,block
           do j=0,block-1
              do k=0,block-1
                 call addm(getm(i+j,i+k,a),i+j,i+k,m_hash)
              enddo
           enddo
        enddo
        m_ija=m_hash
        call reset(m_hash)
        call fspak90('invert',m_ija)
        call fspak90('reset',m_ija)
     endif

   allocate (r(a%n),p(a%n),u(a%n),q(a%n),pqh(a%n),vuh(a%n))
   sol=0
   r=rhs
   call cpu_time(t1)
   do k=1,default_rounds
      tau=dot_product(rhs,r)
      if (tau == 0) then 
         print *, 'CGS Fail '
         sol = 0
         exit
      endif
      if (k==1) then
          u =  r
          p = u
      else
          beta = tau/oldtau
          u = r + beta * q
          p = u + beta * ( q + beta*p )
      endif

      if (block==1) then 
         pqh = m*p
      else
         pqh = mult(m_ija,p)
      endif

      vuh =  mult(A,pqh)
      alpha = tau/dot_product(rhs,vuh)
      q = u - alpha * vuh
      
      if (block==1) then
         vuh = m * (u+q)
      else
         vuh = mult(m_ija,(u+q))
      endif     
 
      sol = sol + alpha*vuh

      pqh = mult(A,vuh)

      r = r - alpha*pqh
     
      conv=dot_product(r,r)/dot_product(rhs,rhs)
      print*,'round ',k,'   convergence=',conv
      if (conv < default_conv) exit
      oldtau=tau
   enddo   
   call cpu_time(t2)
   if (block == 1) then
         deallocate(m)
      else
         call reset(m_ija)
   endif
   deallocate (r,p,u,q,pqh,vuh)
   print*,k,' iterations,   convergence criterion=',conv,' sec/round=',(t2-t1)/k
   end subroutine



   subroutine solve_sparse_hashm_cgsq(A,rhs,sol,blksize)
   ! solve A sol = rhs by Conjugate Gradient Squared for sparse_hashm A
   !   Preconditioner can be block-diagonal if blksize>1
   ! memory share phat qhat -> pqh
   !              vhat uhat -> vuh
   ! TODO reuse more variables
   ! 01/29/09 IA
   type (sparse_hashm)::A
   real (rh)::rhs(:),sol(:)
   integer,optional::blksize
   real (r8),allocatable:: m(:),r(:),p(:),u(:),q(:),pqh(:),vuh(:)
   real (r8)::alpha,beta,tau,oldtau,conv,val
   integer::i,j,k,block
   real :: t1,t2
   !variables for blockdiagonal if any
   type (sparse_hashm) :: m_hash
   type (sparse_ija) :: m_ija

   print*,'Solving by SOLVE_SPARSE_HASHM_CGS'
   if (present(blksize)) then
         block=blksize
      else
         block=1
    endif

    if (block == 1) then
        allocate(m(a%n))
      else
        call zerom(m_hash,a%n) !,int(a%n*block*(block+1)/2.0*1.3))
   endif
   ! find the preconditioner
   if (block == 1 ) then
        !diagonal preconditioner
        m=0
        do i=1,a%n    
            val=getm(i,i,a)
            if (val /= 0) m(i)=1/val
        enddo
     else
        ! block preconditioner; not efficiently implemented but simple
        do i=1,a%n,block
           do j=0,block-1
              do k=0,block-1
                 call addm(getm(i+j,i+k,a),i+j,i+k,m_hash)
              enddo
           enddo
        enddo
        m_ija=m_hash
        call reset(m_hash)
        call fspak90('invert',m_ija)
        call fspak90('reset',m_ija)
     endif

   allocate (r(a%n),p(a%n),u(a%n),q(a%n),pqh(a%n),vuh(a%n))
   sol=0
   r=rhs
   call cpu_time(t1)
   do k=1,default_rounds
      tau=dot_product(rhs,r)
      if (tau == 0) then 
         print *, 'CGS Fail '
         sol = 0
         exit
      endif
      if (k==1) then
          u =  r
          p = u
      else
          beta = tau/oldtau
          u = r + beta * q
          p = u + beta * ( q + beta*p )
      endif

      if (block==1) then 
         pqh = m*p
      else
         pqh = mult(m_ija,p)
      endif

      vuh =  mult(A,pqh)
      alpha = tau/dot_product(rhs,vuh)
      q = u - alpha * vuh
      
      if (block==1) then
         vuh = m * (u+q)
      else
         vuh = mult(m_ija,(u+q))
      endif     
 
      sol = sol + alpha*vuh

      pqh = mult(A,vuh)

      r = r - alpha*pqh
     
      conv=dot_product(r,r)/dot_product(rhs,rhs)
      print*,'round ',k,'   convergence=',conv
      if (conv < default_conv) exit
      oldtau=tau
   enddo   
   call cpu_time(t2)
   if (block == 1) then
         deallocate(m)
      else
         call reset(m_ija)
   endif
   deallocate (r,p,u,q,pqh,vuh)
   print*,k,' iterations,   convergence criterion=',conv,' sec/round=',(t2-t1)/k
   end subroutine


   subroutine solve_densem_bicgstab(A,rhs,sol,blksize)
   ! solve A sol = rhs by BiConjugate Gradient Stabilized for densem A
   !   Preconditioner can be block-diagonal if blksize>1
   ! 
   ! TODO reuse  variables
   ! 01/29/09 IAG
   type (densem)::A
   real (rh)::rhs(:),sol(:)
   integer,optional::blksize
   real (r8),allocatable::m(:),r(:),p(:),v(:),ph(:),sh(:),s(:),t(:)
   real (r8)::alpha,beta,tau,oldtau,conv,val,w
   integer::i,j,k,block
   real :: t1,t2

   !variables for blockdiagonal if any
   type (sparse_hashm) :: m_hash
   type (sparse_ija) :: m_ija

   print*,'Solving by SOLVE_DENSEM_BiCGSTAB'
   if (present(blksize)) then
         block=blksize
      else
         block=1
    endif

    if (block == 1) then
        allocate(m(a%n))
      else
        call zerom(m_hash,a%n) !,int(a%n*block*(block+1)/2.0*1.3))
   endif
   ! find the preconditioner
   if (block == 1 ) then
        !diagonal preconditioner
        m=0
        do i=1,a%n
            val=getm(i,i,a)
            if (val /= 0) m(i)=1/val
        enddo
   else
      ! block preconditioner; not efficiently implemented but simple
      do i=1,a%n,block
         do j=0,block-1
            do k=0,block-1
               call addm(getm(i+j,i+k,a),i+j,i+k,m_hash)
            enddo
         enddo
      enddo
      m_ija=m_hash
      call reset(m_hash)
      call fspak90('invert',m_ija)
      call fspak90('reset',m_ija)
   endif

   allocate (r(a%n),p(a%n),v(a%n),ph(a%n),sh(a%n),s(a%n),t(a%n))

   sol=0
   r=rhs
   call cpu_time(t1)
   do k=1,default_rounds
      tau=dot_product(rhs,r)
      if (tau == 0) then
         print *, 'BiCGSAB Fail '
         sol = 0
         exit
      endif
      if (k==1) then
          p =  r
      else
          beta = (tau/oldtau) * (alpha/w)
          p = r + beta * ( p - w*v )
      endif

      if (block==1) then
         ph = m*p
      else
         ph = mult(m_ija,p)
      endif

      v =  mult(A,ph)
      alpha = tau/dot_product(rhs,v)
      s = r - alpha*v
      conv = dot_product(s,s)/dot_product(rhs,rhs)
     
      if (conv < default_conv ) then
         print*,'round ',k,'   convergence by norms=',conv
         sol = sol + alpha*ph
         exit
      endif 

      if (block==1) then
         sh = m * s
      else
         sh = mult(m_ija,s)
      endif      

      t = mult(A,sh)

      w = dot_product(t,s)/dot_product(t,t)

      if (w == 0) then 
         print *, 'BiCGSAB Fail '
         sol = 0
         exit
      endif

      sol = sol + (alpha * ph) + (w * sh)
      
      if (mod(k,50) /= 0) then
           r = s -  w*t
      else
           r = rhs - mult(A,sol)
      endif
      

      conv=dot_product(r,r)/dot_product(rhs,rhs)
      print*,'round ',k,'   convergence=',conv
      if (conv < default_conv) exit
      oldtau=tau
   enddo
   call cpu_time(t2)
   if (block == 1) then
         deallocate(m)
      else
         call reset(m_ija)
   endif
   deallocate (r,p,v,ph,sh,s,t)
   print*,k,' iterations,   convergence criterion=',conv,' sec/round=',(t2-t1)/k
      
   end subroutine

   subroutine solve_sparse_hashm_bicgstab(A,rhs,sol,blksize)
   ! solve A sol = rhs by BiConjugate Gradient Stabilized for sparse_hashm A
   !   Preconditioner can be block-diagonal if blksize>1
   !
   ! TODO reuse  variables
   ! 01/29/09 IAG
   type (sparse_hashm)::A
   real (rh)::rhs(:),sol(:)
   integer,optional::blksize
   real (r8),allocatable::m(:),r(:),p(:),v(:),ph(:),sh(:),s(:),t(:)
   real (r8)::alpha,beta,tau,oldtau,conv,val,w
   integer::i,j,k,block
   real :: t1,t2

   !variables for blockdiagonal if any
   type (sparse_hashm) :: m_hash
   type (sparse_ija) :: m_ija

   print*,'Solving by SOLVE_SPARSE_HASHM_BiCGSTAB'
   if (present(blksize)) then
         block=blksize
      else
         block=1
    endif

    if (block == 1) then
        allocate(m(a%n))
      else
        call zerom(m_hash,a%n) !,int(a%n*block*(block+1)/2.0*1.3))
   endif
   ! find the preconditioner
   if (block == 1 ) then
        !diagonal preconditioner
        m=0
        do i=1,a%n
            val=getm(i,i,a)
            if (val /= 0) m(i)=1/val
        enddo
   else
      ! block preconditioner; not efficiently implemented but simple
      do i=1,a%n,block
         do j=0,block-1
            do k=0,block-1
               call addm(getm(i+j,i+k,a),i+j,i+k,m_hash)
            enddo
         enddo
      enddo
      m_ija=m_hash
      call reset(m_hash)
      call fspak90('invert',m_ija)
      call fspak90('reset',m_ija)
   endif

   allocate (r(a%n),p(a%n),v(a%n),ph(a%n),sh(a%n),s(a%n),t(a%n))

   sol=0
   r=rhs
   call cpu_time(t1)
   do k=1,default_rounds
      tau=dot_product(rhs,r)
      if (tau == 0) then
         print *, 'BiCGSAB Fail '
         sol = 0
         exit
      endif
      if (k==1) then
          p =  r
      else
          beta = (tau/oldtau) * (alpha/w)
          p = r + beta * ( p - w*v )
      endif

      if (block==1) then
         ph = m*p
      else
         ph = mult(m_ija,p)
      endif

      v =  mult(A,ph)
      alpha = tau/dot_product(rhs,v)
      s = r - alpha*v
      conv = dot_product(s,s)/dot_product(rhs,rhs)
     
      if (conv < default_conv ) then
         print*,'round ',k,'   convergence by norms=',conv
         sol = sol + alpha*ph
         exit
      endif 

      if (block==1) then
         sh = m * s
      else
         sh = mult(m_ija,s)
      endif      

      t = mult(A,sh)

      w = dot_product(t,s)/dot_product(t,t)

      if (w == 0) then 
         print *, 'BiCGSAB Fail '
         sol = 0
         exit
      endif

      sol = sol + (alpha * ph) + (w * sh)

     
      if (mod(k,50) /= 0) then
           r = s -  w*t
      else
           r = rhs - mult(A,sol)
      endif

      conv=dot_product(r,r)/dot_product(rhs,rhs)
      print*,'round ',k,'   convergence=',conv
      if (conv < default_conv) exit
      oldtau=tau
   enddo
   call cpu_time(t2)
   if (block == 1) then
         deallocate(m)
      else
         call reset(m_ija)
   endif
   deallocate (r,p,v,ph,sh,s,t)
   print*,k,' iterations,   convergence criterion=',conv,' sec/round=',(t2-t1)/k
   
   end subroutine
  end module
 
