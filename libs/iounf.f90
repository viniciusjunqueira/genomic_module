module iounf
! unformatted I/O operations with a large buffer
!Ignacy Misztal, 9/98-2/2005
!
implicit none


type unit_info
   integer::status,&	!status is 0 initially, 1 if opened
            bufsize,&	!buffer size	    
	    written,&	! elements written into the current buffer
	    filled,&	! elements filled in the current buffer	
	    read	! elements read from the current buffer    
   real, pointer::buf(:)	!buffer
end type  

! operation codes
integer,parameter:: iob_open=1,&
                    iob_write=2,&
		    iob_read=3,&
		    iob_rewind=4,&
		    iob_close=5,&
		    iob_delete=6


integer, parameter,private::last_unit=90000,&	!largest unit number
		         def_bufsize=1024	!default buffer size
		    	    
type (unit_info),private::xx(last_unit)	! contains information on all units



interface iobuf
   module procedure iobuf_r4, iobuf_int,iobuf_r8
end interface


   contains
   
   subroutine iobuf_r4(op,un,x,stat,bufsize,fname)
   ! works with real numbers only
   integer::op,un
   integer,optional::stat,bufsize
   real,optional::x(:)
   character (*),optional::fname
   integer::i,j,k,n,l,statio,start
   logical::isopen
   !
   
  
  if (un >= last_unit .or. un < 1) then
      print*,'iobuf: un= ',un, ' must be > 0 and < ',last_unit
  endif	   
  
  if (op /=iob_open .and. xx(un)%status == 0) then
     print*,'iobuf: unit ',un,' not opened'
  endif
  
  if ((op == iob_read .or. op == iob_write) .and. .not. present(x)) then
     print*,'iobuff: data variable missing'
  endif   
	      
   select case(op) 
   	case(iob_open)
              inquire(un,opened=isopen)
              if (isopen) close(un)	!close unit if preopened by system
	      if (present(fname)) then
	         open(un,file=fname,form='unformatted')
	        else
  	         open(un,form='unformatted')
	      endif 
	      xx(un)%status=1
	      xx(un)%filled=0
	      xx(un)%written=0
	      xx(un)%read=0
	      if (present(bufsize)) then
	         xx(un)%bufsize=bufsize
		else
		 xx(un)%bufsize=def_bufsize 
	      endif
              allocate(xx(un)%buf(xx(un)%bufsize))
	      
	case (iob_write)
	      if (xx(un)%status /=1) then
	         print*,'iobuf: unit not opened:',un
		 stop
	      endif	 
	      n=size(x)
	      start=0
	      do	!loop if x greater than buf
	         i=xx(un)%bufsize;   j=xx(un)%written	
	         k=i-j		!space left in buffer
		 l=min(n,k)
		 xx(un)%buf(j+1:j+l)=x(start+1:start+l)
		 xx(un)%written=j+l
		 if (j+l == i) then
		    write(un)i
		    !write(8,*)'b',i
		    write(un)xx(un)%buf
		    !write(8,*)'b',xx(un)%buf
		    !print*,'xx(un)%buf',xx(un)%buf
		    xx(un)%written=0
		  endif  
	          if (n <= k) then
		     exit
		    else
		     !x(1:n-l)=x(l+1:n)   
		     start=start+l
		     n=n-l
		   endif
	        enddo     
		
	case(iob_rewind)
		if (xx(un)%written > 0) then
		   write(un)xx(un)%written
		   write(un)xx(un)%buf(1:xx(un)%written)
		   
		   !write(8,*)'c',xx(un)%written
		   !write(8,*)'c',xx(un)%buf(1:xx(un)%written)
		   
		   !print*,'rew',xx(un)%buf(1:xx(un)%written)
		   xx(un)%written=0
		endif
		xx(un)%filled=0
		xx(un)%read=0
		rewind un
		
	case(iob_read)
		n=size(x)
		statio=0
		start=0
	        do	!loop if not enough data in buffer
		   if (xx(un)%filled-xx(un)%read == 0) then
		      read(un,iostat=statio)xx(un)%filled
		      if (statio /= 0) exit
		      read(un,iostat=statio)xx(un)%buf(1:xx(un)%filled)
		      !print*,'xx(un)%filled',xx(un)%filled
		      !print*,'xx(un)%buf(1:xx(un)%filled)',xx(un)%buf(1:xx(un)%filled)
		      if (statio /= 0) exit
		      xx(un)%read=0		      
		   endif   
		   l=min(xx(un)%filled-xx(un)%read,n-start)
		   !print*,'l=',l
 		   x(start+1:start+l)=xx(un)%buf(xx(un)%read+1:xx(un)%read+l)
		   !print*,'x(start+1:start+l)',x(start+1:start+l)
		   start=start+l
		   xx(un)%read=xx(un)%read+l
		   !print*,'xx(un)%filled,xx(un)%read',xx(un)%filled,xx(un)%read
		   !print*,'start,l',start,l
		   if (start == n) exit
		 enddo  
		   
	case(iob_close)
		close(un)
		deallocate(xx(un)%buf)
		xx(un)%status=0
		
	case(iob_delete)
		close(un,status='delete')
		deallocate(xx(un)%buf)
		xx(un)%status=0

	case default
		print*,'iobuf: parameter should be one of: ',&
		       'iob_open, iob_write, iob_rewind, iob_read'
		print*,'        iob_close or iob_delete'       
		stop
      end select		
      if (present(stat)) stat=statio		   
    end subroutine
    
    

   subroutine iobuf_int(op,un,x,stat)
   ! works with integers numbers only
   integer::op,un
   integer,optional::stat
   integer::x(:)
   real::x_r4(size(x))
   !    
   select case(op) 
   	case(iob_write)
	   call iobuf_r4(op,un,transfer(x,x_r4),stat)
	case (iob_read)
	   call iobuf_r4(op,un,x_r4,stat)
	   x=transfer(x_r4,x)
	case default
	   print*,"IOBUF with present integer vector limited to",&
	           " iob_write and iob_read operations"
   end select
   end subroutine iobuf_int
   
   
    subroutine iobuf_r8(op,un,x,stat)
   ! works with integers numbers only
   integer::op,un
   integer,parameter::r8=selected_real_kind(15)
   integer,optional::stat
   real (r8)::x(:)
   real::x_r4(2*size(x))
   !    
   select case(op) 
   	case(iob_write)
	   call iobuf_r4(op,un,transfer(x,x_r4),stat)
	case (iob_read)
	   call iobuf_r4(op,un,x_r4,stat)
	   x=transfer(x_r4,x)
	case default
	   print*,"IOBUF with present real8 vector limited to",&
	           " iob_write and iob_read operations"
   end select
   end subroutine iobuf_r8
  
   		   
end module    
			
		   
 
