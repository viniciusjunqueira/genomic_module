C Legacy subroutines in Fortran 77
C
C Data manipulation
C Residual for inverse of parental dominance matrix


      subroutine nums(a,x,xc,m,n)
c separates line a into items delimited by blanks. character elements 
c kept in a character vector xc, decoded numeric values in vector x.
c version corrected for negative numbers 4/14/93
      integer m,x(m),decode,d,i,j,e,b,lena,n
      character a*(*),xc(m)*(*),buf*1,sp*1,sp12*20
      data sp/' '/,sp12/'                    '/
      d=0
      n=0
      lena=len(a)
      do 5 i=1,m
5       xc(m)=sp12
      j=1
10    buf=a(j:j)
      if (buf.eq.sp) then
        if (d.ne.0) then
          e=j-1
          n=n+1
           if (n>m) then
              print*, "Error in nums splitting string starting: "
              print*,a(1:80)
              print*,'into m:',m,'items'
              stop
           endif
           xc(n)=a(b:e)
          if (e-b.lt.10)  x(n)=decode(xc(n),e-b+1)
          d=0
        endif
        else
        if (d.eq.0) then
          b=j
          d=1
        endif
      endif
      j=j+1
      if (j.le.lena) goto 10
      return
      end
 
      integer function decode(buf,b)
      integer i,j,dd,b0,b,sign
      character*1 d0,buf*(*),minus
      data d0/'0'/,minus/'-'/
      decode=0
      b0=ichar(d0)
      if (buf(1:1).eq.minus) then
           sign=-1
           j=2
         else
           sign=1
           j=1
      endif
      do 10 i=j,b
      dd=ichar(buf(i:i))-b0
      if (dd.ge.0 .and. dd.le.9) then
        decode=10*decode+dd
      else
        goto 20
      endif
10    continue
20    decode=decode*sign
      return
      end
 
      subroutine pack(a,m,x,n)
c packs an x(n) vector of integers into a*80 character variable. Variables
c there are separated by spaces
      character*80 a
      integer n,x(n),y(10),xi,iz,i,iy,j,m
      iz=ichar('0')
      m=0
      do 10 i=1,n
        xi=x(i)
        if (xi.lt.0) then
          xi=-xi
          m=m+1
          a(m:m)='-'
        endif
        iy=0
20      iy=iy+1
        y(iy)=mod(xi,10)
        xi=xi/10
        if (xi.ne.0) goto 20
        do 30 j=iy,1,-1
          m=m+1
30        a(m:m)=char(y(j)+iz)
        if (i.ne.n) then
          m=m+1
          a(m:m)=' '
        endif
10      continue                    
      return
      end


      subroutine creatres(r)
c compute residuals for different classes of missing s-d combinations
      use kinds
      real (rh):: r(0:255),varf(8,8),b(8),b1(8)
      integer i,j,k
c define  b and varf
      data b/4*.5,4*-.25/
      data varf/ 1.0,  .0,  .25, .25, .5,  .0,  .5,  .0,
     +            .0, 1.0,  .25, .25, .0,  .5,  .0,  .5,
     +            .25, .25,1.0,  .0,  .5,  .5,  .0,  .0,
     +            .25, .25, .0, 1.0,  .0,  .0,  .5,  .5,
     +            .5,  .0,  .5,  .0, 1.0,  .0,  .0,  .0,
     +            .0,  .5,  .5,  .0,  .0, 1.0,  .0,  .0,
     +            .5,  .0,  .0,  .5,  .0,  .0, 1.0,  .0,
     +            .0,  .5,  .0,  .5,  .0,  .0,  .0, 1.0/
c
       do i=255,0,-1
          do j=1,8
             if (btest(i,j-1)) then
                 b1(j)=b(j)
                else
                 b1(j)=0
             endif
          enddo
          r(i)=0
          do j=1,8
             do k=1,8
                r(i)=r(i)+b1(j)*varf(j,k)*b1(k)
             enddo
          enddo
          if (r(i).ge..9999) then
              r(i)=0
            else
              r(i)=1/(1-r(i))
          endif
      enddo
c print
c     print '('' r''/,100(15f5.2/))',(r(i),i=1,255)
      end
