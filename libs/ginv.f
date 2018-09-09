      subroutine ginv1(a,n,m,tol,irank)
c returns generalized inverse of matrix x of size n x n declared
c as m x m. tol is working zero and irank returns the rank of
c the matrix. w is a work vector of size m,
c by rohan fernando, slightly structured by i. misztal 05/05/87
      use kinds
      real (rh)::a(m,m),w(m),re,sum,tol
      irank=n
      do 10 i=1,n
         do 20 j=1,i-1
              re=a(i,j)
              do 20 ii=i,n
20                 a(ii,i)=a(ii,i)-re*a(ii,j)
         if (a(i,i).lt.tol) then
              a(i,i)=0.0
              do 45 ii=i+1,n
45                 a(ii,i)=0.
           irank=irank-1
           else
              a(i,i)=sqrt(a(i,i))
              do 40 ii=i+1,n
40                a(ii,i)=a(ii,i)/a(i,i)
         endif
10    continue
 
      do 100 i=1,n
         if (a(i,i).eq.0.) then
              do 150 ii=i+1,n
150                a(ii,i)=0.
           else
              a(i,i)=1.0/ a(i,i)
              do 200 ii=i+1,n
200               w(ii)=0.0
              do 300 ii=i+1,n
                  iim1=ii-1
                  re=a(iim1,i)
                  do 400 iii=ii,n
400                   w(iii)=w(iii)-a(iii,iim1)*re
                  if (a(ii,ii).eq.0.) then
                      a(ii,i)=0.
                    else
                      a(ii,i)=w(ii)/a(ii,ii)
                  endif
300           continue
          endif
100     continue
 
      do 110 j=1,n
         do 110 i=j,n
              sum=0
              do 130 ii=i,n
130                sum=sum+a(ii,j)*a(ii,i)
110           a(i,j)=sum
      do 600 i=1,n
          do 600 j=i,n
600           a(i,j)=a(j,i)
      return
      end

