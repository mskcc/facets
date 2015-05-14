c     Calculate Mann-Whitney statistic for x and y already sorted
      subroutine mwstat(x, n, y, m, ustat)
      integer n, m
      double precision x(n), y(m), ustat

      double precision currenty, dm, dn, numlt, numeq
      integer i, j

      dm = dfloat(m)
      dn = dfloat(n)
c     initialize Mann-Whitney U-statistic
      ustat = 0.0d0

c     set current y to be first observation
      currenty = y(1)
c     number less than currenty is 0
      numlt = 0.0
c     count the number equal to currenty
      j = 0
      numeq = 1.0
      do 10 while ((j .lt. m) .and. (y(j+1) .eq. currenty))
         j = j+1
         numeq = numeq + 1.0
 10   continue

c     loop through the xs
      do 50 i = 1,n
c     if x is less than or equal to currenty
         if (x(i) .le. currenty) then
            ustat = ustat + numlt
         else
c     increment currenty until it equals or exceeds x(i) or end of y
            do 20 while ((j .lt. m) .and. (currenty .lt. x(i)))
               j = j+1
               currenty = y(j)
 20         continue
c     when the loop stops either currenty >= x(i) or j = m (or both)
c     first j < m; so currenty >= x(i)
            if (j .lt. m) then
               numlt = dfloat(j-1)
c     check if there are other y tied at currenty
               numeq = 1.0
               do 30 while ((j .lt. m) .and. (y(j+1) .eq. currenty))
                  j = j+1
                  numeq = numeq + 1.0
 30            continue
            else
c     now j = m; so check if y(m)=currenty < x(i) or equal
               if (currenty .lt. x(i)) then
                  numlt = dm
c     change currenty to be larger than all xs
                  currenty = x(n) + 1.0
               else
                  numlt = dm-1.0
                  numeq = 1.0
               endif
            endif
c     now increment the ustat for x(i)
            ustat = ustat + numlt
         endif
c     finally add the bit for ties
         if (x(i) .eq. currenty) ustat = ustat + 0.5*numeq
c         call dblepr("numlt",5,numlt,1)
 50   continue
      ustat = (ustat - dm*dn/2.0)**2/(dm*dn*(dm+dn+1)/12.0)

      return
      end

c     merge two sorted vectors x and y to get sorted result
      subroutine mergexy(x, n, y, m, xy, mn)
      integer m, n, mn
      double precision x(n), y(m), xy(mn)

      integer i, j, l
      double precision currenty

      l = 0
      j = 1
      do 20 i = 1, n
c     first put all the ys that are <= x(i) into xy vector
         do 10 while ((j .le. m) .and. (y(j) .le. x(i)))
c     y(j) <= x(i) and j <= m so add the y(j) to xy
            l = l+1
            xy(l) = y(j)
c     increment j
            j = j+1
 10      continue
c     now either j > m or y(j) > x(i), so add x(i) to xy
         l = l+1
         xy(l) = x(i)
 20   continue
c     if j <= m add any left over ys to xy
      do 30 while (j .le. m)
         l = l+1
         xy(l) = y(j)
c     increment j
         j = j+1
 30   continue

      return
      end
