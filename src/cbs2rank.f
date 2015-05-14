c     CBS for bivariate (CN, BAF) vector. 
c     using rank transformation for the data 
c     ranks centered at 0 (i.e subtract (n+1)/2 from ranks)
c     using z statistic i.e. rank-sum/sqrt((j-i)*(n-j+i)*(n+1)/12)
c
c     sx is the bivariate vector - centered and scaled
c        will be used to store partial sums (computed locally)
c     ihet is the cumsum of het indicator
c
      subroutine t2maxo(n, sx, ihet, iseg, ostat)
      integer n, iseg(2)
c     the two columns of sx are cnlr, vaflor
      double precision sx(n,2), ihet(n), ostat

c     local variables we need
      integer nb, i, j, k, l, nb2, bi, bj, ilo, ihi, jlo, jhi
      double precision psx1, psx2, psmin1, psmax1, psdiff1, psmin2,
     1     psmax2, psdiff2, rn, zstat, hijinv, dnhet, alen1, alen2,
     1     xm1, xm2
c     
c     use local arrays for working within blocks
c     block absolute max for variables
      double precision, allocatable :: xmax1(:), xmax2(:)
c     block partial sum max and min for variable 1
      double precision, allocatable :: bpsmax1(:), bpsmin1(:)
c     block partial sum max and min for variable 2
      double precision, allocatable :: bpsmax2(:), bpsmin2(:)
c     location of the max and min
      integer, allocatable :: bb(:)

c     t statistic bound for block i,j stored in a vector
      double precision, allocatable :: btmax(:)
c     row, column and order vector for reordering bssbij
      integer, allocatable :: bloci(:), blocj(:), loc(:)

c     number of heterozygous snps
      dnhet = ihet(n)
c     calculate number of blocks (nb) and block boundaries (vector bb)
      rn = dfloat(n)
c     use blocks only if n is at least 100 
      if (n .ge. 100) then
         nb = nint(sqrt(rn))
      else
         nb = 1
      endif

c     the number of paiwise block comparison
      nb2 = nb*(nb+1)/2
c     allocate memory
      allocate(xmax1(nb), xmax2(nb))
      allocate(bpsmax1(nb), bpsmin1(nb),bpsmax2(nb), bpsmin2(nb))
      allocate(bb(nb))
      allocate(btmax(nb2))
      allocate(bloci(nb2), blocj(nb2), loc(nb2))

c     block boundaries
      do 10 i = 1, nb
         bb(i) = nint(rn*(dfloat(i)/dfloat(nb)))
 10   continue

c     find the max, min of partial sums within blocks
      ilo = 1
c     block counter
      l = 0
c     partial sum initialization
      psx1 = 0.0
      psx2 = 0.0
c     loop over blocks
      do 40 j = 1, nb
c     initialize local max of |x1| and |x2|
         xmax1(j) = abs(sx(ilo,1))
         xmax2(j) = abs(sx(ilo,2))
c     update partial sums and store in place
         psx1 = psx1 + sx(ilo,1)
         sx(ilo,1) = psx1
         psx2 = psx2 + sx(ilo,2)
         sx(ilo,2) = psx2
c     initialize block specific min and max at first observation
         psmin1 = sx(ilo,1)
         psmax1 = sx(ilo,1)
         psmin2 = sx(ilo,2)
         psmax2 = sx(ilo,2)
c     loop over observations in block
         do 20 i = ilo+1, bb(j)
c     update local max of |x1| and |x2|
            xmax1(j) = max(xmax1(j), abs(sx(i,1)))
            xmax2(j) = max(xmax2(j), abs(sx(i,2)))
c     update partial sums and store in place
            psx1 =  psx1 + sx(i,1)
            sx(i,1) = psx1
            psx2 =  psx2 + sx(i,2)
            sx(i,2) = psx2
c     update for first variable
            if (sx(i,1) .lt. psmin1) psmin1 = sx(i,1)
            if (sx(i,1) .gt. psmax1) psmax1 = sx(i,1)
c     update for second variable
            if (sx(i,2) .lt. psmin2) psmin2 = sx(i,2)
            if (sx(i,2) .gt. psmax2) psmax2 = sx(i,2)
 20      continue
c     store the block min, max and locations
         bpsmin1(j) = psmin1
         bpsmax1(j) = psmax1
         bpsmin2(j) = psmin2
         bpsmax2(j) = psmax2
c     t-statistic bound for this block matched with all previous blocks
         do 30 i = 1,j
c     increment block counter
            l = l+1
            loc(l) = l
            bloci(l) = i
            blocj(l) = j
c     max abs difference in partial sums (bound of numerator)
            psdiff1 = max(abs(bpsmax1(j)-bpsmin1(i)), 
     1           abs(bpsmax1(i)-bpsmin1(j)))
            psdiff2 = max(abs(bpsmax2(j)-bpsmin2(i)), 
     1           abs(bpsmax2(i)-bpsmin2(j)))
c     mimimum arc lengths (bound for denominator)
            if (i .eq. j) then
               alen1 = 1
               alen2 = min(1.0d0, dnhet-1.0d0)
            else if (i .eq. 1) then
               alen1 = dfloat(min(ilo-bb(i), n-bb(j)+1))
               alen2 = min(ihet(ilo)-ihet(bb(i)), dnhet-ihet(bb(j))+
     1              ihet(1))
            else
               alen1 = dfloat(min(ilo-bb(i), n-bb(j)+bb(i-1)+1))
               alen2 = min(ihet(ilo)-ihet(bb(i)), dnhet-ihet(bb(j))+
     1              ihet(bb(i-1)+1))
            endif
c     minimum number of observations needed to get to psdiff
            if (alen1 .eq. 1) then
               xm1 = max(xmax1(i), xmax1(j))
               if(xm1 .gt. 0) alen1 = dfloat(ceiling(psdiff1/xm1))
            endif
            if (alen2 .le. 1) then
               xm2 = max(xmax2(i), xmax2(j))
               if(xm2 .gt. 0) alen2 = dfloat(ceiling(psdiff2/xm2))
            endif
c     t statistic bounds
c     bound for variable 1
            btmax(l) = rn*psdiff1**2/(alen1*(rn-alen1))
c     bound for variable 2 (dnhet >=2 needed for data to be non zero)
            if (alen2 .gt. 0.5) btmax(l) = btmax(l) +
     1           dnhet*psdiff2**2/(alen2*(dnhet-alen2))
 30      continue
c     reset ilo to be the block boundary + 1
         ilo = bb(j) + 1         
 40   continue

c     Now sort the t-statistics by their magnitude
      call qsort4(btmax, loc, 1, nb2)

c     now go through the blocks in reverse order (largest down)
      ostat = 0.0d0
      l = nb2
      do 100 while ((btmax(l) .gt. ostat) .and. (l .gt. 0))
c     which blocks
         k = loc(l)
         bi = bloci(k)
         bj = blocj(k)
c     start and end indices corresponding to block
         if (bi .eq. 1) then
            ilo = 1
         else
            ilo = bb(bi-1) + 1
         endif
         ihi = bb(bi)
c     if bi and bj are the same indices in a triangle
         if (bi .eq. bj) then
c     loop over indices
            do 70 i = ilo, ihi-1
               do 60 j = i+1, ihi
c     copy number is available for all SNPs
                  zstat = rn*(sx(j,1) - sx(i,1))**2/
     1                 (dfloat(j-i)*dfloat(n-j+i))
c     product of number of het in (i,j] and outside it
                  hijinv = (ihet(j)-ihet(i))*(dnhet-(ihet(j)-ihet(i)))
c     if product is 0 (<0.001 to allow for machine precision) scale is 0
                  if (hijinv .lt. 0.001) then
                     hijinv = 0.0
                  else
                     hijinv = 1.0/hijinv
                  endif
c     baflor is available only for 
                  zstat = zstat + dnhet*hijinv*(sx(j,2) - sx(i,2))**2
                  if (zstat .gt. ostat) then
                     ostat = zstat
                     iseg(1) = i
                     iseg(2) = j
                  endif
 60            continue
 70         continue
         else
            if (bj .eq. 1) then
               jlo = 1
            else
               jlo = bb(bj-1) + 1
            endif
            jhi = bb(bj)
c     loop over indices
            do 90 i = ilo, ihi
               do 80 j = jlo, jhi
c     copy number is available for all SNPs
                  zstat = rn*(sx(j,1) - sx(i,1))**2/
     1                 (dfloat(j-i)*dfloat(n-j+i))
c     product of number of het in (i,j] and outside it
                  hijinv = (ihet(j)-ihet(i))*(dnhet-(ihet(j)-ihet(i)))
c     if product is 0 (<0.001 to allow for machine precision) scale is 0
                  if (hijinv .lt. 0.001) then
                     hijinv = 0.0
                  else
                     hijinv = 1.0/hijinv
                  endif
c     baflor is available only for 
                  zstat = zstat + dnhet*hijinv*(sx(j,2) - sx(i,2))**2
                  if (zstat .gt. ostat) then
                     ostat = zstat
                     iseg(1) = i
                     iseg(2) = j
                  endif
 80            continue
 90         continue
         endif
c     block l is finished step down
         l = l-1
 100  continue
      deallocate(xmax1, xmax2)
      deallocate(bpsmax1, bpsmin1, bpsmax2, bpsmin2, bb)
      deallocate(btmax, bloci, blocj, loc)

      return
      end
