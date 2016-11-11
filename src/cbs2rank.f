c     CBS for bivariate (CN, BAF) vector. 
c     using rank transformation for the data 
c     ranks centered at 0 (i.e subtract (n+1)/2 from ranks)
c     using z statistic i.e. rank-sum/sqrt((j-i)*(n-j+i)*(n+1)/12)
c
c     sx is the bivariate vector - centered and scaled
c        will be used to store partial sums (computed locally)
c     ihet is the cumsum of het indicator
c     nhet is the number of hets in the segment
c     rnij stores  sqrt(n/((j-i)*(n-(j-i))))
c     rhij stores  nhet/((ihet(j)-ihet(i))*(nhet-(ihet(j)-ihet(i))))
c     delij stores  delta*sqrt{((j-i)*(n-(j-i)))/n} - min effect size
c     
      subroutine t2maxo(n, sx, ihet, iseg, ostat, nhet, rnij, rhij,
     1     delij)
      integer n, ihet(n), iseg(2), nhet
c     the two columns of sx are cnlr, vaflor
      double precision sx(n,2), ostat, rnij(n), rhij(nhet), delij(n)

c     local variables we need
      integer nb, i, j, k, l, nb2, bi, bj, ilo, ihi, jlo, jhi, ij, hij,
     1     alen1, alen2
      double precision psx1, psx2, psmin1, psmax1, psdiff1, psmin2,
     1     psmax2, psdiff2, rn, zstat, hijinv, dnhet, xm1, xm2
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

c     when indices are binned into blocks index i(j) is in block bi(bj)
c     ilo-ihi give the range of index i and jlo-jhi for index j

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
      jlo = 1
c     block counter
      l = 0
c     partial sum initialization
      psx1 = 0.0
      psx2 = 0.0
c     loop over block bj for index j
      do 40 bj = 1, nb
c     jlo has been set already; set jhi
         jhi = bb(bj)
c     initialize local max of |x1| and |x2|
         xmax1(bj) = abs(sx(jlo,1))
         xmax2(bj) = abs(sx(jlo,2))
c     update partial sums and store in place
         psx1 = psx1 + sx(jlo,1)
         sx(jlo,1) = psx1
         psx2 = psx2 + sx(jlo,2)
         sx(jlo,2) = psx2
c     initialize block specific min and max at first observation
         psmin1 = sx(jlo,1)
         psmax1 = sx(jlo,1)
         psmin2 = sx(jlo,2)
         psmax2 = sx(jlo,2)
c     loop over rest of the observations in block bj
         do 20 j = jlo+1, jhi
c     update local max of |x1| and |x2|
            xmax1(bj) = max(xmax1(bj), abs(sx(j,1)))
            xmax2(bj) = max(xmax2(bj), abs(sx(j,2)))
c     update partial sums and store in place
            psx1 =  psx1 + sx(j,1)
            sx(j,1) = psx1
            psx2 =  psx2 + sx(j,2)
            sx(j,2) = psx2
c     update for first variable
            if (sx(j,1) .lt. psmin1) psmin1 = sx(j,1)
            if (sx(j,1) .gt. psmax1) psmax1 = sx(j,1)
c     update for second variable
            if (sx(j,2) .lt. psmin2) psmin2 = sx(j,2)
            if (sx(j,2) .gt. psmax2) psmax2 = sx(j,2)
 20      continue
c     store the block min, max and locations
         bpsmin1(bj) = psmin1
         bpsmax1(bj) = psmax1
         bpsmin2(bj) = psmin2
         bpsmax2(bj) = psmax2
c     t-statistic bound for this block matched with all previous blocks
c     loop over block bi for index i
         do 30 bi = 1,bj
c     set the range of i index (ilo-ihi) for block bi
            if (bi .eq. 1) then
               ilo = 1
            else
               ilo = bb(bi-1)+1
            endif
            ihi = bb(bi)
c     increment block counter
            l = l+1
            loc(l) = l
            bloci(l) = bi
            blocj(l) = bj
c     max abs difference in partial sums (bound of numerator)
            psdiff1 = max(abs(bpsmax1(bj)-bpsmin1(bi)), 
     1           abs(bpsmax1(bi)-bpsmin1(bj)))
            psdiff2 = max(abs(bpsmax2(bj)-bpsmin2(bi)), 
     1           abs(bpsmax2(bi)-bpsmin2(bj)))
c     mimimum arc lengths (bound for denominator)
            if (bi .eq. bj) then
c     if blocks bi and bj are the same min arc length is 1
               alen1 = 1
               alen2 = 1
            else
c     o/w min of {end of i to start if j} & {end of j to start of i} 
               alen1 = min(jlo-ihi, n-(jhi-ilo))
               alen2 = min(ihet(jlo)-ihet(ihi),
     1              nhet-(ihet(jhi)-ihet(ilo)))
            endif
c     if alen2 if zero set it to be 1
            if (alen2 .le. 0) alen2 = 1
c     t statistic bounds
c     bound for variable 1
            zstat = rnij(alen1)*abs(psdiff1) - delij(alen1)
            if (zstat .lt. 0.0d0) zstat = 0.0d0
c     bound for variable 2 (if nhet =1 rhij(1) = 0 so just adds zero)
            btmax(l) = zstat**2 + rhij(alen2)*psdiff2**2
 30      continue
c     reset jlo to be the block boundary + 1
         jlo = bb(bj) + 1
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
                  ij = j-i
c     copy number is available for all SNPs
                  zstat = rnij(ij)*abs(sx(j,1) - sx(i,1)) - delij(ij)
                  if (zstat .lt. 0.0d0) zstat = 0.0d0
c     product of number of het in (i,j] and outside it
                  hij = (ihet(j)-ihet(i)) 
                  if (hij .eq. 0) hij = nhet
c     baflor is available only for 
                  zstat = zstat**2 + rhij(hij)*(sx(j,2) - sx(i,2))**2
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
                  ij = j-i
c     copy number is available for all SNPs
                  zstat = rnij(ij)*abs(sx(j,1) - sx(i,1)) - delij(ij)
                  if (zstat .lt. 0.0d0) zstat = 0.0d0
c     product of number of het in (i,j] and outside it
                  hij = (ihet(j)-ihet(i)) 
                  if (hij .eq. 0) hij = nhet
c     baflor is available only for 
                  zstat = zstat**2 + rhij(hij)*(sx(j,2) - sx(i,2))**2
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
