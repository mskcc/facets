      subroutine scansnp(n, maploc, het, keep, nbhd)
      integer n
      double precision maploc(n), het(n), keep(n), nbhd
c     uniform random number from R (declared in rshared.c)
      double precision dunif
      external dunif
c     local variables; i=start of window; j=snps in nbhd; isel=selected snp
      integer i, j, isel
      double precision nsnp, nhet, usnp
c     get random number generator state
      call rndstart()

c     start with the first snp and set keep
      i=1
      keep(i) = 1
c     nsnp, nhet and isel are initialized
      nsnp = 1.0
      nhet = het(1)
      isel = i
c     loop through snps
      do 100 j = 2,n
c     if the next snp is within nbhd
         if (abs(maploc(j) - maploc(i)) .lt. nbhd) then
c     set nsnp and nhet
            nsnp = nsnp + 1.0
            nhet = nhet + het(j)
c     generate a random number to decide which one is kept
            usnp = dunif()
c     if nhet is > 0 choose only among the het snps
            if (nhet > 0) then
c     check if current snp is heterozygous
               if (het(j) .eq. 1) then
c     if random number is 1/nhet switch the selection
                  if (usnp .lt. 1.0/nhet) then
                     keep(isel) =  0
                     keep(j) = 1.0
                     isel = j
c     not selected so don't switch
                  else
                     keep(j) = 0.0
                  endif
c     current snp not a het but there is already at least one het
               else
                  keep(j) = 0.0
               endif
c     no het yet so choose from all snps
            else
c     if random number is 1/nsnp switch the selection
               if (usnp .lt. 1.0/nsnp) then
                  keep(isel) =  0
                  keep(j) = 1.0
                  isel = j
c     not selected so don't switch
               else
                  keep(j) = 0.0
               endif
            endif
c     this snp is at least nbhd bases away from previous snp
         else
c     reset the beginning of the window and select the snp
            i = j
            isel = i
            keep(i) = 1
            nsnp = 1.0
            nhet = het(i)
         endif
 100  continue
c     set random number generator to current state
      call rndend()

      return
      end
