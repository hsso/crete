      SUBROUTINE gauline(v,sigma,gauss)

c (c) Michiel Hogerheijde / Floris van der Tak 2000
c     michiel@strw.leidenuniv.nl, vdtak@sron.nl
c     http://www.strw.leidenuniv.nl/~michiel
c     http://www.sron.rug.nl/~vdtak
c
c     This file is part of the 'ratran' molecular excitation and
c     radiative transfer code. The one-dimensional version of this code
c     is publicly available; the two-dimensional version is available on
c     collaborative basis. Although the code has been thoroughly tested,
c     the authors do not claim it is free of errors or that it gives
c     correct results in all situations. Any publication making use of
c     this code should include a reference to Hogerheijde & van der Tak,
c     2000, A&A, 362, 697.

c     Calculates a gaussian line profile first time called and there
c     after uses a look-up table to return vaule of gaussian input x
c     =distance from center of gaussian profile expressed in sigma.

c     This version: gaussian normalized to have max=1 (not area=1)

c     input:
c        v -- v (or velocity) coordinate
c        sigma --- fwhm/(2*sqrt(2ln2))   --> gauss=exp(-(v/sgima)^2)
c
c     output:
c        gauss --- value of gaussian at x
c
c     maxgau   -- number of steps
c     maxsig   -- number of sigma's to follow profile out to
c     fac      -- maxgau-1/maxsig: divide maxgauu over required sigma's
c     gauvals  -- values of gaussian at step i
c     val      -- steplength*step_number
c     ival     -- corresponding place of v given sigma
c
c     Follow profile out to 4 sigma (exp(-4)**2 = 1e-7)

      IMPLICIT NONE
      DOUBLE PRECISION pi
      PARAMETER (pi=3.14159265d0)

      INTEGER maxgau,maxsig
      DOUBLE PRECISION fac
      PARAMETER (maxgau = 401)
      PARAMETER (maxsig = 4)
      PARAMETER (fac = (maxgau-1)/maxsig )
      INTEGER i,ival
      DOUBLE PRECISION v,sigma,gauss
    1 DOUBLE PRECISION gauvals(maxgau),val
      SAVE gauvals


c     Return gauss=0 if x is larger then the number of sigma in
c     calculated profile

      ival=nint(fac*(dabs(v)/sigma))+1
      if ((ival-1).ge.maxgau) then
        gauss=0.d0
        RETURN
      endif

c     If the gaussian profile is not yet calculated... do it

      if (gauvals(1).lt.1.d-6) then
        do i=1,maxgau
          val = (dble(i)-1.d0)/fac
          gauvals(i) = dexp(-(val**2.d0))
        enddo
      endif

c     Get requested value of gaussian from calculated array

      gauss = gauvals(ival)

      RETURN
      END

c     ------------------------------------------------------------

      FUNCTION length(str)
c     Returns the lengths of a string
      INTEGER length,maxl,i
      PARAMETER(maxl=200)
      CHARACTER*200 str

      do i=1,maxl
         if (str(i:i).eq.' ') then
            length=i-1
            RETURN
         endif
      enddo
      STOP 'Error: File name too long'
      END

c     ------------------------------------------------------------

      FUNCTION planck(iline,t)
c     Returns the value of the planck function at the frequency of
c     iline and for a temperature t.
      IMPLICIT NONE
      INCLUDE 'amccommon.inc'
      INTEGER iline
      DOUBLE PRECISION planck,t

      if (t.le.eps) then
        planck=0.d0
      else
c     Explicitely take Wien approximation for h*nu>>k*T:
        if (hplanck*freq(iline).gt.100.d0*kboltz*t) then
          planck=2.d0*hplanck*((freq(iline)/clight)**2.d0)*freq(iline)
     $      *dexp(-hplanck*freq(iline)/kboltz/t)
        else
          planck=2.d0*hplanck*((freq(iline)/clight)**2.d0)*freq(iline)
     $      /(dexp(hplanck*freq(iline)/kboltz/t)-1.d0)
        endif
      endif

      RETURN
      END

