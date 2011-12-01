c     Use these routines to "translate" numrep routines into
c     slatec routines.

c (c) Michiel Hogerheijde / Floris van der Tak 2000
c     michiel@strw.leidenuniv.nl, vdtak@sron.rug.nl
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
c     2000, A&A 362, 697.

      SUBROUTINE polint(xa,ya,n,x,y,dy)
      IMPLICIT NONE
      INTEGER n,nder,IERR
      DOUBLE PRECISION dy,x,y,xa(n),ya(n),c(n),WORK(2*n)

      nder=1

      call polints(n,xa,ya,c)
      call polyvl(nder,x,y,dy,n,xa,c,WORK,IERR)

      return
      END

c     ------------------------------------------------------------

c     This routines "borrowed" from Numerical Recipes...
      SUBROUTINE locate(xx,n,x,j)
      IMPLICIT NONE
      INTEGER j,n
      DOUBLE PRECISION x,xx(n)
      INTEGER jl,jm,ju
      jl=0
      ju=n+1
10    if(ju-jl.gt.1)then
        jm=(ju+jl)/2
        if((xx(n).gt.xx(1)).eqv.(x.gt.xx(jm)))then
          jl=jm
        else
          ju=jm
        endif
      goto 10
      endif
      j=jl
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..

