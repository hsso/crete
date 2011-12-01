c     Use these routines to "translate" numrep routines into
c     slatec routines.

      FUNCTION ran1(idum)

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
c     2000, A&A, 362, 697.

c     mrh / 27sep01
c     This routines translates the utility of the ran1 routine from
c     Numerical Recipes using SLATEC routines RUNIF and RAND.
c
c     ran1(-1) intializes the routine: at the next call it will take the
c              argument idum as the random seed of the sequence.
c     ran1(seed) will return a random number and changes seed to 0. 
c              If the call is repeated with the same value for seed,
c              the random sequence is repeated.
c     seed/1000 is used as REAL random seed. 1/1000 seemed to work well,
c              given that it is based on the process ID in amc and sky.

      IMPLICIT NONE
      INTEGER idum,N,seed
      PARAMETER (N=32)
      DOUBLE PRECISION ran1
      DOUBLE PRECISION runif,rand,T(N+1),Z(N)
      EXTERNAL runif,rand
      SAVE seed,T

      if (idum.lt.0) then
c     Set seed=0; next call will give determine seed
        seed=0
      else
c     This call determines seed (previous call had idum.lt.0)
        if (seed.eq.0) seed=idum
        if (idum.eq.seed) then
c     idum=seed means (re)initialize sequence (at seed/1000)
c     Initialize rand; initialize runif by calling with N-1 table
          ran1=rand(dble(abs(idum))/1000.)
          ran1=runif(Z,N-1)
        endif
        ran1=runif(T,N)
        idum=0
      endif

      RETURN
      END

c     ------------------------------------------------------------

c     LUDCMP can be skipped entirely in SLATEC (SGEIR does both)

      SUBROUTINE ludcmp(a,n,np,indx,d)
      IMPLICIT NONE
      INTEGER n,np,indx(n)
      DOUBLE PRECISION d,a(np,np)
      continue
      return
      END

c     ------------------------------------------------------------

      SUBROUTINE lubksb(a,n,np,indx,b)
      IMPLICIT NONE
      INTEGER n,np,indx(np)
      DOUBLE PRECISION a(np,np),b(np)

      DOUBLE PRECISION ra(n-1,n-1),rb(n-1)
      INTEGER itask,iwork(n-1),ind,i,j
      DOUBLE PRECISION work((n-1)*(n-1+1))

      itask=1

c     Requires "smaller" matrix: nlev*nlev only, with last equation
c     replaced by sum_of_pops=1
c     This is a bit of a waste of computing time...

      do j=1,n-1
        do i=1,n-2
          ra(i,j)=a(i,j)
        enddo
        ra(n-1,j)=1.
      enddo
      do i=1,n-2
        rb(i)=0.
      enddo
      rb(n-1)=1.

      call sgeir(ra,n-1,n-1,rb,itask,ind,work,iwork)

      do i=1,n-1
        b(i)=rb(i)
      enddo

      return
      END

c     ------------------------------------------------------------

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

