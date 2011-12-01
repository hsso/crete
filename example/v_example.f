c     Velocity field of Shu model. Need to know about a and t...

      SUBROUTINE velo(ix,x,v)
      INTEGER ix
      DOUBLE PRECISION x,v(1)
      DOUBLE PRECISION m0,asound,tsec,xx,vshu
      PARAMETER (m0=0.975)
      EXTERNAL vshu

      asound=0.24*1.d3
      tsec=1.d5*3600.*24.*365.25

      xx=x/(asound*tsec)

      v(1)=-asound*vshu(xx)
c     test
      v(1)=-1.0d0*v(1)

      RETURN
      END
 
c     ------------------------------------------------------------

      FUNCTION vshu(x)

c     Shu's v function

      IMPLICIT NONE
      INTEGER j,n,m,k
      PARAMETER (n=20,m=5)
      DOUBLE PRECISION vshu,x,dd,xv(n),yv(n)
C.....Values of v(x), with x=r/r_s, as given in Shu 1977
      DATA (xv(j),j=1,20)/
     $  0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.55,0.6,0.65,
     $  0.7,0.75,0.8,0.85,0.9,0.95,1.0/
      DATA (yv(j),j=1,20)/
     $    5.44,3.47,2.58,2.05,1.68,1.40,1.18,1.01,0.861,0.735,0.625,0
     $    .528,0.442,0.363,0.291,0.225,0.163,0.106,0.051,0.0/

      if (x.ge.1.) then
        vshu=0.
      else
        if (x.lt.0.05) then
          vshu=dsqrt(2.*0.975/x)
        else
          call locate_v(xv,n,x,j)
          k=min(max(j-(m-1)/2,1),n+1-m)
          call polint_v(xv(k),yv(k),m,x,vshu,dd)
        endif
      endif

      return
      END

c-----------------------------------------------------------------

      SUBROUTINE locate_v(xx,n,x,j)
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

c----------------------------------------------------------------------

      SUBROUTINE polint_v(xa,ya,n,x,y,dy)
      INTEGER n,NMAX
      DOUBLE PRECISION dy,x,y,xa(n),ya(n)
      PARAMETER (NMAX=10)
      INTEGER i,m,ns
      DOUBLE PRECISION den,dif,dift,ho,hp,w,c(NMAX),d(NMAX)
      ns=1
      dif=dabs(x-xa(1))
      do 11 i=1,n
        dift=dabs(x-xa(i))
        if (dift.lt.dif) then
          ns=i
          dif=dift
        endif
        c(i)=ya(i)
        d(i)=ya(i)
11    continue
      y=ya(ns)
      ns=ns-1
      do 13 m=1,n-1
        do 12 i=1,n-m
          ho=xa(i)-x
          hp=xa(i+m)-x
          w=c(i+1)-d(i)
          den=ho-hp
          if(den.eq.0.)pause 'failure in polint_v'
          den=w/den
          d(i)=hp*den
          c(i)=ho*den
12      continue
        if (2*ns.lt.n-m)then
          dy=c(ns+1)
        else
          dy=d(ns)
          ns=ns-1
        endif
        y=y+dy
13    continue
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software =v1.9"217..
