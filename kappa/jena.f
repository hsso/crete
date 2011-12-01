      FUNCTION kappa(id,nu)

c (c) Michiel Hogerheijde / Floris van der Tak 2000
c     michiel@strw.leidenuniv.nl, vdtak@sron.nl
c     http://www.sron.rug.nl/~vdtak/ratran/

c     This file is part of the 'ratran' molecular excitation and
c     radiative transfer code. The one-dimensional version of this code
c     is publicly available; the two-dimensional version is available on
c     collaborative basis. Although the code has been thoroughly tested,
c     the authors do not claim that it is free of errors or that it gives
c     correct results in all situations. Any publication making use of
c     this code should include a reference to Hogerheijde & van der Tak
c     2000, A&A 362, 697.

c     For revision history see http://www.sron.rug.nl/~vdtak/ratran/

c     Returns dust emissivity in m2/kg_dust at frequency freq and
c     in cell id.

c     Ossenkopf & Henning Jena models.
c
c     Useage in amc.inp:  kappa=jena,TYPE,COAG
c     where TYPE = bare / thin / thick
c     and   COAG = no / e5 / e6 / e7 / e8

      IMPLICIT NONE
      INTEGER id,i,n,j,k,m
      PARAMETER (n=67,m=3)
      DOUBLE PRECISION kappa,nu,lamtab(n),kaptab(n),logkap,loglam,dummy
     $  ,clight
      PARAMETER (clight=2.997924562d8) ! speed of light
      SAVE lamtab,kaptab

c     Read table on first call; lamtab(i) might be 0, but only for one i
      if ((lamtab(1).eq.0).and.(lamtab(2).eq.0)) then
        open(22,file='ratranjena.tab',status='old')
        do i=1,n
          read(22,*) lamtab(i),kaptab(i)
          lamtab(i)=dlog10(lamtab(i)/1.d6)
          kaptab(i)=dlog10(kaptab(i))
        enddo
        close(22)
      endif

c     Interpolate/extrapolate table (in log units)

      loglam=dlog10(clight/nu)
      call locate(lamtab,n,loglam,j)
      k=min(max(j-(m-1)/2,1),n+1-m)
      call polint(lamtab(k),kaptab(k),m,loglam,logkap,dummy)

c     ...in m2/kg_dust: (0.1 converts cm2/g_dust to m2/kg_dust)
      kappa=0.1*10.**logkap


      RETURN
      END

