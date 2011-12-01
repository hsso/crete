      FUNCTION kappa(id,nu)

c (c) Michiel Hogerheijde / Floris van der Tak 2000
c     michiel@strw.leidenuniv.nl, vdtak@sron.nl
c     http://www.sron.rug.nl/~vdtak/ratran/
c
c     This file is part of the 'ratran' molecular excitation and
c     radiative transfer code. The one-dimensional version of this code
c     is publicly available; the two-dimensional version is available on
c     collaborative basis. Although the code has been thoroughly tested,
c     the authors do not claim that it is free of errors or that it gives
c     correct results in all situations. Any publication making use of
c     this code should include a reference to Hogerheijde & van der Tak,
c     2000, A&A 362, 697.

c     For revision history see http://www.sron.rug.nl/~vdtak/ratran/

c     Returns dust emissivity in m2/kg_dust at frequency freq and
c     in cell id.

c     Simple power law behavior; nu0/kappa0/beta are defined through
c     common.inc
c
c     Useage in amc.inp:  kappa=powerlaw,NU0,KAPPA0,BETA
c     where NU0 is in Hz, KAPPA0 in cm2/g_dust, and BETA is freq.index.

      IMPLICIT NONE
      INCLUDE 'kappacommon.inc'
      INTEGER id
      DOUBLE PRECISION kappa,nu

c     ...in m2/kg_dust: (0.1 converts cm2/g to m2/kg)
      kappa=0.1*kappa0*(nu/nu0)**beta

      RETURN
      END
