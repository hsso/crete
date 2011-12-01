      FUNCTION kappa(id,freq)

c (c) Michiel Hogerheijde / Floris van der Tak 2000
c     michiel@strw.leidenuniv.nl, vdtak@sron.nl
c     http://www.sron.rug.nl/~vdtak/ratran/

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

c     kappa=0 default

      IMPLICIT NONE
      INTEGER id
      DOUBLE PRECISION kappa,freq

c     ...in m2/kg_dust:
      kappa=0.

      RETURN
      END
