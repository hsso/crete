      SUBROUTINE velo(id,x,v)

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


c     This the default velocity routine for sky and AMC 1D, with
c     the systemic velocity vector set to 0.
c     Do not tamper with this file.
c     Copy v_1d.example.f to your working directory for an example of
c     how to create your own velocity routine.

c     velo: user supplied subroutine returning velocity vector as
c     function of physical coordinates and grid position. The latter
c     ensures that possible discontinuities always fall along grid cell
c     boundaries.

      IMPLICIT NONE
      INTEGER id
      DOUBLE PRECISION x,v(1)

c     ix      cell coordinate
c     x       physical coordinate
c     v       v(1)=v_R  in m/s

      v(1)=0.

      RETURN
      END

