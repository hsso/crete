      FUNCTION vfunc(s)

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
c     2000, A&A, 362, 697.

c     For revision history see http://www.sron.rug.nl/~vdtak/ratran/

c     Determines projected l.o.s. velocity at position s on l.o.s.
c     It needs to know about start positions l.o.s. segment, which
c     are passed on through vproj* common blocks.

      IMPLICIT NONE
      DOUBLE PRECISION vfunc,s
      
      INTEGER id
      COMMON /vproji/ id
      DOUBLE PRECISION rpos,phi,cosphi,sinphi,vphot
      COMMON /vprojd/ rpos,phi,cosphi,sinphi,vphot
      DOUBLE PRECISION r,phis,psis,v(1)

c     vfunc  -- this function
c     s      -- position along line-of-sight
c     id     -- grid cell coordinate
c     rpos,phi,cosphi,sinphi -- position, direction at start of l.o.s. (s=0)
c     psis   -- angle between phi and direction at s
c     r,phis -- physical location, direction at s
c     vphot  -- velocity offset of photon
c     v(1)   -- velocity (vector) of gas at s


c     Get direction and position at location s along l.o.s. 

      psis=datan2(s*sinphi,rpos+s*cosphi)
      phis=phi-psis
      r=dsqrt(rpos**2.d0+s**2.d0+2.d0*rpos*s*cosphi)

c     Get velocity vector of the gas at this position
c     velo is the user-provided velocity field.
      
      call velo(id,r,v)

c     vfunc is velocity difference between betweeb photon and gas
c     projected on l.o.s.

      vfunc=vphot-dcos(phis)*v(1)

      RETURN
      END
