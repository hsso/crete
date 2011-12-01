      SUBROUTINE openim(lout,imname,nxdim,nydim,nzdim)
	
c (c) Michiel Hogerheijde / Floris van der Tak 2000-2006
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


c     Opens miriad image and writes header; requires `mirlib` in compile.
c      
c     inputs: 
c     imname -- name of file to contain image
c     nxdim  -- number of elements along x axis of total image
c     nydim  -- number of elements along y axis of total image
c     nzdim  -- number of elements along z axis of total image
c     angres -- angular resolution of grid
c     velres -- velocity resolution of grid
c      
c     output: 
c     lout -- miriad generated handle for file -- use this
c             to write to file later and close file
c
c     The magic combination to force the correct angular scale is:
c     RA---SIN, DEC--SIN and angular_resolution(in radians)
c
c     Note that quantities in SI are transformed back in more manageable ones.
      
      IMPLICIT NONE
      INCLUDE 'skycommon.inc'
      
      INTEGER lout,nxdim,nydim,nzdim
      CHARACTER imname*60
      INTEGER naxis,nsize(3)
      COMMON /image/ naxis,nsize
      
      naxis = 3
      nsize(1) = nxdim
      nsize(2) = nydim
      nsize(3) = nzdim
      call xyopen(lout,imname,'new',naxis,nsize)
      call wrhdi(lout,'naxis',naxis)
      call wrhdi(lout,'naxis1',nsize(1))
      call wrhdi(lout,'naxis2',nsize(2))
      call wrhdi(lout,'naxis3',nsize(3))
      call wrhda(lout,'ctype1','RA---SIN')
      call wrhda(lout,'ctype2','DEC--SIN')
      call wrhda(lout,'ctype3','VELO-LSR')
      call wrhdd(lout,'crpix1',dble(nxdim/2+0.5))
      call wrhdd(lout,'crpix2',dble(nydim/2+0.5))
      call wrhdd(lout,'crpix3',dble(nzdim/2+1))
      call wrhdd(lout,'cdelt1',dble(-angres))
      call wrhdd(lout,'cdelt2',dble(angres))
      call wrhdd(lout,'cdelt3',dble(velres*1.e-3))
      call wrhdd(lout,'crval1',dble(0.))
      call wrhdd(lout,'crval2',dble(0.))
      call wrhdd(lout,'crval3',dble(0.))
      call wrhda(lout,'object',imname)
      if ((units(1:1).eq.'K').or.(units(1:1).eq.'k'))
     $  call wrhda(lout,'bunit','K')
      if ((units(1:4).eq.'Jypx').or.(units(1:4).eq.'jypx'))
     $  call wrhda(lout,'bunit','JY/PIXEL')
      if ((units(1:7).eq.'Wm2Hzsr').or.(units(1:7).eq.'wm2hzsr'))
     $  call wrhda(lout,'bunit','Wm2Hzsr')
      if ((units(1:3).eq.'lnu').or.(units(1:3).eq.'Lnu'))
     $  call wrhda(lout,'bunit','Lsun/PIXEL')
      
      RETURN
      END

