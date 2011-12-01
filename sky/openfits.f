      SUBROUTINE openfits(lout,outname,naxis,naxes,status)
      IMPLICIT none

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
c     this code should include a reference to Hogerheijde & van der Tak
c     2000, A&A 362, 697.

c     fvdt  20nov00  initial version 
c     fvdt  02jul09  added status=0 on line 2 (some compilers do not initialize)

      INCLUDE 'skycommon.inc'

      INTEGER lout,blocksize,bitpix,status,naxis,naxes(3)
      LOGICAL simple,extend
      CHARACTER outname*60

      blocksize = 1
      status    = 0
      call ftinit(lout,outname,blocksize,status)
      if (status.ne.0) then
         write(*,'(A)') 'SKY: WARNING: overwriting existing FITS file'
         write(*,'(A)') 'SKY:'
         status=0
         call ftnopn(lout,outname,'',status)
         call ftdelt(lout,status)
         call ftinit(lout,outname,blocksize,status)
      endif

C     Write FITS headers
      simple=.true.
      bitpix=-32
      naxes(1)=nsky
      naxes(2)=nsky
      if (naxis.eq.3) naxes(3)=nchan
      extend=.true.
      call ftphps(lout,bitpix,naxis,naxes,status)

      call ftpkys(lout,'CTYPE1','RA---SIN','',status)
      call ftpkys(lout,'CTYPE2','DEC--SIN','',status)
      call ftpkyd(lout,'CDELT1',(-1.8d2*angres/pi),-8,'',status)
      call ftpkyd(lout,'CDELT2',(angres*1.8d2/pi),-8,'',status)
      call ftpkyd(lout,'CRPIX1',dble(nsky/2+0.5),-8,'',status)
      call ftpkyd(lout,'CRPIX2',dble(nsky/2+0.5),-8,'',status)
      call ftpkyd(lout,'CRVAL1',0.d0,-8,'',status)
      call ftpkyd(lout,'CRVAL2',0.d0,-8,'',status)
      call ftpkys(lout,'CTYPE3','VELO-LSR','',status)
      call ftpkyd(lout,'CDELT3',dble(velres),-8,'',status)
      call ftpkyd(lout,'CRPIX3',dble(nchan/2+1),-8,'',status)
      call ftpkyd(lout,'CRVAL3',0.d0,-8,'',status)
c     note: default velocity unit is km/s in Fits format.

      if ((units(1:1).eq.'K').or.(units(1:1).eq.'k'))
     $  call ftpkys(lout,'bunit','K','',status)
      if ((units(1:4).eq.'Jypx').or.(units(1:4).eq.'jypx'))
     $  call ftpkys(lout,'bunit','JY/PIXEL','',status)
      if ((units(1:7).eq.'Wm2Hzsr').or.(units(1:7).eq.'wm2hzsr'))
     $  call ftpkys(lout,'bunit','Wm2Hzsr','',status)
      if ((units(1:3).eq.'lnu').or.(units(1:3).eq.'Lnu'))
     $  call ftpkys(lout,'bunit','Lsun/PIXEL','',status)

      call ftpkys(lout,'OBJECT',outname,'',status)

      return
      end
