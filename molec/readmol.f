      PROGRAM readmol

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

c     A utility to read the molecular data files, and return
c     maxlev, maxline, maxtrans,maxtrans2,maxtemp,maxtemp2

      IMPLICIT NONE
      CHARACTER*80 molfile
      INTEGER length
      EXTERNAL length
      REAL realy
      INTEGER ilev,nlev,iline,nline,itrans,ntrans,itemp,ntemp,ipart
     $     ,npart,ntrans2,ntemp2,inty

      read(*,'(A)') molfile
      molfile=molfile(1:length(molfile))

      open(unit=11,file=molfile,status='old',err=1999)

c     Preamble
      read(11,*) 
      read(11,*) !specref
      read(11,*) 
      read(11,*) !amass
      read(11,*) 
      read(11,*) nlev

c     Term energies and statistical weights
      read(11,*)
      do ilev=1,nlev
         read(11,*) inty,realy,realy !dummy,eterm(ilev),gstat(ilev)
      enddo

c     Radiative upper & lower levels and Einstein coefficients
      read(11,*) 
      read(11,*) nline
      read(11,*)
      do iline=1,nline
         read(11,*) inty,inty,inty,realy,realy  !dummy,lau(iline),lal(iline),aeinst(iline),freq(iline),eup(iline)
      enddo

      read(11,*)
      read(11,*) npart
      if (npart.eq.1) then
         ntrans2=1
         ntemp2=1
      endif
      if (npart.gt.2) goto 2001

      do ipart=1,npart
         read(11,*)
         read(11,*) !collref 
         read(11,*)
         if (ipart.eq.1) then
            read(11,*) ntrans
            read(11,*)
            read(11,*) ntemp
            read(11,*)
            read(11,*) (realy,itemp=1,ntemp)  !(temp(itemp),itemp=1,ntemp)
            read(11,*)
            do itrans=1,ntrans
               read(11,*) inty,inty,inty,(realy,itemp=1,ntemp)
c               read(11,*) dummy,lcu(itrans),lcl(itrans),
c     $              (colld(ipart,itrans,itemp),itemp=1,ntemp)
            enddo
         else
            read(11,*) ntrans2
            read(11,*)
            read(11,*) ntemp2
            read(11,*)
            read(11,*) (realy,itemp=1,ntemp2)  !(temp2(itemp),itemp=1,ntemp2)
            read(11,*)
            do itrans=1,ntrans2
               read(11,*) inty,inty,inty,(realy,itemp=1,ntemp2)
c               read(11,*) dummy,lcu(itrans),lcl(itrans),
c     $              (colld(ipart,itrans,itemp),itemp=1,ntemp2)
            enddo
         endif
      enddo

      close(11)

      write(*,*) nlev
      write(*,*) nline
      write(*,*) ntrans
      write(*,*) ntrans2
      write(*,*) ntemp
      write(*,*) ntemp2

      STOP
 1999 write(*,*) 'ERROR: cannot open molecular data file'
 2001 write(*,*) 'Error: Maximum 2 collision partners'
      stop
      END

c     ----------------------------------------------------------------------

      FUNCTION length(str)
      INTEGER length,maxl,i
      PARAMETER(maxl=200)
      CHARACTER*200 str

      do i=1,maxl
         if (str(i:i).eq.' ') then
            length=i-1
            RETURN
         endif
      enddo
      STOP 'Error: File name too long'
      END

