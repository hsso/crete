      SUBROUTINE molinit(molec,molec2)

c (c) Michiel Hogerheijde / Floris van der Tak 2000
c     michiel@strw.leidenuniv.nl, vdtak@sron.nl
c     http://www.sron.rug.nl/~vdtak/ratran
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

c     Reads molecular data from file, calculates up/down collision
c     rates at all grid positions, and writes them out to an array or a
c     disk file for later reference.

      IMPLICIT NONE
      INCLUDE 'amccommon.inc'
      INTEGER id,ilev,iline,itrans,itemp,ntemp,ntemp2,tnint,length
      INTEGER ipart,npart,maxpart,idum
      PARAMETER(maxpart=2)
      DOUBLE PRECISION fac,uprate,downrate,eterm(maxlev),hckb
     $     ,gstat(maxlev),colld(maxpart,maxtrans,max(maxtemp,maxtemp2))
     $     ,temp(maxtemp),temp2(maxtemp2),kappa,planck
      PARAMETER (hckb = 100.d0 * hplanck * clight / kboltz )
      CHARACTER*80 molec,molec2,specref,collref
      EXTERNAL length,kappa,planck
      LOGICAL debug
      PARAMETER(debug=.false.)

c     id,i,ix,itrans,itemp: counters
c     ntemp:      number of different temperatures for which collision
c                 rates are present
c     ntemp2:     idem, for second collision partner
c     tnint:      nearest (lower) Tkin value for wich collision rates in file
c     length:     function that returns length of character string
c     idum:       dummy
c     fac:        unit conversion factor (various uses)
c     tkin:       kinetic temperature at grid point
c     te:         temperature of 2nd collision partner (electrons)
c     uprate,downrate: rates at local temperature
c     eterm:     molecular rotational level energies (cm^-1)
c     hckb:       100.*hplanck*clight / k_boltz: converts energy in cm-1
c                 to J
c     gstat:      statistical weights of rotational levels.
c     colld:      downward rates from molecular data file
c     temp,temp2: temperatures for which downward rates are given
c     molec,collfile: input files, output collision rates file
c     specref     text with reference to spectroscopic data
c     collref     text with reference to collisional data
c     kappa:      (user provided) function that returns dust emissivity
c     debug: turns debugging output on/off

c     ------------------------------------------------------------


      open(unit=11,file=molec,status='old',err=911) ! Open molecular data file
      if (debug) print*,'[debug] opened molecular data file'

c     Preamble
      read(11,*) 
      read(11,*) specref
      read(11,*) 
      read(11,*) amass
      read(11,*) 
      read(11,*) nlev
      if (debug) print*,'[debug] read header'

c     Term energies and statistical weights
      read(11,*)
      do ilev=1,nlev
         read(11,*) idum,eterm(ilev),gstat(ilev)
      enddo
      if (debug) print*,'[debug] read level energies'

c     Radiative upper & lower levels and Einstein coefficients
      read(11,*) 
      read(11,*) nline
      read(11,*)
      do iline=1,nline
         read(11,*) idum,lau(iline),lal(iline),aeinst(iline)
      enddo
      if (debug) print*,'[debug] read Einstein coeffs'

c     Compute line frequencies and Einstein B coeffs (energy in cm-1).
c     Also add thermal broadening to line widths.

      amass=amass*amu
      do iline=1,nline
         freq(iline)=(eterm(lau(iline))-eterm(lal(iline)))*100.d0*clight
         beinstu(iline)=aeinst(iline)*(clight/freq(iline))*(clight
     $        /freq(iline))/(hplanck*freq(iline))/2.d0
         beinstl(iline)=gstat(lau(iline))/gstat(lal(iline))
     $        *beinstu(iline)
      enddo

      do id=1,ncell
        doppb(id)=dsqrt(doppb(id)**2.d0+2.d0*kboltz/amass*tkin(id))
      enddo

c     Collisional data
      second=.false.
      read(11,*)
      read(11,*) npart
      if (npart.gt.1) second=.true.

c     The ugly if-statement should be superseded by loops over ipart
C     throughout the code, multiplying with respective (up to 6) densities. 
c     That's major surgery, left for later times ...

      do ipart=1,npart
         read(11,*)
         read(11,*) collref 
         read(11,*)
         if (ipart.eq.1) then
            read(11,*) ntrans
            read(11,*)
            read(11,*) ntemp
            read(11,*)
            read(11,*) (temp(itemp),itemp=1,ntemp)
            read(11,*)
            do itrans=1,ntrans
               read(11,*) idum,lcu(itrans),lcl(itrans),
     $              (colld(ipart,itrans,itemp),itemp=1,ntemp)
               do itemp=1,ntemp
                  colld(ipart,itrans,itemp)=colld(ipart,itrans,itemp)/1 ! [cm3/s] -> [m3/s]
     $                 .d6 
               enddo
            enddo
         else
            read(11,*) ntrans2
            read(11,*)
            read(11,*) ntemp2
            read(11,*)
            read(11,*) (temp2(itemp),itemp=1,ntemp2)
            read(11,*)
            do itrans=1,ntrans2
               read(11,*) idum,lcu2(itrans),lcl2(itrans),
     $              (colld(ipart,itrans,itemp),itemp=1,ntemp2)
               do itemp=1,ntemp2
                  colld(ipart,itrans,itemp)=colld(ipart,itrans,itemp)/1 ! [cm3/s] -> [m3/s]
     $                 .d6 
               enddo
            enddo
         endif
      enddo
      if (debug) print*,'[debug] read collision data'

      close(11)

c     Finished reading molecular data file.
      
      write(*,'(A)') 'AMC: Read molecular data file ' !//comm_s3(1:46)
      write(*,'(''AMC: Including '',I3,
     $     '' transitions among '',I3,'' levels.'')') nline,nlev

c     ------------------------------------------------------------

c     Calculate upward/downward rates in all non-empty cells.
c     Interpolate downward rates, but do not extrapolate.

      do id=1,ncell

        if (nmol(id).gt.eps) then
          
          do itrans=1,ntrans    ! First partner
            
            if (((tkin(id).gt.temp(1)).and.(tkin(id).lt.temp(ntemp))))
     $        then
              do itemp=1,ntemp-1
                if (tkin(id).gt.temp(itemp).and.(tkin(id).le.temp(itemp
     $            +1)))tnint=itemp
              enddo
              fac=(tkin(id)-temp(tnint))/(temp(tnint+1)-temp(tnint))
              downrate=colld(1,itrans,tnint) +
     $          fac*(colld(1,itrans,tnint+1)-colld(1,itrans,tnint))
            else
              if (tkin(id).le.temp(1)) downrate=colld(1,itrans,1)
              if (tkin(id).ge.temp(ntemp))
     $             downrate=colld(1,itrans,ntemp)
            endif
            
            uprate=gstat(lcu(itrans))/gstat(lcl(itrans))*downrate*
     $        dexp(-hckb*(eterm(lcu(itrans))-eterm(lcl(itrans)))
     $        /tkin(id))
            
            if (disk) then
              write(12,'(2(1X,E10.3))') uprate,downrate
            else
              up(itrans,id)=uprate
              down(itrans,id)=downrate
            endif

          enddo

          if (npart.eq.2) then
            do itrans=1,ntrans2 ! Second partner
              
              if ((te(id).gt.temp2(1)).and.(te(id).lt.temp2(ntemp2)))
     $          then
                do itemp=1,ntemp2-1
                  if ((te(id).gt.temp2(itemp)).and.
     $              (te(id).le.temp2(itemp+1))) 
     $              tnint=itemp
                enddo
                fac=(te(id)-temp2(tnint))/
     $            (temp2(tnint+1)-temp2(tnint))
                downrate=colld(2,itrans,tnint) + fac*
     $            (colld(2,itrans,tnint+1)-colld(2,itrans,tnint))
              else
                if (te(id).le.temp2(1)) downrate=colld(2,itrans,1)
                if (te(id).ge.temp2(ntemp2)) downrate=colld(2,itrans
     $            ,ntemp2)
              endif
              
              uprate=gstat(lcu2(itrans))/gstat(lcl2(itrans))*downrate*
     $          dexp(-hckb*(eterm(lcu2(itrans))-eterm(lcl2(itrans)))
     $          /te(id))
              
c     Write rates to arrays or file (optional).

              if (disk) then
                write(12,'(2(1X,E10.3))') uprate,downrate
              else
                up2(itrans,id)=uprate
                down2(itrans,id)=downrate
              endif

            enddo
          endif

        endif

      enddo
      if (debug) print*,'[debug] interpolated rate coeffs'


c     -----------------------------------------------------------------

c     Also initialize dust emissivity, converting from m2/kg_dust
c     to "m-1 per H2/cm3" (so that tau_dust=knu)

      do iline=1,nline
        norm(iline)=planck(1,tnorm)
        if (tcmb.gt.0.) then
          cmb(iline)=planck(iline,tcmb)/norm(iline)
        else
          cmb(iline)=0.
        endif
      enddo

c     ! Do not normalize dust; will be done in photon !
      do iline=1,nline
        do id=1,ncell
          knu(iline,id)=kappa(id,freq(iline))*2.4d0*amu/gas2dust*nh2(id)
          dust(iline,id)=planck(iline,tdust(id))
        enddo
      enddo

      RETURN

  911 stop 'AMC: error opening molecular data file...abort'
      END

