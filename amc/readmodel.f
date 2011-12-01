      SUBROUTINE readmodel(modelfile)

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

c     Reads the input model, in "flexible" format.

      IMPLICIT NONE
      INCLUDE 'amccommon.inc'
      INTEGER i,j,k,length,maxcol,ncol,start
      PARAMETER (maxcol=16)
      INTEGER colno(maxcol)
      DOUBLE PRECISION coldum(maxcol+maxlev-1)
      CHARACTER*32 lcol,ccol
      CHARACTER*80 modelfile,line,keyw,valu,columns
      EXTERNAL length
      LOGICAL debug
      PARAMETER(debug=.false.)

c     i,j,k:     counters
c     length:    function that returns the length of a string
c     maxcol:    maximum number of columns; set to 15
c     ncol:      actual number of columns in input file
c     start:     start position of columns
c     colno:     translate column to location in lcol/ccol
c     coldum:    stores column value before assigning to proper array
c     lcol,ccol: shorthand list of 2-letter column codes
c     modelfile: model input file
c     line:      reads line from input file
c     keyw,valu: used in peeling of header keywords and their values
c     columns:   column description header
c     debug:     turns debugging output on/off

      nlev=maxlev


      lcol='idrarbzazbnhtknmnetevrvzvadbtdlp' ! Currently known columns ids
      ccol='IDRARBZAZBNHTKNMNETEVRVZVADBTDLP'

      rmax=0                    ! Initial settings
      zmax=0
      ncell=0
      tcmb=2.735d0
      gas2dust=100.d0
      columns=' '
      
      open(unit=11,file=modelfile,status='old',err=911)
      if (debug) print*,'[debug] opened model file'

c     Read header in free format

   10 read(11,'(A)') line
      if (line(1:1).eq.'#') goto 10           ! # signals comment line
      if (line(1:1).eq.'@') goto 100          ! @ signals start of grid
 
      do i=1,200                              ! Split line in keyw and valu
        if (line(i:i).eq.'=') then
          do j=i+1,200
            if (line(j:j).eq.' ') goto 20
          enddo
          j=201
          goto 20
        endif
      enddo
      write(*,'(A)') 'AMC: cannot understand input'
      write(*,'(A)') line
      write(*,'(A)') 'AMC: skipping...'
      goto 10
   20 keyw=line(1:i-1)
      valu=line(i+1:j-1)
      if (valu(1:1).eq.' ') then
        write(*,'(A)') 'AMC: cannot understand input'
        write(*,'(A)') line
        write(*,'(A)') 'AMC: skipping...'
        goto 10
      endif
      if (debug) print*,'[debug] read header'

c     Search keywords (case sensitive)

      if (keyw(1:4).eq.'rmax') read(valu,*) rmax
      if (keyw(1:4).eq.'zmax') read(valu,*) zmax
      if (keyw(1:5).eq.'ncell') read(valu,*) ncell
      if (keyw(1:4).eq.'tcmb') read(valu,*) tcmb
      if (keyw(1:8).eq.'gas:dust') read(valu,*) gas2dust
      if (keyw(1:7).eq.'columns') columns=valu(1:length(valu))

      goto 10                   ! Next line


c     End of header reached: check on validity

  100 if (rmax.le.0.) stop 'AMC: <rmax> missing or 0...abort'
      if (ncell.le.0) stop 'AMC: <ncell> missing or 0...abort'
      if (columns(1:1).eq.' ') 
     $  stop 'AMC: <columns> must be defined...abort'

      if (debug) print*,'[debug] validated header'

c     Interpret columns

      do j=1,maxcol
        colno(j)=0
      enddo
      ncol=0
      start=1
      k=length(columns)
      do i=1,k
        if ((columns(i:i).eq.',').or.(i.eq.k)) then
          ncol=ncol+1
          keyw=columns(start:i-1)
          if (i.eq.k) keyw=columns(start:k)
          do j=1,maxcol         ! Retrieve column number
            if ((keyw.eq.lcol((j-1)*2+1:j*2)).or.
     $        (keyw.eq.ccol((j-1)*2+1:j*2))) then
              if (colno(j).eq.0) then
                colno(j)=ncol
                goto 101
              else
                stop 'AMC: ambiguous <columns>...abort'
              endif
            endif
          enddo
          write(*,'(A,A,A)') 
     $      'AMC: unknown column identifier:',keyw,'...abort'
          stop 
  101     start=i+1
        endif
      enddo

      if (debug) print*,'[debug] read column IDs'

c     Check for missing columns

      if (colno(1).eq.0) stop 'AMC: column ID missing...abort'
      if (colno(2).eq.0) stop 'AMC: column RA missing...abort'
      if (colno(3).eq.0) stop 'AMC: column RB missing...abort'
      if (zmax.gt.0) then
        if (colno(4).eq.0) stop 'AMC: column ZA missing...abort'
        if (colno(5).eq.0) stop 'AMC: column ZB missing...abort'
      endif
      if (colno(6).eq.0) stop 'AMC: column NH missing...abort'
      if (colno(7).eq.0) stop 'AMC: column TK missing...abort'
      if (colno(8).eq.0) stop 'AMC: column NM missing...abort'
      if (colno(14).eq.0) stop 'AMC: column DB missing...abort'

      if (debug) print*,'[debug] validated column IDs'

c     Start reading grid in format defined by COLUMNS
c     Convert from 'astro' units to SI
c     Level populations (lp) always start at the last defined column (ncol)

      do i=1,ncell
        if (colno(16).gt.0) then ! Are there populations defined?
          read(11,*) (coldum(j),j=1,ncol-1+nlev)
        else
          read(11,*) (coldum(j),j=1,ncol)
        endif
        if (debug) print*,'[debug] icell = ',i,' coldum = ',coldum
        k=coldum(colno(1))
        ra(k)=coldum(colno(2))                            ! [m]
        rb(k)=coldum(colno(3))                            ! [m]
        if (colno(4).gt.0) za(k)=coldum(colno(4))         ! [m]
        if (colno(5).gt.0) zb(k)=coldum(colno(5))         ! [m]
        nh2(k)=coldum(colno(6))*1.d6                      ! [cm-3] -> [m-3]
        tkin(k)=coldum(colno(7))                          ! [K]
        nmol(k)=coldum(colno(8))*1.d6                     ! [cm-3] -> [m-3]
        if (colno(9).gt.0) ne(k)=coldum(colno(9))*1.d6    ! [cm-3] -> [m-3]
        if (colno(10).gt.0) te(k)=coldum(colno(10))       ! [K]
        if (colno(11).gt.0) vr(k)=coldum(colno(11))*1.d3  ! [km/s] -> [m/s]
        if (colno(12).gt.0) vz(k)=coldum(colno(12))*1.d3  ! [km/s] -> [m/s]
        if (colno(13).gt.0) va(k)=coldum(colno(13))*1.d3  ! [km/s] -> [m/s]
        doppb(k)=coldum(colno(14))*1.d3                   ! [km/s] -> [m/s]
        if (colno(15).gt.0) tdust(k)=coldum(colno(15))    ! [K]
        if (colno(16).gt.0) then
          do j=1,nlev
            pops(j,k)=coldum(ncol-1+j)
          enddo
        endif
      enddo

      close(11)                 ! Close model file

      RETURN

  911 stop 'AMC: Error opening <source>...abort'

      END

