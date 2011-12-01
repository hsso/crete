      SUBROUTINE readmodel(modelfile,molfile,tkin)

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

c     Reads the source model with or without populations

      IMPLICIT NONE
      INCLUDE 'skycommon.inc'
      INTEGER i,j,k,length,maxcol,ncol,start
      PARAMETER (maxcol=16)
      INTEGER colno(maxcol)
      DOUBLE PRECISION coldum(maxcol+maxlev-1)
      DOUBLE PRECISION tkin(maxcell),totpop
      CHARACTER*32 lcol,ccol
      CHARACTER*80 modelfile,molfile,line,keyw,valu,columns
      EXTERNAL length
      LOGICAL debug
      PARAMETER (debug=.false.)
      
c     modelfile :        name of input file
c     molfile :          name of molecular data file
c     tkin:              kinetic temperature
c     i,j,k:             counters
c     maxcol,ncol,start: helps in reading input file
c     colno:             column number of variable in modelfile
c     columns,coldum:    helps in reading columns
c     lcol,ccol:         string of column identifyers
c     line,keyw,valu:    reads header of modelfile
c     totpop:            checks if populations sum to unity
c     debug:             turns debugging on/off

      lcol='idrarbzazbnhtknmnetevrvzvadbtdlp' ! Currently known column IDs
      ccol='IDRARBZAZBNHTKNMNETEVRVZVADBTDlp'

      rmax=0                    ! Initial settings
      zmax=0
      ncell=0
      tbg=2.735d0
      gas2dust=100.d0
      columns=' '
      molfile='continuum'

      nlev=maxlev               ! sky+ c-script makes sure that maxlev=nlev
      
      modelfile=modelfile(1:length(modelfile))

      open(unit=11,file=modelfile,status='old',err=911)
      if (debug) print*,'[debug] opened modelfile'

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
      write(*,'(A)') 'SKY: cannot understand input'
      write(*,'(A)') line
      write(*,'(A)') 'SKY: skipping...'
      goto 10
   20 keyw=line(1:i-1)
      valu=line(i+1:j-1)
      if (valu(1:1).eq.' ') then
        write(*,'(A)') 'SKY: cannot understand input'
        write(*,'(A)') line
        write(*,'(A)') 'SKY: skipping...'
        goto 10
      endif


c     Search keywords (case sensitive)

      if (keyw(1:4).eq.'rmax') read(valu,*) rmax
      if (keyw(1:4).eq.'zmax') read(valu,*) zmax
      if (keyw(1:5).eq.'ncell') read(valu,*) ncell
      if (keyw(1:4).eq.'tcmb') read(valu,*) tbg
      if (keyw(1:8).eq.'gas:dust') read(valu,*) gas2dust
      if (keyw(1:7).eq.'molfile') molfile=valu(1:length(valu))
      if (keyw(1:7).eq.'columns') columns=valu(1:length(valu))

      goto 10                   ! Next line


c     End of header reached: check on validity

  100 if (rmax.le.0.) stop 'SKY: <rmax> missing or 0...abort'
      if (ncell.le.0) stop 'SKY: <ncell> missing or 0...abort'
      if (columns(1:1).eq.' ') 
     $  stop 'SKY: <columns> must be defined...abort'

      if (debug) print*,'[debug] done reading header'

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
                stop 'SKY: ambiguous <columns>...abort'
              endif
            endif
          enddo
          write(*,'(A,A,A)') 
     $      'SKY: unknown column identifier:',keyw,'...abort'
          stop 
  101     start=i+1
        endif
      enddo

      if (debug) print*,'[debug] ncell,rmax: ',ncell,rmax

c     Check for missing columns

      if (colno(1).eq.0) stop 'SKY: column ID missing...abort'
      if (colno(2).eq.0) stop 'SKY: column RA missing...abort'
      if (colno(3).eq.0) stop 'SKY: column RB missing...abort'
      if (zmax.gt.0) then
        if (colno(4).eq.0) stop 'SKY: column ZA missing...abort'
        if (colno(5).eq.0) stop 'SKY: column ZB missing...abort'
      endif
      if (colno(6).eq.0) stop 'SKY: column NH missing...abort'
      if (colno(7).eq.0) stop 'SKY: column TK missing...abort'
      if (colno(8).eq.0) stop 'SKY: column NM missing...abort'
      if (colno(14).eq.0) stop 'SKY: column DB missing...abort'


c     Start reading grid in format defined by COLUMNS
c     Convert from 'natural' units to SI
c     Level populations (lp) always start at the last defined column (maxcol)

      do i=1,ncell
c        if (debug) print*,'[debug] icell=',i
        if (colno(16).gt.0) then ! Are there populations defined?
          read(11,*) (coldum(j),j=1,ncol-1+nlev)
        else
          read(11,*) (coldum(j),j=1,ncol)
        endif
        if (debug) print*,'[debug] coldum = ',coldum
        k=coldum(colno(1))
        ra(k)=coldum(colno(2))                            ! [m]
        rb(k)=coldum(colno(3))                            ! [m]
        if (colno(4).gt.0) za(k)=coldum(colno(4))         ! [m]
        if (colno(5).gt.0) zb(k)=coldum(colno(5))         ! [m]
        nh2(k)=coldum(colno(6))*1.d6                      ! [cm-3] -> [m-3]
        tkin(k)=coldum(colno(7))                          ! [K]
        nmol(k)=coldum(colno(8))*1.d6                     ! [cm-3] -> [m-3]
        if (colno(11).gt.0) vr(k)=coldum(colno(11))*1.d3  ! [km/s] -> [m/s]
        if (colno(12).gt.0) vz(k)=coldum(colno(12))*1.d3  ! [km/s] -> [m/s]
        if (colno(13).gt.0) va(k)=coldum(colno(13))*1.d3  ! [km/s] -> [m/s]
        doppb(k)=coldum(colno(14))*1.d3                   ! [km/s] -> [m/s]
        if (colno(15).gt.0) tdust(k)=coldum(colno(15))    ! [K]
        if (debug) print*,'[debug] reading pops'
        if (colno(16).gt.0) then
           totpop = 0.d0
          do j=1,nlev
            pops(j,k)=coldum(ncol-1+j)
            totpop = totpop + pops(j,k)
          enddo
c     Check if populations sum roughly to unity in non-empty cells
          if (dabs(totpop-1.d0).gt.1.d-3.and.nmol(k).gt.eps) then
             print*,'ERROR: Populations do not sum to unity in cell ',k
ccc             print*,'nlev, totpop = ',nlev, totpop
             stop
          endif
        endif
      enddo

      close(11)                 ! Close model file

      RETURN


  911 stop 'SKY: Error opening <source>...abort'

      END
