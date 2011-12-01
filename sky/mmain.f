      PROGRAM sky

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

c     SKY generates MIRIAD images/cubes of the sky brightness, based
c     on a model with level populations and/or dust emissivity values.

c     Main program of sky.

      IMPLICIT NONE
      INCLUDE 'skycommon.inc'
      CHARACTER filein*60,fileout*60,outname*60
      INTEGER ix,iy,ixprime,ichan,iline,lout(maxrx)
      INTEGER length
      EXTERNAL length
      REAL row(maxsky)          ! REAL!!
      LOGICAL debug
      PARAMETER (debug=.false.)

c     filein           : name input file
c     fileout          : (prefix) of output file/image
c     outname          : completed output name; should be CHAR*60 !!!
c     ix,iy,ichan,iline: used to loop through sky grid, passband and
c                         transitions
c     lout(iline)      : miriad file handle number for image(iline)
c     ixprime          : de-rotated ix correspodong to (ix,iy) (in 1D)
c     row              : contains image row prior to writing to miriad file
c     debug            : turns debugging on/off

      write(*,'(A)') 'SKY:'
      write(*,'(A)') 'SKY: Starting calculations'
      write(*,'(A)') 'SKY:'


c     Read inputs, model, molecular data, fg model...

      if (debug) write(*,*) '[debug] calling getinputs'
      call getinputs (filein,fileout)
      if (debug) write(*,*) '[debug] back from getinputs'

c     Open output images "outfile_filter(iline)" or "outfile_nu(iline)"

      do iline=1,nrx
        if (filter(1).eq.0) then
          write(outname,'(A,''_'',1p,E9.3)') 
     $      fileout(1:length(fileout)),nu(iline)
        else
          write(outname,'(A,''_'',I3.3)') 
     $      fileout(1:length(fileout)),filter(iline)
        endif
        call openim(lout(iline),outname,nsky,nsky,nchan)
      enddo
      if (debug) write(*,*) '[debug] opened output files'

c     Do radiative transfer

      write(*,'(A)') 'SKY: Starting line of sight integration'

      do ix=1,nsky              ! Initialize sky frame row
        do ichan=1,nchan
          do iline=1,nrx
            intens(iline,ichan,ix)=0.d0
          enddo
        enddo
      enddo

c     Do l.o.s. integration.

      if (twodee) then

        do iy=1,nsky            ! 2D: loop through full image
          do ix=1,nsky
            call losintegr(ix,iy)
            call progress(ix,iy)
          enddo                 ! ix=1,nsky
          do ichan=1,nchan
            do iline=1,nrx
              do ix=1,nsky
                row(ix)=real(
     $            (ucon(iline)*norm(iline))*intens(iline,ichan,ix))
                if (row(ix).lt.real(eps)) row(ix)=real(eps)
              enddo             ! ix=1,nsky
              call xysetpl(lout(iline),1,ichan)
              call xywrite(lout(iline),iy,row)
            enddo               ! iline=1,nrx
          enddo                 ! ichan=1,nchan
        enddo                   ! iy=1,nsky
        write(*,'(A)') 'SKY: Completed 100%'

      else
        
        do ix=int(xycen+1.),nsky
          call losintegr(ix,ix)
          call progress(ix,ix)
        enddo                   ! ix=xycen,nsky
        write(*,'(A)') 'SKY: Completed 100%'
        write(*,'(A)') 'SKY: Expanding 1D to full image plane'
        do iy=1,nsky            ! 1D: ...and rotate to get full image
          do iline=1,nrx
            do ichan=1,nchan
              do ix=1,nsky
                if ((ix.eq.iy).and.(ix.ge.int(xycen+1.))) then
                  row(ix)=real((ucon(iline)*norm(iline))*
     $              intens(iline,ichan,ix))
                  if (row(ix).lt.eps) row(ix)=eps
                else
                  ixprime=nint(xycen+0.5+dsqrt(
     $              (dble(ix)-xycen-0.5)**2.+(dble(iy)-xycen-0.5)**2.)
     $              /dsqrt(2.d0))
                  if ((ixprime.ge.int(xycen+1.)).and.(ixprime.le.nsky))
     $              then
                    row(ix)=real((ucon(iline)*norm(iline))*
     $                intens(iline,ichan,ixprime))
                    if (row(ix).lt.eps) row(ix)=eps
                  else
                    row(ix)=eps
                  endif
                endif
              enddo             ! ix=1,nsky
              call xysetpl(lout(iline),1,ichan)
              call xywrite(lout(iline),iy,row)
            enddo               ! iline=1,nrx
          enddo                 ! ichan=1,nchan
        enddo                   ! iy=1,nsky

      endif                     ! (twodee)


      do iline=1,nrx           ! Close output images
        call closeim(lout(iline))
      enddo
      write(*,'(A)') 'SKY:'
      write(*,'(A,A)') 'SKY: Written output to ',
     $  fileout(1:length(fileout))
      write(*,'(A)') 'SKY:'


      STOP
      END

c      ---------------------------------------------------------------------

      SUBROUTINE progress(ix,iy)

c     Keeps track of progress

      IMPLICIT NONE
      INCLUDE 'skycommon.inc'
      INTEGER k,done,i,total,ix,iy,length
      EXTERNAL length
      SAVE k,done

      k=k+1

      if (twodee) then
        total=nsky*nsky
      else
        total=nsky-int(xycen+1.)
      endif
      
      do i=9,done,-1
        if (dble(k)/dble(total).gt.dble(i)/10.d0) then
          if (done.gt.0) write(*,'(A,I3,A)') 'SKY: Completed ',10*i,'%'
          done=done+1
          goto 100
        endif
      enddo
     

  100 if ((ix.eq.nint(xycen+0.5)).and.(iy.eq.nint(xycen+0.5))) then
        write(*,'(A)') 'SKY:'
        write(*,'(A)')
     $    'SKY: Statistics at line and source center:'
        if (filter(1).ne.0) then
          write(*,'(A,A,A)') 
     $          'SKY: line |   opacity  |  intensity ('
     $          ,units(1:length(units)),')'
          do i=1,nrx
            write(*
     $        ,'(''SKY:  '',I3.3,1p,'' | '',E10.3,'' | '',E10.3)')
     $            filter(i),taucen(i,vcen),
     $            (ucon(i)*norm(i))*intens(i,vcen,ix) 
          enddo
        else
          write(*,'(A,A,A)') 
     $      'SKY:  frequency |   opacity  |  intensity ('
     $      ,units(1:length(units)),')'
          do i=1,nrx
            write(*
     $        ,'(''SKY: '',1p,E10.3,'' | '',E10.3,'' | '',E10.3)')
     $        nu(i),taucen(i,vcen),
     $        (ucon(i)*norm(i))*intens(i,vcen,ix)
          enddo
        endif
        write(*,'(A)') 'SKY:'
      endif


      RETURN
      END
