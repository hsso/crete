      PROGRAM amc

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

c     Main program of amc.


      IMPLICIT NONE
      INCLUDE 'amccommon.inc'
      INTEGER id,ilev,iter,niter,dummyi,newphot,conv,exceed,totphot,
     $  totphot2,stage,fixseed
      DOUBLE PRECISION var,goal,snr,minsnr,dummy,
     $  opops(maxlev,maxcell),oopops(maxlev,maxcell),percent,avepops
     $  ,ran1,fixset,staterr
      CHARACTER*80 molfile,modelfile,outfile,collfile
      LOGICAL trace
      EXTERNAL ran1
      LOGICAL debug
      PARAMETER(debug=.false.)

c     id:            cell counter
c     ilev:          line counter
c     iter:          iteration counter
c     niter:         required number of iterations
c     dummyi:        integer dummy variable
c     newphot:       number of photons for next iteration
c     conv:          number of converged cells
c     exceed:        number of cells where required photons exceeds max_phot
c     totphot:       total number of photons in this round
c     totphot2:      total number of photons in next round
c     var:           maximum difference between last 3 solutions and
c     _                their average
c     goal:          requested s/n ratio
c     snr:           lowest s/n ratio (=1/var) in cell
c     minsnr:        lowest s/n ratio of all cells
c     dummy:         a dble dummy
c     opops:         previous solution
c     oopops:        the solution before that
c     percent:       percentage of converged cells
c     avepops:       the average of the last 3 solutions
c     ran1:          the NumRep ran1 function
c     molfile:         molecular data file
c     modelfile:     the input model file
c     outfile:       the output file
c     collfile:      the collision rate file (if used)
c     stage:         =1 FIXSET, =2 RANDOM
c     fixseed:       remembers the random seed, for the FIXSET stage
c     fixset:        relative accuracy of fixset stage
c     staterr:       keeps track of non-convergence in stateq
c     trace:         logical to signal if convergence history should be written
c     debug: turns debugging output on/off


      write(*,'(A)') 'AMC: '
      write(*,'(A)') 'AMC: Starting calculations'
      write(*,'(A)') 'AMC:'


c     Read inputs, model, collision rates (if disk true), molecular data

      do id=1,maxcell
        do ilev=1,maxlev
          pops(ilev,id)=0.d0     ! Initialize populations
        enddo
      enddo
      if (debug) print*,'[debug] calling getinputs'
      call getinputs 
     $  (molfile,modelfile,outfile,goal,newphot,fixset,trace)
      if (debug) print*,'[debug] calling readmodel'
      call readmodel(modelfile)
      if (disk) open(12,file=collfile,status='unknown')
      if (debug) print*,'[debug] calling molinit'
      call molinit(molfile)


c     Set-up for Monte Carlo simulation

      do id=1,ncell
        nphot(id)=newphot       ! Set nphot to initial number
      enddo
      niter=ncell               ! Estimated crossing time
      dummyi=-1                 ! Initialize random number generator
      dummy=ran1(dummyi)

      write(*,'(A)') 'AMC:'
      write(*,'(A,1p,E11.5)') 
     $  'AMC: Starting with FIXSET convergence; limit=',fixset

      stage=1                   ! 1=initial phase with fixed photon paths=FIXSET
      fixseed=seed
      percent=0

   50 conv=0
      exceed=0
      totphot=0
      totphot2=0
      minsnr=1./fixset          ! fixset is smallest number to be counted
      staterr=0.

      do id=1,ncell             ! Loop over all cells          
        do iter=1,3             ! always do sets of 3 iterations to build snr

          if (stage.eq.1) then  ! Stage 1=FIXSET -> re-initialize ran1 each time
            dummyi=-1
            dummy=ran1(dummyi)
            seed=fixseed
          endif

          do ilev=1,nlev        ! Keep pops statistics
            oopops(ilev,id)=opops(ilev,id)
            opops(ilev,id)=pops(ilev,id)
          enddo

          if (nh2(id).ge.eps) then ! Propagate photons and solve excitation
             if (debug) print*,'[debug] calling photon for cell ',id
            call photon(id)
            if (debug) print*,'[debug] calling stateq for cell ',id
            call stateq(id,collfile,stage,staterr)
          endif

        enddo                   ! iter=1,3

        if (debug) print*,'[debug] calculate s/n for cell ',id

        snr=fixset               ! Determine snr in cell
        var=0.d0
        totphot=totphot+nphot(id)
        do ilev=1,nlev
          avepops=(pops(ilev,id)+opops(ilev,id)+oopops(ilev,id))/3.d0
          if (avepops.ge.minpop) then
            var=dmax1(dabs(pops(ilev,id)-avepops)/avepops,
     $        dabs(opops(ilev,id)-avepops)/avepops,
     $        dabs(oopops(ilev,id)-avepops)/avepops)
            snr=dmax1(snr,var)
          endif
        enddo
        snr=1./snr
        minsnr=dmin1(snr,minsnr)

        if (stage.eq.1) then
          if (snr.ge.1./fixset) conv=conv+1 ! Stage 1=FIXSET 
        else
          if (snr.ge.goal) then
            conv=conv+1
          else
            newphot=nphot(id)*2 ! Double photons if cell not converged
            if (newphot.gt.max_phot) then
              newphot=max_phot
              exceed=exceed+1
              write(*,'(A,I4)') 'AMC: *** Limiting nphot in cell ',id
            endif
            nphot(id)=newphot
          endif
        endif

        totphot2=totphot2+nphot(id)

      enddo                     ! id=1,ncell


c     Report any convergence problems if they occurred

      if (staterr.gt.0.d0) write(*,'(A,1p,E11.5,A)')
     $  '### WARNING: stateq did not converge everywhere (err='
     $  ,staterr,')'

      if (stage.eq.1) then

        percent=dble(conv)/dble(ncell)*100.d0
        call blowpops(outfile,molfile,goal,minsnr,percent,stage,fixset
     $    ,trace)
        write(*,'(A,1p,E11.5,A,0p,F6.2,A)') 
     $    'AMC: FIXSET fractional error ',1./minsnr,', ',percent
     $    ,'% converged'
        if (conv.eq.ncell) then
          stage=2
          write(*,'(A)') 'AMC:'
          write(*,'(A)') 
     $      'AMC: FIXSET convergence reached...starting RANDOM'
          write(*,'(A)') 'AMC:'
          write(*,'(A,A)') 
     $      'AMC: minimum S/N  |  converged  |     photons  |  ',
     $      'increase to'
          write(*,'(A,A)')
     $      'AMC: -------------|-------------|--------------|--',
     $      '-----------'
        endif
        goto 50

      else

        if (conv.eq.ncell) then
          percent=100.d0
        else
          if (exceed.lt.ncell) then
            percent=dble(conv)/dble(ncell)*100.d0
            call blowpops(outfile,molfile,goal,minsnr,percent,stage
     $        ,fixset,trace)
            write(*,1000) minsnr,percent,totphot,totphot2
 1000       format('AMC:',1p,E12.5,0p,'  |  ',F9.2,'% |  ',I10,'  |  '
     $        ,I10)
            goto 50             ! Next iteration
          else
            write(*,'(A)')
     $        '### WARNING: Insufficient photons. Not converged.'
          endif                 ! if (exceed.lt.ncell)
        endif                   ! if (conv.eq.ncell)

      endif                     ! if (stage.eq.1)


c     Convergence reached (or bailed out)
      
   60 write(*,1001) minsnr,percent,totphot
 1001 format('AMC:',1p,E12.5,0p,'  |  ',F9.2,'% |  ',I10
     $  ,'  |   converged')
   70 write(*,'(A)') 'AMC:'

      call blowpops(outfile,molfile,goal,minsnr,percent,stage,fixset
     $  ,trace)
      write(*,'(A,A)') 'AMC: Written output to ',outfile
      if (disk) close(12)

      END
