      SUBROUTINE stateq(id,collfile,stage,staterr)

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

c     Solves statistical equilibrium in cell id subject to the incident
c     radiation field determined in "photon" + the local contribution 
c     over "ds" and "vfac". Uses LU decomposition and back substitution
c     to solve the level population equations. Calls getjbar to get
c     the radiation field in the cell, and iterates until they have
c     converged.

c     calls: getmatrix, ludcmp, lubksb


      IMPLICIT NONE
      INCLUDE 'amccommon.inc'
      INTEGER id,itrans,s,indx(maxlev+1),iter,miniter,maxiter,stage,numb
      DOUBLE PRECISION ratem(maxlev+1,maxlev+1),newpop(maxlev+1),d
     $     ,up_loc(maxtrans),down_loc(maxtrans),up2_loc(maxtrans2)
     $  ,down2_loc(maxtrans2),opop(maxlev),oopop(maxlev),diff,tol
     $  ,staterr
      CHARACTER*80 collfile
      PARAMETER (tol=1.d-6)
      PARAMETER (miniter=10,maxiter=100)
      LOGICAL debug
      PARAMETER(debug=.false.)

c     id,itrans,iter,s:  counters
c     indx:              keeps track of row permutation (a NumRep thing)
c     miniter,maxiter:   minimum & maximum number of iterations
c     collfile:          transition rates file
c     ratem:             rate matrix
c     newpop:            solution population vector
c     d:   +/-1          odd/even row interchanges (a NumRep thing)
c     up_loc,down_loc,
c     up2_loc,down2_loc: up/down rates at local position, 
c                        either copied from rate arrays or read from disk.
c     opop,oopop:        previous 2 solutions to check for convergence
c     diff:              max difference between subsequent solutions
c     tol:               max allowed difference
c     numb:              no. levels above minpop
c     staterr:           keeps track of non-convergence
c     debug:             turns debugging output on/off


c     First, get the collision rates from array or file

      if (disk) then            ! Get collision rates
        do itrans=1,ntrans
          read(12,*) up_loc(itrans),down_loc(itrans)
        enddo
        if (second) then
          do itrans=1,ntrans2
            read(12,*) up2_loc(itrans),down2_loc(itrans)
          enddo
        endif
      else
        do itrans=1,ntrans
          up_loc(itrans)=up(itrans,id)
          down_loc(itrans)=down(itrans,id)
        enddo
        if (second) then
          do itrans=1,ntrans2
            up2_loc(itrans)=up2(itrans,id)
            down2_loc(itrans)=down2(itrans,id)
          enddo
        endif
      endif


c     Iterate for convergence between radiation field and excitation

      do s=1,nlev               ! reset previous solutions
        opop(s)=0.d0
        oopop(s)=0.d0
      enddo

      do iter=1,maxiter         ! iterate until convergence reached
            
         if (debug) print*,'[debug] calling getjbar, iter= ',iter

        call getjbar(id)        ! get updated jbar (depends on pops!)

        do s=1,nlev             ! set newpop to 0 (null vector)
          newpop(s)=0.d0
        enddo
        newpop(nlev+1)=1.d0     ! nlev+1 is RHS of equation
              
                                ! fill collision rate matrix            
        if (debug) print*,'[debug] calling getmatrix, iter= ',iter
        call getmatrix(id,ratem,up_loc,down_loc,up2_loc,down2_loc)

                                ! solve with LU-decomposition
        if (debug) print*,'[debug] calling ludcmp, iter= ',iter
        call ludcmp(ratem,nlev+1,maxlev+1,indx,d)
        if (debug) print*,'[debug] calling lubksb, iter= ',iter
        call lubksb(ratem,nlev+1,maxlev+1,indx,newpop)
                                ! newpop now contains solution

        numb=0
        diff=0.
        do s=1,nlev             ! get statistics; converged?
          newpop(s)=dmax1(newpop(s),eps)
          oopop(s)=opop(s)
          opop(s)=pops(s,id)
          pops(s,id)=newpop(s)
          if (dmin1(newpop(s),opop(s),oopop(s)).gt.minpop) then
             numb=numb+1
             diff=dmax1(dabs(newpop(s)-opop(s))/newpop(s),
     $            dabs(newpop(s)-oopop(s))/newpop(s),diff)
          endif
        enddo

        if ((iter.gt.miniter).and.(diff.lt.tol)) goto 91

      enddo                     ! iter=1,maxiter

 91   continue

c     If not all converged, save diff in staterr
      if (diff.gt.tol) staterr=diff

      RETURN
      END

