      SUBROUTINE getmatrix(id,mtrx,up_loc,down_loc,up2_loc,down2_loc)

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
      
c     getmatrix: builds matrix of (local) transition rates + equation
c     of conservation in row nlev+1
      
      IMPLICIT NONE
      INCLUDE 'amccommon.inc'      
      INTEGER id,s,t,k,l
      DOUBLE PRECISION up_loc(maxtrans),down_loc(maxtrans)
     $     ,up2_loc(maxtrans2),down2_loc(maxtrans2),mtrx(maxlev+1,maxlev
     $     +1) ,colli(maxlev,maxlev),ctot(maxlev)
      DOUBLE PRECISION colli2(maxlev,maxlev),ctot2(maxlev)
      
c     id:   current grid position
c     s,t:     counters
c     k,l:     upper, lower level of current transition
c     up/down_loc: local up/downward collision rates (1st/2nd partner).
c     mtrx:    rate matrix to be determined
c     colli:   off-axis elements mtrx
c     ctot:    axis elements mtrx

c     ------------------------------------------------------------
      
      do t=1,nlev+1             ! initialization
        do s=1,nlev+1
          mtrx(s,t)=0.d0
        enddo
      enddo
      do t=1,nlev
        do s=1,nlev               
          colli(s,t)=0.d0
        enddo
      enddo
      do t=1,nlev
        do s=1,nlev
          colli2(s,t)=0.d0
        enddo
      enddo
      
      do t=1,nline              ! radiative transitions (beinstl = Blu)
        k=lau(t)
        l=lal(t)
        mtrx(k,k)=mtrx(k,k)+beinstu(t)*jbar(t)+aeinst(t)
        mtrx(l,l)=mtrx(l,l)+beinstl(t)*jbar(t)
        mtrx(k,l)=mtrx(k,l)-beinstl(t)*jbar(t)
        mtrx(l,k)=mtrx(l,k)-beinstu(t)*jbar(t)-aeinst(t)
      enddo

      do t=1,ntrans             ! create collision rate matrix
        colli(lcu(t),lcl(t))=down_loc(t)
        colli(lcl(t),lcu(t))=up_loc(t)
      enddo
      do t=1,ntrans2
        colli2(lcu2(t),lcl2(t))=down2_loc(t)
        colli2(lcl2(t),lcu2(t))=up2_loc(t)
      enddo
            
c     Sum by rows and add the conservation equation to the
c     left-hand side (at position nlev+1).
      
      do s=1,nlev
        ctot(s)=0.d0
        do t=1,nlev
          ctot(s)=ctot(s)+colli(s,t)
        enddo
      enddo
      do s=1,nlev
        ctot2(s)=0.d0
        do t=1,nlev
          ctot2(s)=ctot2(s)+colli2(s,t)
        enddo
      enddo

c     Fill the rate matrix with the collisional contributions.
      
      do s=1,nlev
        mtrx(s,s)=mtrx(s,s)+nh2(id)*ctot(s)+ne(id)*ctot2(s)
        do t=1,nlev
          if (s.ne.t) 
     $      mtrx(s,t)=mtrx(s,t)-nh2(id)*colli(t,s)-ne(id)*colli2(t,s)
       enddo
        mtrx(nlev+1,s)=1.d0
        mtrx(s,nlev+1)=0.d0
      enddo    

      RETURN
      END
