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
      DOUBLE PRECISION colli2(maxlev,maxlev),ctot2(maxlev), 
     $ g_ir(maxlev,maxlev),ctot_ir(maxlev)
      
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
          g_ir(s,t)=0.d0
        enddo
      enddo
      do t=1,nlev
        do s=1,nlev
          colli2(s,t)=0.d0
        enddo
      enddo

! IR pumping ortho-water
!      g_ir = reshape(
!     $ (/0., 1.654e-5,2.464e-5, 8.486e-5, 1.471e-4, 1.359e-5, 2.905e-5,
!     $ 1.423e-5, 0., 1.882e-4, 1.696e-5, 1.118e-5, 9.143e-5, 1.294e-5,
!     $ 1.323e-5, 1.154e-4, 0., 1.315e-5, 1.492e-5, 1.341e-4, 1.421e-5,
!     $ 5.619e-5, 1.602e-5, 1.531e-5, 0., 2.846e-5, 1.694e-5, 1.415e-4,
!     $ 6.620e-5, 5.189e-6, 1.113e-5, 1.982e-5, 0., 1.013e-5, 5.302e-5,
!     $ 4.428e-6, 4.090e-5, 1.006e-4, 5.843e-6, 7.287e-6, 0., 9.177e-6,
!     $ 1.448e-5, 7.867e-6, 9.250e-6, 1.038e-4, 5.593e-5, 1.166e-5, 0./),
!     $ (/ maxlev,maxlev /) )

! IR pumping para-water
!      g_ir = reshape(
!     $ (/0., 2.1944012E-05, 2.4522215E-04, 1.5621181E-05, 1.0322739E-04,
!     $   1.4237256E-05, 6.7647539E-07,
!     $ 6.1202104E-06, 0., 1.7902237E-05, 1.6408926E-04, 1.5259558E-05,
!     $   1.0526829E-04, 1.1883916E-05,
!     $ 5.0827490E-05, 1.3236418E-05, 0., 2.2540904E-05, 5.0350245E-05,
!     $   2.0026466E-05, 8.6539898E-05,
!     $ 3.0782348E-06, 1.0119543E-04, 1.9672434E-05, 0., 1.4988847E-05,
!     $   1.2663435E-04, 1.7797252E-05,
!     $ 2.3122320E-05, 1.3558578E-05, 5.3738728E-05, 2.4144989E-05, 0.,
!     $   9.3817425E-06, 1.4817977E-04,
!     $ 1.9010073E-06, 4.7502926E-05, 1.1559701E-05, 9.3329596E-05, 
!     $   6.0966895E-06, 0., 8.9717369E-06,
!     $ 3.4442709E-07, 7.2605194E-06, 6.8272238E-05, 1.5878764E-05,
!     $   1.0893742E-04, 7.1137065E-06, 0./),
!     $ (/ maxlev,maxlev /) )

! scale by rh**2 (in AU) input from comet_mdl
!       g_ir = g_ir/rh**2
      
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
      do s=1,nlev
        ctot_ir(s)=0.d0
        do t=1,nlev
          ctot_ir(s)=ctot_ir(s) + g_ir(t,s) ! IR pumping
        enddo
      enddo

c     Fill the rate matrix with the collisional contributions.
      
      do s=1,nlev
        mtrx(s,s)=mtrx(s,s)+nh2(id)*ctot(s)+ne(id)*ctot2(s)+ctot_ir(s)
        do t=1,nlev
          if (s.ne.t) 
     $      mtrx(s,t)=mtrx(s,t)-nh2(id)*colli(t,s)-ne(id)*colli2(t,s)
     $      -g_ir(s,t) ! IR pumping
       enddo
        mtrx(nlev+1,s)=1.d0
        mtrx(s,nlev+1)=0.d0
      enddo    

      RETURN
      END
