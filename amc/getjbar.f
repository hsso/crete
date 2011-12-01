      SUBROUTINE getjbar(id)

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

c     Determines direction-averaged radiation field (jbar) in cell id
c     according to the local population and the incident radiation field
c     from all other cells.

      IMPLICIT NONE
      INCLUDE 'amccommon.inc'
      
      INTEGER iline,iphot,id
      DOUBLE PRECISION ds,tau,snu,vsum,vfac,jnu,alpha,jnu_dust(maxline),
     $  alpha_dust(maxline),planck
      EXTERNAL planck
      LOGICAL debug
      PARAMETER(debug=.false.)
      
c     iline,iphot,id: counters
c     ds:             length of path of photon through cell
c     jnu,alpha:      local emission and absorption coeffs
c     jnu_dust,alpha_dust: jnu, alpha associated with dust only
c     snu:            local source function
c     tau:            local opacity
c     vfac:           "velocity" factor of photon
c     vsum:           sum of vfac = normalization factor
c     planck:         Planck function, from numrep.f
c     debug: turns debugging output on/off

      vsum=0.d0                   ! Initialization
      do iline=1,nline
        jbar(iline)=0.d0
        jnu_dust(iline)=dust(iline,id)*knu(iline,id)
        alpha_dust(iline)=knu(iline,id)
      enddo

      do iphot=1,nphot(id)      ! loop through photons

         if (debug) print*,'[debug] iphot = ',iphot

        ds=phot(1,iphot)
        vfac=phot(2,iphot)

        do iline=1,nline

          jnu=jnu_dust(iline)+
     $      vfac/doppb(id)*hpip*nmol(id)*
     $      pops(lau(iline),id)*aeinst(iline)

          alpha=alpha_dust(iline)+
     $      vfac/doppb(id)*hpip*nmol(id)*
     $      (pops(lal(iline),id)*beinstl(iline)-
     $      pops(lau(iline),id)*beinstu(iline))

          if (dabs(alpha).lt.eps) then   !added dabs 08jan02
            snu=0.d0
          else
            snu=jnu/alpha/norm(iline)
          endif

          tau=alpha*ds
          if (tau.lt.negtaulim) then ! Limit negative opacity
            tau=negtaulim
          endif
           
c     Add intensity along line segment 
          
          jbar(iline)=jbar(iline)+
     $      vfac*(dexp(-tau)*phot(iline+2,iphot)+(1.d0-dexp(-tau))*snu)
            
        enddo

        vsum=vsum+vfac

      enddo                     ! next photon



      if (vsum.gt.0.d0) then      ! scale by norm and vsum
        do iline=1,nline
          jbar(iline)=norm(iline)*jbar(iline)/vsum ! Normalize and scale
        enddo
      endif

c     ------------------------------------------------------------------

      return
      end


