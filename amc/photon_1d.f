      SUBROUTINE photon(id)
      
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

c     For nphot(id) photons, this routines chooses a random start
c     position and direction, propagates it to the nearest cell edge,
c     and  calculates the distance ds and "velocity factor" vfac. These
c     are stored in the array "phot", together with the incident
c     intensity that that photons represents: calculated by stepping
c     from cell edge to cell edge, and solving radiative transfer (and
c     finally adding the cmb field).

c     Together with vfunc_?d.f, this is the only routine where the
c     dimensionality of the source is apparent.


      IMPLICIT NONE
      INCLUDE 'amccommon.inc'
      INTEGER iphot,id,iline,ispline,nspline,naver,iaver,dpos,posn
      DOUBLE PRECISION rpos,phi,deltav,cosphi,sinphi,dummy,ran1
      DOUBLE PRECISION bac,ds,b,dsplus,dsminn,v1,v2,vfunc,s1,s2
     $     ,vfac,s,v,vfacsub,jnu,alpha
      EXTERNAL ran1,vfunc
      DOUBLE PRECISION snu,tau(maxline),dtau,psi,phif,rnext
      LOGICAL firststep
      DOUBLE PRECISION vel(3)
      COMMON /vproji/ posn
      COMMON /vprojd/ rpos,phi,cosphi,sinphi,deltav
      LOGICAL debug
      PARAMETER(debug=.false.)

c     iphot:     no. of current photon
c     id:        no. of current "origin" cell
c     iline:     no. of current radiative transition
c     nspline, ispline: total + current no. of sub-steps through cell.
c     naver,iaver:      total + current no. of velocity steps
c     dpos:      +/-1 when photon moves out/in
c     posn:      running position of photon (cell no.)
c     rpos, phi, deltav:  location and velocity of photon
c     cosphi, sinphi:     cosine and sine of phi
c     dummy, seed, ran1:   stuff for the random number generator
c     ds:                 distance to cell edge
c     bac,dsplus,dsmin:   discriminant + roots of eqn. for ds
c     b:                  local Doppler width
c     v1,v2               local gas velocity at begin/end of photon path
c     s1,s2               local gas velocity at begin/end of substeps
c     vfunc               determines gas velocities
c     vfac:               exponential velocity factor
c     s,v:                running values of actual & projected velocity
c     vfacsub:            -- dummy for average vfac
c     snu, tau, dtau, i0: source function, optical depth, intensity
c     jnu,alpha:          combined emission and absorption coeffs lines&dust
c     psi, phif, rnext:   angle and position in new cell
c     firststep:          .true. if this step is from origin to first edge.
c     debug: turns debugging output on/off


      do iphot=1,nphot(id)      ! Loop through photons


        posn=id                 ! Initialize posn,intensity, tau
        firststep=.true.
        do iline=1,nline
          phot(iline,iphot)=0.d0
          tau(iline)=0.d0
        enddo
        phot(nline+1,iphot)=0.d0
        phot(nline+2,iphot)=0.d0        

c     Assign random position within cell id, direction and velocity, and
c     determine distance ds to edge of cell. Choose position so that the
c     cell volume is equally sampled in volume, the direction so that
c     all solid angles are equally sampled, and the velocity offset
c     homogeneously distributed over +/-2.15 * local Doppler b from
c     local velocity.

        dummy=ran1(seed)
        if (ra(id).gt.0.d0) then
          rpos=ra(id)*(1.d0+dummy*((rb(id)/ra(id))**3.d0-1.d0))**(1.d0/3
     $      .d0)
        else
          rpos=rb(id)*dummy**(1.d0/3.d0)
        endif

        dummy=2.d0*ran1(seed)-1.d0
        phi=dasin(dummy)+pi/2.d0

        if (debug) print*,'[debug] calling velo, iphot = ',iphot
        call velo(id,rpos,vel)
        deltav=(ran1(seed)-0.5d0)*4.3d0*doppb(id)+dcos(phi)*vel(1)


c     Propagate to edge of cloud by moving from cell edge to cell edge.
c     After the first step (i.e., after moving to the edge of cell id),
c     store ds and vfac in phot(1,iphot) and phot(2,iphot). After
c     leaving the source, add CMB and store intensities in 
c     phot(iline+2,iphot)


  100   cosphi=dcos(phi)        ! Return here for moving thru next cell
        sinphi=dsin(phi)


c     Find distance to nearest cell edge
            
        if ((dabs(phi).le.(pi/2.d0)).or.(posn.eq.1)) then
          dpos=+1
          rnext=rb(posn)*(1.+delta)
          bac=4.d0*((rpos*cosphi)**2.d0-rpos**2.d0+rnext**2.d0)
        else
          dpos=-1
          rnext=ra(posn)*(1.d0-delta)
          bac=4.d0*((rpos*cosphi)**2.d0-rpos**2.d0+rnext**2.d0)
          if (bac.lt.0.d0) then
            dpos=+1
            rnext=rb(posn)*(1.d0+delta)
            bac=4.d0*((rpos*cosphi)**2.d0-rpos**2.d0+rnext**2.d0)
          endif
        endif
        dsplus=-0.5d0*(2.d0*rpos*cosphi+dsqrt(bac))
        dsminn=-0.5d0*(2.d0*rpos*cosphi-dsqrt(bac))
        if (dsminn*dsplus.eq.0.d0) then
          ds=0.d0
        else
          if (dsplus.lt.0.d0) ds=dsminn
          if (dsminn.lt.0.d0) ds=dsplus
          if (dsminn*dsplus.gt.0.d0) ds=dmin1(dsplus,dsminn)
        endif


c     Find "vfac", the velocity line profile factor
c     Number of splines nspline=los_delta_v/local_line_width
c     Number of averaging steps naver=local_delta_v/local_line_width
         
        if (nmol(posn).gt.eps) then

          b=doppb(posn)
          v1=vfunc(0.d0)
          v2=vfunc(ds)
          nspline=max(1,int((dabs(v1-v2))/b))
          vfac=0.d0
          do ispline=1,nspline
            s1=ds*dble(ispline-1)/dble(nspline)
            s2=ds*dble(ispline)/dble(nspline)
            v1=vfunc(s1)
            v2=vfunc(s2)
            naver=max(1,int(dabs(v1-v2)/b))
            do iaver=1,naver
              s=s1+(s2-s1)*(dble(iaver)-0.5d0)/dble(naver)
              v=vfunc(s)
              call gauline(v,b,vfacsub)
              vfac=vfac+vfacsub/dble(naver)
            enddo
          enddo
          vfac=vfac/dble(nspline)


          do iline=1,nline      ! backward integration of dI/ds
              
            jnu=dust(iline,posn)*knu(iline,posn)+
     $        vfac*hpip/b*nmol(posn)*
     $        pops(lau(iline),posn)*aeinst(iline)

            alpha=knu(iline,posn)+
     $        vfac*hpip/b*nmol(posn)*
     $        (pops(lal(iline),posn)*beinstl(iline)-
     $        pops(lau(iline),posn)*beinstu(iline))
            
            if (dabs(alpha).lt.eps) then  !added dabs 08jan02
              snu=0.d0
            else
              snu=jnu/alpha/norm(iline)
            endif

            dtau=alpha*ds
            if (dtau.lt.negtaulim) then ! Limit negative opacity
              dtau=negtaulim
ccc              write(*,'(A,I4)')
ccc     $          'AMC: ***WARNING*** negtaulim exceeded in cell ',posn
            endif
            
            if (firststep.neqv..true.) then   !new 08jan02
               phot(iline+2,iphot)=phot(iline+2,iphot)+
     $              dexp(-tau(iline))*(1.d0-dexp(-dtau))*snu
               tau(iline)=tau(iline)+dtau ! update cumulative tau
               if (tau(iline).lt.negtaulim) then ! Limit negative opacity
                  tau(iline)=negtaulim
ccc   write(*,'(A,I4)')
ccc   $          'AMC: ***WARNING*** negtaulim exceeded in cell ',posn
               endif
            endif   !new 08jan02

          enddo                 ! iline=1,nline


          if (firststep) then   ! memorize ds and vfac after first step
            phot(1,iphot)=ds
            phot(2,iphot)=vfac
            firststep=.false.
          endif

        endif                   ! if (nmol(posn).gt.eps)

c     Update photon position, direction; check if escaped

        posn=posn+dpos
        if (posn.gt.ncell) goto 200 ! edge of cloud reached

        psi=datan2(ds*sinphi,rpos+ds*cosphi)
        phif=phi-psi
        phi=dmod(phif,pi)       ! phi = [0,pi]
        rpos=rnext
        
        
        goto 100                ! next step


c     Finally, add cmb to memorized i0 incident on cell id
        
  200   if (tcmb.gt.0.) then
          do iline=1,nline
            phot(iline+2,iphot)=phot(iline+2,iphot)+
     $        dexp(-tau(iline))*cmb(iline)
          enddo
        endif

      enddo                     ! next photon

c     ------------------------------------------------------------------

      return
      end
