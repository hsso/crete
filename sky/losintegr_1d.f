      SUBROUTINE losintegr(ix,iy)

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

      IMPLICIT NONE
      INCLUDE 'skycommon.inc'
      INTEGER ix,iy,ichan,iline,ju,jl,subi,subj,sup,posn
     $  ,nspline,ispline,naver,iaver,dpos
      DOUBLE PRECISION dx,dy,dxdy,dv,vfac,vfacsub
     $  ,subix,subiy,dtau,tau(maxrx,maxchan),subintens(maxrx,maxchan)
     $  ,subtau(maxrx,maxchan),plus,minn,rpos,dlos,theta,phi,sintheta
     $  ,costheta,sinphi,cosphi,psi,bac,ds,b
     $  ,s1,s2,v1,v2,v,s,rnext,dsplus,dsminn,jnu,alpha,snu
      LOGICAL doline
      DOUBLE PRECISION vfunc
      EXTERNAL vfunc
      COMMON /vproji/ posn
      COMMON /vprojd/ rpos,phi,cosphi,sinphi,dv
      LOGICAL debug
      PARAMETER(debug=.false.)

c     ix,iy: sky-grid position
c     ichan,iline,subi: counters
c     ju,jl: upper, lower level numbers of transition
c     sup: number of super-resolution integrations per sky-pixel
c     posn: current source cell position number
c     ispline,nspline,iaver,naver: chop up los-segment in subunits
c     _    according to local velocity field for proper integration
c     dpos: change in (1d!) cell position 
c     dx,dy,dxdy: physical offsets associated with sky position
c     dv: velocity offset of channel
c     vfac: associated velocity factor of line profile
c     vfacsub: idem, for los-segment
c     subix,subiy: physical offsets of super-resolution position
c     dtau: opacity of los step
c     tau: cumulative tau of total los
c     subintens,subtau: intensity, tau of super-resolution los
c     plus,minn: roots of transfer solution
c     rpos: current radial position in source
c     dlos: length of initial los distance from sky-plane (where los
c     _     enters source)
c     theta, phi: direction of los in source grid
c     sintheta,costheta,sinphi,cosphi: their sine cosine values
c     psi: difference between old and new phi (before/after propagation)
c     bac: helps in solving transfer equation
c     ds: length of los step
c     b: Doppler b linewidth
c     s,v: position, velocity along los step
c     s1,s2,v1,v2: help in solving velocity along los step
c     rnext: next r position
c     dsplus,dsminn: roots of transfer equation
c     jnu, alpha: emission/absorption coefficients
c     snu: source function
c     doline: .true. if line transfer required
c     debug: turns debugging output on/off

      do ichan=1,nchan          ! Clear cumulative opacity
        do iline=1,nrx
          tau(iline,ichan)=0.d0
          intens(iline,ichan,ix)=0.d0
        enddo
      enddo
      if (filter(1).eq.0) then
        doline=.false.
      else
        doline=.true.
      endif
      

      sup=1                     ! Super resolution, if requested
      if (abs(ix-int(xycen)).le.zoom) sup=super
      

      do subi=1,sup             ! Loop over sub-l.o.s.
        do subj=1,sup
        
          subix=dble(ix)+(dble(subi-0.5))/dble(sup)-0.5
          subiy=dble(ix)+(dble(subj-0.5))/dble(sup)-0.5
        
        
          do ichan=1,nchan      ! Initialize this l.o.s.
            do iline=1,nrx
              subtau(iline,ichan)=0.d0
              subintens(iline,ichan)=0.d0
            enddo
          enddo
          

c     Convert pixels and channels in m and m/s w.r.t. center

          dx=(dble(subix)-xycen-0.5)*angres*distance
          dy=(dble(subiy)-xycen-0.5)*angres*distance
          dxdy=dsqrt(dx**2.+dy**2.)


c     Find rpos,zpos,theta,phi of los where it enters front of source

          if (dabs(dxdy).ge.rmax) goto 200 ! los misses source...
          theta=0.
          rpos=rmax*(1.d0-delta)
          plus=+dsqrt(rpos**2.d0-dxdy**2.d0)
          minn=-dsqrt(rpos**2.d0-dxdy**2.d0)
          dlos=dmin1(plus,minn)
          phi=pi-datan2(dxdy,-dlos)
          phi=dmod(phi+pi,2.d0*pi)-pi

          posn=ncell            ! Grid position of start is always ncell in 1D
          
          costheta=1.d0
          sintheta=0.d0
          cosphi=dcos(phi)
          sinphi=dsin(phi)
          

  100     continue              ! Return here for moving thru next cell

c     Find distance to nearest cell edge
 
      if (debug) write(*,*) '[debug] finding closest cell edge'

          if ((dabs(phi).le.pi/2.d0).or.(posn.eq.1)) then
            goto 200
            dpos=+1
            rnext=rb(posn)*(1.d0+delta)
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
          if (dsminn*dsplus.eq.0.) then
            ds=0.d0
          else
            if (dsplus.lt.0.d0) ds=dsminn
            if (dsminn.lt.0.d0) ds=dsplus
            if (dsminn*dsplus.gt.0.d0) ds=dmin1(dsplus,dsminn)
          endif
        
 
c     If there is material in this cell:
c     Find "vfac", the velocity line profile factor
c     Number of splines nspline=los_delta_v/local_line_width
c     Number of averaging steps naver=local_delta_v/local_line_width

          if (nh2(posn).gt.0.d0) then

            do ichan=1,nchan
            
      if (debug) write(*,*) '[debug] finding vfac for channel ',ichan

              if (doline) then            
                dv=(ichan-vcen)*velres
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
              endif
          

              do iline=1,nrx    ! Backward integration of dI/ds

                jnu=dust(iline,posn)*knu(iline,posn)
                alpha=knu(iline,posn)

c     (hpip=hplanck*clight/4./pi/dsqrt(pi))
                if (doline) then
                  ju=lau(iline)
                  jl=lal(iline)
                  jnu=jnu+
     $              vfac/b*hpip*nmol(posn)*
     $              pops(lau(iline),posn)*aeinst(iline)
                  alpha=alpha+
     $              vfac/b*hpip*nmol(posn)*
     $              (pops(lal(iline),posn)*beinstl(iline)-
     $              pops(lau(iline),posn)*beinstu(iline))
                endif
              
                if (dabs(alpha).lt.eps) then !added dabs 08jan02
                  snu=0.d0
                else
                  snu=jnu/alpha/norm(iline)
                endif
              
                dtau=alpha*ds
                if (dtau.lt.negtaulim) then ! Limit negative opacity
                  dtau=negtaulim
                  write(*,'(A,I4)')
     $              'SKY: ***WARNING*** negtaulim exceeded in cell '
     $              ,posn
                endif
            
                subintens(iline,ichan)=subintens(iline,ichan)+
     $            dexp(-subtau(iline,ichan))*(1.d0-dexp(-dtau))*snu
                subtau(iline,ichan)=subtau(iline,ichan)+dtau
            
c     Remember tau toward center of source, when
c     -> at center of sky grid
c     -> at first "super resolution" integration
c     -> heading inward
c     -> correct for back half of inner cell

                if ((ix.eq.int(xycen+1.0)).and.(subi.eq.1)
     $               .and.(dpos.lt.0))
     $            taucen(iline,ichan)=subtau(iline,ichan)-dtau/2.d0

              enddo             ! iline=1,nrx 
            enddo               ! ichan=1,nchan
          endif                 ! nh2(posn).gt.0.


c     Update photon position & direction; first check if escaped

      if (debug) write(*,*) '[debug] updating photon position'

          psi=datan2(ds*sinphi,rpos+ds*cosphi)
          phi=phi-psi
          phi=dmod(phi+pi,2.d0*pi)-pi
          sinphi=dsin(phi)
          cosphi=dcos(phi)
        
          posn=posn+dpos
          if (posn.gt.ncell) goto 200 ! end of integration
          if (dpos.eq.+1) rpos=ra(posn)*(1.d0+delta)
          if (dpos.eq.-1) rpos=rb(posn)*(1.d0-delta)
          goto 100              ! next step
          


  200     do ichan=1,nchan      ! add and subtract cmb
            do iline=1,nrx
              if (subtau(iline,ichan).le.15.d0) then
                subintens(iline,ichan)=subintens(iline,ichan)+
     $            (dexp(-subtau(iline,ichan))-1.d0)*cmb(iline)
              else
                subintens(iline,ichan)=subintens(iline,ichan)-cmb(iline)
              endif
            enddo
          enddo
        

          do ichan=1,nchan      ! Add sub-integration to total
            do iline=1,nrx
              intens(iline,ichan,ix)=intens(iline,ichan,ix)+
     $          subintens(iline,ichan)/sup/sup
              tau(iline,ichan)=tau(iline,ichan)+subtau(iline,ichan)/sup
     $          /sup
            enddo
          enddo
          
        enddo                   ! subj=1,sup
      enddo                     ! subi=1,sup
      

      if ((tcen.gt.0.d0).and.(ix.eq.nint(xycen))) then ! Add central source, if any
        do ichan=1,nchan
          do iline=1,nrx
            intens(iline,ichan,ix)=intens(iline,ichan,ix)+
     $        dexp(-taucen(iline,ichan))*cen(iline)
          enddo
        enddo
      endif

      if (fgdv.gt.0.d0) then      ! Add fg opacity/emission, if defined
        do ichan=1,nchan
          dv=(ichan-vcen)*velres-fvel
          vfac=dexp(-((dv/fgdv)**2.d0))
          do iline=1,nrx
            intens(iline,ichan,ix)=
     $        dexp(-fgtau(iline)*vfac)*intens(iline,ichan,ix)+
     $        (1.d0-dexp(-fgtau(iline)*vfac))*fgtr(iline)
          enddo
        enddo
      endif


      RETURN
      END



