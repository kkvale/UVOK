C
C Metos3D: A Marine Ecosystem Toolkit for Optimization and Simulation in 3-D
C Copyright (C) 2019  Jaroslaw Piwonski, CAU, jpi@informatik.uni-kiel.de
C
C This program is free software: you can redistribute it and/or modify
C it under the terms of the GNU General Public License as published by
C the Free Software Foundation, either version 3 of the License, or
C (at your option) any later version.
C
C This program is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C GNU General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with this program.  If not, see <http://www.gnu.org/licenses/>.
C

#include "UVOK_TMM_OPTIONS.h"

      subroutine uvok_ini(
     & nzmax, z, drF, DeltaT, S_surf_glob, TR_surf_glob, debugFlag)

      implicit none
#include "size.h"
#include "npzd.h"
#include "calendar.h"
#include "coord.h"
#include "grdvar.h"
#include "pconst.h"
#include "scalar.h"
#include "csbc.h"
#include "mw.h"
#include "atm.h"
#include "cembm.h"
#include "levind.h"
#include "switch.h"

      integer nzmax
      real z(nzmax), drF(nzmax), DeltaT
      real S_surf_glob, TR_surf_glob(nsrc)
      integer debugFlag

      print *, nzmax, nsrc
      print *, z
      print *, drF
      print *, DeltaT
      print *, S_surf_glob
      print *, TR_surf_glob
      print *, debugFlag

      ! set some constants
      pi = 4.0 * atan(1.0)
      radian = 360./(2.*pi)
      rho0 = 1.035
      rho0r = c1/rho0
      socn = 0.03475

      ! levind.h
      kmt(1,1) = nzmax ! this will be set to the local value in uvok_calc

      call sbc_init
      call tracer_init
      call read_namelist

      call exit(1)

      end subroutine uvok_ini
C
C     integer m, n, k, itr
C     real dmsk(imt,jmt)  ! SPK: setmom.F
C
C
CCSPK  overwrite some values
C     dtts=DeltaT  ! make sure this matches the value in control.in
C     daylen=86400.0
C     c2dtts  = dtts
CC
C
C!-----------------------------------------------------------------------
C!     read grid
C!-----------------------------------------------------------------------
C
CC     call grids
C
C!-----------------------------------------------------------------------
C!     read topography
C!-----------------------------------------------------------------------
C
CC      call topog (kmt, kmu, map, xt, yt, zt, xu, yu, zw, imt, jmt, km
CC     &,           sg_bathy)
C
CC      call isleperim (kmt, map, iperm, jperm, iofs, nippts, nisle, imt
CC     &,               jmt, km, mnisle, maxipp, xu, yu, zw)
C
C
CCSPK copy over/set grid variables (based on grids.F)
C     dzt(:) = c0
C     dzw(:) = c0
C     zt(:) = c0
C     zw(:) = c0
C
C     do k=1,km
C       dzt(k) = drF(k)
C       zt(k) = z(k)
C     enddo
C
C     zw(1) = drF(1)
C     do k=2,km
C       zw(k) = zw(k-1) + drF(k)
C     enddo
C
C     do k=1,km-1
C       dzw(k) = zt(k+1) - zt(k)
C     enddo
C     dzw(0) = zt(1)
C     dzw(km) = zw(km) - zt(km)
C
C     do k=1,km
C       c2dzt(k) = c2*dzt(k)
C       dzt2r(k) = c1/c2dzt(k)
C     enddo
C
C     dzwr(km)  = c1/dzw(km)
C     dzw2r(km) = p5/dzw(km)
C
C     do k=1,km
C       dzwr(k-1)    = c1/dzw(k-1)
C       dzw2r(k-1)   = p5/dzw(k-1)
C     enddo
C
C     do k=1,km
C       dztr(k)  = c1/dzt(k)
C     enddo
C
CCKK from S/R setmom.F
C#if defined O_npzd
C!     taken from setmom.F
C!     convert units of NPZD parameters to MOM units
C     redctn = redctn*1.e-3
C     redotn = redotn*1.e-3
C     redotp = redotn/redptn
C     redctp = redctn/redptn
C     redntp = 1./redptn
C     k1p   = k1n*redptn
C     kw = kw*1.e-2
C     kc = kc*1.e-2
C     ki = ki*1.e-2
C     wd0 = wd0*1.e2
C     alpha = alpha/daylen
C     abio = abio/daylen
C     nup = nup/daylen
C     nupt0 = nupt0/daylen
C     gbio = gbio/daylen
C     nuz = nuz/daylen
C     nud0 = nud0/daylen
C#endif
C
C     do k=1,km
C!     linear increase wd0-200m with depth
C        wd(k) = (wd0+6.0e-2*zt(k))/daylen/dzt(k) ! [s-1]
C        rkwz(k) = 1./(kw*dzt(k))
C     enddo
C     ztt(1)=0.0
C     do k=1,km-1
C        ztt(k+1)=(-1)*zw(k)
C     enddo
C     ztt(km)=(-1)*zw(km-1)
C# if defined O_carbon || defined O_npzd_alk
C
C!---------------------------------------------------------------------
C!     calculate variables used in calcite remineralization
C!---------------------------------------------------------------------
C
C     rcak(1) = -(exp(-zw(1)/dcaco3)-1.0)/dzt(1)
C     rcab(1) = 1./dzt(1)
C     do k=2,km
C       rcak(k) = -(exp(-zw(k)/dcaco3))/dzt(k)
C    &          + (exp(-zw(k-1)/dcaco3))/dzt(k)
C       rcab(k) = (exp(-zw(k-1)/dcaco3))/dzt(k)
C     enddo
C# endif
C!-----------------------------------------------------------------------
C!     common setup
C!-----------------------------------------------------------------------
C
CC      call setcom (1, imt, 1, jmt)
C
C!-----------------------------------------------------------------------
C!     ocean setup
C!-----------------------------------------------------------------------
C
CC      call setmom (1, imt, 1, jmt)
C
CCSPK we set these to all point to the same place
C     tau   = 0
C     taum1 = 0
C     taup1 = 0
CC
C
CCSPK copy initial surface values
C     gaost(:) = 0.
C     if (isalt .ne. 0 .and. isss .ne. 0) then
C       gaost(isalt) = S_surf_glob
CC        write(*,*)'Salt: ',gaost(isalt)
C     endif
C     itr = 0
C     do n=3,nt
C       if (itrc(n) .ne. 0) then
C        itr = itr + 1
C	 gaost(n) = TR_surf_glob(itr)
CC		 write(*,*)mapt(n),gaost(n)
C       endif
C     enddo
CC     check
C     if (itr.ne.nsrc) then
C       print*,'ERROR in UVOK_INI: number of tracers incorrect!'
C       stop
C     endif
CC
C
C     if (debugFlag.eq.1) then
CCSPK taken from setmom.F, S/R rowi
C!-----------------------------------------------------------------------
C!     find inital average surface references
C!-----------------------------------------------------------------------
C     print*, " "
C     print*, "inital average surface references: "
C     dmsk(:,:) = 1.
C     where (kmt(:,:) .eq. 0) dmsk(:,:) = 0.
CC     gaost(:) = 0.
C     if (isalt .ne. 0 .and. isss .ne. 0) then
C       call areaavg (sbc(1,1,isss), dmsk, gaost(isalt))
C       gaost(isalt) = gaost(isalt) + 0.035
C       socn = gaost(isalt)
C       print*, "global average sea surface salinity (psu) = "
C    &,   gaost(isalt)*1000.
C     endif
C     if (idic .ne. 0 .and. issdic .ne. 0) then
C       call areaavg (sbc(1,1,issdic), dmsk, gaost(idic))
C       print*, "global average sea surface dic (mol m-3) = "
C    &,   gaost(idic)
C     endif
CC      if (idic13 .ne. 0 .and. issdic13 .ne. 0) then
CC        call areaavg (sbc(1,1,issdic13), dmsk, gaost(idic13))
CC        print*, "global average sea surface dic 13 (mol m-3) = "
CC     &,   gaost(idic13)
CC      endif
C     if (io2 .ne. 0 .and. isso2 .ne. 0) then
C       call areaavg (sbc(1,1,isso2), dmsk, gaost(io2))
C       print*, "global average sea surface oxygen (mol m-3) = "
C    &,   gaost(io2)
C     endif
C     if (ialk .ne. 0 .and. issalk .ne. 0) then
C       call areaavg (sbc(1,1,issalk), dmsk, gaost(ialk))
C       print*, "global average sea surface alkalinity (mol m-3) = "
C    &,   gaost(ialk)
C     endif
C     if (ipo4 .ne. 0 .and. isspo4 .ne. 0) then
C       call areaavg (sbc(1,1,isspo4), dmsk, gaost(ipo4))
C       print*, "global average sea surface phosphate (mol m-3) = "
C    &,   gaost(ipo4)*0.001
C     endif
CC      if (idop .ne. 0 .and. issdop .ne. 0) then
CC        call areaavg (sbc(1,1,issdop), dmsk, gaost(idop))
CC        print*, "global average sea surface DOP (mol m-3) = "
CC     &,   gaost(idop)*0.001
CC      endif
C# if !defined O_npzd_no_vflux
C     if (iphyt .ne. 0 .and. issphyt .ne. 0) then
C       call areaavg (sbc(1,1,issphyt), dmsk, gaost(iphyt))
C       print*, "global average sea surface phytoplankton (mol m-3) = "
C    &,   gaost(iphyt)*0.001
C     endif
C     if (izoop .ne. 0 .and. isszoop .ne. 0) then
C       call areaavg (sbc(1,1,isszoop), dmsk, gaost(izoop))
C       print*, "global average sea surface zooplankton (mol m-3) = "
C    &,   gaost(izoop)*0.001
C     endif
C     if (idetr .ne. 0 .and. issdetr .ne. 0) then
C       call areaavg (sbc(1,1,issdetr), dmsk, gaost(idetr))
C       print*, "global average sea surface detritus (mol m-3) = "
C    &,   gaost(idetr)*0.001
C     endif
C# endif
C     if (ino3 .ne. 0 .and. issno3 .ne. 0) then
C       call areaavg (sbc(1,1,issno3), dmsk, gaost(ino3))
C       print*, "global average sea surface nitrate (mol m-3) = "
C    &,   gaost(ino3)*0.001
C     endif
CC      if (idon .ne. 0 .and. issdon .ne. 0) then
CC        call areaavg (sbc(1,1,issdon), dmsk, gaost(idon))
CC        print*, "global average sea surface DON (mol m-3) = "
CC     &,   gaost(idon)*0.001
CC      endif
C# if !defined O_npzd_no_vflux
C     if (idiaz .ne. 0 .and. issdiaz .ne. 0) then
C       call areaavg (sbc(1,1,issdiaz), dmsk, gaost(idiaz))
C       print*, "global average sea surface diazotrophs (mol m-3) = "
C    &,   gaost(idiaz)*0.001
C     endif
C# endif
CC      if (idin15 .ne. 0 .and. issdin15 .ne. 0) then
CC        call areaavg (sbc(1,1,issdin15), dmsk, gaost(idin15))
CC        print*, "global average sea surface nitrate 15 (mol m-3) = "
CC     &,   gaost(idin15)*0.001
CC      endif
CC      if (idon15 .ne. 0 .and. issdon15 .ne. 0) then
CC        call areaavg (sbc(1,1,issdon15), dmsk, gaost(idon15))
CC        print*, "global average sea surface DON15 (mol m-3) = "
CC     &,   gaost(idon15)*0.001
CC      endif
C# if !defined O_npzd_no_vflux
CC      if (iphytn15 .ne. 0 .and. issphytn15 .ne. 0) then
CC        call areaavg (sbc(1,1,issphytn15), dmsk, gaost(iphytn15))
CC        print*, "global average sea surface phytoplankton n15 "
CC     &,         "(mol m-3) = ", gaost(iphytn15)*0.001
CC      endif
CC      if (izoopn15 .ne. 0 .and. isszoopn15 .ne. 0) then
CC        call areaavg (sbc(1,1,isszoopn15), dmsk, gaost(izoopn15))
CC        print*, "global average sea surface zooplankton n15 (mol m-3)
CC     &           = ", gaost(izoopn15)*0.001
CC      endif
CC      if (idetrn15 .ne. 0 .and. issdetrn15 .ne. 0) then
CC        call areaavg (sbc(1,1,issdetrn15), dmsk, gaost(idetrn15))
CC        print*, "global average sea surface detritus n15 (mol m-3)
CC     &           = ", gaost(idetrn15)*0.001
CC      endif
CC      if (idiazn15 .ne. 0 .and. issdiazn15 .ne. 0) then
CC        call areaavg (sbc(1,1,issdiazn15), dmsk, gaost(idiazn15))
CC        print*, "global average sea surface diazotrophs n15 (mol m-3)
CC     &           = ", gaost(idiazn15)*0.001
CC      endif
C# endif
CC      if (idoc13 .ne. 0 .and. issdoc13 .ne. 0) then
CC        call areaavg (sbc(1,1,issdoc13), dmsk, gaost(idoc13))
CC        print*, "global average sea surface DOC13"
CC     &,         " (mol m-3) = ", gaost(idoc13)*0.001
CC      endif
C# if !defined O_npzd_no_vflux
CC      if (iphytc13 .ne. 0 .and. issphytc13 .ne. 0) then
CC        call areaavg (sbc(1,1,issphytc13), dmsk, gaost(iphytc13))
CC        print*, "global average sea surface phytoplankton C13"
CC     &,         " (mol m-3) = ", gaost(iphytc13)*0.001
CC      endif
CC      if (izoopc13 .ne. 0 .and. isszoopc13 .ne. 0) then
CC        call areaavg (sbc(1,1,isszoopc13), dmsk, gaost(izoopc13))
CC        print*, "global average sea surface zooplankton C13"
CC     &,         " (mol m-3) = ", gaost(izoopc13)*0.001
CC      endif
CC      if (idetrc13 .ne. 0 .and. issdetrc13 .ne. 0) then
CC        call areaavg (sbc(1,1,issdetrc13), dmsk, gaost(idetrc13))
CC        print*, "global average sea surface detritus c13"
CC     &,         " (mol m-3) = ", gaost(idetrc13)*0.001
CC      endif
CC      if (idiazc13 .ne. 0 .and. issdiazc13 .ne. 0) then
CC        call areaavg (sbc(1,1,issdiazc13), dmsk, gaost(idiazc13))
CC        print*, "global average sea surface diazotrophs c13"
CC     &,         " (mol m-3) = ", gaost(idiazc13)*0.001
CC      endif
C# endif
CC      if (idfe .ne. 0 .and. issdfe .ne. 0) then
CC        call areaavg (sbc(1,1,issdfe), dmsk, gaost(idfe))
CC        print*, "global average sea surface iron (mol m-3) = "
CC     &,   gaost(idfe)*0.001
CC      endif
C# if !defined O_npzd_no_vflux
CC      if (idetrfe .ne. 0 .and. issdetrfe .ne. 0) then
CC        call areaavg (sbc(1,1,issdetrfe), dmsk, gaost(idetrfe))
CC        print*, "global average sea surface particulate iron (mol m-3)="
CC     &,   gaost(idetrfe)*0.001
CC      endif
C# endif
C     if (ic14 .ne. 0 .and. issc14 .ne. 0) then
C       call areaavg (sbc(1,1,issc14), dmsk, gaost(ic14))
C       print*, "global average sea surface carbon 14 (mol m-3) = "
C    &,   gaost(ic14)
C     endif
C     if (icfc11 .ne. 0 .and. isscfc11 .ne. 0) then
C       call areaavg (sbc(1,1,isscfc11), dmsk, gaost(icfc11))
C       print*, "global average sea surface cfc 11 (mol m-3) = "
C    &,   gaost(icfc11)
C     endif
C     if (icfc12 .ne. 0 .and. isscfc12 .ne. 0) then
C       call areaavg (sbc(1,1,isscfc12), dmsk, gaost(icfc12))
C       print*, "global average sea surface cfc 12 (mol m-3) = "
C    &,   gaost(icfc12)
C     endif
C     print*, " "
CCSPK
C     endif ! debugFlag
C
CCSPK from setmom.F, S/R setmom
CC
C
C!-----------------------------------------------------------------------
C!     atmosphere setup
C!-----------------------------------------------------------------------
CC      call setembm (1, imt, 1, jmt)
C
CCSPK based on setembm.F
C     tmsk(:,:) = 1.0 ! land/ocean mask
CC
C
CCSPK this is generally set in S/R loadmw called from S/R mom
C     tmask(:,:,:) = c1
CC
C
CCSPK set some flags so that UVic doesn't start accumulating time averages
C     euler2 = .false.
C     timavgperts = .false.
C     tsiperts = .false.
C
C# if defined O_npzd_fe_limitation
C     do k=km,1,-1
C        if (zw(k) .gt. 24000) kmfe = k-1
C     enddo
C# endif
C
CCSPK  Check
C     if (dtts.ne.DeltaT) then
C        print*,'ERROR: UVOK_INI: dtts does not match DeltaT'
C        stop
C     endif
C
C     if (debugFlag.eq.1) then
C       print*,'UVOK time steps: dtts=',dtts,' c2dtts=',c2dtts,
C    &         ' dtnpzd=',dtnpzd
C#ifdef O_npzd_fe_limitation
C       print*,'Iron limitation depth index kmfe=',kmfe
C#endif
C     endif
C
C     END
C#endif
