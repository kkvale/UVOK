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
C   This routine is originally part of the official UVOK model:
C       https://raw.githubusercontent.com/samarkhatiwala/tmm/master/models/current/uvok1.0/src/uvok_ini.F
C

#include "UVOK_TMM_OPTIONS.h"

      subroutine uvok_ini(
     &      nzmax,
     &      z,
     &      drF,
     &      DeltaT,
     &      S_surf_glob,
     &      TR_surf_glob,
     &      debugFlag)

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

        integer k, n, itr

        ! set some constants
        pi = 4.0 * atan(1.0)
        radian = 360./(2.*pi)
        rho0 = 1.035
        rho0r = c1/rho0
        socn = 0.03475

        kmt(1,1) = nzmax      ! this will be set to the local value in uvok_calc

        call sbc_init
        call tracer_init
        call read_namelist

        ! SPK overwrite some values
        dtts = DeltaT         ! make sure this matches the value in control.in
        daylen = 86400.0
        c2dtts = dtts

        ! SPK copy over/set grid variables (based on grids.F)
        dzt(:) = c0
        dzw(:) = c0
        zt(:) = c0
        zw(:) = c0

!        do k=1,km
        do k=1,nzmax
            dzt(k) = drF(k)
            zt(k) = z(k)
        enddo

        zw(1) = drF(1)
!        do k=2,km
        do k=2,nzmax
            zw(k) = zw(k-1) + drF(k)
        enddo

!        do k=1,km-1
        do k=1,nzmax-1
            dzw(k) = zt(k+1) - zt(k)
        enddo
        dzw(0) = zt(1)
!        dzw(km) = zw(km) - zt(km)
        dzw(nzmax) = zw(nzmax) - zt(nzmax)

!        do k=1,km
        do k=1,nzmax
            c2dzt(k) = c2*dzt(k)
            dzt2r(k) = c1/c2dzt(k)
        enddo

!        dzwr(km)  = c1/dzw(km)
!        dzw2r(km) = p5/dzw(km)
        dzwr(nzmax)  = c1/dzw(nzmax)
        dzw2r(nzmax) = p5/dzw(nzmax)

!        do k=1,km
        do k=1,nzmax
            dzwr(k-1)    = c1/dzw(k-1)
            dzw2r(k-1)   = p5/dzw(k-1)
        enddo

!        do k=1,km
        do k=1,nzmax
            dztr(k)  = c1/dzt(k)
        enddo

        ! taken from setmom.F
        ! convert units of NPZD parameters to MOM units
        redctn = redctn*1.e-3
        redotn = redotn*1.e-3
        redotp = redotn/redptn
        redctp = redctn/redptn
        redntp = 1./redptn
        k1p   = k1n*redptn
        kw = kw*1.e-2
        kc = kc*1.e-2
        ki = ki*1.e-2
        wd0 = wd0*1.e2
        alpha = alpha/daylen
        abio = abio/daylen
        nup = nup/daylen
        nupt0 = nupt0/daylen
        gbio = gbio/daylen
        nuz = nuz/daylen
        nud0 = nud0/daylen

!        do k=1,km
        do k=1,nzmax
            ! linear increase wd0-200m with depth
            wd(k) = (wd0+6.0e-2*zt(k))/daylen/dzt(k) ! [s-1]
            rkwz(k) = 1./(kw*dzt(k))
        enddo
        ztt(1)=0.0
!        do k=1,km-1
        do k=1,nzmax-1
            ztt(k+1)=(-1)*zw(k)
        enddo
!        ztt(km)=(-1)*zw(km-1)
        ztt(nzmax)=(-1)*zw(nzmax-1)

        ! calculate variables used in calcite remineralization
        rcak(1) = -(exp(-zw(1)/dcaco3)-1.0)/dzt(1)
        rcab(1) = 1./dzt(1)
!        do k=2,km
        do k=2,nzmax
            rcak(k) = -(exp(-zw(k)/dcaco3))/dzt(k)
     &         + (exp(-zw(k-1)/dcaco3))/dzt(k)
            rcab(k) = (exp(-zw(k-1)/dcaco3))/dzt(k)
        enddo

        ! SPK we set these to all point to the same place
        tau   = 0
        taum1 = 0
        taup1 = 0

        ! SPK copy initial surface values
        ! csbc.h
        ! gaost = global average ocean surface tracer
        gaost(:) = 0.
        if (isalt .ne. 0 .and. isss .ne. 0) then
            gaost(isalt) = S_surf_glob
            socn = gaost(isalt) + 0.035
        endif
        itr = 0
        do n=3,nt
            if (itrc(n) .ne. 0) then
                itr = itr + 1
                gaost(n) = TR_surf_glob(itr)
            endif
        enddo

        ! SPK from setmom.F, S/R setmom
        ! atmosphere setup
        ! SPK based on setembm.F
        tmsk(:,:) = 1.0                       ! land/ocean mask
        ! SPK this is generally set in S/R loadmw called from S/R mom
        tmask(:,:,:) = c1
        ! SPK set some flags so that UVic doesn't start accumulating time averages
        euler2 = .false.
        timavgperts = .false.
        tsiperts = .false.

!        do k=km,1,-1
!            if (zw(k) .gt. 24000) kmfe = k-1
!        enddo

      end subroutine uvok_ini


