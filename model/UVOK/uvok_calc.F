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
C       https://raw.githubusercontent.com/samarkhatiwala/tmm/master/models/current/uvok1.0/src/uvok_calc.F
C

#include "UVOK_TMM_OPTIONS.h"

      subroutine uvok_calc(
     &      kmt_loc, tlat_loc, day_loc, relyr_loc,
     &      TEMP, SALT, TR_surf_glob,dz_loc,z_loc,
     &      winds_loc,
     &      fe_dissolved_loc,
     &      swr_loc,
     &      aice_loc, hice_loc, hsno_loc,
     &      emp_loc, emp_glob,
     &      gasexfluxloc, totfluxloc,
     &      debugFlag)

        implicit none

#include "size.h"

        integer kmt_loc
        real tlat_loc, day_loc, relyr_loc
        real TEMP(kmt_loc), SALT(kmt_loc)
        real TR_surf_glob(nsrc)
        real dz_loc(kmt_loc)
        real z_loc(kmt_loc)
!        real z(km)
        real winds_loc
        real fe_dissolved_loc(kmt_loc)            ! only the first min(kmfe,kmt_loc) elements are used
        real swr_loc
        real aice_loc, hice_loc, hsno_loc
        real emp_loc, emp_glob
        real gasexfluxloc, totfluxloc
        integer debugFlag

        integer i, j, js, je, is, ie, joff, itr, k, n

#include "param.h"
#include "pconst.h"
#include "stdunits.h"
#include "coord.h"
#include "csbc.h"
#include "grdvar.h"
#include "levind.h"
#include "mw.h"
#include "scalar.h"
#include "tmngr.h"
#include "npzd.h"
#include "calendar.h"
#include "diaga.h"
        real rrc14std
#include "ice.h"
#include "atm.h"

        real src(imt,km,jsmw:jemw,nsrc)
        common/uvokcomm/src

        ! SPK Set index of deepest point and other parameters
        kmt(1,1) = kmt_loc
        tlat(1,1) = tlat_loc
        dayoyr = day_loc        ! relative day number referenced to the beginning of the current year
        relyr = relyr_loc       ! number of years (and fractional years) of model

        ! copy metos3d variables to UVOK variables
        t(1,1:kmt_loc,1,itemp,taum1) = TEMP(1:kmt_loc)      ! degree Celsius
        t(1,1:kmt_loc,1,isalt,taum1) = SALT(1:kmt_loc)

        print *
        print *, imt, km, jmw, jsmw, jemw, nsrc, nt
        print *, itemp, isalt, idic, ic14, ialk, io2
        print *, ipo4, iphyt, izoop, idetr, ino3, idiaz
        print *, itrc

        print *, 't'
        print *, t(1,1:kmt_loc,1,itemp,taum1)
        print *, t(1,1:kmt_loc,1,isalt,taum1)
        print *, t(1,1:kmt_loc,1,idic,taum1)
        print *, t(1,1:kmt_loc,1,ic14,taum1)
        print *, t(1,1:kmt_loc,1,ialk,taum1)
        print *, t(1,1:kmt_loc,1,io2,taum1)
        print *, t(1,1:kmt_loc,1,ipo4,taum1)
        print *, t(1,1:kmt_loc,1,iphyt,taum1)
        print *, t(1,1:kmt_loc,1,izoop,taum1)
        print *, t(1,1:kmt_loc,1,idetr,taum1)
        print *, t(1,1:kmt_loc,1,ino3,taum1)
        print *, t(1,1:kmt_loc,1,idiaz,taum1)

        print *, 'src'
        print *, src(1,1:kmt_loc,1,itrc(idic))
        print *, src(1,1:kmt_loc,1,itrc(ic14))
        print *, src(1,1:kmt_loc,1,itrc(ialk))
        print *, src(1,1:kmt_loc,1,itrc(io2))
        print *, src(1,1:kmt_loc,1,itrc(ipo4))
        print *, src(1,1:kmt_loc,1,itrc(iphyt))
        print *, src(1,1:kmt_loc,1,itrc(izoop))
        print *, src(1,1:kmt_loc,1,itrc(idetr))
        print *, src(1,1:kmt_loc,1,itrc(ino3))
        print *, src(1,1:kmt_loc,1,itrc(idiaz))
!        call exit(1)

        ! SPK set SBC arrays: code based on UVic_ESCM.F
        ! SPK note: we set tau, taum1 and taup1 to all point to the same location in S/R uvok_ini
        i = 1
        j = 1
        tau = 0
        taum1 = 0
        taup1 = 0

        ! TEMP, SALT
        if (isst .ne. 0)    sbc(i,j,isst)       = t(i,1,j,itemp,taup1)
        if (isss .ne. 0)    sbc(i,j,isss)       = t(i,1,j,isalt,taup1)
        ! DIC,C14,ALK,O2,PO4,PHY,ZOO,DET,NO3,DIAZ
        if (issdic .ne. 0)  sbc(i,j,issdic)     = t(i,1,j,idic,taup1)
        if (issc14 .ne. 0)  sbc(i,j,issc14)     = t(i,1,j,ic14,taup1)
        if (issalk .ne. 0)  sbc(i,j,issalk)     = t(i,1,j,ialk,taup1)
        if (isso2 .ne. 0)   sbc(i,j,isso2)      = t(i,1,j,io2,taup1)
        if (isspo4 .ne. 0)  sbc(i,j,isspo4)     = t(i,1,j,ipo4,taup1)
        if (issphyt .ne. 0) sbc(i,j,issphyt)    = t(i,1,j,iphyt,taup1)
        if (isszoop .ne. 0) sbc(i,j,isszoop)    = t(i,1,j,izoop,taup1)
        if (issdetr .ne. 0) sbc(i,j,issdetr)    = t(i,1,j,idetr,taup1)
        if (issno3 .ne. 0)  sbc(i,j,issno3)     = t(i,1,j,ino3,taup1)
        if (issdiaz .ne. 0) sbc(i,j,issdiaz)    = t(i,1,j,idiaz,taup1)

        print *, 'sbc'
        print *, numsbc
        print *, isst, isss
        print *, issdic, issc14, issalk, isso2, isspo4
        print *, issphyt, isszoop, issdetr, issno3, issdiaz
        print *
        if (isst .ne. 0)    print *, sbc(i,j,isst)
        if (isss .ne. 0)    print *, sbc(i,j,isss)
        if (issdic .ne. 0)  print *, sbc(i,j,issdic)
        if (issc14 .ne. 0)  print *, sbc(i,j,issc14)
        if (issalk .ne. 0)  print *, sbc(i,j,issalk)
        if (isso2 .ne. 0)   print *, sbc(i,j,isso2)
        if (isspo4 .ne. 0)  print *, sbc(i,j,isspo4)
        if (issphyt .ne. 0) print *, sbc(i,j,issphyt)
        if (isszoop .ne. 0) print *, sbc(i,j,isszoop)
        if (issdetr .ne. 0) print *, sbc(i,j,issdetr)
        if (issno3 .ne. 0)  print *, sbc(i,j,issno3)
        if (issdiaz .ne. 0) print *, sbc(i,j,issdiaz)
!        call exit(1)

        ! SPK The sequence here is based on the time stepping in UVic_ESCM.F
        ! SPK copy over winds and other data needed subsequently by S/R gasbc and S/R tracer
        sbc(i,j,iws) = winds_loc
        dnswr(1,1) = swr_loc
        aice(1,1,2) = aice_loc
        hice(1,1,2) = hice_loc
        hsno(1,1,2) = hsno_loc

        ! SPK set atmospheric CO2 here as the call to co2ccndata in gasbc is to a stub routine where
        ! it would otherwise be read from file
!        co2ccn = co2atm
!        print *, 'co2ccn'
!        print *, co2ccn

        ! SPK set global averages here if necessary as the calls in gasbc are to stub routines
        itr = 0
        do n=3,nt
            itr = itr + 1
 	        gaost(n) = TR_surf_glob(itr)
        enddo

        print *, 'gasbc, pre'
        print *, 'sbc'
        print *, sbc
        print *, 'stf'
        print *, stf

        ! get the atmospheric S.B.C.
        call gasbc (1, imt, 1, jmt)

        print *, 'gasbc, post'
        print *, 'sbc'
        print *, sbc
        print *, 'stf'
        print *, stf
!        call exit(1)

        ! SPK S/R gasbc computes the air-sea gas exchange flux in umol cm-2 s-1. We convert it here
        ! to mol m-2 s-1.
        gasexfluxloc = sbc(1,1,idicflx)*1.e-2

        do k=1,kmt_loc
            if (k.le.kmfe) then
                fe_dissolved(1,1,k,1) = fe_dissolved_loc(k)*1.e9  ! input is in nM; convert here as per code in npzd_src.F where it is otherwise done when read from file
            endif
        enddo

        ! Salt flux for virtual flux calculation in gosbc; emp_loc should be in cm/s
        ! Note: emp_glob is now not used and global mean freshwater flux is set to zero in gosbc.
        sbc(1,1,isflx) = emp_loc*socn

        print *, 'emp_loc, socn', emp_loc, socn

        print *, 'gosbc, pre'
        print *, 'sbc'
        print *, sbc
        print *, 'stf'
        print *, stf

        ! get ocean S.B.C.s
        call gosbc (1, imt, 1, jmt)

        print *, 'gosbc, post'
        print *, 'sbc'
        print *, sbc
        print *, 'stf'
        print *, stf
!        call exit(1)

        ! call the ocean model once for each time step until one
        ! segment of "segtim" days is complete. hold ocean S.B.C. fixed
        ! during each segment and predict average S.B.C. for atmos

        ! SPK S/R mom calls S/R setvbc and S/R tracer
        js = 1
        je = 1
        is = 1
        ie = 1
        joff = 0

        print *, 'setvbc, pre'
        print *, 't'
        print *, t(1,1:kmt_loc,1,itemp,taum1)
        print *, t(1,1:kmt_loc,1,isalt,taum1)
        print *, t(1,1:kmt_loc,1,idic,taum1)
        print *, t(1,1:kmt_loc,1,ic14,taum1)
        print *, t(1,1:kmt_loc,1,ialk,taum1)
        print *, t(1,1:kmt_loc,1,io2,taum1)
        print *, t(1,1:kmt_loc,1,ipo4,taum1)
        print *, t(1,1:kmt_loc,1,iphyt,taum1)
        print *, t(1,1:kmt_loc,1,izoop,taum1)
        print *, t(1,1:kmt_loc,1,idetr,taum1)
        print *, t(1,1:kmt_loc,1,ino3,taum1)
        print *, t(1,1:kmt_loc,1,idiaz,taum1)
        print *, 'src'
        print *, src(1,1:kmt_loc,1,itrc(idic))
        print *, src(1,1:kmt_loc,1,itrc(ic14))
        print *, src(1,1:kmt_loc,1,itrc(ialk))
        print *, src(1,1:kmt_loc,1,itrc(io2))
        print *, src(1,1:kmt_loc,1,itrc(ipo4))
        print *, src(1,1:kmt_loc,1,itrc(iphyt))
        print *, src(1,1:kmt_loc,1,itrc(izoop))
        print *, src(1,1:kmt_loc,1,itrc(idetr))
        print *, src(1,1:kmt_loc,1,itrc(ino3))
        print *, src(1,1:kmt_loc,1,itrc(idiaz))
        print *, 'sbc'
        if (isst .ne. 0)    print *, sbc(i,j,isst)
        if (isss .ne. 0)    print *, sbc(i,j,isss)
        if (issdic .ne. 0)  print *, sbc(i,j,issdic)
        if (issc14 .ne. 0)  print *, sbc(i,j,issc14)
        if (issalk .ne. 0)  print *, sbc(i,j,issalk)
        if (isso2 .ne. 0)   print *, sbc(i,j,isso2)
        if (isspo4 .ne. 0)  print *, sbc(i,j,isspo4)
        if (issphyt .ne. 0) print *, sbc(i,j,issphyt)
        if (isszoop .ne. 0) print *, sbc(i,j,isszoop)
        if (issdetr .ne. 0) print *, sbc(i,j,issdetr)
        if (issno3 .ne. 0)  print *, sbc(i,j,issno3)
        if (issdiaz .ne. 0) print *, sbc(i,j,issdiaz)
        print *, 'sbc'
        print *, sbc
        print *, 'stf'
        print *, stf

        call setvbc(joff, js, je, is, ie)

        print *, 'setvbc, post'
        print *, 't'
        print *, t(1,1:kmt_loc,1,itemp,taum1)
        print *, t(1,1:kmt_loc,1,isalt,taum1)
        print *, t(1,1:kmt_loc,1,idic,taum1)
        print *, t(1,1:kmt_loc,1,ic14,taum1)
        print *, t(1,1:kmt_loc,1,ialk,taum1)
        print *, t(1,1:kmt_loc,1,io2,taum1)
        print *, t(1,1:kmt_loc,1,ipo4,taum1)
        print *, t(1,1:kmt_loc,1,iphyt,taum1)
        print *, t(1,1:kmt_loc,1,izoop,taum1)
        print *, t(1,1:kmt_loc,1,idetr,taum1)
        print *, t(1,1:kmt_loc,1,ino3,taum1)
        print *, t(1,1:kmt_loc,1,idiaz,taum1)
        print *, 'src'
        print *, src(1,1:kmt_loc,1,itrc(idic))
        print *, src(1,1:kmt_loc,1,itrc(ic14))
        print *, src(1,1:kmt_loc,1,itrc(ialk))
        print *, src(1,1:kmt_loc,1,itrc(io2))
        print *, src(1,1:kmt_loc,1,itrc(ipo4))
        print *, src(1,1:kmt_loc,1,itrc(iphyt))
        print *, src(1,1:kmt_loc,1,itrc(izoop))
        print *, src(1,1:kmt_loc,1,itrc(idetr))
        print *, src(1,1:kmt_loc,1,itrc(ino3))
        print *, src(1,1:kmt_loc,1,itrc(idiaz))
        print *, 'sbc'
        if (isst .ne. 0)    print *, sbc(i,j,isst)
        if (isss .ne. 0)    print *, sbc(i,j,isss)
        if (issdic .ne. 0)  print *, sbc(i,j,issdic)
        if (issc14 .ne. 0)  print *, sbc(i,j,issc14)
        if (issalk .ne. 0)  print *, sbc(i,j,issalk)
        if (isso2 .ne. 0)   print *, sbc(i,j,isso2)
        if (isspo4 .ne. 0)  print *, sbc(i,j,isspo4)
        if (issphyt .ne. 0) print *, sbc(i,j,issphyt)
        if (isszoop .ne. 0) print *, sbc(i,j,isszoop)
        if (issdetr .ne. 0) print *, sbc(i,j,issdetr)
        if (issno3 .ne. 0)  print *, sbc(i,j,issno3)
        if (issdiaz .ne. 0) print *, sbc(i,j,issdiaz)
        print *, 'sbc'
        print *, sbc
        print *, 'stf', isalt
        print *, stf
!        call exit(1)

        print *, 'tracer, pre'
        call tracer(joff, js, je, is, ie)
        print *, 'tracer, post'
        print *, 't'
        print *, t(1,1:kmt_loc,1,itemp,taum1)
        print *, t(1,1:kmt_loc,1,isalt,taum1)
        print *, t(1,1:kmt_loc,1,idic,taum1)
        print *, t(1,1:kmt_loc,1,ic14,taum1)
        print *, t(1,1:kmt_loc,1,ialk,taum1)
        print *, t(1,1:kmt_loc,1,io2,taum1)
        print *, t(1,1:kmt_loc,1,ipo4,taum1)
        print *, t(1,1:kmt_loc,1,iphyt,taum1)
        print *, t(1,1:kmt_loc,1,izoop,taum1)
        print *, t(1,1:kmt_loc,1,idetr,taum1)
        print *, t(1,1:kmt_loc,1,ino3,taum1)
        print *, t(1,1:kmt_loc,1,idiaz,taum1)
        print *, 'src'
        print *, src(1,1:kmt_loc,1,itrc(idic))
        print *, src(1,1:kmt_loc,1,itrc(ic14))
        print *, src(1,1:kmt_loc,1,itrc(ialk))
        print *, src(1,1:kmt_loc,1,itrc(io2))
        print *, src(1,1:kmt_loc,1,itrc(ipo4))
        print *, src(1,1:kmt_loc,1,itrc(iphyt))
        print *, src(1,1:kmt_loc,1,itrc(izoop))
        print *, src(1,1:kmt_loc,1,itrc(idetr))
        print *, src(1,1:kmt_loc,1,itrc(ino3))
        print *, src(1,1:kmt_loc,1,itrc(idiaz))
        print *, 'sbc'
        print *, sbc
        print *, 'stf', isalt
        print *, stf
        call exit(1)


        ! SPK Add surface fluxes (calculated in setvbc.F) as a body force term in surface layer. Be very careful of tracer indexing!
        do j=js,je
            do i=is,ie
                do n=3,nt
                    src(i,1,j,itrc(n)) =
     &                  src(i,1,j,itrc(n))+stf(i,j,n)*dztr(1)
                enddo
            enddo
        enddo

        ! SPK Virtual flux is added above to the air-sea gas exchange flux computed in S/R gasbc. Units
        ! are umol cm-2 s-1. We convert here to mol m-2 s-1.
        totfluxloc = sbc(1,1,idicflx)*1.e-2

      end subroutine uvok_calc


