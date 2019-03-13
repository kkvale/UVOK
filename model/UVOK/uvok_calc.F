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

        if (kmt_loc .eq. km) then
            print *, kmt_loc
            print *, tlat_loc, day_loc, relyr_loc
            print *, TEMP(kmt_loc), SALT(kmt_loc)
            print *, TR_surf_glob(nsrc)
            print *, dz_loc(kmt_loc)
            print *, z_loc(kmt_loc)

            print *, winds_loc
            print *, fe_dissolved_loc(kmt_loc)
            print *, swr_loc
            print *, aice_loc, hice_loc, hsno_loc
            print *, emp_loc, emp_glob
            print *, gasexfluxloc, totfluxloc
            print *, debugFlag
            call exit(1)
        endif

      end subroutine uvok_calc



C     integer i, k, j, ip, n
C     real tarea, tsflx, rsocn, tmp
C!       joff = offset relating "j" in the MW to latitude "jrow"
C!       js   = starting row in the MW
C!       je   = ending row in the MW
C!       is   = starting longitude index in the MW
C!       ie   = ending longitude index in the MW
C
C     integer joff, is, ie, js, je, itr
C
C#include "param.h"
C#include "pconst.h"
C#include "stdunits.h"
C#include "coord.h"
C#include "csbc.h"
C#include "grdvar.h"
C#include "levind.h"
C#include "mw.h"
C#include "scalar.h"
C#include "tmngr.h"
C# if defined O_npzd || defined O_carbon_14
C#include "npzd.h"
C#include "calendar.h"
C# endif
C# if defined O_save_convection || defined O_carbon_14
C#include "diaga.h"
C     real rrc14std
C# endif
C# if defined O_ice
C#  if defined O_ice_cpts
C#include "cpts.h"
C#  endif
C#include "ice.h"
C# endif
C# if defined O_npzd
CC     real t_in(km), po4_in(km)
CC#  if defined O_npzd_o2
CC     real o2_in(km)
CC#  endif
CC#  if defined O_carbon_13
CC     real s_in(km), dic_in(km), alk_in(km), co2_in, dic13_in(km)
CC#  endif
CC#  if defined O_npzd_nitrogen
CC     real no3_in(km), sgb_in(km)
CC#   if defined O_npzd_nitrogen_15
CC     real din15_in(km)
CC#   endif
CC#  endif
CC#  if defined O_npzd_fe_limitation
CC     real felimit_in(km), felimit_D_in(km)
CC#  endif
CC     real expo, tnpzd(km,ntnpzd)
C#  if defined O_embm
C#include "atm.h"
C#   if defined O_carbon_13
C#include "cembm.h"
C#   endif
C#  endif
CC# endif
CC# if defined O_carbon_fnpzd
CC#include "calendar.h"
C# endif
C
C# if defined O_npzd || defined O_carbon_14
C     real src(imt,km,jsmw:jemw,nsrc)
C     common/uvokcomm/src
C# endif
C
C#ifdef O_TMM_partial_step_topo
CCSPK Recompute local grid-dependent data
C#if defined O_mom && defined O_npzd
C     dzt(:) = c0
C     dzw(:) = c0
C     zt(:) = c0
C     zw(:) = c0
C
C     do k=1,kmt_loc
C       dzt(k) = dz_loc(k)
C       zt(k) = z(k)
C     enddo
C
C     zw(1) = dz_loc(1)
C     do k=2,kmt_loc
C       zw(k) = zw(k-1) + dz_loc(k)
C     enddo
C
C     do k=1,kmt_loc-1
C       dzw(k) = zt(k+1) - zt(k)
C     enddo
C     dzw(0) = zt(1)
C     dzw(kmt_loc) = zw(kmt_loc) - zt(kmt_loc)
C
C     do k=1,kmt_loc
C       c2dzt(k) = c2*dzt(k)
C       dzt2r(k) = c1/c2dzt(k)
C     enddo
C
C     dzwr(kmt_loc)  = c1/dzw(kmt_loc)
C     dzw2r(kmt_loc) = p5/dzw(kmt_loc)
C
C     do k=1,kmt_loc
C       dzwr(k-1)    = c1/dzw(k-1)
C       dzw2r(k-1)   = p5/dzw(k-1)
C     enddo
C
C     do k=1,kmt_loc
C       dztr(k)  = c1/dzt(k)
C     enddo
C#endif
CCCCCC
C#endif !O_TMM_partial_step_topo
C
CCSPK  Set index of deepest point and other parameters
C  kmt(1,1) = kmt_loc
C  tlat(1,1) = tlat_loc
C  dayoyr = day_loc ! relative day number referenced to the beginning of the current year
C  relyr = relyr_loc ! number of years (and fractional years) of model
C
CCSPK Copy T/S from TMM to UVOK
CCSPK Other tracers should already have been copied by a previous call to S/R uvok_copy_data
C  t(1,1:kmt_loc,1,itemp,taum1) = TEMP(1:kmt_loc) ! !degree Celsius
C  t(1,1:kmt_loc,1,isalt,taum1) = SALT(1:kmt_loc)
C
CCSPK set SBC arrays: code based on UVic_ESCM.F
CCSPK note: we set tau, taum1 and taup1 to all point to the same location in S/R uvok_ini
C      i = 1
C      j = 1
CC       do j=1,jmt
CC         do i=1,imt
C# if !defined O_replacst
C           if (isst .ne. 0) sbc(i,j,isst) = t(i,1,j,itemp,taup1)
C           if (isss .ne. 0) sbc(i,j,isss) = t(i,1,j,isalt,taup1)
C# endif
C# if defined O_carbon
C           if (issdic .ne. 0) sbc(i,j,issdic) = t(i,1,j,idic,taup1)
C#  if defined O_carbon_13
C           if (issdic13 .ne. 0) sbc(i,j,issdic13) =
C    &           t(i,1,j,idic13,taup1)
C#  endif
C#  if defined O_carbon_14
C           if (issc14 .ne. 0) sbc(i,j,issc14) = t(i,1,j,ic14,taup1)
C#  endif
C# endif
C# if defined O_npzd_alk
C           if (issalk .ne. 0) sbc(i,j,issalk) = t(i,1,j,ialk,taup1)
C# endif
C# if defined O_npzd_o2
C           if (isso2 .ne. 0) sbc(i,j,isso2) = t(i,1,j,io2,taup1)
C# endif
C# if defined O_npzd
C           if (isspo4 .ne. 0) sbc(i,j,isspo4) = t(i,1,j,ipo4,taup1)
CC            if (issdop .ne. 0) sbc(i,j,issdon) = t(i,1,j,idop,taup1)
C           if (issphyt .ne. 0) sbc(i,j,issphyt) = t(i,1,j,iphyt,taup1)
C           if (isszoop .ne. 0) sbc(i,j,isszoop) = t(i,1,j,izoop,taup1)
C           if (issdetr .ne. 0) sbc(i,j,issdetr) = t(i,1,j,idetr,taup1)
C#  if defined O_npzd_nitrogen
C           if (issno3 .ne. 0) sbc(i,j,issno3) = t(i,1,j,ino3,taup1)
CC            if (issdon .ne. 0) sbc(i,j,issdon) = t(i,1,j,idon,taup1)
C           if (issdiaz .ne. 0) sbc(i,j,issdiaz) = t(i,1,j,idiaz,taup1)
C#   if defined O_npzd_nitrogen_15
C           if (issdin15 .ne. 0) sbc(i,j,issdin15)
C    &           = t(i,1,j,idin15,taup1)
C           if (issdon15 .ne. 0) sbc(i,j,issdon15)
C    &           = t(i,1,j,idon15,taup1)
C           if (issphytn15 .ne. 0) sbc(i,j,issphytn15)
C    &           = t(i,1,j,iphytn15,taup1)
C           if (isszoopn15 .ne. 0) sbc(i,j,isszoopn15)
C    &           = t(i,1,j,izoopn15,taup1)
C           if (issdetrn15 .ne. 0) sbc(i,j,issdetrn15)
C    &           = t(i,1,j,idetrn15,taup1)
C           if (issdiazn15 .ne. 0) sbc(i,j,issdiazn15)
C    &           = t(i,1,j,idiazn15,taup1)
C#   endif
C#  endif
C#  if defined O_npzd_iron
C      if (issdfe .ne. 0) sbc(i,j,issdfe) = t(i,1,j,idfe,taup1)
C      if (issdetrfe .ne. 0) sbc(i,j,issdetrfe) = t(i,1,j,idetrfe,taup1)
C#  endif
C#  if defined O_carbon_13
C           if (issdoc13 .ne. 0) sbc(i,j,issdoc13) =
C    &           t(i,1,j,idoc13,taup1)
C           if (issphytc13 .ne. 0) sbc(i,j,issphytc13) =
C    &           t(i,1,j,iphytc13,taup1)
C           if (isszoopc13 .ne. 0) sbc(i,j,isszoopc13) =
C    &           t(i,1,j,izoopc13,taup1)
C           if (issdetrc13 .ne. 0) sbc(i,j,issdetrc13) =
C    &           t(i,1,j,idetrc13,taup1)
C#   if defined O_npzd_nitrogen
C           if (issdiazc13 .ne. 0) sbc(i,j,issdiazc13) =
C    &           t(i,1,j,idiazc13,taup1)
C#   endif
C#  endif
C# endif
C# if defined O_cfcs_data || defined O_cfcs_data_transient
C           if (isscfc11 .ne. 0) sbc(i,j,isscfc11) =
C    &        t(i,1,j,icfc11,taup1)
C           if (isscfc12 .ne. 0) sbc(i,j,isscfc12) =
C    &        t(i,1,j,icfc12,taup1)
C# endif
CC         enddo
CC       enddo
C
CCSPK The sequence here is based on the time stepping in UVic_ESCM.F
C
CCSPK copy over winds and other data needed subsequently by S/R gasbc and S/R tracer
C#ifdef O_carbon
C     sbc(i,j,iws) = winds_loc
C#endif
C
C#  if defined O_npzd_nitrogen_15
C  do k=1,kmt_loc
C    sg_bathy(1,1,k) = sg_bathy_loc(k)
C  enddo
C#  endif
C
C#  if defined O_embm
C  dnswr(1,1) = swr_loc
C#  endif
C
C#  if defined O_ice
C#   if !defined O_ice_cpts
C  aice(1,1,2) = aice_loc
C  hice(1,1,2) = hice_loc
C  hsno(1,1,2) = hsno_loc
C#   endif
C#  endif
C
CCSPK set atmospheric CO2 here as the call to co2ccndata in gasbc is to a stub routine where
CC    it would otherwise be read from file
C#if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
C#   if defined O_carbon_co2_2d
C       at(1,1,2,ico2) = co2atm
C#   else
C       co2ccn = co2atm
C#   endif
C#endif
C
CCSPK set global averages here if necessary as the calls in gasbc are to stub routines
C#  if !defined O_constant_flux_reference
C     itr = 0
C     do n=3,nt
CC        if (itrc(n) .ne. 0) then
C        itr = itr + 1
C	 gaost(n) = TR_surf_glob(itr)
CC        endif
C     enddo
C# endif
C
CCSPK set atmospheric dc14 here as the call to c14data in gasbc is to a stub routine where
CC    it would otherwise be read from file
C# if defined O_c14ccn_data
C     dc14ccnn = dc14ccnnatm
C     dc14ccne = dc14ccneatm
C     dc14ccns = dc14ccnsatm
C#endif
C
CCSPK These arrays are first set to zero in S/R glsbc and then updated in S/R mtlm. They
CC    are used in S/R gasbc but as we're not using the land model they should not be
CC    referenced anywhere. Still, we set them to zero anyway to be safe.
C#if defined O_mtlm && defined O_carbon
C     sbc(i,j,inpp) = 0.0
C     sbc(i,j,isr) = 0.0
C     sbc(i,j,iburn) = 0.0
C#endif
C
C!-----------------------------------------------------------------------
C!       get the atmospheric S.B.C.
C!-----------------------------------------------------------------------
C
C       call gasbc (1, imt, 1, jmt)
C
CCSPK S/R gasbc modifies co2ccn. This should not have any effect if co2emit and global mean
CC    dic flux calculated by a call to stub routine areaavg are exactly zero as we want them
CC    to be. But to be safe, we reset its value again here.
C#if defined O_co2ccn_data || defined O_TMM_interactive_atmosphere
C#   if defined O_carbon_co2_2d
C       at(1,1,2,ico2) = co2atm
C#   else
C       co2ccn = co2atm
C#   endif
C#endif
C
CCSPK S/R gasbc computes the air-sea gas exchange flux in umol cm-2 s-1. We convert it here
CC    to mol m-2 s-1.
C#if defined O_carbon
C     gasexfluxloc = sbc(1,1,idicflx)*1.e-2
C#endif
C
C#  if defined O_npzd_fe_limitation
C  do k=1,kmt_loc
C    if (k.le.kmfe) then
C      fe_dissolved(1,1,k,1) = fe_dissolved_loc(k)*1.e9  ! input is in nM; convert here as per code in npzd_src.F where it is otherwise done when read from file
C    endif
C  enddo
C#  endif
C
CC Note: we do this after call to S/R gasbc because sbc(1,1,idfeadep) is set to zero there (before being read from file)
C#if defined O_npzd_iron
C     sbc(1,1,idfeadep)=fe_adep_loc ! note input is in molFe/m^2/s; it is converted to mmolFe/m^2/s in tracer.F so we don't do it here
C  do k=1,kmt_loc
C       fe_hydr(1,1,k) = fe_hydr_loc(k)*1000.0 ! input is molFe/m^3/s; we convert it to mmolFe/m^3/s as per code in npzd_src.F where it is otherwise done when read from file
C     enddo
C#endif
C
CC Salt flux for virtual flux calculation in gosbc; emp_loc should be in cm/s
CC Note: emp_glob is now not used and global mean freshwater flux is set to zero in gosbc.
C#if defined O_embm
C     sbc(1,1,isflx) = emp_loc*socn
C#endif
C
CC#if defined O_mtlm
CC!-----------------------------------------------------------------------
CC!       get land S.B.C.s
CC!-----------------------------------------------------------------------
CC
CC       call glsbc (1, imt, 1, jmt)
CC
CC!----------------------------------------------------------------------
CC!       call the land-surface and vegetation  model once for each time
CC!       step until one segment of "segtim" days is complete.
CC!-----------------------------------------------------------------------
CC
CC       do loop=1,ntspls
CC         call mtlm (1, imt, 1, jmt)
CC       enddo
CC
CC#endif
C
CC
CC!-----------------------------------------------------------------------
CC!       get ocean S.B.C.s
CC!-----------------------------------------------------------------------
C
C     call gosbc (1, imt, 1, jmt)
C
CC
CC!-----------------------------------------------------------------------
CC!       call the ocean model once for each time step until one
CC!       segment of "segtim" days is complete. hold ocean S.B.C. fixed
CC!       during each segment and predict average S.B.C. for atmos
CC!-----------------------------------------------------------------------
CC
CC       do loop=1,ntspos
CC         call mom
C
CCSPK S/R mom calls S/R setvbc and S/R tracer
C
C      js = 1
C      je = 1
C      is = 1
C      ie = 1
C      joff = 0
C
C      call setvbc(joff, js, je, is, ie)
C      call tracer(joff, js, je, is, ie)
C
CC      if (debugFlag.eq.1) then
C!      n=ialk
C!      print*,'TEST:',idic,isdic,itrc(n),src(1,1,1,itrc(n)),stf(1,1,n),
CC    &    sbc(1,1,idicflx),dztr(1),tmask(1,1,1)
CC      endif
CCSPK Add surface fluxes (calculated in setvbc.F) as a body force term in surface layer. Be very careful of tracer indexing!
C     do j=js,je
C       do i=is,ie
C	  do n=3,nt
CC	  if (itrc(n) .ne. 0) then
C	   src(i,1,j,itrc(n)) = src(i,1,j,itrc(n))+stf(i,j,n)*dztr(1)
CC	  endif
C	  enddo
C	enddo
C  enddo
CC      if (debugFlag.eq.1) then
CC      n=idic
CC      print*,'TEST:',idic,isdic,itrc(n),src(1,1,1,itrc(n)),stf(1,1,n),
CC    &    sbc(1,1,idicflx),dztr(1),tmask(1,1,1)
CC      endif
C
CCSPK Virtual flux is added above to the air-sea gas exchange flux computed in S/R gasbc. Units
CC    are umol cm-2 s-1. We convert here to mol m-2 s-1.
C#if defined O_carbon
C     totfluxloc = sbc(1,1,idicflx)*1.e-2
C#endif
C
CC      if (flag.eq.1) then
CC      do n=3,nt
CC        if (itrc(n) .ne. 0) then
CC         print*,n-2,mapt(n),src(1,1:kmt_loc,1,itrc(n))
CC        endif
CC      enddo
CC      endif
C
C     END
C#endif
