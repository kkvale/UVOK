!
! Metos3D: A Marine Ecosystem Toolkit for Optimization and Simulation in 3-D
! Copyright (C) 2014  Jaroslaw Piwonski, CAU, jpi@informatik.uni-kiel.de
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!

!
!   metos3dbgcinit
!
subroutine metos3dbgcinit(ny, nx, nu, nb, nd, dt, q, t, y, u, b, d, ndiag, diag)
    implicit none
    ! input variables
    integer :: ny, nx, nu, nb, nd, ndiag
    real(8) :: dt, q(nx, ny), t, y(nx, ny), u(nu), b(nb), d(nx, nd), diag(nx, ndiag)



!    print *,
!    call uvok_ini(nx, )

!SUBROUTINE UVOK_INI(
!nzmax,
!z,
!drF,
!DeltaT,
!S_surf_glob,
!TR_surf_glob,
!debugFlag)
!
!integer nzmax
!real z(km), drF(km), DeltaT
!real S_surf_glob, TR_surf_glob(nsrc)
!integer debugFlag

!uvok_ini_(
!&nzmax,
!zt,
!drF,
!&DeltaT,
!&Sglobavg,
!TRglobavg,
!&debugFlag);

end subroutine

!
!   metos3dbgcbegin
!
subroutine metos3dbgcbegin(ny, nx, nu, nb, nd, dt, q, t, y, u, b, d, ndiag, diag)
    implicit none
    ! input variables
    integer :: ny, nx, nu, nb, nd, ndiag
    real(8) :: dt, q(nx, ny), t, y(nx, ny), u(nu), b(nb), d(nx, nd), diag(nx, ndiag)
end subroutine

!
!   metos3dbgc
!
subroutine metos3dbgc(ny, nx, nu, nb, nd, dt, q, t, y, u, b, d, ndiag, diag)
    implicit none
    ! input variables
    integer :: ny, nx, nu, nb, nd, ndiag
    real(8) :: dt, q(nx, ny), t, y(nx, ny), u(nu), b(nb), d(nx, nd), diag(nx, ndiag)

!SUBROUTINE UVOK_CALC(
!kmt_loc,
!tlat_loc,
!day_loc,
!relyr_loc,
!&     TEMP, SALT, TR_surf_glob,dz_loc,z,
!&     winds_loc,
!&     fe_dissolved_loc,
!&     swr_loc,
!&     aice_loc, hice_loc, hsno_loc,
!& emp_loc, emp_glob,
!& gasexfluxloc, totfluxloc,
!& debugFlag)

!uvok_calc_(&nzloc,&locallatitude[ip],&day,&relyr,
!&localTs[kl],&localSs[kl],&TRglobavg[0],&localdz[kl],zt,
!&localwind[ip],
!&localFe_dissolved[kl],
!&localswrad[ip],
!&localaice[ip], &localhice[ip], &localhsno[ip],
!&localEmP[ip], &EmPglobavg,
!&localgasexflux, &localtotflux,
!&debugFlag);

end subroutine

!
!   metos3dbgcend
!
subroutine metos3dbgcend(ny, nx, nu, nb, nd, dt, q, t, y, u, b, d, ndiag, diag)
    implicit none
    ! input variables
    integer :: ny, nx, nu, nb, nd, ndiag
    real*8  :: dt, q(nx, ny), t, y(nx, ny), u(nu), b(nb), d(nx, nd), diag(nx, ndiag)
end subroutine

!
!   metos3dbgcfinal
!
subroutine metos3dbgcfinal(ny, nx, nu, nb, nd, dt, q, t, y, u, b, d, ndiag, diag)
    implicit none
    ! input variables
    integer :: ny, nx, nu, nb, nd, ndiag
    real*8  :: dt, q(nx, ny), t, y(nx, ny), u(nu), b(nb), d(nx, nd), diag(nx, ndiag)
end subroutine



! ------------------------------------------------------------------------------------------
! (un)used subroutine stubs
! ------------------------------------------------------------------------------------------

subroutine areaavg (data, dmsk, avg)
!    implicit none
end

subroutine setbcx (a, imt, jmtorkm)
!    implicit none
!    dimension a(imt,jmtorkm)
end

subroutine data (is, ie, js, je)
!    implicit none
!    integer is, ie, js, je
end

subroutine co2forc
!    implicit none
end

subroutine c14data
!    implicit none
end

subroutine co2ccndata
!    implicit none
end

subroutine defvar (name, ncid, nd, id, rmin, rmax, axis, type, lname, sname, units)
!    implicit none
!    character(*) name, axis, lname, sname, type, units
!    integer nd, id(nd), ncid
!    real rmax, rmin
end

subroutine putvaramsk (name, ncid, ln, is, ic, din, dm, s, o)
!    implicit none
!    character(*) name
!    integer ic(10), is(10), ln, ncid
!    real din(ln)
!    real dm(ln), o, s
end

subroutine getunit (iounit, oldfilename, optionlist)
    implicit none
    integer iounit
    character(*) :: oldfilename, optionlist
    iounit = 1
end

subroutine relunit (iounit)
end

character(120) function new_file_name(filename)
    implicit none
    character(*) filename
    new_file_name = filename
end
