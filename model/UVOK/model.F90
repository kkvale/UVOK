!
! Metos3D: A Marine Ecosystem Toolkit for Optimization and Simulation in 3-D
! Copyright (C) 2019  Jaroslaw Piwonski, CAU, jpi@informatik.uni-kiel.de
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

    real(8) :: salt_avg = -3.451231022787175d-004
    real(8) :: trace_avg(10) = (/ 2.315d0, 2.3140740d-012, 2.429d0, 0.1692d0, 0.543d0, 0.14d0, 1.4d-002, 1.0d-004, 5.3d0, 1.4d-002/)
    real(8) :: zt(nx), z(nx)
    real(8) :: dtbgc
    logical :: debug

    zt = d(:, 1)                ! layers mids
    z = d(:, 2)                 ! layers depths
    dtbgc = 28800.0d0
    debug = .true.

    print *, nx, ny
    print *, zt
    print *, z
    print *, dtbgc
    print *, salt_avg
    print *, trace_avg
    print *, debug

    call uvok_ini(nx, zt, z, dtbgc, salt_avg, trace_avg, debug)

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

!19
!1750.00000000000        8250.00000000000        17750.0000000000
!30250.0000000000        45750.0000000000        64250.0000000000
!85750.0000000000        110250.000000000        137750.000000000
!168250.000000000        201750.000000000        238250.000000000
!277750.000000000        320250.000000000        365750.000000000
!414250.000000000        465750.000000000        520250.000000000
!577750.000000000
!5000.00000000000        8000.00000000000        11000.0000000000
!14000.0000000000        17000.0000000000        20000.0000000000
!23000.0000000000        26000.0000000000        29000.0000000000
!32000.0000000000        35000.0000000000        38000.0000000000
!41000.0000000000        44000.0000000000        47000.0000000000
!50000.0000000000        53000.0000000000        56000.0000000000
!59000.0000000000
!28800.0000000000
!-3.451231022787175E-004
!2.31500000000000       2.314074000000000E-012   2.42900000000000
!0.169200000000000       0.543000000000000       0.140000000000000
!1.400000000000000E-002  9.999999999999999E-005   5.30000000000000
!1.400000000000000E-002
!1

!19
!1750.0000000000000        8250.0000000000000        17750.000000000000        30250.000000000000        45750.000000000000        64250.000000000000        85750.000000000000        110250.00000000000        137750.00000000000        168250.00000000000        201750.00000000000        238250.00000000000        277750.00000000000        320250.00000000000        365750.00000000000        414250.00000000000        465750.00000000000        520250.00000000000        577750.00000000000
!5000.0000000000000        8000.0000000000000        11000.000000000000        14000.000000000000        17000.000000000000        20000.000000000000        23000.000000000000        26000.000000000000        29000.000000000000        32000.000000000000        35000.000000000000        38000.000000000000        41000.000000000000        44000.000000000000        47000.000000000000        50000.000000000000        53000.000000000000        56000.000000000000        59000.000000000000
!28800.000000000000
!-3.4512310227871747E-004
!2.3149999999999999        2.3140739999999999E-012   2.4289999999999998       0.16919999999999999       0.54300000000000004       0.14000000000000001        1.4000000000000000E-002   1.0000000000000000E-004   5.2999999999999998        1.4000000000000000E-002
!T

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
    open(unit=iounit, file=oldfilename)
end

subroutine relunit (iounit)
    implicit none
    integer iounit
    close(iounit)
end

character(120) function new_file_name(filename)
    implicit none
    character(*) filename
    new_file_name = filename
end
