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

    real(8) :: salt_avg = -3.451231022787175e-004
    real(8) :: trace_avg(10) = (/ &
        2.315e0, 2.3140740e-012, 2.429e0, 0.1692e0, 0.543e0, 0.14e0, 1.4e-002, 1.0e-004, 5.3e0, 1.4e-002 /)
    real(8) :: dtbgc = 28800.0d0
    integer :: debug = 1

    call uvok_ini( &
        nx,         &   ! nzmax
        d(1,1),     &   ! z
        d(1,2),     &   ! drF
        dtbgc,      &   ! DeltaT
        salt_avg,   &   ! S_surf_glob
        trace_avg,  &   ! TR_surf_glob
        debug       &   ! debugFlag
    )

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

    real(8) :: trace_avg(10) = (/ &
        2.315e0, 2.3140740e-012, 2.429e0, 0.1692e0, 0.543e0, 0.14e0, 1.4e-002, 1.0e-004, 5.3e0, 1.4e-002 /)
    integer :: debug = 1
    real(8) :: emp_glob = 0.d0
    real(8) :: gasexfluxloc = 0.d0
    real(8) :: totfluxloc = 0.d0
    real(8) :: day_loc

    integer :: istep
    real(8) :: day_frac
    real(8) :: dtbgc = 28800.0d0

    ! hard coded number of time steps
    istep = nint(t*1095.0)
    day_frac = 365.0/1095.0
    day_loc = istep*day_frac

    call uvok_copy_to(ny, nx, y)
    call uvok_calc( &
        nx,             &   ! kmt_loc
        b(1),           &   ! tlat_loc
        day_loc,        &   ! day_loc
        t,              &   ! relyr_loc
        d(1,3),         &   ! TEMP
        d(1,4),         &   ! SALT
        trace_avg,      &   ! TR_surf_glob
        d(1,2),         &   ! dz_loc
        d(1,1),         &   ! z
        b(2),           &   ! winds_loc
        d(1,5),         &   ! fe_dissolved_loc
        b(3),           &   ! swr_loc
        b(4),           &   ! aice_loc
        b(5),           &   ! hice_loc
        b(6),           &   ! hsno_loc
        b(7),           &   ! emp_loc
        emp_glob,       &   ! emp_glob
        gasexfluxloc,   &   ! gasexfluxloc
        totfluxloc,     &   ! totfluxloc
        debug           &   ! debugFlag
    )
    call uvok_copy_from(ny, nx, q)

    ! scale
    q = dtbgc * q

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
    implicit none
    integer imt, jmtorkm
    real a(imt,jmtorkm)
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



