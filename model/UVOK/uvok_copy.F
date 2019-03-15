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
C       https://raw.githubusercontent.com/samarkhatiwala/tmm/master/models/current/uvok1.0/src/uvok_copy_data.F
C

#include "UVOK_TMM_OPTIONS.h"

      !
      ! uvok_copy_to(ny, nx, y)
      !     copy tracers from metos3d variable `y` to uvok variable `t`
      !
      subroutine uvok_copy_to(ny, nx, y)
        implicit none
        integer ny, nx
        real y(nx, ny)

#include "size.h"
#include "mw.h"

        integer i, j

        ! DIC, C14, ALK, O2, PO4, PHY, ZOO, DET, NO3, DIAZ
        i = 1
        j = 1
        taup1 = 0
        t(i, 1:nx, j, idic, taup1)     = y(1:nx, 1)
        t(i, 1:nx, j, ic14, taup1)     = y(1:nx, 2)
        t(i, 1:nx, j, ialk, taup1)     = y(1:nx, 3)
        t(i, 1:nx, j, io2, taup1)      = y(1:nx, 4)
        t(i, 1:nx, j, ipo4, taup1)     = y(1:nx, 5)
        t(i, 1:nx, j, iphyt, taup1)    = y(1:nx, 6)
        t(i, 1:nx, j, izoop, taup1)    = y(1:nx, 7)
        t(i, 1:nx, j, idetr, taup1)    = y(1:nx, 8)
        t(i, 1:nx, j, ino3, taup1)     = y(1:nx, 9)
        t(i, 1:nx, j, idiaz, taup1)    = y(1:nx, 10)

      end subroutine uvok_copy_to

      !
      ! uvok_copy_from(ny, nx, q)
      !     copy tracers from uvok variable `src` to metos3d variable `q`
      !
      subroutine uvok_copy_from(ny, nx, q)
        implicit none
        integer ny, nx
        real q(nx, ny)

#include "size.h"
#include "mw.h"

        real src(imt,km,jsmw:jemw,nsrc)
        common/uvokcomm/src

        q(1:nx, 1)  = src(1, 1:nx, 1, idic)
        q(1:nx, 2)  = src(1, 1:nx, 1, ic14)
        q(1:nx, 3)  = src(1, 1:nx, 1, ialk)
        q(1:nx, 4)  = src(1, 1:nx, 1, io2)
        q(1:nx, 5)  = src(1, 1:nx, 1, ipo4)
        q(1:nx, 6)  = src(1, 1:nx, 1, iphyt)
        q(1:nx, 7)  = src(1, 1:nx, 1, izoop)
        q(1:nx, 8)  = src(1, 1:nx, 1, idetr)
        q(1:nx, 9)  = src(1, 1:nx, 1, ino3)
        q(1:nx, 10) = src(1, 1:nx, 1, idiaz)

      end subroutine uvok_copy_from



