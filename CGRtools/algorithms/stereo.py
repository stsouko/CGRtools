# -*- coding: utf-8 -*-
#
#  Copyright 2017 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from importlib.util import find_spec


def pyramid_volume(n, u, v, w):
    nx, ny, nz = n[0], n[1], n[2]
    ux, uy, uz = u[0] - nx, u[1] - ny, u[2] - nz
    vx, vy, vz = v[0] - nx, v[1] - ny, v[2] - nz
    wx, wy, wz = w[0] - nx, w[1] - ny, w[2] - nz

    return ux * (vy * wz - vz * wy) + uy * (vz * wx - vx * wz) + uz * (vx * wy - vy * wx)


if find_spec('numba'):  # jitting if available
    from numba import njit
    pyramid_volume = njit(pyramid_volume)
