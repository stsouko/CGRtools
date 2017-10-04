# -*- coding: utf-8 -*-
#
#  Copyright 2017 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
from importlib.util import find_spec


def fix_stereo(g):
    g.reset_query_marks()
    weights = g.get_morgan()

    for atom, attr in g.nodes(data=True):
        if attr['element'] in ('C', 'Si'):
            if attr['s_neighbors'] in (3, 4) and attr['s_hyb'] == 1:
                pass

    return g


def pyramid_volume(n, u, v, w):
    res = {}
    for d, kx, ky, kz in (('s', 's_x', 's_y', 's_z'), ('p', 'p_x', 'p_y', 'p_z')):
        zx, zy, zz = n[kx], n[ky], n[kz]

        ux, uy, uz = u['s_x'] - zx, u['s_y'] - zy, u['s_z'] - zz
        vx, vy, vz = v['s_x'] - zx, v['s_y'] - zy, v['s_z'] - zz
        wx, wy, wz = w['s_x'] - zx, w['s_y'] - zy, w['s_z'] - zz

        res[d] = ux * (vy * wz - vz * wy) + \
            uy * (vz * wx - vx * wz) + \
            uz * (vx * wy - vy * wx)  # 14 operations / det

    return res


if find_spec('numba'):  # jitting if available
    from numba import njit
    pyramid_volume = njit(pyramid_volume)
