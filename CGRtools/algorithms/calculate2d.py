# -*- coding: utf-8 -*-
#
#  Copyright 2019 Ramil Nugmanov <stsouko@live.ru>
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
# from networkx.drawing.layout import kamada_kawai_layout
from math import hypot, pi, atan2, cos, sin
from numpy import linalg, dot
from random import uniform


def rotate_vector(x1, y1, x2, y2, alpha=.0):
    """
    rotate x,y vector over x2-x1, y2-y1 angle
    """
    if alpha:
        alpha = alpha / 2
    angle = atan2(y2, x2) + alpha
    cos_rad = cos(angle)
    sin_rad = sin(angle)
    return cos_rad * x1 + sin_rad * y1, -sin_rad * x1 + cos_rad * y1


def calculate_angle(v1, v2):
    '''
    :params v1, v2: input vectors
    :return: angle between vectors v1, v2 in radians
    '''
    return atan2(linalg.det([v1, v2]), dot(v1, v2))


class Calculate2D:
    def calculate2d(self, force=False, scale=1):
        """
        recalculate 2d coordinates. currently rings can be calculated badly.

        :param scale: rescale calculated positions.
        :param force: ignore existing coordinates of atoms
        """
        plane = self._plane
        atoms = self._atoms
        bonds = self._bonds

        force_fields = {k: (0, 0) for k in atoms}
        dist = {}

        # distance forces
        for n, m, bond in self.bonds():
            nx, ny = plane[n]
            mx, my = plane[m]
            x, y = mx - nx, my - ny
            distance = hypot(x, y)
            difference = normal_distance - distance
            d = normal_distance/distance - 1
            dx, dy = force_fields[n]
            force_fields[n] = (dx + mx*d, dy + my*d)
            dx, dy = force_fields[m]
            force_fields[m] = (dx - mx*d, dy - my*d)

        # angle forces
        for n, m_bond in bonds.items():
            angles = []
            neighbors = list(m_bond)
            if len(m_bond) == 2:
                angles.append((neighbors[0], n, neighbors[1], False))
            elif len(m_bond) == 4:
                angles.append((neighbors[0], n, neighbors[2], 2*pi))
                angles.append((neighbors[1], n, neighbors[3], 2*pi))
                angles.append((neighbors[-1], n, neighbors[0], pi))
                for v, w in zip(neighbors, neighbors[1:]):
                    angles.append((v, n, w, pi))
            elif len(m_bond) in (3, 5, 6, 7):
                ll = len(m_bond)
                angles.append((neighbors[-1], n, neighbors[0], 4*pi/ll))
                for v, w in zip(neighbors, neighbors[1:]):
                    angles.append((v, n, w, 4*pi/ll))

            for a, b, c, ang in angles:
                ax, ay = plane[a]
                bx, by = plane[b]
                cx, cy = plane[c]
                abx, aby = bx - ax, by - ay
                bcx, bcy = cx - bx, cy - by
                this_angle = calculate_angle((abx, aby), (bcx, bcy))
                if ang:
                    optimal_angle = ang
                else:
                    optimal_angle = normal_angles[bonds[a][b].order][bonds[b][c].order]

                force = 1 * (this_angle - optimal_angle) ** 2
                dx, dy = rotate_vector(force, 0, abx, aby, alpha=this_angle)
                kx, ky = force_fields[b]
                force_fields[b] = (dx + kx, dy + ky)

        self.flush_cache()


normal_distance = .825
normal_angles = {1: {1: 2/3*pi, 2: 2/3*pi, 3: 2*pi},
                 2: {1: 2/3*pi, 2: 2*pi, 3: 2*pi},
                 3: {1: 2*pi, 2: 2*pi, 3: 2*pi}}

__all__ = ['Calculate2D']
