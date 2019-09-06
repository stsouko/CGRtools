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
from random import uniform
from itertools import combinations, chain
from numpy import dot, linalg, zeros
from math import hypot, pi, atan2, cos, sin


def rotate_vector(x1, y1, x2, y2, for_angle=None):
    """
    rotate x,y vector over x2-x1, y2-y1 angle
    """
    alpha = 0
    if for_angle:
        alpha = for_angle / 2
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

        s = len(atoms)
        matrix_of_force_fields = zeros((2, s))
        dist = {}
        # length forces
        for n, m_bond in bonds.items():
            dist[n] = {}
            for m in m_bond:
                dist[n][m] = .825

        # angle forces
        for n, m_bond in bonds.items():
            if len(m_bond) == 2:  # single-single or single-double bonds has angle = 120, other 180
                (m1, b1), (m2, b2) = m_bond.items()
                dist[m1][m2] = dist[m2][m1] = 1.43 if b1.order + b2.order in (2, 3) else 1.7  # +.05
            elif len(m_bond) == 3:
                m1, m2, m3 = m_bond
                dist[m1][m2] = dist[m1][m3] = dist[m2][m3] = dist[m3][m2] = dist[m2][m1] = dist[m3][m1] = 1.43
            elif len(m_bond) == 4:
                #    1
                #
                # 2  X  4
                #
                #    3
                m1, m2, m3, m4 = m_bond
                dist[m1][m2] = dist[m1][m4] = dist[m2][m1] = dist[m2][m3] = 1.17
                dist[m3][m2] = dist[m3][m4] = dist[m4][m1] = dist[m4][m3] = 1.17
                dist[m1][m3] = dist[m3][m1] = dist[m2][m4] = dist[m4][m2] = 1.7  # +.05

        # cycle forces
        for r in self.sssr:
            if len(r) == 6:
                #    6
                #
                # 1     5
                #
                # 2     4
                #
                #    3
                m1, m2, m3, m4, m5, m6 = r
                dist[m1][m4] = dist[m4][m1] = dist[m2][m5] = dist[m5][m2] = dist[m3][m6] = dist[m6][m3] = 1.7  # +.05

        # distance forces
        for n, m, bond in self.bonds():
            nx, ny = plane[n]
            mx, my = plane[m]
            x, y = mx - nx, my - ny
            distance = hypot(x, y)
            difference = dist[n][m] - distance
            matrix_of_force_fields[0][n - 1], matrix_of_force_fields[1][n - 1] = \
                rotate_vector(0, difference, x, y)

        m, n = n, m
        nx, ny = plane[n]
        mx, my = plane[m]
        x, y = mx - nx, my - ny
        distance = hypot(x, y)
        difference = dist[n][m] - distance
        matrix_of_force_fields[0][n - 1], matrix_of_force_fields[1][n - 1] = \
            rotate_vector(0, -difference, x, y)

        # angle forces
        for n, m_bond in bonds.items():
            angles = []
            neighbors = list(m_bond.keys())
            if len(m_bond) == 2:
                for i, p in enumerate(neighbors):
                    try:
                        angles.append((p, n, neighbors[i + 1]))
                    except:
                        pass
            elif len(m_bond) > 2:
                for i, p in enumerate(neighbors):
                    try:
                        angles.append((p, n, neighbors[i + 1]))
                    except:
                        angles.append((neighbors[-1], n, neighbors[0]))

            for a, b, c in angles:
                ax, ay = plane[a]
                bx, by = plane[b]
                cx, cy = plane[c]
                abx, aby = bx - ax, by - ay
                bcx, bcy = cx - bx, cy - by
                this_angle = calculate_angle((abx, aby), (bcx, bcy))

                Vangle = 1 * (this_angle - normal_angles[bonds[a][b].order][bonds[b][c].order]) ** 2
                dx, dy = rotate_vector(Vangle, 0, abx, aby, for_angle=this_angle)
                kx, ky = matrix_of_force_fields[0][n - 1], matrix_of_force_fields[1][n - 1]
                matrix_of_force_fields[0][n - 1], matrix_of_force_fields[1][n - 1] = dx + kx, dy + ky

        if force:
            pos = None
        else:
            pos = {n: (plane[n][0] or uniform(0, .01), plane[n][1] or uniform(0, .01)) for n, atom in atoms()}

        for n, xy in kamada_kawai_layout(self, dist=dict(dist), pos=pos, scale=scale).items():
            x, y = plane[n]

        self.flush_cache()


normal_angles = {1: {1: 2/3*pi, 2: 2/3*pi, 3: 2*pi},
                 2: {1: 2/3*pi, 2: 2*pi, 3: 2*pi},
                 3: {1: 2*pi, 2: 2*pi, 3: 2*pi}}

__all__ = ['Calculate2D']
