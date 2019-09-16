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


def dist_force(n, m):
    nx, ny = n
    mx, my = m
    x, y = mx - nx, my - ny
    distance = hypot(x, y)
    return (normal_distance / distance - 1, *m)


def ang_force(a, b, c, optimal_ang):
    ax, ay = a
    bx, by = b
    cx, cy = c
    abx, aby = bx - ax, by - ay
    bcx, bcy = cx - bx, cy - by
    this_angle = pi - atan2(bcx, bcy)
    force = 1 * (this_angle - optimal_ang) ** 2
    return rotate_vector(force, 0, abx, aby, alpha=this_angle)


class Calculate2D:

    def calculate2d(self, initial=True):
        """
        recalculate 2d coordinates. currently rings can be calculated badly.

        :param scale: rescale calculated positions.
        :param force: ignore existing coordinates of atoms
        :param initial: force fields are free
        """
        plane = self._plane
        atoms = self._atoms
        bonds = self._bonds

        if not isinstance(initial, dict):
            force_fields = {k: (0, 0) for k in atoms}
        else:
            force_fields = initial

        # distance forces
        for n, m, bond in self.bonds():
            d, mx, my = dist_force(plane[n], plane[m])
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
                if ang:
                    optimal_angle = ang
                else:
                    optimal_angle = normal_angles[bonds[a][b].order][bonds[b][c].order]

                dx, dy = ang_force(plane[a], plane[b], plane[c], optimal_angle)
                kx, ky = force_fields[b]
                force_fields[b] = (dx + kx, dy + ky)
        self.flush_cache()
        return force_fields

    def clean2d(self):
        """
        steps for calculate 2d coordinates
        :return: None
        """
        # original force fields
        forces = self.calculate2d()

        # next steps for calculate 2d
        plane = self._plane
        t = .2
        stack = 1
        while stack:
            print('step')
            for atom, v in forces.items():
                x2, y2 = v
                x2, y2 = x2 * t, y2 * t
                x1, y1 = plane[atom]
                plane[atom] = (x1 + x2, y1 + y2)
            self._plane = plane
            stack = 0

            forces = self.calculate2d(initial=forces)
            for atom, v in forces.items():
                force = hypot(*v)
                if force > .05:
                    stack = 1


normal_distance = .825
normal_angles = {1: {1: 2/3*pi, 2: 2/3*pi, 3: 2*pi},
                 2: {1: 2/3*pi, 2: 2*pi, 3: 2*pi},
                 3: {1: 2*pi, 2: 2*pi, 3: 2*pi}}

__all__ = ['Calculate2D']
