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
    """
    :param n, m: atoms n and m
    :return: force of distance
    """
    nx, ny = n
    mx, my = m
    x, y = nx - mx, ny - my
    elong = normal_distance / hypot(x, y) - 1
    return x * elong, y * elong


def ang_force(x, y, current_ang, optimal_ang):
    """
    :param x, y: vector
    :return: force of angle
    """
    force = 1 * (current_ang - optimal_ang) ** 2
    return rotate_vector(force, 0, x, y, alpha=current_ang)


class Calculate2D:
    __slots__ = ()

    def __get_forces(self):
        plane = self._plane
        atoms = self._atoms
        bonds = self._bonds

        force_fields = {k: (0, 0) for k in atoms}

        # distance forces
        for n, m, bond in self.bonds():
            dx, dy = dist_force(plane[n], plane[m])
            x, y = force_fields[n]
            force_fields[n] = (x + dx, y + dy)
            x, y = force_fields[m]
            force_fields[m] = (x - dx, y - dy)

        # angle forces
        for n, m_bond in bonds.items():
            angles = []
            neighbors = list(m_bond)
            ln = len(m_bond)
            if ln == 2:
                angles.append((neighbors[0], n, neighbors[1], False))
            elif ln == 4:
                angles.append((neighbors[0], n, neighbors[2], 2*pi))
                angles.append((neighbors[1], n, neighbors[3], 2*pi))
                angles.append((neighbors[-1], n, neighbors[0], pi))
                for v, w in zip(neighbors, neighbors[1:]):
                    angles.append((v, n, w, pi))
            elif ln != 1:
                ln = len(m_bond)
                angles.append((neighbors[-1], n, neighbors[0], 4*pi/ln))
                for v, w in zip(neighbors, neighbors[1:]):
                    angles.append((v, n, w, 4*pi/ln))

            for a, b, c, ang in angles:
                ax, ay = plane[a]
                bx, by = plane[b]
                cx, cy = plane[c]
                abx, aby = bx - ax, by - ay
                bcx, bcy = cx - bx, cy - by
                angle = pi - calculate_angle((abx, aby), (bcx, bcy))
                if ang:
                    optimal_angle = ang
                else:
                    optimal_angle = normal_angles[bonds[a][b].order][bonds[b][c].order]

                if angle < .01:
                    dx, dy = rotate_vector(0.2, .0, abx, aby)
                    fax, fay = force_fields[a]
                    force_fields[a] = fax + dx, fay + dy
                    fcx, fcy = force_fields[c]
                    force_fields[a] = fcx - dx, fcy - dy
                    # cx, cy = cx + dx, cy + dy
                    # plane[c] = (cx, cy)
                    # bcx, bcy = cx - bx, cy - by
                    # angle = pi - calculate_angle((abx, aby), (bcx, bcy))

                dx, dy = ang_force(abx, aby, angle, optimal_angle)
                kx, ky = force_fields[b]
                force_fields[b] = (dx + kx, dy + ky)
        self.flush_cache()
        return force_fields

    def clean2d(self):
        """
        steps for calculate 2d coordinates
        :return: None
        """
        plane = self._plane
        t = .2
        stack = 1
        while stack:
            stack = 0
            forces = self.__get_forces()
            for atom, v in forces.items():
                x2, y2 = v
                if hypot(x2, y2) > .05:
                    stack = 1
                x2, y2 = x2 * t, y2 * t
                x1, y1 = plane[atom]
                plane[atom] = (x1 + x2, y1 + y2)
            self._plane = plane


normal_distance = .825
normal_angles = {1: {1: 2/3*pi, 2: 2/3*pi, 3: 2*pi},
                 2: {1: 2/3*pi, 2: 2*pi, 3: 2*pi},
                 3: {1: 2*pi, 2: 2*pi, 3: 2*pi}}

__all__ = ['Calculate2D']
