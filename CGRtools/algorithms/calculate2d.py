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
from math import hypot, pi, atan2, cos, sin, acos


def rotate_vector(x1, y1, angle):
    cos_rad = cos(angle)
    sin_rad = sin(angle)
    return cos_rad * x1 + sin_rad * y1, sin_rad * x1 - cos_rad * y1


def calculate_angle(v1, v2):
    """
    :params v1, v2: input vectors
    :return: angle between vectors v1, v2 in radians
    """
    x1, y1 = v1
    x2, y2 = v2
    return abs(atan2(y1, x1) - atan2(y2, x2))


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


class Calculate2D:
    __slots__ = ()

    def __get_forces(self):
        plane = self._plane
        atoms = self._atoms
        bonds = self._bonds
        sssr = self.sssr

        force_fields = {k: (0, 0) for k in atoms}

        # distance forces
        for n, m, bond in self.bonds():
            dx, dy = dist_force(plane[n], plane[m])
            x, y = force_fields[n]
            force_fields[n] = (x + dx, y + dy)
            x, y = force_fields[m]
            force_fields[m] = (x - dx, y - dy)

        # angle forces
        angles = []
        for n, m_bond in bonds.items():
            neighbors = list(m_bond)
            ln = len(m_bond)
            cycle_angle = [pi / 2 - pi / len(cycle) for cycle in sssr if n in cycle]
            if ln == 2:
                a = neighbors[0]
                c = neighbors[1]
                ange = normal_angles[bonds[a][n].order][bonds[n][c].order]
                if cycle_angle:
                    ange = cycle_angle[0]
                angles.append((a, n, c, ange))
            elif ln == 4:
                angles.append((neighbors[0], n, neighbors[2], 2*pi))
                angles.append((neighbors[1], n, neighbors[3], 2*pi))
                angles.append((neighbors[-1], n, neighbors[0], pi))
                for v, w in zip(neighbors, neighbors[1:]):
                    angles.append((v, n, w, pi))
            elif ln != 1:
                ange = 4*pi/ln
                if cycle_angle:
                    ange = cycle_angle[0]
                angles.append((neighbors[-1], n, neighbors[0], ange))
                for v, w in zip(neighbors, neighbors[1:]):
                    angles.append((v, n, w, ange))

        #
        # a
        # |
        # b---c
        #
        for a, b, c, opt in angles:
            ax, ay = plane[a]
            bx, by = plane[b]
            cx, cy = plane[c]
            bax, bay = ax - bx, ay - by
            bcx, bcy = cx - bx, cy - by
            angle_ba = atan2(bay, bax)
            angle_bc = atan2(bcy, bcx)
            angle = acos((bax*bcx + bay*bcy) / (hypot(bay, bax) * hypot(bcy, bcx)))
            angle_biss = (angle_bc + angle_ba) / 2
            print('angle', angle, angle_biss, angle_ba, angle_bc)

            # if angle < .01:
            #     dx, dy = rotate_vector(0.2, .0, abx, aby)
            #     fax, fay = force_fields[a]
            #     force_fields[a] = fax + dx, fay + dy
            #     fcx, fcy = force_fields[c]
            #     force_fields[a] = fcx - dx, fcy - dy

            dx, dy = rotate_vector(1 * (angle - opt) ** 2, 0, angle_biss)
            print(dx, dy, 1 * (angle - opt) ** 2)
            kx, ky = force_fields[b]
            # if b == 2:
            #     print(dx + kx, dy + ky)
            force_fields[b] = (kx - dx, ky - dy)
        return force_fields

    def clean2d(self):
        """
        steps for calculate 2d coordinates
        :return: None
        """
        plane = self._plane
        for x in range(1):
            forces = self.__get_forces()
            print(forces, '\n')
            for atom, v in forces.items():
                x1, y1 = plane[atom]
                x2, y2 = v
                force = hypot(x2, y2)
                ratio = .2 / force
                if ratio < 1:
                    x2, y2 = x2 * ratio, y2 * ratio
                plane[atom] = (x1 + x2, y1 + y2)
            # if steps >= 100:
            #     break


normal_distance = .825
normal_angles = {1: {1: 2/3*pi, 2: 2/3*pi, 3: 2*pi},
                 2: {1: 2/3*pi, 2: 2*pi, 3: 2*pi},
                 3: {1: 2*pi, 2: 2*pi, 3: 2*pi}}

__all__ = ['Calculate2D']
