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
from CachedMethods import cached_property
from math import hypot, pi, acos
from ..algorithms.depict import rotate_vector


def hypot2(x ,y):
    """
    x and y are zeroes so hypot(x, y) is zeroes
    :param x, y: vector x, y
    :return: vectors len
    """
    lnm = hypot(x, y)
    if not lnm:
        lnm = .01
    return lnm


class Calculate2D:
    __slots__ = ()

    @cached_property
    def __angles(self):
        bonds = self._bonds
        sssr = self.sssr
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
                angles.append((neighbors[0], n, neighbors[2], 2 * pi))
                angles.append((neighbors[1], n, neighbors[3], 2 * pi))
                angles.append((neighbors[-1], n, neighbors[0], pi))
                for v, w in zip(neighbors, neighbors[1:]):
                    angles.append((v, n, w, pi))
            elif ln != 1:
                ange = 4 * pi / ln
                if cycle_angle:
                    ange = cycle_angle[0]
                angles.append((neighbors[-1], n, neighbors[0], ange))
                for v, w in zip(neighbors, neighbors[1:]):
                    angles.append((v, n, w, ange))
        return angles

    def __get_forces(self):
        plane = self._plane
        atoms = self._atoms

        force_fields = {k: (0, 0) for k in atoms}

        # distance forces
        for n, m, bond in self.bonds():
            nx, ny = plane[n]
            mx, my = plane[m]
            x, y = nx - mx, ny - my
            # x,y = 0,0 then hypot(y, x) = 0
            elong = normal_distance / hypot2(y, x) / 2 - .5
            dx, dy = x * elong, y * elong
            x, y = force_fields[n]
            force_fields[n] = (x + dx, y + dy)
            x, y = force_fields[m]
            force_fields[m] = (x - dx, y - dy)

        # angle forces
        #
        # a
        # |
        # b---c
        #
        for a, b, c, opt in self.__angles:
            ax, ay = plane[a]
            bx, by = plane[b]
            cx, cy = plane[c]
            bax, bay = ax - bx, ay - by
            bcx, bcy = cx - bx, cy - by

            l_ba, l_bc = hypot2(bay, bax), hypot2(bcy, bcx)
            bax /= l_ba
            bay /= l_ba
            bcx /= l_bc
            bcy /= l_bc

            bis_x = bax + bcx
            bis_y = bay + bcy
            bis_l = hypot2(bis_y, bis_x)

            angle = acos(bax * bcx + bay * bcy)
            if angle < .01:
                dx, dy = rotate_vector(0, .2, bcx, bcy)
                fax, fay = force_fields[a]
                force_fields[a] = fax + dx, fay + dy
                fcx, fcy = force_fields[c]
                force_fields[c] = fcx - dx, fcy - dy
            force = 2 * (opt - angle) * abs(angle - opt)

            ratio = force / bis_l
            bis_x *= ratio
            bis_y *= ratio

            kx, ky = force_fields[b]
            force_fields[b] = (kx + bis_x, ky + bis_y)
        return force_fields

    def clean2d(self):
        """
        steps for calculate 2d coordinates
        :return: None
        """
        plane = self._plane
        stack = 1
        steps = 1
        while stack:
            stack = 0
            forces = self.__get_forces()
            force = max(hypot2(x, y) for x, y in forces.values())
            ratio = .2 / force
            if ratio > 1:
                ratio = 1
            for atom, (x2, y2) in forces.items():
                x1, y1 = plane[atom]
                x2, y2 = x2 * ratio, y2 * ratio
                plane[atom] = (x1 + x2, y1 + y2)
            steps += 1
            if force > .05:
                stack = 1
            if steps >= 100:
                break


normal_distance = .825
normal_angles = {1: {1: 2/3*pi, 2: 2/3*pi, 3: 2*pi},
                 2: {1: 2/3*pi, 2: 2*pi, 3: 2*pi},
                 3: {1: 2*pi, 2: 2*pi, 3: 2*pi}}

__all__ = ['Calculate2D']
