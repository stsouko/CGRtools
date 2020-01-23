# -*- coding: utf-8 -*-
#
#  Copyright 2019 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from array import array
from CachedMethods import cached_property
from collections import defaultdict
from math import hypot, pi, acos, cos, sin, sqrt
from ..algorithms.depict import rotate_vector
from random import uniform


def rotate_vector2(x1, y1, angle):
    """
    rotate x,y vector over angle
    """
    cos_rad = cos(angle)
    sin_rad = sin(angle)
    return cos_rad * x1 + sin_rad * y1, -sin_rad * x1 + cos_rad * y1


class Calculate2D:
    __slots__ = ()

    def _update(self):
        atoms = self._atoms
        bonds = list(self.bonds())
        mapping = {}
        dd = defaultdict(dict)
        update_atoms = list(atoms)
        xyz_matrix = []

        # add virtual atoms
        for i, (n, m, bond) in enumerate(bonds):
            order1 = bond.order
            dd[n][m] = dd[m][n] = 8.0 if order1 == 8 else 5.0
            for _n, _m, _bond in bonds[i + 1:]:
                order2 = _bond.order
                dd[_n][_m] = dd[_m][_n] = 8.0 if order2 == 8 else 5.0
                common = {n, m}.intersection({_n, _m})
                if common:
                    atom_number = list(common)[0]
                    atom = atoms[atom_number]
                    if atom.neighbors == 2 and \
                            atom.atomic_number in (5, 6, 7, 8, 14, 15, 16, 17, 33, 34, 35, 52, 53, 85) and not \
                            (order1 == 1 and order2 == 3 or order1 == 3 and order2 == 1 or order1 == 2 and order2 == 2):
                        dd[atom_number][0] = 5.0
                        update_atoms.append(0)

        # create matrix of coordinates
        n = len(update_atoms)
        n = n / 2
        ind = 0
        for atm in update_atoms:
            if atm:
                mapping[atm] = ind
            xyz_matrix.append(array('f', [uniform(-n, n), uniform(-n, n), uniform(-n, n)]))
            ind += 1
        return update_atoms, mapping, dd, tuple(xyz_matrix)

    @staticmethod
    def _distances(atoms, matrix):
        sum_sqr = 0
        dists = []
        for i, a in enumerate(atoms):
            for _a in atoms[i + 1:]:
                for g, j in zip(matrix[a], matrix[_a]):
                    sum_sqr += (int(g) - int(j)) ** 2
                dists.append([a, _a, sqrt(sum_sqr)])
        return dists

    def _repulsive_force(self):
        pass

    def clean2d(self):
        """
        steps for calculate 2d coordinates
        :return: None
        """
        up_atoms, mapping, adj, xyz = self._update()
        dist = self._distances(up_atoms, xyz)
        # stack = 1
        # steps = 1
        # primary_forces, dif_dist = self.__get_dist_forces()
        # secondary_forces, dif_ang = self.__get_angle_forces()
        # print(dif_dist, dif_ang)
        # if dif_ang > dif_dist:
        #     primary_forces, secondary_forces = secondary_forces, primary_forces
        #
        # while stack:
        #     # for x in range(1):
        #     stack = 0
        #     self._changes(primary_forces)
        #     self._changes(secondary_forces)
        #
        #     primary_forces, dif_dist = self.__get_dist_forces()
        #     secondary_forces, dif_ang = self.__get_angle_forces()
        #     if dif_ang > dif_dist:
        #         primary_forces, secondary_forces = secondary_forces, primary_forces
        #     force_p = max(hypot(x, y) for x, y in primary_forces.values())
        #     force_s = max(hypot(x, y) for x, y in secondary_forces.values())
        #     print(force_p, force_s)
        #     if force_p < .05 and force_s < 0.5:
        #         stack = 0
        #
        #     if steps >= 200:
        #         break
        #     steps += 1


__all__ = ['Calculate2D']
