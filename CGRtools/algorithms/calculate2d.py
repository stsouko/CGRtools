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
from itertools import combinations
from math import hypot, pi, acos, cos, sin, sqrt
from numba import njit, f8, int64, jit
from ..algorithms.depict import rotate_vector
from random import uniform


def rotate_vector2(x1, y1, angle):
    """
    rotate x,y vector over angle
    """
    cos_rad = cos(angle)
    sin_rad = sin(angle)
    return cos_rad * x1 + sin_rad * y1, -sin_rad * x1 + cos_rad * y1


# @jit(locals={'sum_sqr': f8, 'i': int64, 'g': f8, 'j': f8, 'a': int64, '_a': int64})
def distances(atoms, matrix):
    sum_sqr = 0.0
    dists = []
    i = 0
    for a in atoms:
        for _a in atoms[i + 1:]:
            for g, j in zip(matrix[a], matrix[_a]):
                sum_sqr += (int(g) - int(j)) ** 2
            dists.append([a, _a, sqrt(sum_sqr)])
    return dists


class Calculate2D:
    __slots__ = ()

    def _update(self):
        atoms = self._atoms
        bonds = self._bonds
        mapping = []
        dd = {k: {i: 1.32 if b.order == 8 or len(v) > 4 else .825 for i, b in v.items()} for k, v in bonds.items()}
        xyz_matrix = []

        # add virtual atoms
        for atm, neighbors in bonds.items():
            if len(neighbors) == 2:
                atom = atoms[atm]
                if atom.atomic_number in {5, 6, 7, 8, 14, 15, 16, 17, 33, 34, 35, 52, 53, 85}:
                    (a1, bond1), (a2, bond2) = neighbors.items()
                    order1, order2 = bond1.order, bond2.order
                    if not (order1 == order2 == 2 or order1 == 3 or order2 == 3 or order1 == 8 or order2 == 8):
                        dd[atm][-atm] = .825
                        dd[-atm] = {atm: .825}

        # add virtual bonds
        for atom, neighbors in dd.items():
            if len(neighbors) == 3:
                for n, m in combinations(neighbors, 2):
                    dd[n][m] = dd[m][n] = 1.43

        # create matrix of coordinates
        n = len(dd)
        for i, atm in enumerate(dd):
            if atm > 0:
                mapping.append((atm, i))
            xyz_matrix.append(array('f', [uniform(-n, n), uniform(-n, n), uniform(-n, n)]))
        return tuple(mapping), dd, tuple(xyz_matrix)

    @staticmethod
    # @njit(locals={'n': int64, 'm': int64, 'distance': f8, 'c_rep': f8})
    def _repulsive_force(d_matrix, adj, c_rep):
        rep_forces = {k: 0.0 for k in adj}
        for n, m, distance in d_matrix:
            n_item = rep_forces[n]
            m_item = rep_forces[m]
            diff = c_rep / distance
            rep_forces[n] = n_item + diff
            rep_forces[m] = m_item - diff
        return rep_forces

    def clean2d(self):
        """
        steps for calculate 2d coordinates
        :return: None
        """
        mapping, adj, xyz = self._update()
        dist = tuple(distances(tuple(adj), xyz))
        rep = self._repulsive_force(dist, adj, .04)
        e = 7
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
