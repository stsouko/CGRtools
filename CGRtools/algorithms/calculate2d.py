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
from CachedMethods import cached_property
from collections import defaultdict
from itertools import combinations
from math import hypot, pi, acos, cos, sin, sqrt
from numba import njit, f8, i8, jit, u2
from numpy import array, zeros
from ..algorithms.depict import rotate_vector
from random import uniform


def rotate_vector2(x1, y1, angle):
    """
    rotate x,y vector over angle
    """
    cos_rad = cos(angle)
    sin_rad = sin(angle)
    return cos_rad * x1 + sin_rad * y1, -sin_rad * x1 + cos_rad * y1


@njit(f8[:, :](f8[:, :], u2, f8),
      {'i': u2, 'j': u2, 'n': u2, 'xi': f8, 'yi': f8, 'zi': f8, 'xj': f8, 'yj': f8, 'zj': f8, 'c': f8, 'dx': f8,
       'dy': f8, 'dz': f8, 'ni': f8, 'mi': f8, 'oi': f8, 'nj': f8, 'mj': f8, 'oj': f8})
def repulsive_force(matrix, n, c_rep):
    forces = zeros((n, 3))
    for i in range(n - 1):
        for j in range(i + 1, n):
            xi, yi, zi = matrix[i]
            xj, yj, zj = matrix[j]
            dx, dy, dz = xi - xj, yi - yj, zi - zj
            c = c_rep / (dx ** 2 + dy ** 2 + dz ** 2)

            # calculate repulsive force for each dimension
            rfx, rfy, rfz = dx * c, dy * c, dz * c
            ni, mi, oi = forces[i]
            nj, mj, oj = forces[j]
            forces[i] = ni + rfx, i + rfy, oi + rfz
            forces[j] = nj - rfx, mj - rfy, oj - rfz
    return forces


@njit(locals={'n': i8, 'm': i8, 'distance': f8, 'c_rep': f8})
def spring_force(adj, xyz, mapping, c_bond):
    forces = zeros((len(adj), 3))
    for i, v in adj.items():
        for j, r in v.items():
            n, m = mapping[i], mapping[j]
            nx, ny, nz = xyz[n]
            mx, my, mz = xyz[m]
            f = c_bond * (sqrt((nx - mx) ** 2 + (ny - my) ** 2 + (nz - mz) ** 2) - r)

            xi, yi, zi = forces[n]
            xj, yj, zj = forces[m]
            forces[n] = xi + f * nx, yi + f * ny, zi + f * nz
            forces[m] = xj + f * mx, yj + f * my, zj + f * mz
    return forces


class Calculate2D:
    __slots__ = ()

    def __prepare(self):
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
        springs = []
        for atom, neighbors in dd.items():
            if len(neighbors) == 3:
                for n, m in combinations(neighbors, 2):
                    springs.append((n, m))
                    springs.append((m, n))
        for i, j in springs:
            dd[i][j] = dd[j][i] = 1.43

        # create matrix of coordinates
        n = len(dd) / 2
        for i, atm in enumerate(dd):
            mapping.append((atm, i))
            xyz_matrix.append([uniform(-n, n), uniform(-n, n), uniform(-n, n)])
        return tuple(mapping), dd, array(xyz_matrix)

    def clean2d(self):
        """
        steps for calculate 2d coordinates
        :return: None
        """
        mapping, adj, xyz = self.__prepare()
        forces = repulsive_force(xyz, len(self._atoms), .04)
        rep = spring_force(adj, xyz, mapping, .1)
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
