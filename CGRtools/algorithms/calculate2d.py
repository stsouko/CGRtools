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
from numba import njit, f8, i8, i2, u2
from numpy import array, zeros, uint16
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
      {'i': u2, 'j': u2, 'xi': f8, 'yi': f8, 'zi': f8, 'xj': f8, 'yj': f8, 'zj': f8, 'c': f8, 'dx': f8,
       'dy': f8, 'dz': f8, 'fxi': f8, 'fyi': f8, 'fzi': f8, 'fxj': f8, 'fyj': f8, 'fzj': f8})
def repulsive_force(xyz, n, c_rep):
    forces = zeros((len(xyz), 3))
    for i in range(n - 1):
        xi, yi, zi = xyz[i]
        fxi, fyi, fzi = forces[i]
        for j in range(i + 1, n):
            xj, yj, zj = xyz[j]
            dx, dy, dz = xi - xj, yi - yj, zi - zj
            c = c_rep / (dx ** 2 + dy ** 2 + dz ** 2)

            # calculate repulsive force for each dimension
            dx, dy, dz = dx * c, dy * c, dz * c

            fxi += dx
            fyi += dy
            fzi += dz

            fxj, fyj, fzj = forces[j]
            forces[j] = fxj - dx, fyj - dy, fzj - dz
        forces[i] = fxi, fyi, fzi
    return forces


@njit(f8[:, :](f8[:, :], u2[:, :], f8[:], u2, f8),
      {'i': u2, 'n': u2, 'm': u2, 'r': f8, 'nx': f8, 'ny': f8, 'nz': f8, 'mx': f8, 'my': f8, 'mz': f8, 'dx': f8,
       'dy': f8, 'dz': f8, 'distance': f8, 'f': f8, 'fdx': f8, 'fdy': f8, 'fdz': f8, 'fxi': f8, 'fyi': f8, 'fzi': f8,
       'fxj': f8, 'fyj': f8, 'fzj': f8})
def spring_force(xyz, springs, springs_distances, c, c_bond):
    forces = zeros((len(xyz), 3))
    for i in range(c):
        n, m = springs[i]
        r = springs_distances[i]
        nx, ny, nz = xyz[n]
        mx, my, mz = xyz[m]

        dx, dy, dz = nx - mx, ny - my, nz - mz
        distance = sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        f = c_bond * (distance - r) / distance

        fdx, fdy, fdz = f * dx, f * dy, f * dz
        fxi, fyi, fzi = forces[n]
        fxj, fyj, fzj = forces[m]
        forces[n] = fxi - fdx, fyi - fdy, fzi - fdz
        forces[m] = fxj + fdx, fyj + fdy, fzj + fdz

    return forces


@njit(f8[:, :](f8[:, :]), {'fc': f8, 'n': u2, 'i': u2, 'x': f8, 'y': f8, 'z': f8})
def flatting(forces):
    fc = .2
    n = len(forces)
    ff = zeros((n, 3))
    for i in range(n):
        x, y, z = forces[i]
        if z > 0:
            z -= fc
        else:
            z += fc
        ff[i] = x, y, z

    return ff


class Calculate2D:
    __slots__ = ()

    def __prepare(self):
        atoms = self._atoms
        bonds = self._bonds

        mapping = {}
        springs_distances = []
        springs = []
        xyz_matrix = []

        dd = {n: {m: 1.32 for m, b in bs.items()} if len(bs) > 4 else
                 {m: 1.32 if b.order == 8 else .825 for m, b in bs.items()}
              for n, bs in bonds.items()}

        # add virtual atoms
        for n, m_bond in bonds.items():
            if len(m_bond) == 2:
                atom = atoms[n]
                if atom.atomic_number in {5, 6, 7, 8, 14, 15, 16, 17, 33, 34, 35, 52, 53, 85}:
                    (a1, bond1), (a2, bond2) = m_bond.items()
                    order1, order2 = bond1.order, bond2.order
                    if not (order1 == order2 == 2 or order1 == 3 or order2 == 3 or order1 == 8 or order2 == 8):
                        dd[n][-n] = .825
                        dd[-n] = {n: .825}

        # add virtual bonds
        v_bonds = []
        for n, m_bond in dd.items():
            if len(m_bond) == 3:
                for m1, m2 in combinations(m_bond, 2):
                    if m2 not in dd[m1]:
                        v_bonds.append((m1, m2))
        for n, m in v_bonds:
            dd[n][m] = dd[m][n] = 1.43

        # create matrix of coordinates
        cube = len(dd) / 2
        for i, n in enumerate(dd):
            mapping[n] = i
            xyz_matrix.append([uniform(-cube, cube), uniform(-cube, cube), uniform(-cube, cube)])

        # create matrices of bonds
        seen = set()
        for n, m_bond in dd.items():
            seen.add(n)
            n = mapping[n]
            for m, bond in m_bond.items():
                if m not in seen:
                    springs.append((n, mapping[m]))
                    springs_distances.append(bond)

        return mapping, array(xyz_matrix), array(springs, dtype=uint16), array(springs_distances)

    def clean2d(self):
        """
        reference for using article:
        Frączek, T. (2016). Simulation-Based Algorithm for Two-Dimensional Chemical Structure Diagram Generation of
        Complex Molecules and Ligand–Protein Interactions. Journal of Chemical Information and Modeling, 56(12),
        2320–2335. doi:10.1021/acs.jcim.6b00391 (https://doi.org/10.1021/acs.jcim.6b00391)
        """
        plane = self._plane
        count = len(self._atoms)
        mapping, xyz, springs, springs_distances = self.__prepare()
        for _ in range(5000):
            forces = repulsive_force(xyz, count, .04) + spring_force(xyz, springs, springs_distances, count, .1)
            xyz = forces + xyz

        for n in self._atoms:
            m = mapping[n]
            plane[n] = (xyz[m, 0], xyz[m, 1])


__all__ = ['Calculate2D']
