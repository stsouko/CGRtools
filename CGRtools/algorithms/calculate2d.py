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
from collections import defaultdict
from itertools import combinations
from math import sqrt
from numba import njit, f8, u2
from numpy import array, zeros, uint16, zeros_like
from random import uniform


@njit(f8[:, :](f8[:, :], f8, f8, f8),
      {'n': u2, 'i': u2, 'j': u2, 'xi': f8, 'yi': f8, 'zi': f8, 'xj': f8, 'yj': f8, 'zj': f8, 'c': f8, 'dx': f8,
       'dy': f8, 'dz': f8, 'fxi': f8, 'fyi': f8, 'fzi': f8, 'fxj': f8, 'fyj': f8, 'fzj': f8, 'distance': f8})
def repulsive_force(xyz, coef, cut, inf):
    forces = zeros_like(xyz)
    n = len(xyz)
    inf **= 2
    cut **= 2

    for i in range(n - 1):
        xi, yi, zi = xyz[i]
        fxi, fyi, fzi = forces[i]
        for j in range(i + 1, n):
            xj, yj, zj = xyz[j]
            dx, dy, dz = xi - xj, yi - yj, zi - zj

            distance = dx ** 2 + dy ** 2 + dz ** 2
            if distance < cut:
                c = coef / cut
            elif distance > inf:
                c = 0.
            else:
                c = coef / distance

            # calculate repulsive force for each dimension
            dx, dy, dz = dx * c, dy * c, dz * c

            fxi += dx
            fyi += dy
            fzi += dz

            fxj, fyj, fzj = forces[j]
            forces[j] = fxj - dx, fyj - dy, fzj - dz
        forces[i] = fxi, fyi, fzi
    return forces


@njit(f8[:, :](f8[:, :], u2[:, :], f8[:], f8),
      {'i': u2, 'n': u2, 'm': u2, 'r': f8, 'nx': f8, 'ny': f8, 'nz': f8, 'mx': f8, 'my': f8, 'mz': f8, 'dx': f8,
       'dy': f8, 'dz': f8, 'distance': f8, 'f': f8, 'fxi': f8, 'fyi': f8, 'fzi': f8, 'fdx': f8, 'fdy': f8, 'fdz': f8,
       'fxj': f8, 'fyj': f8, 'fzj': f8})
def spring_force(xyz, springs, springs_distances, coef):
    forces = zeros((len(xyz), 3))
    for i in range(len(springs)):
        n, m = springs[i]
        r = springs_distances[i]
        nx, ny, nz = xyz[n]
        mx, my, mz = xyz[m]

        dx, dy, dz = nx - mx, ny - my, nz - mz
        distance = sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        f = coef * (distance - r)

        if -.001 < distance < .001:
            fdx, fdy, fdz = f, .0, .0
        else:
            f /= distance
            fdx, fdy, fdz = f * dx, f * dy, f * dz

        fxi, fyi, fzi = forces[n]
        fxj, fyj, fzj = forces[m]
        forces[n] = fxi - fdx, fyi - fdy, fzi - fdz
        forces[m] = fxj + fdx, fyj + fdy, fzj + fdz

    return forces


@njit(f8[:, :](f8[:, :], f8[:, :], f8),
      {'n': u2, 'ff': f8[:, :], 'i': u2, 'x': f8, 'y': f8, 'z': f8, 'fx': f8, 'fy': f8, 'fz': f8})
def flattening(forces, xyz, fc):
    n = len(forces)
    ff = zeros_like(xyz)
    for i in range(n):
        x, y, z = xyz[i]
        fx, fy, fz = forces[i]
        if -fc <= z <= fc:
            fz = -z
        else:
            if fz > 0 < z:
                fz -= fc
            elif fz < 0 > z:
                fz += fc
        ff[i] = fx, fy, fz
    return ff


@njit(f8[:, :](f8[:, :], f8), {'n': u2, 'dx': f8, 'dy': f8, 'dz': f8, 'distance': f8, 'f': f8})
def cutoff(forces, cut):
    for n in range(len(forces)):
        dx, dy, dz = forces[n]
        distance = sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        if distance > cut:
            f = cut / distance
            forces[n] = f * dx, f * dy, f * dz
    return forces


@njit(f8[:, :](f8[:, :], u2[:, :], f8[:]))
def steps(xyz, springs, springs_distances):
    # step 1
    for _ in range(1000):
        r_forces = repulsive_force(xyz, .05, .2, 10)
        s_forces = spring_force(xyz, springs, springs_distances, .2)
        forces = r_forces + s_forces
        forces = cutoff(forces, .3)
        xyz = forces + xyz

    # step 2
    for _ in range(1000):
        r_forces = repulsive_force(xyz, .02, .3, 5)
        s_forces = spring_force(xyz, springs, springs_distances, .2)
        forces = r_forces + s_forces
        forces = flattening(forces, xyz, .2)
        forces = cutoff(forces, .3)
        xyz = forces + xyz

    # step 3
    for _ in range(1000):
        r_forces = repulsive_force(xyz, .005, .4, 2)
        s_forces = spring_force(xyz, springs, springs_distances, .2)
        forces = r_forces + s_forces
        forces = flattening(forces, xyz, .2)
        forces = cutoff(forces, .3)
        xyz = forces + xyz

    return xyz


class Calculate2D:
    __slots__ = ()

    def __prepare(self):
        atoms = self._atoms
        bonds = self._bonds

        mapping = {}
        springs_distances = []
        springs = []
        xyz_matrix = []
        dd = defaultdict(list)

        # create matrix of coordinates
        cube = len(atoms)
        for i, n in enumerate(atoms):
            mapping[n] = i
            xyz_matrix.append([uniform(-cube, cube), uniform(-cube, cube), uniform(-cube, cube)])

        # create matrix of connecting, springs and springs distances
        for n, m_bond in bonds.items():
            i = mapping[n]
            for m, b in m_bond.items():
                dd[n].append(m)
                j = mapping[m]
                if (j, i) not in springs:
                    springs.append((i, j))
                    if b.order == 8 or len(m_bond) > 4:
                        springs_distances.append(1.32)
                    else:
                        springs_distances.append(.825)

        # add virtual atoms and complement matrices of bonds and coordinates
        end = len(atoms)
        for n, m_bond in bonds.items():
            if len(m_bond) == 2:
                atom = atoms[n]
                if atom.atomic_number in {5, 6, 7, 8, 14, 15, 16, 17, 33, 34, 35, 52, 53, 85}:
                    (a1, bond1), (a2, bond2) = m_bond.items()
                    if self._is_angle(bond1, bond2):
                        dd[n].append(-n)
                        dd[-n].append(n)
                        mapping[-n] = end
                        xyz_matrix.append([uniform(-cube, cube), uniform(-cube, cube), uniform(-cube, cube)])
                        springs.append((mapping[n], end))
                        springs_distances.append(.825)
                        end += 1

        # add virtual bonds and complement matrices of bonds and coordinates
        for n, m_bond in dd.items():
            if len(m_bond) == 3:
                for m1, m2 in combinations(m_bond, 2):
                    if m2 not in dd[m1]:
                        i, j = mapping[m1], mapping[m2]
                        if (j, i) not in springs or (i, j) not in springs:
                            springs.append((i, j))
                            springs_distances.append(1.43)

        return array(xyz_matrix), array(list(springs), dtype=uint16), array(springs_distances)

    def _is_angle(self, bond1, bond2):
        pass

    def clean2d(self):
        """
        reference for using article:
        Frączek, T. (2016). Simulation-Based Algorithm for Two-Dimensional Chemical Structure Diagram Generation of
        Complex Molecules and Ligand–Protein Interactions. Journal of Chemical Information and Modeling, 56(12),
        2320–2335. doi:10.1021/acs.jcim.6b00391 (https://doi.org/10.1021/acs.jcim.6b00391)
        """
        plane = self._plane
        atoms = self._atoms

        xyz, springs, springs_distances = self.__prepare()
        xyz = steps(xyz, springs, springs_distances)
        for m, n in enumerate(atoms):
            plane[n] = (xyz[m, 0], xyz[m, 1])
        self.__dict__.pop('__cached_method__repr_svg_', None)


class Calculate2DMolecule(Calculate2D):
    __slots__ = ()

    def _is_angle(self, bond1, bond2):
        order1, order2 = bond1.order, bond2.order
        return not (order1 == order2 == 2 or order1 == 3 or order2 == 3 or order1 == 8 or order2 == 8)


class Calculate2DCGR(Calculate2D):
    __slots__ = ()

    @staticmethod
    def __primary(bond):
        w_order, p_order = bond.order, bond.p_order
        if w_order is None:
            order = p_order
        elif p_order is None:
            order = w_order
        else:
            if w_order == 8:
                order = p_order
            elif p_order == 8:
                order = w_order
            else:
                if w_order > p_order:
                    order = w_order
                else:
                    order = p_order
        return order

    def _is_angle(self, bond1, bond2):
        order1 = self.__primary(bond1)
        order2 = self.__primary(bond2)
        return not (order1 == order2 == 2 or order1 == 3 or order2 == 3 or order1 == 8 or order2 == 8)


__all__ = ['Calculate2DMolecule', 'Calculate2DCGR']
