# -*- coding: utf-8 -*-
#
#  Copyright 2019, 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019, 2020 Dinar Batyrshin <batyrshin-dinar@mail.ru>
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
from importlib.util import find_spec
from itertools import combinations
from math import sqrt, pi, atan2, cos, sin
from random import uniform

if find_spec('numpy') and find_spec('numba'):  # try to load numba jit
    from numpy import array, zeros, uint16, zeros_like, empty
    from numba import b1, njit, f8, u2
else:
    def njit(*args, **kwargs):
        def wrapper(f):
            return f
        return wrapper

    class NumbaType:
        def __getitem__(self, item):
            return self

        def __call__(self, *args, **kwargs):
            return self

    b1 = f8 = u2 = NumbaType()


@njit(f8[:, :](f8[:, :], f8, f8, f8), cache=True)
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


@njit(f8[:, :](f8[:, :], u2[:, :], f8, f8, f8), cache=True)
def straight_repulsive(xyz, straights, coef, cut, inf):
    forces = zeros_like(xyz)
    inf **= 2
    cut **= 2

    for n in range(len(straights)):
        i, j = straights[n]
        xi, yi, zi = xyz[i]
        xj, yj, zj = xyz[j]
        dx, dy, dz = xi - xj, yi - yj, zi - zj

        distance = dx ** 2 + dy ** 2 + dz ** 2
        if distance < cut:
            c = coef / cut
        elif distance > inf:
            c = 0.
        else:
            c = coef / distance

        dx, dy, dz = dx * c, dy * c, dz * c
        fxi, fyi, fzi = forces[i]
        fxj, fyj, fzj = forces[j]
        forces[i] = fxi + dx, fyi + dy, fzi + dz
        forces[j] = fxj - dx, fyj - dy, fzj - dz
    return forces


@njit(f8[:, :](f8[:, :], u2[:, :], f8[:, :]), cache=True)
def spring_force(xyz, springs, springs_distances):
    forces = zeros_like(xyz)
    for i in range(len(springs)):
        n, m = springs[i]
        r, stiff = springs_distances[i]
        nx, ny, nz = xyz[n]
        mx, my, mz = xyz[m]

        dx, dy, dz = nx - mx, ny - my, nz - mz
        distance = sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        f = stiff * (distance - r)

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


@njit(f8[:, :](f8[:, :], f8[:, :], f8), cache=True)
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


@njit(f8[:, :](f8[:, :], f8), cache=True)
def cutoff(forces, cut):
    for n in range(len(forces)):
        dx, dy, dz = forces[n]
        distance = sqrt(dx ** 2 + dy ** 2 + dz ** 2)
        if distance > cut:
            f = cut / distance
            forces[n] = f * dx, f * dy, f * dz
    return forces


@njit((f8[:, :], b1[:, :], u2), cache=True)
def calculate_center(xyz, sssr_matrix, start):
    for n in range(len(sssr_matrix)):
        k = 0
        line = sssr_matrix[n]
        center_x, center_y, center_z = .0, .0, .0
        for i, b in enumerate(line):
            if b:
                k += 1
                x, y, z = xyz[i]
                center_x += x
                center_y += y
                center_z += z
        xyz[n + start] = center_x / k, center_y / k, center_z / k


@njit(f8[:, :](f8[:, :], u2[:, :], u2[:, :], f8[:, :], b1[:, :], u2), cache=True)
def steps(xyz, springs, straights, distances_stiffness, sssr_matrix, start_centers):
    # step 1
    for _ in range(2000):
        r_forces = repulsive_force(xyz, .05, .2, 10)
        s_forces = spring_force(xyz, springs, distances_stiffness)
        forces = r_forces + s_forces
        forces = cutoff(forces, .3)
        xyz = forces + xyz
        calculate_center(xyz, sssr_matrix, start_centers)

    # step 2
    for _ in range(1000):
        r_forces = repulsive_force(xyz, .02, .3, 3)
        s_forces = spring_force(xyz, springs, distances_stiffness)
        forces = r_forces + s_forces
        forces = flattening(forces, xyz, .1)
        forces = cutoff(forces, .3)
        xyz = forces + xyz
        calculate_center(xyz, sssr_matrix, start_centers)

    # step 3
    for _ in range(1000):
        r_forces = repulsive_force(xyz, .005, .4, 1.17) + straight_repulsive(xyz, straights, .05, .4, 2.)
        s_forces = spring_force(xyz, springs, distances_stiffness)
        forces = r_forces + s_forces
        forces = flattening(forces, xyz, .1)
        forces = cutoff(forces, .3)
        xyz = forces + xyz
        calculate_center(xyz, sssr_matrix, start_centers)

    return xyz


@njit(f8[:](f8[:, :], u2[:, :], u2), cache=True)
def get_angles(xyz, springs, bonds_count):
    angles = zeros(bonds_count)
    for i in range(bonds_count):
        n, m = springs[i]
        nx, ny = xyz[n, :2]
        mx, my = xyz[m, :2]
        nm = atan2(mx - nx, my - ny)
        if nm < 0:
            nm += pi
        angles[i] = nm
    return angles


@njit(f8[:, :](f8[:, :], u2, f8, f8), {'p': u2}, cache=True)
def rotate(xyz, atoms_count, shift_x, angle):
    cos_rad = cos(angle)
    sin_rad = sin(angle)
    xy = zeros((atoms_count, 2))

    dx, dy = xyz[0, :2]
    for p in range(1, atoms_count):
        px, py = xyz[p, :2]
        px, py = px - dx, py - dy
        xy[p] = cos_rad * px - sin_rad * py, sin_rad * px + cos_rad * py

    shift_y = xy[:, 1].mean()
    shift_r = shift_x - xy[:, 0].min()
    for p in range(atoms_count):
        px, py = xy[p]
        xy[p] = px + shift_r, py - shift_y
    return xy


class Calculate2D:
    __slots__ = ()

    def __prepare(self, component, randomize, c_stiff, r_stiff):
        atoms = self._atoms
        bonds = self._bonds
        plane = self._plane
        rings = [x for x in self.sssr if all(y in component for y in x)]

        mapping = {}
        springs = []
        straights = []
        xyz_matrix = []
        distances_stiffness = []
        dd = defaultdict(dict)

        # for cycles of 4 and 6 atoms
        c6 = []
        ac4 = set()
        ac6 = set()
        cycles_atoms = set()
        for ccl in rings:
            k = len(ccl)
            cycles_atoms.update(set(ccl))
            if k == 4:
                for a in ccl:
                    ac4.add(a)
            elif k == 6:
                c6.append(ccl)
                for a in ccl:
                    ac6.add(a)

        # create matrix of coordinates
        cube = len(component)
        if randomize:
            for i, n in enumerate(component):
                mapping[n] = i
                xyz_matrix.append([uniform(-cube, cube), uniform(-cube, cube), uniform(-cube, cube)])
        else:
            for i, n in enumerate(component):
                mapping[n] = i
                xyz_matrix.append([*plane[n], .1])

        # create matrix of connecting, springs and springs distances
        for n, m_bond in sorted(bonds.items(), reverse=True, key=lambda x: len(x[1])):
            if n not in component:
                continue
            for m, b in m_bond.items():
                if n not in dd[m]:
                    b = b.order
                    if b == 8 or len(m_bond) > 4:
                        dd[n][m] = dd[m][n] = True
                    else:
                        dd[n][m] = dd[m][n] = False

            # for angle = 180'
            if len(m_bond) == 2:
                (a1, bond1), (a2, bond2) = m_bond.items()
                order1, order2 = bond1.order, bond2.order
                if order1 == order2 == 2 or order1 == 3 or order2 == 3:
                    straights.append((mapping[a1], mapping[a2]))

        for n, m_bond in dd.items():
            i = mapping[n]
            for m, dist in m_bond.items():
                j = mapping[m]
                if (j, i) not in springs:
                    springs.append((i, j))
                    if n in cycles_atoms and m in cycles_atoms:
                        stiff = c_stiff
                    else:
                        stiff = r_stiff
                    distances_stiffness.append([1.32, stiff] if dist else [.825, stiff])

        bonds_count = len(distances_stiffness)
        # add virtual atoms and complement matrices of bonds and coordinates
        end = cube
        for n, m_bond in bonds.items():
            if n not in component:
                continue
            if len(m_bond) == 2:
                atom = atoms[n]
                if atom.atomic_number in {5, 6, 7, 8, 14, 15, 16, 17, 33, 34, 35, 52, 53, 85}:
                    (a1, bond1), (a2, bond2) = m_bond.items()
                    if self._is_angle(bond1, bond2):
                        mapping[-n] = end
                        dd[n][-n] = dd[-n][n] = False
                        if randomize:
                            xyz_matrix.append([uniform(-cube, cube), uniform(-cube, cube), uniform(-cube, cube)])
                        else:
                            xyz_matrix.append([*plane[n], -.1])
                        springs.append((mapping[n], end))
                        distances_stiffness.append([.825, r_stiff])
                        end += 1

        # add cycle centers
        start_centers = end
        sssr_matrix = []
        for ring in rings:
            k = len(ring)
            if k >= 11:  # skip big rings
                continue
            mapping[-end] = end
            sme = [False] * start_centers
            sssr_matrix.append(sme)
            for a, m in zip(ring, [mapping[x] for x in ring]):
                ns = bonds[a]
                sme[m] = True
                if len(ns) == 3:
                    for n in ns:
                        if n not in ring:
                            springs.append((mapping[n], end))
                            distances_stiffness.append([.825 + c_long[k], r_stiff])
            xyz_matrix.append([.0, .0, .0])
            end += 1
        xyz_matrix = array(xyz_matrix)
        if sssr_matrix:
            sssr_matrix = array(sssr_matrix, dtype=bool)
        else:
            sssr_matrix = zeros((0, start_centers), dtype=bool)
        calculate_center(xyz_matrix, sssr_matrix, start_centers)

        # add springs between cycles
        ini = start_centers
        for i, ring in enumerate(rings):
            ln = len(ring)
            if ln >= 11:
                continue
            clc = ini
            for ccl in rings[i+1:]:
                k = len(ccl)
                if k >= 11:
                    continue
                clc += 1
                if not set(ring).isdisjoint(set(ccl)):
                    springs.append((ini, clc))
                    distances_stiffness.append([c_short[ln] + c_short[k], c_stiff])
            ini += 1

        # add springs for cycles of six atoms
        for cycle in c6:
            for i, j in zip(cycle[:3], cycle[3:]):
                springs.append((mapping[i], mapping[j]))
                distances_stiffness.append([1.65, c_stiff])
        # add virtual bonds and complement matrices of bonds and coordinates
        for n, m_bond in dd.items():
            if len(m_bond) == 3:
                if n in cycles_atoms:
                    stiff = c_stiff
                else:
                    stiff = r_stiff
                for m1, m2 in combinations(m_bond, 2):
                    if m2 not in dd[m1]:
                        i, j = mapping[m1], mapping[m2]
                        if (j, i) not in springs or (i, j) not in springs:
                            springs.append((i, j))
                            long1, long2 = m_bond[m1], m_bond[m2]
                            if m1 in ac4 and m2 in ac4:
                                distances_stiffness.append([1.17, stiff])
                            elif not long1 and not long2:
                                distances_stiffness.append([1.43, stiff])
                            elif long1 or long2:
                                distances_stiffness.append([1.87, stiff])
                            else:
                                distances_stiffness.append([2.29, stiff])

        if not straights:
            straights = empty(shape=(0, 2), dtype=uint16)
        else:
            straights = array(straights, dtype=uint16)
        return xyz_matrix, array(springs, dtype=uint16), straights, array(distances_stiffness), cube, bonds_count, \
            sssr_matrix, start_centers

    def _is_angle(self, bond1, bond2):
        pass

    @staticmethod
    def __finish_xyz(xyz, springs, atoms_count, bonds_count, shift_x):
        angles = get_angles(xyz, springs, bonds_count)

        clusters = {}
        bonds = list(range(bonds_count))
        d = pi / 60
        while bonds:
            n = bonds.pop(0)
            nm = angles[n]
            values = [nm]
            for m in range(1, len(bonds)):
                m = angles[-m]
                if nm - d <= m <= nm + d:
                    values.append(m)
            clusters[n] = values

        angles = max((x for x in clusters.items()), key=lambda x: len(x[1]))[1]
        angle = sum(angles) / len(angles)
        xy = rotate(xyz, atoms_count, shift_x, angle)
        return xy, xy[:, 0].max() + .825

    def clean2d(self, *, randomize=False, cycle_stiff=.1, bond_stiff=.05):
        """
        Calculate 2d layout of graph.

        Idea from article:
        Frączek, T. (2016). Simulation-Based Algorithm for Two-Dimensional Chemical Structure Diagram Generation of
        Complex Molecules and Ligand–Protein Interactions. Journal of Chemical Information and Modeling, 56(12),
        2320–2335. doi:10.1021/acs.jcim.6b00391 (https://doi.org/10.1021/acs.jcim.6b00391)

        :param randomize: if True generating random coordinates for molecule
        :param cycle_stiff: stiffness for springs in cycles
        :param bond_stiff: stiffness for other springs
        """
        plane = self._plane

        shift_x = .0
        for component in self.connected_components:
            if len(component) == 1:
                plane[component[0]] = (shift_x, .0)
                shift_x += .825
            elif len(component) == 2:
                plane[component[0]] = (shift_x, .0)
                plane[component[1]] = (shift_x, .825)
                shift_x += .825
            else:
                if not randomize and all(-.0001 < x[0] < .0001 and -.0001 < x[1] < .0001 for x in plane.values()):
                    randomize = True
                xyz, springs, straights, distances_stiffness, atoms_count, bonds_count, sssr_matrix, start_centers = \
                    self.__prepare(component, randomize, cycle_stiff, bond_stiff)
                xyz = steps(xyz, springs, straights, distances_stiffness, sssr_matrix, start_centers)
                xy, shift_x = self.__finish_xyz(xyz, springs, atoms_count, bonds_count, shift_x)
                for i, n in enumerate(component):
                    plane[n] = tuple(xy[i])
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


if find_spec('numpy'):  # load numpy not yet loaded with numba
    if not find_spec('numba'):
        from numpy import array, zeros, uint16, zeros_like, empty
else:  # disable clean2d support
    class Calculate2DMolecule:
        __slots__ = ()

        def clean2d(self, **kwargs):
            raise NotImplemented('numpy required for clean2d')


    class Calculate2DCGR:
        __slots__ = ()

        def clean2d(self, **kwargs):
            raise NotImplemented('numpy required for clean2d')

c_long = {3: 1.65 * cos(2 * pi / 3), 4: 1.65 * cos(pi / 2), 5: 1.65 * cos(2 * pi / 5), 6: 1.65 * cos(pi / 3),
          7: 1.65 * cos(2 * pi / 7), 8: 1.65 * cos(pi / 4), 9: 1.65 * cos(2 * pi / 9), 10: 1.65 * cos(pi / 5)}

c_short = {3: c_long[3] * cos(pi / 3), 4: c_long[4] * cos(pi / 4), 5: c_long[5] * cos(pi / 5),
           6: c_long[6] * cos(pi / 6), 7: c_long[7] * cos(pi / 7), 8: c_long[8] * cos(pi / 8),
           9: c_long[9] * cos(pi / 9), 10: c_long[10] * cos(pi / 10)}

__all__ = ['Calculate2DMolecule', 'Calculate2DCGR']
