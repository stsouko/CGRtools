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
from collections import defaultdict
from itertools import combinations, product
#from ..exceptions import InvalidAtomNumber, InvalidWedgeMark, InvalidStereoCenter


def _pyramid_sign(n, u, v, w):
    #
    #  |   n /
    #  |   |\
    #  |   | \
    #  |  /|  \
    #  | / u---v
    #  |/___\_/___
    #        w
    #
    nx, ny, nz = n
    ux, uy, uz = u
    vx, vy, vz = v
    wx, wy, wz = w

    q1x = ux - nx
    q1y = uy - ny
    q1z = uz - nz
    q2x = vx - nx
    q2y = vy - ny
    q2z = vz - nz
    q3x = wx - nx
    q3y = wy - ny
    q3z = wz - nz

    vol = q1x * (q2y * q3z - q2z * q3y) + q1y * (q2z * q3x - q2x * q3z) + q1z * (q2x * q3y - q2y * q3x)
    if vol > 0:
        return 1
    elif vol < 0:
        return -1
    return 0


def _dihedral_sign(n, u, v, w):
    # n    w
    # |   /
    # u--v
    nx, ny, nz = n
    ux, uy, uz = u
    vx, vy, vz = v
    wx, wy, wz = w

    q1x = ux - nx
    q1y = uy - ny
    q1z = uz - nz
    q2x = vx - ux
    q2y = vy - uy
    q2z = vz - uz
    q3x = wx - vx
    q3y = wy - vy
    q3z = wz - vz

    # cross vectors
    q1q2x = q1y * q2z - q1z * q2y
    q1q2y = q1z * q2x - q1x * q2z
    q1q2z = q1x * q2y - q1y * q2x
    q2q3x = q2y * q3z - q2z * q3y
    q2q3y = q2z * q3x - q2x * q3z
    q2q3z = q2x * q3y - q2y * q3x

    # angle calculation
    # len_q1q2 = sqrt(q1q2x ** 2 + q1q2y ** 2 + q1q2z ** 2)
    # n1x = q1q2x / len_q1q2
    # n1y = q1q2y / len_q1q2
    # n1z = q1q2z / len_q1q2
    # len_q2q3 = sqrt(q2q3x ** 2 + q2q3y ** 2 + q2q3z ** 2)
    # u1x = q2q3x / len_q2q3
    # u1y = q2q3y / len_q2q3
    # u1z = q2q3z / len_q2q3
    # len_q2 = sqrt(q2x ** 2 + q2y ** 2 + q2z ** 2)
    # u3x = q2x / len_q2
    # u3y = q2y / len_q2
    # u3z = q2z / len_q2
    # u2x = u3y * u1z - u3z * u1y
    # u2y = u3z * u1x - u3x * u1z
    # u2z = u3x * u1y - u3y * u1x
    # cos_theta = n1x * u1x + n1y * u1y + n1z * u1z
    # sin_theta = n1x * u2x + n1y * u2y + n1z * u2z
    # return -atan2(sin_theta, cos_theta)

    dot = q1q2x * q2q3x + q1q2y * q2q3y + q1q2z * q2q3z
    if dot > 0:
        return 1
    elif dot < 0:
        return -1
    return 0


class Stereo:
    def add_wedge(self, n, m, mark):
        try:
            if self._bonds[n][m].order not in (1, 4):
                raise InvalidWedgeMark((n, m))
        except KeyError:
            raise InvalidAtomNumber((n, m))

        if n in self.__tetrahedrons:
            env = self.__tetrahedrons[n]

        elif n in self.__allenes:
            ...
        else:
            raise InvalidStereoCenter(n)

    @cached_property
    def chiral_atoms(self):
        bonds = self._bonds
        morgan = self.atoms_order
        atoms_stereo = self._atoms_stereo
        tetra = self.__tetrahedrons

        grouped_stereo = defaultdict(list)
        for n, s in atoms_stereo.items():
            grouped_stereo[morgan[n]].append(n)  # collect equal stereo atoms

        morgan_update = {}
        for group in grouped_stereo.values():
            if len(group) % 2 == 0:  # only even number of equal stereo atoms give new stereo center
                if group[0] in tetra:  # tetrahedron stereo
                    s = [n for n in group if self._translate_tetrahedron_stereo(n, sorted(bonds[n], key=morgan.get))]
                    if 0 < len(s) < len(group):  # RS pair required
                        for n in s:
                            morgan_update[n] = -morgan[n]

        morgan = self._morgan(self._sorted_primed({**morgan, **morgan_update}))

        chiral_atoms = {n for n, env in tetra.items()
                        if n not in atoms_stereo and len(set(morgan[x] for x in env)) == len(env)}

        #for (n, m), (n1, m1, n2, m2) in self.__cumulenes.items():
        #    if n2 and morgan[n1] == morgan[n2]:
        #        continue
        #    if m2 and morgan[m1] == morgan[m2]:
        #        continue
        #    chiral_atoms.add(n)
        #    chiral_atoms.add(m)

        #    morgan = self._morgan(self._sorted_primed({**morgan,
        #                                               **{m: -n for n, m in enumerate(chiral_atoms, start=1)}}))
        return chiral_atoms

    def _translate_tetrahedron_stereo(self, n, env):
        order = self.__tetrahedrons[n]
        if len(order) != len(env):
            raise ValueError('invalid atoms list')
        s = self._atoms_stereo[n]

        translate = tuple(order.index(x) + 1 for x in env[:3])
        if _tetrahedron_translate[translate]:
            return not s
        return s

    @cached_property
    def __tetrahedrons(self):
        #    2
        #    |
        # 1--K--3
        #    |
        #    4?
        atoms = self._atoms
        bonds = self._bonds
        tetrahedrons = {}
        for n in self.tetrahedrons:
            env = tuple(x for x in bonds[n] if atoms[x].atomic_number != 1)
            if len(env) in (3, 4):
                tetrahedrons[n] = env
        return tetrahedrons

    @cached_property
    def __cumulenes(self):
        # 5       4
        #  \     /
        #   2---3
        #  /     \
        # 1       6
        bonds = self._bonds
        atoms = self._atoms
        cumulenes = {}
        for path in self.cumulenes:
            n, m = path[0], path[-1]
            n1, m1 = path[1], path[-2]

            nn = [x for x in bonds[n] if x != n1 and atoms[x].atomic_number != 1]
            mn = [x for x in bonds[m] if x != m1 and atoms[x].atomic_number != 1]
            if nn and mn:
                sn = nn[1] if len(nn) == 2 else None
                sm = mn[1] if len(mn) == 2 else None
                cumulenes[(n, m)] = (nn[0], mn[0], sn, sm)
        return cumulenes

    def __atropoisomers(self):
        #     ___
        #    |   |
        # 7--3   5--9
        #     \ /
        #      1[K]
        #      |
        #      2[K]
        #     / \
        # 8--4   6--10
        #    |___|
        #
        adj = self._bonds
        if len(self.sssr) < 2:
            return {}
        aromatic = self.aromatic_rings
        if len(aromatic) < 2:
            return {}

        reduced_aromatic = [[x for x in x if len(adj[x]) == 3] for x in aromatic]  # remove :[CH]: atoms
        connections = {}
        for rings in combinations(range(len(aromatic)), 2):
            ring1, ring2 = rings
            for n, m in product(reduced_aromatic[ring1], reduced_aromatic[ring2]):
                if n in adj[m]:
                    if rings in connections:  # remove condensed rings or twice-bonded rings
                        del connections[rings]
                        break  # skip rings
                    connections[rings] = (n, m)

        atropos = {}
        for (ring1, ring2), (n, m) in connections.items():
            # neighbors of connection atoms in rings
            r1n, r1m = (x for x in adj[n] if x != m)
            r2n, r2m = (x for x in adj[m] if x != n)
            # substituents of neighbors
            nr1n = next((x for x in adj[r1n] if x not in aromatic[ring1]), None)
            nr1m = next((x for x in adj[r1m] if x not in aromatic[ring1]), None)
            nr2n = next((x for x in adj[r2n] if x not in aromatic[ring2]), None)
            nr2m = next((x for x in adj[r2m] if x not in aromatic[ring2]), None)

            # skip rings without substituents
            # todo: rings bounded with chain
            if nr1n is None and (nr1m is None or nr2n is None or nr2m is None):
                continue
            elif nr1m is None and (nr2n is None or nr2m is None):
                continue
            elif nr2n is None and nr2m is None:
                continue
            atropos[(n, m)] = (n, m, r1n, r2n, r1m, r2m, nr1n, nr2n, nr1m, nr2m)
        return atropos


# 2  3
#  \ |
#   \|
#    0---4
#   /
#  /
# 1
_tetrahedron_translate = {(1, 2, 3): False, (2, 3, 1): False, (3, 1, 2): False,
                          (1, 3, 2): True, (2, 1, 3): True, (3, 2, 1): True,
                          (1, 4, 2): False, (4, 2, 1): False, (2, 1, 4): False,
                          (1, 2, 4): True, (2, 4, 1): True, (4, 1, 2): True,
                          (1, 3, 4): False, (3, 4, 1): False, (4, 1, 3): False,
                          (1, 4, 3): True, (4, 3, 1): True, (3, 1, 4): True,
                          (2, 4, 3): False, (4, 3, 2): False, (3, 2, 4): False,
                          (2, 3, 4): True, (3, 4, 2): True, (4, 2, 3): True}
# 5       4
#  \     /
#   2---3
#  /     \
# 1       6
_alkene_translate = {(1, 2, 3, 4): False, (4, 3, 2, 1): False, (1, 2, 3, 6): True, (6, 3, 2, 1): True,
                     (5, 2, 3, 6): False, (6, 3, 2, 5): False, (5, 2, 3, 4): True, (4, 3, 2, 5): True}


__all__ = ['Stereo']
