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
from collections import defaultdict
from itertools import combinations, product
from typing import List
from ..cache import cached_property


def _pyramid_sign(n, u, v, w):
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
    @cached_property
    def tetrahedrons(self):
        # tetrahedral should be single bonded and contain zero or one H
        atoms = self._node
        return {n: tuple(env) for n, env in self._adj.items()
                if atoms[n].element == 'C' and all(x.order == 1 for x in env.values())}

    def __tetrahedrons(self):
        atoms = self._node
        tetrahedrons = {}
        for n, env in self.tetrahedrons.items():
            env = tuple(x for x in env if atoms[x].element != 'H')
            if len(env) >= 3:
                tetrahedrons[n] = env
        return tetrahedrons

    @cached_property
    def cumulenes(self) -> List[List[int]]:
        """
        list of alkenes, allenes and cumulenes atoms lists
        """
        atoms = self._node
        adj = defaultdict(set)  # carbon double bonds adjacency matrix
        for n, m, bond in self.bonds():
            if bond.order == 2 and atoms[n].element == atoms[m].element == 'C':
                adj[n].add(m)
                adj[m].add(n)

        if not adj:
            return []

        terminals = {x for x, y in adj.items() if len(y) == 1}
        cumulenes = []
        while terminals:
            m = terminals.pop()
            path = [m]
            cumulenes.append(path)
            while m not in terminals:
                n, m = m, adj[m].pop()
                adj[m].discard(n)
                path.append(m)
            terminals.discard(m)
        return cumulenes

    def __allenes(self):
        # structure of allene:
        # (central atom, first neighbor atom of enter atom, enter atom, exit atom, first neighbor atom of exit atom,
        #  second enter or None, second exit or None)
        allenes = []
        for path in self.cumulenes:
            lp = len(path)
            if lp % 2:
                try:
                    allenes.append((path[lp // 2], *self.__cumulene_filter(path)))
                except TypeError:
                    pass
        return allenes

    def __alkenes(self):
        # structure of alkene:
        # (first neighbor atom of enter atom, enter atom, exit atom, first neighbor atom of exit atom,
        #  second enter or None, second exit or None)
        alkenes = []
        for path in self.cumulenes:
            if len(path) == 2:
                sp = set(path)
                if any(len(ring) <= 7 and sp.issubset(ring) for ring in self.sssr):
                    continue
                env = self.__cumulene_filter(path)
                if env:
                    alkenes.append(env)
        return alkenes

    def __cumulenes(self):
        # structure of cis_trans:
        # (first central atom, second central atom, first neighbor atom of enter atom, enter atom, exit atom,
        #  first neighbor atom of exit atom, second enter or None, second exit or None)
        cis_trans = []
        for path in self.cumulenes:
            lp = len(path)
            if lp > 3 and not lp % 2:
                try:
                    cis_trans.append((path[lp // 2 - 1], path[lp // 2], *self.__cumulene_filter(path)))
                except TypeError:
                    pass
        return cis_trans

    def __cumulene_filter(self, path):
        adj = self._adj
        atoms = self._node
        n, m = path[0], path[-1]
        n1, m1 = path[1], path[-2]

        nn = [x for x in adj[n] if x != n1 and atoms[x].element != 'H']
        mn = [x for x in adj[m] if x != m1 and atoms[x].element != 'H']
        if nn and mn:
            if len(nn) == 2:
                sn = nn[1]
            else:
                sn = None
            if len(mn) == 2:
                sm = mn[1]
            else:
                sm = None
            return nn[0], n, m, mn[0], sn, sm

    @cached_property
    def chiral_atoms(self):
        """
        dict of chiral atoms valued with ordered neighbors
        """
        atoms = {}
        bonds = {}
        order = self.atoms_order

        chance_tetrahedron = {}
        for a, n in self.__potentially_tetrahedron.items():
            # chiral atom all times contain unique neighbors
            if len(n) == len({order[x] for x in n}):
                atoms[a] = n
            else:
                chance_tetrahedron[a] = n

        cis_trans, allene = self.__potentially_alkene
        # k         n
        #  \       /
        #   l==a==m
        #  /       \
        # i         j
        chance_allene = []
        for k, l, m, n, a, i, j in allene:
            if i and order[i] == order[k] or j and order[j] == order[n]:
                chance_allene.append((k, l, m, n, a, i, j))
            else:
                atoms[a] = (k, l, m, n)
        # k            n
        #  \          /
        #   l==a==b==m
        #  /          \
        # i            j
        chance_cis_trans = []
        for k, l, m, n, a, b, i, j in cis_trans:
            if i and order[i] == order[k] or j and order[j] == order[n]:
                chance_cis_trans.append((k, l, m, n, a, b, i, j))
            else:
                bonds[(a, b)] = (k, l, m, n)

        #    ___
        #   |   |
        #   k   i
        #  / \ / \
        # x   l   y
        #     |
        # v   m   w
        #  \ / \ /
        #   n   j
        #   |___|
        #
        for k, l, m, n, i, j, x, y, v, w in self.__potentially_atropisomer:
            pass
        return

    @cached_property
    def __potentially_atropisomer(self):
        if len(self.sssr) < 2:
            return []
        aromatic = [ring for ring in self.sssr if len(ring) in (5, 6) and self._adj[ring[0]][ring[-1]].order == 4
                    and all(self._adj[n][m].order == 4 for n, m in zip(ring, ring[1:]))]
        if len(aromatic) < 2:
            return []

        reduced_aromatic = [[x for x in x if len(self._adj[x]) == 3] for x in aromatic]  # remove :[CH]: atoms
        connections = {}
        for rings in combinations(range(len(reduced_aromatic)), 2):
            ring1, ring2 = rings
            for n, m in product(reduced_aromatic[ring1], reduced_aromatic[ring2]):
                if n in self._adj[m]:
                    if rings in connections:  # remove condensed rings or twice-bonded rings
                        del connections[rings]
                        break  # skip rings
                    connections[rings] = (n, m)

        # structure of atropisomer:
        # (first neighbor atom of ring1 atom, ring1 atom, ring2 atom, first neighbor atom of ring2 atom,
        #  second neighbor atom of ring1 atom, second neighbor atom of ring2 atom,
        #  nonring neighbor of first neighbor atom of ring1 atom or None,
        #  nonring neighbor of second neighbor atom of ring1 atom or None,
        #  nonring neighbor of first neighbor atom of ring2 atom or None,
        #  nonring neighbor of second neighbor atom of ring2 atom or None)

        atropos = []
        for (ring1, ring2), (n, m) in connections.items():
            r1n, r1m = (x for x in self._adj[n] if x != m)
            r2n, r2m = (x for x in self._adj[m] if x != n)

            nr1n = next((x for x in self._adj[r1n] if x not in aromatic[ring1]), None)
            nr1m = next((x for x in self._adj[r1m] if x not in aromatic[ring1]), None)
            nr2n = next((x for x in self._adj[r2n] if x not in aromatic[ring2]), None)
            nr2m = next((x for x in self._adj[r2m] if x not in aromatic[ring2]), None)

            # skip rings without substituents
            if nr1n is None and (nr1m is None or nr2n is None or nr2m is None):
                continue
            elif nr1m is None and (nr2n is None or nr2m is None):
                continue
            elif nr2n is nr2m:
                continue
            atropos.append((r1n, n, m, r2n, r1m, r2m, nr1n, nr1m, nr2n, nr2m))
        return atropos

    @cached_property
    def __potentially_plane(self):
        if not self.sssr:  # no rings
            return []

        order = self.atoms_order
        if len(set(order.values())) == len(self):  # not symmetric
            return []

        not_aromatic = [ring for ring in self.sssr if any(self._adj[n][m].order != 4 for n, m in zip(ring, ring[1:]))]
        if not not_aromatic:
            return []

        target = []
        tetrahedron = self.__potentially_tetrahedron
        for ring in not_aromatic:
            sring = set(ring)
            if sring.isdisjoint(tetrahedron):
                continue  # hasn't candidates
            for x in sring.intersection(tetrahedron):
                if len(self._adj[x]) != len({order[x] for x in self._adj[x]}):
                    target.append(x)


__all__ = ['Stereo']
