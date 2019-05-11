# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
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
from itertools import repeat
from typing import List
from ..cache import cached_property
from ..exceptions import InvalidAromaticRing
from ..periodictable import C


class Aromatize:
    @cached_property
    def aromatic_rings(self) -> List[List[int]]:
        adj = self._adj
        return [ring for ring in self.sssr if len(ring) in (5, 6, 7) and adj[ring[0]][ring[-1]].order == 4
                and all(adj[n][m].order == 4 for n, m in zip(ring, ring[1:]))]

    def dearomatize(self):
        adj = defaultdict(set)  # aromatic skeleton
        for n, m_bond in self._adj.items():
            for m, bond in m_bond.items():
                if bond.order == 4:
                    adj[n].add(m)

    def dummy_aromatize(self):
        """
        convert structure to aromatic form (dummy algorithm. don't detect quinones)

        :return: number of processed rings
        """
        adj = self._adj
        atom = self._node
        total = 0
        unsaturated = {n for n, m_bond in adj.items() if any(bond.order in (2, 4) for bond in m_bond.values())}

        for ring in self.sssr:
            lr = len(ring)
            if lr in (5, 6, 7) and unsaturated.issuperset(ring):
                for n, m in zip(ring, ring[1:]):
                    b = adj[n][m]
                    if b.order != 4:
                        b.order = 4
                b = adj[ring[0]][ring[-1]]
                if b.order != 4:
                    b.order = 4
                total += 1
            elif lr == 5:
                sr = set(ring)
                if len(unsaturated & sr) == 4 and atom[(sr - unsaturated).pop()]._atom in _pyrole_atoms:
                    for n, m in zip(ring, ring[1:]):
                        b = adj[n][m]
                        if b.order != 4:
                            b.order = 4
                    b = adj[ring[0]][ring[-1]]
                    if b.order != 4:
                        b.order = 4
                    total += 1
        if total:
            self.flush_cache()
        return total

    def aromatize(self) -> int:
        """
        convert structure to aromatic form

        :return: number of processed rings
        """
        self.dummy_aromatize()
        adj = self._adj
        patch = set()
        total = 0
        double_bonded = {n for n, m_bond in adj.items() if any(bond.order == 2 for bond in m_bond.values())}

        quinones = []
        condensed_rings = defaultdict(lambda: defaultdict(list))
        for ring in self.aromatic_rings:
            ring = tuple(ring)
            if not double_bonded.isdisjoint(ring):  # search quinones
                quinones.append(ring)

            for n, m in zip(ring, ring[1:]):  # condensed rings graph
                condensed_rings[n][m].append(ring)
                condensed_rings[m][n].append(ring)
            n, *_, m = ring
            condensed_rings[n][m].append(ring)
            condensed_rings[m][n].append(ring)

        while quinones:
            total += 1
            ring = quinones.pop()
            for n, m in zip(ring, ring[1:]):  # remove from condensed rings graph
                condensed_rings[n][m].remove(ring)
                condensed_rings[m][n].remove(ring)
            n, *_, m = ring
            condensed_rings[n][m].remove(ring)
            condensed_rings[m][n].remove(ring)

            start = next(n for n, m in enumerate(ring) if m in double_bonded)
            if start:  # reorder double bonded to starting position
                ordered_ring = ring[start:] + ring[:start]
            else:
                ordered_ring = ring

            bond = 1
            n = ordered_ring[0]
            for m in ordered_ring[1:]:
                if bond == 1:
                    if m not in double_bonded:
                        bond = 2
                    if not condensed_rings[n][m]:
                        patch.add((n, m, 1))
                    elif n in double_bonded:  # found new quinone ring (Y)
                        q = condensed_rings[n][m][0]
                        if q not in quinones:
                            quinones.insert(0, q)  # low priority
                else:
                    if m in double_bonded:
                        raise InvalidAromaticRing(ring)
                    bond = 1
                    if not condensed_rings[n][m]:
                        patch.add((n, m, 2))
                        double_bonded.add(n)
                        double_bonded.add(m)
                        if condensed_rings[n][p]:
                            q = condensed_rings[n][p][0]
                            if q in quinones:  # up priority
                                quinones.remove(q)
                                quinones.append(q)
                            else:
                                quinones.insert(0, q)
                p, n = n, m
            else:
                m = ordered_ring[0]
                if bond != 1 and not condensed_rings[n][p]:
                    raise InvalidAromaticRing(ring)
                patch.add((n, m, 1))

        if patch:
            for n, m, b in patch:
                adj[n][m].order = b
            self.flush_cache()
        return total

    def aromatize_(self) -> int:
        """
        convert structure to aromatic form

        :return: number of processed rings
        """
        rings5, rings6, rings10 = self.__prepare_rings()
        if not (rings6 or rings5 or rings10):
            return 0

        init = len(rings6) + len(rings5) + len(rings10)
        total = 0
        while True:
            # quinone rules don't match if condensed aromatics rings has invalid bond orders.
            # this require repeating of quinonize step after aromatize step
            old = 0
            new = init
            rings5_c = rings5.copy()
            rings6_c = rings6.copy()
            rings10_c = rings10.copy()
            while new != old:
                found = []
                for n, r in enumerate(rings6_c):
                    if self.__quinonize_benzene(r):
                        found.insert(0, n)
                for n in found:
                    del rings6_c[n]

                found = []
                for n, r in enumerate(rings5_c):
                    if self.__quinonize_pyrole(r):
                        found.insert(0, n)
                for n in found:
                    del rings5_c[n]

                found = []
                for n, r in enumerate(rings10_c):
                    if self.__quinonize_azulene(r):
                        found.insert(0, n)
                for n in found:
                    del rings10_c[n]

                old, new = new, len(rings6_c) + len(rings5_c) + len(rings10_c)
            c = init - new
            if c:
                total += c
            elif total:
                break

            old = 0
            new = init
            rings5_c = rings5.copy()
            rings6_c = rings6.copy()
            rings10_c = rings10.copy()
            while new != old:
                found = []
                for n, r in enumerate(rings6_c):
                    if self.__aromatize_benzene(r):
                        found.insert(0, n)
                for n in found:
                    del rings6_c[n]

                found = []
                for n, r in enumerate(rings5_c):
                    if self.__aromatize_pyrole(r):
                        found.insert(0, n)
                for n in found:
                    del rings5_c[n]

                found = []
                for n, r in enumerate(rings10_c):
                    if self.__aromatize_azulene(r):
                        found.insert(0, n)
                for n in found:
                    del rings10_c[n]

                old, new = new, len(rings6_c) + len(rings5_c) + len(rings10_c)
            c = init - new
            if not c:
                break
            total += c

        if total:
            self.flush_cache()
        return total

    def __prepare_rings(self):
        rings5 = {}
        rings7 = {}
        rings6 = []
        for ring in self.sssr:
            lr = len(ring)
            if lr == 6:
                rings6.append(ring)
            elif lr == 5:
                rings5[frozenset(ring)] = ring
            elif lr == 7:
                rings7[frozenset(ring)] = ring

        # for azulene only external contour need for detection
        #
        #    7    8
        #  /   \ / \
        # 1     6   9
        # |     |   |
        # 2     5---10
        #  \   /
        #   3-4
        #
        possible_rings10 = defaultdict(list)
        for r7_k, r7 in rings7.items():
            for r5_k, r5 in rings5.items():
                c = r5_k & r7_k
                if len(c) == 2:
                    possible_rings10[r7_k].append((c, r7, r5))  # found r5r7

        rings10 = []
        ambiguous = []
        for r7_k, cr75_list in possible_rings10.items():
            if len(cr75_list) == 1:
                rings10.append(cr75_list[0])

            n, m = c
            # rearrange ring to n(6)----m(5)
            i = r7.index(n)
            ij = i - r7.index(m)
            if ij == 1:  # normal direction. e.g. 1-2-3-4-5-6-7
                r7 = r7[i:] + r7[:i]
            elif ij == -1:  # reverse direction. e.g.  1-7-6-5-4-3-2
                r7 = r7[i::-1] + r7[:i:-1]
            elif i:
                r7 = r7[::-1]

            # rearrange ring to m(5)---n(6)
            r5 = rings5.pop(r5)
            i = r5.index(m)
            j = r5.index(n)
            ij = i - j
            if ij == 1:  # normal direction. e.g. 8-6-5-10-9
                r5 = r5[i + 1:] + r5[:j]  # 10-9-8
            elif ij == -1:  # reverse direction. e.g 10-5-6-8-9
                if i:
                    r5 = r5[i - 1::-1] + r5[:j:-1]
                else:
                    r5 = r5[:1:-1]
            elif i:
                r5 = r5[-2:0:-1]
            else:
                r5 = r5[1:-1]
            rings10.append(r7 + r5)

        return list(rings5.values()), rings6, rings10

    def __quinonize_benzene(self, ring):
        adj = self._adj
        r1, r2, r3, r4, r5, r6 = ring
        r12 = adj[r1][r2]
        r23 = adj[r2][r3]
        r34 = adj[r3][r4]
        r45 = adj[r4][r5]
        r56 = adj[r5][r6]
        r61 = adj[r6][r1]
        key = (r12.order, r23.order, r34.order, r45.order, r56.order, r61.order)
        if 4 not in key:
            return False

        doubles = tuple(y for y, x in enumerate(ring) if len(adj[x]) == 3 and
                        next(b.order for a, b in adj[x].items() if a not in ring) == 2)
        if not doubles:
            return False

        if len(doubles) == 6:
            r12.order = r23.order = r34.order = r45.order = r56.order = r61.order = 1
            return True
        if key in _quinone_pattern.get(doubles, {}):
            r12.order, r23.order, r34.order, r45.order, r56.order, r61.order = _quinone_fix.get(doubles)
            return True
        return False

    def __quinonize_pyrole(self, ring):
        adj = self._adj
        atoms = self._node
        r1, r2, r3, r4, r5 = ring
        r12 = adj[r1][r2]
        r23 = adj[r2][r3]
        r34 = adj[r3][r4]
        r45 = adj[r4][r5]
        r51 = adj[r5][r1]
        key = (r12.order, r23.order, r34.order, r45.order, r51.order)
        if 4 not in key:
            return False

        positions = _pyrole_pattern.get(key)
        if positions is None:
            return False

        for m, pos in enumerate(positions):
            if atoms[ring[pos]]._atom in _pyrole_atoms:
                r12.order, r23.order, r34.order, r45.order, r51.order = _pyrole_fix[key][m]
                return True
        return False

    def __quinonize_azulene(self, ring):
        adj = self._adj
        r1, r2, r3, r4, r5, r6, r7, r8, r9, r0 = ring
        r12 = adj[r1][r2]
        r23 = adj[r2][r3]
        r34 = adj[r3][r4]
        r45 = adj[r4][r5]
        r56 = adj[r5][r6]
        r67 = adj[r6][r7]
        r78 = adj[r7][r8]
        r89 = adj[r8][r9]
        r90 = adj[r9][r0]
        r01 = adj[r0][r1]
        r17 = adj[r1][r7]
        return False

    def __aromatize_benzene(self, ring):
        adj = self._adj
        r1, r2, r3, r4, r5, r6 = ring
        r12 = adj[r1][r2]
        r23 = adj[r2][r3]
        r34 = adj[r3][r4]
        r45 = adj[r4][r5]
        r56 = adj[r5][r6]
        r61 = adj[r6][r1]
        if (r12.order, r23.order, r34.order, r45.order, r56.order, r61.order) in _benzene:
            r12.order = r23.order = r34.order = r45.order = r56.order = r61.order = 4
            return True
        return False

    def __aromatize_pyrole(self, ring):
        adj = self._adj
        atoms = self._node
        r1, r2, r3, r4, r5 = ring
        r12 = adj[r1][r2]
        r23 = adj[r2][r3]
        r34 = adj[r3][r4]
        r45 = adj[r4][r5]
        r51 = adj[r5][r1]
        position = _pyrole.get((r12.order, r23.order, r34.order, r45.order, r51.order))

        if position is not None and atoms[ring[position]]._atom in _pyrole_atoms:
            r12.order = r23.order = r34.order = r45.order = r51.order = 4
            return True
        return False

    def __aromatize_azulene(self, ring):
        adj = self._adj
        r1, r2, r3, r4, r5, r6, r7, r8, r9, r0 = ring
        r12 = adj[r1][r2]
        r23 = adj[r2][r3]
        r34 = adj[r3][r4]
        r45 = adj[r4][r5]
        r56 = adj[r5][r6]
        r67 = adj[r6][r7]
        r78 = adj[r7][r8]
        r89 = adj[r8][r9]
        r90 = adj[r9][r0]
        r01 = adj[r0][r1]
        r17 = adj[r1][r7]
        if (r12.order, r23.order, r34.order, r45.order, r56.order,
            r67.order, r78.order, r89.order, r90.order, r01.order) in _azulene:
            r12.order = r23.order = r34.order = r45.order = r56.order = 4
            r67.order = r78.order = r89.order = r90.order = r01.order = r17.order = 4
            return True
        return False


def _clock(a):
    yield a
    for _ in range(1, len(a)):
        a = a[1:] + a[:1]
        yield a


_pyrole_atoms = ('N', 'O', 'S', 'Se', 'P', C(-1))

_azulene = {(1, 2, 1, 2, 1, 2, 1, 2, 1, 2), (2, 1, 2, 1, 2, 1, 2, 1, 2, 1)}

_benzene = set()
_benzene.update(_clock((1, 2, 1, 2, 1, 2)))
_benzene.update(_clock((1, 2, 1, 2, 1, 4)))
_benzene.update(_clock((1, 2, 1, 2, 4, 2)))
_benzene.update(_clock((1, 2, 1, 2, 4, 4)))
_benzene.update(_clock((1, 2, 1, 4, 4, 2)))
_benzene.update(_clock((1, 2, 1, 4, 1, 4)))
_benzene.update(_clock((1, 2, 4, 2, 4, 2)))
_benzene.update(_clock((1, 2, 4, 1, 4, 2)))
_benzene.update(_clock((1, 2, 4, 2, 1, 4)))
_benzene.update(_clock((1, 2, 1, 4, 4, 4)))
_benzene.update(_clock((1, 2, 4, 4, 4, 2)))
_benzene.update(_clock((1, 4, 1, 4, 1, 4)))
_benzene.update(_clock((4, 2, 4, 2, 4, 2)))
_benzene.update(_clock((1, 2, 4, 2, 4, 4)))
_benzene.update(_clock((1, 4, 4, 2, 4, 2)))
_benzene.update(_clock((1, 2, 4, 4, 1, 4)))
_benzene.update(_clock((1, 4, 4, 2, 1, 4)))
_benzene.update(_clock((1, 2, 4, 4, 4, 4)))
_benzene.update(_clock((1, 4, 4, 4, 4, 2)))
_benzene.update(_clock((1, 4, 4, 4, 4, 4)))
_benzene.update(_clock((4, 2, 4, 4, 4, 4)))

_ind = (0, 4, 3, 2, 1)
_pyrole = {}
_pyrole.update(zip(_clock((1, 2, 1, 2, 1)), _ind))
_pyrole.update(zip(_clock((1, 4, 1, 2, 1)), _ind))
_pyrole.update(zip(_clock((1, 2, 1, 4, 1)), _ind))
_pyrole.update(zip(_clock((1, 2, 4, 2, 1)), _ind))
_pyrole.update(zip(_clock((1, 4, 1, 4, 1)), _ind))
_pyrole.update(zip(_clock((1, 4, 4, 2, 1)), _ind))
_pyrole.update(zip(_clock((1, 2, 4, 4, 1)), _ind))
_pyrole.update(zip(_clock((1, 4, 4, 4, 1)), _ind))
# fixes after quinonize
_pyrole.update(zip(_clock((4, 2, 4, 4, 4)), _ind))
_pyrole.update(zip(_clock((4, 4, 4, 2, 4)), _ind))
_pyrole.update(zip(_clock((4, 2, 4, 2, 4)), _ind))
_pyrole.update(zip(_clock((4, 4, 2, 4, 4)), _ind))

_ind = ((0, 1), (0, 5), (4, 5), (3, 4), (2, 3), (1, 2))  # o-quinones
_quinone_pattern = {}
_quinone_fix = {}
for i, *p in zip(_ind, repeat((4, 4, 4, 4, 4, 4)),
                 _clock((4, 4, 4, 2, 4, 4)),
                 _clock((4, 4, 4, 4, 2, 4)),
                 _clock((4, 4, 2, 4, 4, 4)),
                 _clock((4, 4, 2, 4, 2, 4)),
                 _clock((4, 4, 4, 1, 2, 4)),
                 _clock((4, 4, 2, 1, 4, 4)),
                 _clock((4, 4, 2, 1, 2, 4)),

                 _clock((1, 4, 4, 4, 4, 4)),
                 _clock((1, 4, 4, 2, 4, 4)),
                 _clock((1, 4, 4, 4, 2, 4)),
                 _clock((1, 4, 2, 4, 4, 4)),
                 _clock((1, 4, 2, 4, 2, 4)),
                 _clock((1, 4, 4, 1, 2, 4)),
                 _clock((1, 4, 2, 1, 4, 4)),
                 _clock((1, 4, 2, 1, 2, 4)),

                 _clock((1, 1, 4, 4, 4, 4)),
                 _clock((1, 4, 4, 4, 4, 1))
                 ):
    _quinone_pattern[i] = set(p)

_quinone_fix.update(zip(_ind, _clock((1, 1, 2, 1, 2, 1))))

_ind = ((0, 1, 2, 3), (0, 1, 2, 5), (0, 1, 4, 5), (0, 3, 4, 5), (2, 3, 4, 5), (1, 2, 3, 4))  # 1,2,3,4-quinones
for i, *p in zip(_ind, repeat((4, 4, 4, 4, 4, 4)),
                 _clock((1, 4, 1, 4, 4, 4)),
                 _clock((1, 4, 1, 4, 2, 4))
                 ):
    _quinone_pattern[i] = set(p)

_quinone_fix.update(zip(_ind, _clock((1, 1, 1, 1, 2, 1))))

_ind = ((0, 3), (2, 5), (1, 4))  # p-quinones
for i, *p in zip(_ind, repeat((4, 4, 4, 4, 4, 4)),
                 _clock((4, 2, 4, 4, 4, 4)),
                 _clock((4, 2, 4, 4, 2, 4))
                 ):
    _quinone_pattern[i] = set(p)

_quinone_fix.update(zip(_ind, _clock((1, 2, 1, 1, 2, 1))))

# pyroles condensed with quinones fixes
_ind = (0, 4, 3, 2, 1)
_pyrole_pattern = defaultdict(list)
_pyrole_fix = defaultdict(list)
for i, *p in zip(_ind,
                 zip(_clock((4, 1, 4, 4, 4)), _clock((1, 1, 1, 2, 1))),
                 zip(_clock((4, 4, 4, 1, 4)), _clock((1, 2, 1, 1, 1))),
                 zip(_clock((4, 1, 4, 2, 4)), _clock((1, 1, 1, 2, 1))),
                 zip(_clock((4, 2, 4, 1, 4)), _clock((1, 2, 1, 1, 1))),
                 zip(_clock((4, 1, 4, 1, 4)), repeat((1, 1, 1, 1, 1)))
                 ):
    for x, y in p:
        _pyrole_pattern[x].append(i)
        _pyrole_fix[x].append(y)

_pyrole_pattern = dict(_pyrole_pattern)
_pyrole_fix = dict(_pyrole_fix)


del x, y, i, p, _ind


__all__ = ['Aromatize']
