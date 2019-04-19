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
from ..periodictable import C


class Aromatize:
    @cached_property
    def aromatic_rings(self) -> List[List[int]]:
        adj = self._adj
        return [ring for ring in self.sssr if len(x) in (5, 6, 7) and adj[ring[0]][ring[-1]].order == 4
                and all(adj[n][m].order == 4 for n, m in zip(ring, ring[1:]))]

    def aromatize(self) -> int:
        """
        convert structure to aromatic form

        :return: number of processed rings
        """
        rings5 = {}
        rings7 = {}
        rings6 = []
        rings10 = []
        for ring in self.sssr:
            lr = len(ring)
            if lr == 6:
                rings6.append(ring)
            if lr == 5:
                rings5[frozenset(ring)] = ring
            elif lr == 7:
                rings7[frozenset(ring)] = ring

        #
        #    7    8
        #  /   \ / \
        # 1     6   9
        # |     |   |
        # 2     5---10
        #  \   /
        #   3-4
        #
        for r7_k, r7 in rings7.items():
            for r5 in rings5:
                c = r5 & r7_k
                if len(c) == 2:
                    break  # found r5r7
            else:
                continue  # not found
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

        rings5 = list(rings5.values())

        if not (rings6 or rings5 or rings10):
            return 0

        init = len(rings6) + len(rings5) + len(rings10)

        old = 0
        new = init
        rings5_c = rings5.copy()
        rings6_c = rings6.copy()
        while new != old:
            found = []
            for n, r in enumerate(rings6_c):
                if self.__quinonize_benzene(r):
                    found.insert(0, n)
            for n in found:
                del rings6_c[n]

            old, new = new, len(rings6_c) + len(rings5_c)
        total = init - new

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
        total += init - new

        if total:
            self.flush_cache()
        return total

    def __quinonize_benzene(self, ring):
        r1, r2, r3, r4, r5, r6 = r
        key = (self._adj[r1][r2][bond], self._adj[r2][r3][bond], self._adj[r3][r4][bond],
               self._adj[r4][r5][bond], self._adj[r5][r6][bond], self._adj[r6][r1][bond])
        if 4 not in key:
            continue

        doubles = tuple(y for y, x in enumerate(r) if len(self._adj[x]) == 3 and
                        next(attr[bond] for a, attr in self._adj[x].items() if a not in r) == 2)
        if not doubles:
            continue

        if len(doubles) == 6:
            self._adj[r1][r2][bond] = self._adj[r2][r3][bond] = self._adj[r3][r4][bond] = 1
            self._adj[r4][r5][bond] = self._adj[r5][r6][bond] = self._adj[r6][r1][bond] = 1
            found.append(n)
        else:
            if key in _quinone_pattern.get(doubles, {}):
                dear = _quinone_fix.get(doubles)
                self._adj[r1][r2][bond], self._adj[r2][r3][bond], self._adj[r3][r4][bond], \
                self._adj[r4][r5][bond], self._adj[r5][r6][bond], self._adj[r6][r1][bond] = dear
                found.append(n)

    def _quinonize(self, rings):
        bond = 'order'
        rings = rings.copy()
        init = len(rings)
        old = 0
        new = len(rings)
        while new != old:
            old = new
            found = []
            for n, r in enumerate(rings):
                if len(r) == 6:
                    r1, r2, r3, r4, r5, r6 = r
                    key = (self._adj[r1][r2][bond], self._adj[r2][r3][bond], self._adj[r3][r4][bond],
                           self._adj[r4][r5][bond], self._adj[r5][r6][bond], self._adj[r6][r1][bond])
                    if 4 not in key:
                        continue

                    doubles = tuple(y for y, x in enumerate(r) if len(self._adj[x]) == 3 and
                                    next(attr[bond] for a, attr in self._adj[x].items() if a not in r) == 2)
                    if not doubles:
                        continue

                    if len(doubles) == 6:
                        self._adj[r1][r2][bond] = self._adj[r2][r3][bond] = self._adj[r3][r4][bond] = 1
                        self._adj[r4][r5][bond] = self._adj[r5][r6][bond] = self._adj[r6][r1][bond] = 1
                        found.append(n)
                    else:
                        if key in _quinone_pattern.get(doubles, {}):
                            dear = _quinone_fix.get(doubles)
                            self._adj[r1][r2][bond], self._adj[r2][r3][bond], self._adj[r3][r4][bond], \
                                self._adj[r4][r5][bond], self._adj[r5][r6][bond], self._adj[r6][r1][bond] = dear
                            found.append(n)
                elif len(r) == 5:
                    r1, r2, r3, r4, r5 = r
                    key = (self._adj[r1][r2][bond], self._adj[r2][r3][bond], self._adj[r3][r4][bond],
                           self._adj[r4][r5][bond], self._adj[r5][r1][bond])
                    if 4 not in key:
                        continue

                    positions = _pyrole_pattern.get(key)
                    if positions is None:
                        continue

                    for m, pos in enumerate(positions):
                        if self._node[r[pos]]._atom in _pyrole_atoms:
                            dear = _pyrole_fix[key][m]
                            self._adj[r1][r2][bond], self._adj[r2][r3][bond], self._adj[r3][r4][bond], \
                                self._adj[r4][r5][bond], self._adj[r5][r1][bond] = dear
                            found.append(n)

            for n in found[::-1]:
                del rings[n]
            new = len(rings)
        return init - old

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
        return False


def _clock(a):
    yield a
    for _ in range(1, len(a)):
        a = a[1:] + a[:1]
        yield a


_pyrole_atoms = ('N', 'O', 'S', 'Se', 'P', C(-1))
_benzene = set()
_pyrole_pattern = defaultdict(list)
_pyrole_fix = defaultdict(list)
_pyrole = {}
_quinone_pattern = {}
_quinone_fix = {}


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
