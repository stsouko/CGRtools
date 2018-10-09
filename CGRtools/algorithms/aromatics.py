# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
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
from .sssr import find_sssr


def aromatize(g):
    rings = [x for x in find_sssr(g) if 4 < len(x) < 7]
    if not rings:
        return 0
    total = 0
    while True:
        c = _quinonize(g, rings, 's_bond')
        if c:
            total += c
        elif total:
            break

        c = _aromatize(g, rings, 's_bond')
        if not c:
            break
        total += c

    return total


def aromatize_cgr(g):
    rings = [x for x in find_sssr(g) if 4 < len(x) < 7]
    if not rings:
        return 0, 0
    total_s = total_p = 0
    while True:
        c_s = _quinonize(g, rings, 's_bond')
        c_p = _quinonize(g, rings, 'p_bond')
        if c_s or c_p:
            total_s += c_s
            total_p += c_p
        elif total_s or total_p:
            break

        c_s = _aromatize(g, rings, 's_bond')
        c_p = _aromatize(g, rings, 'p_bond')
        if not (c_s or c_p):
            break
        total_s += c_s
        total_p += c_p

    return total_s, total_p


def _quinonize(g, rings, bond):
    rings = rings.copy()
    init = len(rings)
    old = 0
    while len(rings) != old:
        old = len(rings)
        found = []
        for n, r in enumerate(rings):
            if len(r) == 6:
                r1, r2, r3, r4, r5, r6 = r
                key = (g[r1][r2][bond], g[r2][r3][bond], g[r3][r4][bond], g[r4][r5][bond],
                       g[r5][r6][bond], g[r6][r1][bond])
                if 4 not in key:
                    continue

                doubles = tuple(y for y, x in enumerate(r)
                                if len(g[x]) == 3 and next(attr[bond] for a, attr in g[x].items() if a not in r) == 2)
                if not doubles:
                    continue

                if len(doubles) == 6:
                    g[r1][r2][bond] = g[r2][r3][bond] = g[r3][r4][bond] = 1
                    g[r4][r5][bond] = g[r5][r6][bond] = g[r6][r1][bond] = 1
                    found.append(n)
                else:
                    if key in _quinone_pattern.get(doubles, {}):
                        dear = _quinone_fix.get(doubles)
                        g[r1][r2][bond], g[r2][r3][bond], g[r3][r4][bond], g[r4][r5][bond], \
                        g[r5][r6][bond], g[r6][r1][bond] = dear
                        found.append(n)
            elif len(r) == 5:
                r1, r2, r3, r4, r5 = r
                key = (g[r1][r2][bond], g[r2][r3][bond], g[r3][r4][bond], g[r4][r5][bond],
                       g[r5][r1][bond])
                if 4 not in key:
                    continue

                positions = _pyrole_pattern.get(key)
                if positions is None:
                    continue

                for m, pos in enumerate(positions):
                    if g.nodes[r[pos]]['element'] in ('N', 'O', 'S', 'Se', 'P'):
                        dear = _pyrole_fix[key][m]
                        g[r1][r2][bond], g[r2][r3][bond], g[r3][r4][bond], g[r4][r5][bond], \
                        g[r5][r1][bond] = dear
                        found.append(n)

        for n in found[::-1]:
            del rings[n]
    return init - old


def _aromatize(g, rings, bond):
    rings = rings.copy()
    init = len(rings)
    old = 0
    while len(rings) != old:
        old = len(rings)
        found = []
        for n, r in enumerate(rings):
            if len(r) == 6:
                r1, r2, r3, r4, r5, r6 = r
                if (g[r1][r2][bond], g[r2][r3][bond], g[r3][r4][bond], g[r4][r5][bond],
                    g[r5][r6][bond], g[r6][r1][bond]) in _benzene:
                    g[r1][r2][bond] = g[r2][r3][bond] = g[r3][r4][bond] = 4
                    g[r4][r5][bond] = g[r5][r6][bond] = g[r6][r1][bond] = 4
                    found.append(n)
            elif len(r) == 5:
                r1, r2, r3, r4, r5 = r
                position = _pyrole.get((g[r1][r2][bond], g[r2][r3][bond], g[r3][r4][bond],
                                        g[r4][r5][bond], g[r5][r1][bond]))

                if position is not None and g.nodes[r[position]]['element'] in ('N', 'O', 'S', 'Se', 'P'):
                    g[r1][r2][bond] = g[r2][r3][bond] = g[r3][r4][bond] = 4
                    g[r4][r5][bond] = g[r5][r1][bond] = 4
                    found.append(n)

        for n in found[::-1]:
            del rings[n]
    return init - old


def _clock(a):
    yield a
    for _ in range(1, len(a)):
        a = a[1:] + a[:1]
        yield a


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


_quinone_pattern = {}
_quinone_fix = {}


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
