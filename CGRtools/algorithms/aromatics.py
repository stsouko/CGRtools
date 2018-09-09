# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
from .sssr import find_sssr


def aromatize(g):
    rings = [x for x in find_sssr(g) if 4 < len(x) < 7]
    total = 0
    while True:
        c = _quinonize(g, rings)
        if c:
            total += c
        elif total:
            break

        c = _aromatize(g, rings)
        if not c:
            break
        total += c

    return total


def _quinonize(g, rings):
    rings = rings.copy()
    init = len(rings)
    old = 0
    while len(rings) != old:
        old = len(rings)
        found = []
        for n, r in enumerate(rings):
            if len(r) != 6:
                continue

            r1, r2, r3, r4, r5, r6 = r
            key = (g[r1][r2]['s_bond'], g[r2][r3]['s_bond'], g[r3][r4]['s_bond'], g[r4][r5]['s_bond'],
                   g[r5][r6]['s_bond'], g[r6][r1]['s_bond'])
            if 4 not in key:
                continue

            doubles = tuple(y for y, x in enumerate(r)
                            if len(g[x]) == 3 and next(attr['s_bond'] for a, attr in g[x].items() if a not in r) == 2)
            if not doubles:
                continue

            ld = len(doubles)
            if ld == 6:
                g[r1][r2]['s_bond'] = g[r2][r3]['s_bond'] = g[r3][r4]['s_bond'] = 1
                g[r4][r5]['s_bond'] = g[r5][r6]['s_bond'] = g[r6][r1]['s_bond'] = 1
                found.append(n)
            else:
                if key in _quinone_pattern.get(doubles, {}):
                    dear = _quinone_fix.get(doubles)
                    g[r1][r2]['s_bond'], g[r2][r3]['s_bond'], g[r3][r4]['s_bond'], g[r4][r5]['s_bond'], \
                    g[r5][r6]['s_bond'], g[r6][r1]['s_bond'] = dear
                    found.append(n)

        for n in found[::-1]:
            del rings[n]
    return init - old


def _aromatize(g, rings):
    rings = rings.copy()
    init = len(rings)
    old = 0
    while len(rings) != old:
        old = len(rings)
        found = []
        for n, r in enumerate(rings):
            if len(r) == 6:
                r1, r2, r3, r4, r5, r6 = r
                if (g[r1][r2]['s_bond'], g[r2][r3]['s_bond'], g[r3][r4]['s_bond'], g[r4][r5]['s_bond'],
                    g[r5][r6]['s_bond'], g[r6][r1]['s_bond']) in _benzene:
                    g[r1][r2]['s_bond'] = g[r2][r3]['s_bond'] = g[r3][r4]['s_bond'] = 4
                    g[r4][r5]['s_bond'] = g[r5][r6]['s_bond'] = g[r6][r1]['s_bond'] = 4
                    found.append(n)
            elif len(r) == 5:
                r1, r2, r3, r4, r5 = r
                k = _pyrole.get((g[r1][r2]['s_bond'], g[r2][r3]['s_bond'], g[r3][r4]['s_bond'], g[r4][r5]['s_bond'],
                                 g[r5][r1]['s_bond']))

                if k is not None and g.nodes[r[k]]['element'] in ('N', 'O', 'S', 'Se'):
                    g[r1][r2]['s_bond'] = g[r2][r3]['s_bond'] = g[r3][r4]['s_bond'] = 4
                    g[r4][r5]['s_bond'] = g[r5][r1]['s_bond'] = 4
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


_quinone_pattern = {}
_quinone_fix = {}


_ind = ((0, 1), (0, 5), (4, 5), (3, 4), (2, 3), (1, 2))
for i, *p in zip(_ind,
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
    _quinone_pattern[i] = {(4, 4, 4, 4, 4, 4), *p}

_quinone_fix.update(zip(_ind, _clock((1, 1, 2, 1, 2, 1))))


_ind = ((0, 1, 2, 3), (0, 1, 2, 5), (0, 1, 4, 5), (0, 3, 4, 5), (2, 3, 4, 5), (1, 2, 3, 4))
for i, *p in zip(_ind,
                 _clock((1, 4, 1, 4, 4, 4)),
                 _clock((1, 4, 1, 4, 2, 4))
                 ):
    _quinone_pattern[i] = {(4, 4, 4, 4, 4, 4), *p}

_quinone_fix.update(zip(_ind, _clock((1, 1, 1, 1, 2, 1))))
