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
from itertools import permutations
from ..attributes import QueryAtom, QueryBond


class Standardize:
    def __init__(self, query, patch):
        self.__center = query[0]
        self.__center_patch = patch[0]
        self.__shell_patch = list(zip(permutations(query[1:]), permutations(patch[1:])))

    def __call__(self, structure):
        seen = set()
        for n, c in structure._node.items():
            if n in seen:
                continue
            if self.__center == c:
                shell = tuple((bond, structure._node[m]) for m, bond in structure._adj[n].items())
                for q, p in self.__shell_patch:
                    if q == shell:
                        c.update(self.__center_patch)
                        for (b_p, a_p), (b, a) in zip(p, shell):
                            b.update(b_p)
                            a.update(a_p)
                        seen.add(n)
                        seen.update(structure._adj[n])
                        break


class Nitro(Standardize):
    def __init__(self):
        n = QueryAtom()
        n.update(element='N', neighbors=3, hybridization=3)
        o = QueryAtom()
        o.update(element='O', neighbors=1, hybridization=2)
        a = QueryAtom()
        a.update(element='A')
        b = QueryBond()
        b.order = 2

        super().__init__((n, (b, o), (b, o), (QueryBond(), a)),
                         ({'charge': 1}, ({'order': 1}, {'charge': -1}), ({}, {}), ({}, {})))


__all__ = ['Nitro']
