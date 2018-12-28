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
from itertools import permutations, repeat
from ..attributes import QueryAtom, QueryBond


class Standardize:
    def standardize(self):
        """
        standardize functional groups

        :return: amoun of found groups
        """
        self.aromatize()
        self.reset_query_marks()
        seen = set()
        total = 0
        for k, center in central.items():
            for n, atom in self._node.items():
                if n in seen:
                    continue
                if center == atom:
                    shell = tuple((bond, self._node[m]) for m, bond in self._adj[n].items())
                    for shell_query, shell_patch, atom_patch in query_patch[k]:
                        if shell_query == shell:
                            total += 1
                            atom.update(atom_patch)
                            for (b_p, a_p), (b, a) in zip(shell_patch, shell):
                                b.update(b_p)
                                a.update(a_p)
                            seen.add(n)
                            seen.update(self._adj[n])
                            break
        return total


def _prepare(q, p):
    d = len(q) - len(p) + 1
    if d:
        p.extend([({}, {})] * d)
    return list(zip(permutations(q), permutations(p[1:]), repeat(p[0])))


central = {}
query_patch = defaultdict(list)

# patterns
b1 = QueryBond()
b2 = QueryBond()
b3 = QueryBond()
b2.order = 2
b3.order = 3

# any atom need for full neighbors describing
a = QueryAtom()
a.element = 'A'

o2 = QueryAtom()
o2.update(element='O', neighbors=1, hybridization=2)

# Nitro [A]-N(=O)=O>>[A]-[N+](=O)-[O-]
n = QueryAtom()
n.update(element='N', neighbors=3, hybridization=3)
central[str(n)] = n
query_patch[str(n)].extend(_prepare([(b2, o2), (b2, o2), (b1, a)], [{'charge': 1}, ({'order': 1}, {'charge': -1})]))

# Azide [A]-[N-]-[N+]#N>>[A]-N=[N+]=[N-]
n1 = QueryAtom()
n2 = QueryAtom()
n3 = QueryAtom()
n1.update(element='N', charge=-1, neighbors=2, hybridization=1)
n2.update(element='N', charge=1, neighbors=2, hybridization=3)
n3.update(element='N', neighbors=1, hybridization=3)
central[str(n2)] = n2
query_patch[str(n2)].extend(_prepare([(b1, n1), (b3, n3)],
                                     [{}, ({'order': 2}, {'charge': 0}), ({'order': 2}, {'charge': -1})]))

del b1, b2, b3, a, o2, n, n1, n2, n3


__all__ = ['Standardize']
