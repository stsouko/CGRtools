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
b14 = QueryBond()
b124 = QueryBond()
b2.order = 2
b3.order = 3
b14.order = (1, 4)
b124.order = (1, 2, 4)

# any atom need for full neighbors describing
a = QueryAtom()
a.element = 'A'

o2 = QueryAtom()
o2.update(element='O', neighbors=1, hybridization=2)

# N-Oxide [A]<-:>N(<-=:>[A])=O>>[A]-[N+](=O)-[O-]
n = QueryAtom()
n.update(element='N', neighbors=3, hybridization=(2, 3, 4))
central[str(n)] = n
query_patch[str(n)].extend(_prepare([(b2, o2), (b14, a), (b124, a)], [{'charge': 1}, ({'order': 1}, {'charge': -1})]))


# Azide [N-]-[N+]#N>>N=[N+]=[N-]
n1 = QueryAtom()
n2 = QueryAtom()
n3 = QueryAtom()
n1.update(element='N', charge=-1, neighbors=2, hybridization=1)
n2.update(element='N', charge=1, neighbors=2, hybridization=3)
n3.update(element='N', neighbors=1, hybridization=3)
central[str(n2)] = n2
query_patch[str(n2)].extend(_prepare([(b1, n1), (b3, n3)],
                                     [{}, ({'order': 2}, {'charge': 0}), ({'order': 2}, {'charge': -1})]))

# Diazo [C]-[N+]#N>>[C]-N=[N+]=[N-]
c1 = QueryAtom()
n2 = QueryAtom()
n3 = QueryAtom()
n3.update(element='N', neighbors=1, hybridization=3)
n2.update(element='N', charge=1, neighbors=2, hybridization=3)
c1.update(element='C', charge=-1, hybridization=1)
central[str(n2)] = n2
query_patch[str(n2)].extend(_prepare([(b3, n3), (b1, c1)],
                                     [{}, ({'order': 2}, {'charge': -1}), ({'order': 2}, {'charge': 0})]))

# Diazonium [C]-[N]=[N+]>>C-[N+]#[N]
c1 = QueryAtom()
n2 = QueryAtom()
n3 = QueryAtom()
c1.update(element='C')
n2.update(element='N', neighbors=2, hybridization=2)
n3.update(element='N', charge=1, neighbors=1, hybridization=2)
central[str(n2)] = n2
query_patch[str(n2)].extend(_prepare([(b2, n2), (b2, n3)],
                                     [{}, ({'order': 3}, {'charge': 1}), ({'order': 3}, {'charge': 0})]))

# Iminium [C+]-N>>C=[N+]
c1 = QueryAtom()
n2 = QueryAtom()
c1.update(element='C', charge=1, neighbors=3, hybridization=1)
n2.update(element='N', neighbors=3, hybridization=1)
central[str(n2)] = n2
query_patch[str(n2)].extend(_prepare([(b1, c1)], [{'charge': 1}, ({'order': 2}, {'charge': 0})]))

# Isocyanate [N+]-[C-]=O>>N=C=O
n1 = QueryAtom()
c2 = QueryAtom()
o3 = QueryAtom()
n1.update(element='N', charge=1, neighbors=2, hybridization=1)
c2.update(element='C', charge=-1, neighbors=2, hybridization=2)
o3.update(element='O', hybridization=2)
central[str(c2)] = c2
query_patch[str(c2)].extend(_prepare([(b1, n1)], [{'charge': 0}, ({'order': 2}, {'charge': 0})]))

# Nitrilium [C+]=N>>C#[N+]
c1 = QueryAtom()
n2 = QueryAtom()
c1.update(element='C', charge=1, neighbors=2, hybridization=2)
n2.update(element='N', neighbors=2, hybridization=2)
central[str(n2)] = n2
query_patch[str(n2)].extend(_prepare([(b1, c1)], [{'charge': 1}, ({'order': 3}, {'charge': 0})]))

# Nitrone Nitronate [A]-N(=O)=O>>[A]-[N+](=A)-[O-]
c1 = QueryAtom()
n2 = QueryAtom()
o3 = QueryAtom()
c1.update(element='C', neighbors=3, hybridization=2)
n2.update(element='N', neighbors=2, hybridization=1)
o3.update(element='O', hybridization=2)
central[str(n2)] = n2
query_patch[str(n2)].extend(_prepare([(b2, o3)], [{'charge': 1}, ({'order': 1}, {'charge': -1})]))

# Nitroso
c1 = QueryAtom()
n2 = QueryAtom()
o3 = QueryAtom()
c1.update(element='C', hybridization=(1, 2, 4))
n2.update(element='N', charge=1, neighbors=2, hybridization=1)
o3.update(element='O', charge=-1, hybridization=1)
central[str(n2)] = n2
query_patch[str(n2)].extend(_prepare([(b2, o3)], [{'charge': 0}, ({'order': 2}, {'charge': 0})]))

#Phosphonic
c1 = QueryAtom()
p2 = QueryAtom()
o3 = QueryAtom()
o4 = QueryAtom()
o5 = QueryAtom()
c1.update(element='C',  hybridization=(1, 2, 4))
p2.update(element='P', charge=1, neighbors=4, hybridization=1)
o3.update(element='O', neighbors=2, hybridization=1)
o4.update(element='O', neighbors=2, hybridization=1)
o5.update(element='O', charge=-1, neighbors=1, hybridization=1)
central[str(p2)] = p2
query_patch[str(p2)].extend(_prepare([(b1, o5)], [{'charge': 0}, ({'order': 2}, {'charge': 0})]))

#Phosphonic Ylide
c1 = QueryAtom()
p2 = QueryAtom()
c3 = QueryAtom()
c4 = QueryAtom()
c5 = QueryAtom()
c1.update(element='C')
p2.update(element='P', charge=-1, neighbors=4, hybridization=1)
c3.update(element='C')
c4.update(element='C')
c5.update(element='C', charge=1,  hybridization=1)
central[str(p2)] = p2
query_patch[str(p2)].extend(_prepare([(b1, c5)], [{'charge': 0}, ({'order': 2}, {'charge': 0})]))

# N=N#N >> N=[N+]=[N-]
n1 = QueryAtom()
n2 = QueryAtom()
n3 = QueryAtom()
n1.update(element='N', hybridization=2)
n2.update(element='N', neighbors=2, hybridization=3)
n3.update(element='N', neighbors=1, hybridization=3)
central[str(n2)] = n2
query_patch[str(n2)].extend(_prepare([(b3, n3)],
                                     [{'charge': 1}, ({'order': 2}, {'charge': -1})]))


# N=N(C)=0 >> N=[N+](C)-[0-]
n1 = QueryAtom()
n2 = QueryAtom()
c3 = QueryAtom()
o4 = QueryAtom()
n1.update(element='N', hybridization=2)
n2.update(element='N', neighbors=3, hybridization=2)
c3.update(element='C')
o4.update(element='O', neighbors=1, hybridization=2)
central[str(n2)] = n2
query_patch[str(n2)].extend(_prepare([(b2, o4)],
                                     [{'charge': 1}, ({'order': 1}, {'charge': -1})]))


# [N-]-[N+](C)=0 >> N=[N+](C)-[0-]
n1 = QueryAtom()
n2 = QueryAtom()
c3 = QueryAtom()
o4 = QueryAtom()
n1.update(element='N', hybridization=1, charge=-1)
n2.update(element='N', neighbors=3, hybridization=2, charge=1)
c3.update(element='C')
o4.update(element='O', neighbors=1, hybridization=2)
central[str(n2)] = n2
query_patch[str(n2)].extend(_prepare([(b2, o4)],
                                     [{'charge': 1}, ({'order': 1}, {'charge': -1})]))


del b1, b2, b3, a, o2, n, n1, n2, n3


__all__ = ['Standardize']
