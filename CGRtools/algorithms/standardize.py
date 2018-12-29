# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
#  Copyright 2018 Tagir Akhmetshin <tagirshin@gmail.com>
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
from ..attributes import QueryAtom, Bond


class Standardize:
    def standardize(self):
        """
        standardize functional groups

        :return: number of found groups
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
    doubles = []
    for q, p, c in zip(permutations(q), permutations(p[1:]), repeat(p[0])):
        if q in doubles:
            continue
        doubles.append(q)
        yield q, p, c


central = {}
query_patch = defaultdict(list)

# patterns
b1 = Bond()
b2 = Bond()
b3 = Bond()
b4 = Bond()
b2.order = 2
b3.order = 3
b4.order = 4

# any atom need for full neighbors describing
a = QueryAtom()

# =O
o2 = QueryAtom()
o2.update(element='O', hybridization=2)

# 1. Nitro
#
#       O          O-
#      //         /
#  A - N  >> A - N+
#      \\        \\
#       O         O
#
n33 = QueryAtom()
n33.update(element='N', neighbors=3, hybridization=3)
central['N3;3'] = n33
query_patch['N3;3'].extend(_prepare([(b2, o2), (b2, o2), (b1, a)], [{'charge': 1}, ({'order': 1}, {'charge': -1})]))


# 2. Aromatic N-Oxide
#
#  A : N : A  >>  A : N+ : A
#      \\             |
#       O             O-
#
n = QueryAtom()
n.update(element='N', neighbors=3, hybridization=4)
central['N3;4'] = n
query_patch['N3;4'].extend(_prepare([(b2, o2), (b4, a), (b4, a)], [{'charge': 1}, ({'order': 1}, {'charge': -1})]))


# 3. Azide
#
#  N- - N+ # N  >>  N = N+ = N-
#
n1 = QueryAtom()
n2 = QueryAtom()
n3 = QueryAtom()
n1.update(element='N', charge=-1, neighbors=2, hybridization=1)
n2.update(element='N', charge=1, neighbors=2, hybridization=3)
n3.update(element='N', neighbors=1)
central['N+2;3'] = n2
query_patch['N+2;3'].extend(_prepare([(b1, n1), (b3, n3)],
                                     [{}, ({'order': 2}, {'charge': 0}), ({'order': 2}, {'charge': -1})]))


# 4. Diazo
#
#  C- - N+ # N  >>  C = N+ = N-
#
c = QueryAtom()
n1 = n2  # same as in azide
n2 = QueryAtom()
c.update(element='C', charge=-1, hybridization=1)
n2.update(element='N', neighbors=1)


central['N+2;3'] = n1
query_patch['N+2;3'].extend(_prepare([(b1, c), (b3, n2)],
                                     [{}, ({'order': 2}, {'charge': 0}), ({'order': 2}, {'charge': -1})]))


# 5. Diazonium
#
#  C - N = N+  >>  C - N+ # N
#
c = QueryAtom()
n1 = QueryAtom()
n2 = QueryAtom()
c.update(element='C')
n1.update(element='N', neighbors=2, hybridization=2)
n2.update(element='N', charge=1, neighbors=1)
central['N2;2'] = n1
query_patch['N2;2'].extend(_prepare([(b2, n2), (b1, c)], [{'charge': 1}, ({'order': 3}, {'charge': 0})]))


# 6. Iminium
#
#  C+ - N  >> C = N+
#
c = QueryAtom()
n = QueryAtom()
c.update(element='C', charge=1, hybridization=1)
n.update(element='N', neighbors=3, hybridization=1)
central['N3;1'] = n
query_patch['N3;1'].extend(_prepare([(b1, c), (b1, a), (b1, a)], [{'charge': 1}, ({'order': 2}, {'charge': 0})]))


# 7. Isocyanate
#
#  N+ - C- = O  >>  N = C = O
#
n = QueryAtom()
c = QueryAtom()
n.update(element='N', charge=1, neighbors=(1, 2), hybridization=1)
c.update(element='C', charge=-1, neighbors=2, hybridization=2)
central['C-2;2'] = c
query_patch['C-2;2'].extend(_prepare([(b1, n), (b2, o2)], [{'charge': 0}, ({'order': 2}, {'charge': 0})]))


# 8. Nitrilium
#
#  C+ = N  >>  C # N+
#
c = QueryAtom()
n = QueryAtom()
c.update(element='C', charge=1, neighbors=2, hybridization=2)
n.update(element='N', neighbors=(1, 2), hybridization=2)
central['C+2;2'] = c
query_patch['C+2;2'].extend(_prepare([(b2, n), (b1, a)], [{'charge': 0}, ({'order': 3}, {'charge': 1})]))


# 9. Nitrone Nitronate
#
#      O          O-
#     //         /
# C = N  >> C = N+
#      \         \
#       A         A
#
c = QueryAtom()
c.update(element='C', hybridization=2)
query_patch['N3;3'].extend(_prepare([(b2, o2), (b2, c), (b1, a)], [{'charge': 1}, ({'order': 1}, {'charge': -1})]))


# 10. Nitroso
#
# C - N+ - O-  >>  C - N = O
c = QueryAtom()
n = QueryAtom()
o = QueryAtom()
c.update(element='C')
n.update(element='N', charge=1, neighbors=2, hybridization=1)
o.update(element='O', charge=-1)
central['N+2;1'] = n
query_patch['N+2;1'].extend(_prepare([(b1, o), (b1, c)], [{'charge': 0}, ({'order': 2}, {'charge': 0})]))


del b1, b2, b3, b4, a, c, o2, n, n1, n2, n3


__all__ = ['Standardize']
