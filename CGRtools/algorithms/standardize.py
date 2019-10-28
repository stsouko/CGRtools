# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
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
        self.reset_query_marks()
        seen = set()
        total = 0
        for n, atom in self.atoms():
            if n in seen:
                continue
            for k, center in central.items():
                if center != atom:
                    continue
                shell = tuple((bond, self._node[m]) for m, bond in self._adj[n].items())
                for shell_query, shell_patch, atom_patch in query_patch[k]:
                    if shell_query != shell:
                        continue
                    total += 1
                    for attr_name, attr_value in atom_patch.items():
                        setattr(atom, attr_name, attr_value)
                    for (bond_patch, atom_patch), (bond, atom) in zip(shell_patch, shell):
                        bond.update(bond_patch)
                        for attr_name, attr_value in atom_patch.items():
                            setattr(atom, attr_name, attr_value)
                    seen.add(n)
                    seen.update(self._adj[n])
                    break
                else:
                    continue
                break
        if total:
            self.flush_cache()
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


# 1. Nitro
#
#       O          O-
#      //         /
#  A - N  >> A - N+
#      \\        \\
#       O         O
#
a = QueryAtom()
o = QueryAtom()
n33 = QueryAtom()
o.update(element='O')
n33.update(element='N', neighbors=3, hybridization=3)
central['N3;3'] = n33
query_patch['N3;3'].extend(_prepare([(b2, o), (b2, o), (b1, a)],
                                    [{'charge': 1, '_hybridization': 2},
                                     ({'order': 1}, {'charge': -1, '_hybridization': 1})]))


# 2. Aromatic N-Oxide
#
#  A : N : A  >>  A : N+ : A
#      \\             |
#       O             O-
#
n34 = QueryAtom()
n34.update(element='N', neighbors=3, hybridization=4)
central['N3;4'] = n34
query_patch['N3;4'].extend(_prepare([(b2, o), (b4, a), (b4, a)],
                                    [{'charge': 1},
                                     ({'order': 1}, {'charge': -1, '_hybridization': 1})]))


# 3. Azide
#
#  N- - N+ # N  >>  N = N+ = N-
#
nn21 = QueryAtom()
np23 = QueryAtom()
n1_ = QueryAtom()
nn21.update(element='N', charge=-1, neighbors=2, hybridization=1)
np23.update(element='N', charge=1, neighbors=2, hybridization=3)
n1_.update(element='N', neighbors=1)
central['N+2;3'] = np23
query_patch['N+2;3'].extend(_prepare([(b1, nn21), (b3, n1_)],
                                     [{}, ({'order': 2}, {'charge': 0, '_hybridization': 2}),
                                      ({'order': 2}, {'charge': -1, '_hybridization': 2})]))


# 3.1. Azide ChemAxoned
#
#  N+ # N = N-  >>  N = N+ = N-
#
n23 = QueryAtom()
nn12 = QueryAtom()
n23.update(element='N', neighbors=2, hybridization=3)
nn12.update(element='N', charge=-1, neighbors=1, hybridization=2)
central['N2;3'] = n23
query_patch['N2;3'].extend(_prepare([(b3, np23), (b2, nn12)],
                                    [{'charge': 1}, ({'order': 2}, {'charge': 0, '_hybridization': 2})]))


# 4. Diazo
#
#  C- - N+ # N  >>  C = N+ = N-
#
cn_1 = QueryAtom()
cn_1.update(element='C', charge=-1, hybridization=1)
query_patch['N+2;3'].extend(_prepare([(b1, cn_1), (b3, n1_)],
                                     [{}, ({'order': 2}, {'charge': 0, '_hybridization': 2}),
                                      ({'order': 2}, {'charge': -1, '_hybridization': 2})]))


# 5. Diazonium
#
#  C - N = N+  >>  C - N+ # N
#
c = QueryAtom()
n22 = QueryAtom()
np1_ = QueryAtom()
c.update(element='C')
n22.update(element='N', neighbors=2, hybridization=2)
np1_.update(element='N', charge=1, neighbors=1)
central['N2;2'] = n22
query_patch['N2;2'].extend(_prepare([(b2, np1_), (b1, c)],
                                    [{'charge': 1, '_hybridization': 3},
                                     ({'order': 3}, {'charge': 0, '_hybridization': 3})]))


# 6. Iminium
#
#  C+ - N  >> C = N+
#
cp_1 = QueryAtom()
n31 = QueryAtom()
cp_1.update(element='C', charge=1, hybridization=1)
n31.update(element='N', neighbors=3, hybridization=1)
central['N3;1'] = n31
query_patch['N3;1'].extend(_prepare([(b1, cp_1), (b1, a), (b1, a)],
                                    [{'charge': 1, '_hybridization': 2},
                                     ({'order': 2}, {'charge': 0, '_hybridization': 2})]))


# 7. Isocyanate
#
#  N+ - C- = O  >>  N = C = O
#
np121 = QueryAtom()
cn22 = QueryAtom()
np121.update(element='N', charge=1, neighbors=(1, 2), hybridization=1)
cn22.update(element='C', charge=-1, neighbors=2, hybridization=2)
central['C-2;2'] = cn22
query_patch['C-2;2'].extend(_prepare([(b1, np121), (b2, o)],
                                     [{'charge': 0, '_hybridization': 3},
                                      ({'order': 2}, {'charge': 0, '_hybridization': 2})]))


# 8. Nitrilium
#
#  C+ = N  >>  C # N+
#
cp22 = QueryAtom()
n122 = QueryAtom()
cp22.update(element='C', charge=1, neighbors=2, hybridization=2)
n122.update(element='N', neighbors=(1, 2), hybridization=2)
central['C+2;2'] = cp22
query_patch['C+2;2'].extend(_prepare([(b2, n122), (b1, a)],
                                     [{'charge': 0, '_hybridization': 3},
                                      ({'order': 3}, {'charge': 1, '_hybridization': 3})]))


# 9. Nitrone
#
#      O          O-
#     //         /
# C = N  >> C = N+
#      \         \
#       C         C
#
c_2 = QueryAtom()
c_2.update(element='C', hybridization=2)
query_patch['N3;3'].extend(_prepare([(b2, o), (b2, c_2), (b1, c)],
                                    [{'charge': 1, '_hybridization': 2},
                                     ({'order': 1}, {'charge': -1, '_hybridization': 1})]))


# 10. Nitronate
#
#      O         O
#     //        //
# C = N  >> C - N+
#      \         \
#       OH        O-
#
o1_ = QueryAtom()
o1_.update(element='O', neighbors=1)
query_patch['N3;3'].extend(_prepare([(b2, c_2), (b1, o1_), (b2, o)],
                                    [{'charge': 1, '_hybridization': 2},
                                     ({'order': 1}, {'_hybridization': 1}), ({}, {'charge': -1})]))


# 10.1 Nitronate ChemAxoned
#
#       O-         O-
#      /          /
# C = N+  >>  C - N+
#      \          \\
#       OH         O
#
np32 = QueryAtom()
np32.update(element='N', charge=1, neighbors=3, hybridization=2)
on = QueryAtom()
on.update(element='O', charge=-1)
central['N+3;2'] = np32
query_patch['N+3;2'].extend(_prepare([(b2, c_2), (b1, o1_), (b1, on)],
                                     [{}, ({'order': 1}, {'_hybridization': 1}),
                                      ({'order': 2}, {'_hybridization': 2})]))


# 11. Nitroso
#
# C - N+ - O-  >>  C - N = O
#
np21 = QueryAtom()
on = QueryAtom()
np21.update(element='N', charge=1, neighbors=2, hybridization=1)
on.update(element='O', charge=-1)
central['N+2;1'] = np21
query_patch['N+2;1'].extend(_prepare([(b1, on), (b1, c)], [{'charge': 0, '_hybridization': 2},
                                                           ({'order': 2}, {'charge': 0, '_hybridization': 2})]))

# 12. Tetriary N-oxide
#
#      C              C
#      |              |
#  A - N - A  >>  A - N+ - A
#      \\             |
#       O             O-
#
n42 = QueryAtom()
n42.update(element='N', neighbors=4, hybridization=2)
central['N4;2'] = n42
query_patch['N4;2'].extend(_prepare([(b2, o), (b1, c), (b1, a), (b1, a)],
                                    [{'charge': 1, '_hybridization': 1},
                                     ({'order': 1}, {'charge': -1, '_hybridization': 1})]))


# 13. Phosphonic
#
#      O                O
#      |                |
#  C - P+ - O-  >>  C - P = O
#      |                |
#      O                O
#
pp41 = QueryAtom()
pp41.update(element='P', charge=1, neighbors=4, hybridization=1)
central['P+4;1'] = pp41
query_patch['P+4;1'].extend(_prepare([(b1, on), (b1, c), (b1, o), (b1, o)],
                                     [{'charge': 0, '_hybridization': 2},
                                      ({'order': 2}, {'charge': 0, '_hybridization': 2})]))


# 14. Phosphonium ylide
#
#      C                C
#      |                |
#  C - P- - C+  >>  C - P = C
#      |                |
#      C                C
#
pn41 = QueryAtom()
pn41.update(element='P', charge=-1, neighbors=4, hybridization=1)
central['P-4;1'] = pn41
query_patch['P-4;1'].extend(_prepare([(b1, cp_1), (b1, c), (b1, c), (b1, c)],
                                     [{'charge': 0, '_hybridization': 2},
                                      ({'order': 2}, {'charge': 0, '_hybridization': 2})]))

# 15. Silicate Selenite
#
#        O            O
#       /            /
# O- - Si+  >>  O = Si
#       \            \
#        O            O
#
sesip31 = QueryAtom()
sesip31.update(element=('Se', 'Si'), charge=1, neighbors=3, hybridization=1)
central['SeSi+3;1'] = sesip31
query_patch['SeSi+3;1'].extend(_prepare([(b1, on), (b1, o), (b1, o)],
                                        [{'charge': 0, '_hybridization': 2},
                                         ({'order': 2}, {'charge': 0, '_hybridization': 2})]))


# 16. Sulfine
#
#  O- - S+ = C  >>  O = S = C
#
sp22 = QueryAtom()
sp22.update(element='S', charge=1, neighbors=2, hybridization=2)
central['S+2;2'] = sp22
query_patch['S+2;2'].extend(_prepare([(b1, on), (b2, c_2)],
                                     [{'charge': 0, '_hybridization': 3},
                                      ({'order': 2}, {'charge': 0, '_hybridization': 2})]))


# 17. Sulfon
#
#        C           C
#       /           /
# O- - S+  >>  O = S
#       \           \
#        C           C
#
sp31 = QueryAtom()
sp31.update(element='S', charge=1, neighbors=3, hybridization=1)
central['S+3;1'] = sp31
query_patch['S+3;1'].extend(_prepare([(b1, on), (b1, c), (b1, c)],
                                     [{'charge': 0, '_hybridization': 2},
                                      ({'order': 2}, {'charge': 0, '_hybridization': 2})]))


# 18. Sulfonium ylide
#
#  C - S- - C+  >>  C - S = C
#      |                |
#      C                C
#
sn31 = QueryAtom()
sn31.update(element='S', charge=-1, neighbors=3, hybridization=1)
central['S-3;1'] = sn31
query_patch['S-3;1'].extend(_prepare([(b1, cp_1), (b1, c), (b1, c)],
                                     [{'charge': 0, '_hybridization': 2},
                                      ({'order': 2}, {'charge': 0, '_hybridization': 2})]))


# 19. Sulfoxide
#
#      C                C
#      |                |
#  O = S+ - O-  >>  O = S = O
#      |                |
#      C                C
#
sp42 = QueryAtom()
sp42.update(element='S', charge=1, neighbors=4, hybridization=2)
central['S+4;2'] = sp42
query_patch['S+4;2'].extend(_prepare([(b1, on), (b1, c), (b1, c), (b2, o)],
                                     [{'charge': 0, '_hybridization': 3},
                                      ({'order': 2}, {'charge': 0, '_hybridization': 2})]))


# 20. Sulfoxonium ylide
#
#      C                C
#      |                |
#  C = S+ - O-  >>  C = S = O
#      |                |
#      C                C
#
query_patch['S+4;2'].extend(_prepare([(b1, on), (b1, c), (b1, c), (b2, c_2)],
                                     [{'charge': 0, '_hybridization': 3},
                                      ({'order': 2}, {'charge': 0, '_hybridization': 2})]))


# 21
#
# N = N # N >> N = N+ = N-
#
central['N2;3'] = n23
query_patch['N2;3'].extend(_prepare([(b3, n1_), (b2, n22)],
                                    [{'charge': 1}, ({'order': 2}, {'charge': -1, '_hybridization': 2})]))


# 22
#
# C = N # N >> C = N+ = N-
#
central['N2;3'] = n23
query_patch['N2;3'].extend(_prepare([(b3, n1_), (b2, c_2)],
                                    [{'charge': 1}, ({'order': 2}, {'charge': -1, '_hybridization': 2})]))


# 23
#
# - N = N = N >> - N = N+ = N-
#
central['N2;3'] = n23
query_patch['N2;3'].extend(_prepare([(b2, n1_), (b2, n22)],
                                    [{'charge': 1}, ({}, {'charge': -1, '_hybridization': 2})]))


# 24
#
# N = N = O >> N = N+ - O-
#     |            |
#     C            C
#
o12 = QueryAtom()
o12.update(element='O', neighbors=1, hybridization=2)
central['N3;3'] = n33
query_patch['N3;3'].extend(_prepare([(b2, o12), (b2, n22), (b1, c)],
                                    [{'charge': 1}, ({'order': 1}, {'charge': -1, '_hybridization': 1})]))


# 25
#
# N- - N+ = O >> N = N+ - O-
#      |             |
#      C             C
#
central['N+3;2'] = np32
query_patch['N+3;2'].extend(_prepare([(b2, o12), (b1, nn21), (b1, c)],
                                     [{'charge': 1}, ({'order': 1}, {'charge': -1, '_hybridization': 1}),
                                     ({'order': 2}, {'charge': 0, '_hybridization': 2})]))


# 26
#
# - C # N = O >> - C # N+ - O-
#
с23 = QueryAtom()
с23.update(element='C', neighbors=2, hybridization=3)
central['N2;3'] = n23
query_patch['N2;3'].extend(_prepare([(b2, o12), (b3, с23)],
                                    [{'charge': 1}, ({'order': 1}, {'charge': -1, '_hybridization': 1})]))


# 27
#
# - N - N # N >> - N = N+ = N-
#   |
#   H
#
n21 = QueryAtom()
n21.update(element='N', neighbors=2, hybridization=1)
central['N2;3'] = n23
query_patch['N2;3'].extend(_prepare([(b1, n21), (b3, n1_)],
                                    [{'charge': 1}, ({'order': 2}, {'_hybridization': 2}),
                                     ({'order': 2}, {'charge': -1, '_hybridization': 2})]))

# 29
# Carbon Monoxide
#
# [CX1] = O  >> С- # O+
#
c_co = QueryAtom()
o_co = QueryAtom()
c_co.update(element='C', neighbors=1, hybridization=2, multiplicity=3)
o_co.update(element='O', neighbors=1, hybridization=2)
central['C1;2;3'] = c_co

query_patch['C1;2;3'].extend(_prepare([(b2, o_co)],
                                      [{'charge': -1, '_hybridization': 3, 'multiplicity': None},
                                       ({'order': 3}, {'charge': 1, '_hybridization': 3})]))

__all__ = ['Standardize']
