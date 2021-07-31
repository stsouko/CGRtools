# -*- coding: utf-8 -*-
#
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from lazy_object_proxy import Proxy
from ...periodictable import ListElement, AnyMetal


def _rules():
    from ...containers import QueryContainer
    rules = []

    #
    #        R - N - C                 R - N - C
    #           /    ||                   /    ||
    #  A - M = C     || >>    A - M .. [C-]    ||
    #           \    ||                  \\    ||
    #        R - N - C               R - [N+]- C
    #
    q = QueryContainer()
    q.add_atom(AnyMetal())
    q.add_atom('C')
    q.add_atom('N', hybridization=1, neighbors=3, heteroatoms=0)
    q.add_atom('N', hybridization=1, neighbors=3, heteroatoms=0)
    q.add_atom('C', hybridization=2)
    q.add_atom('C', hybridization=2)
    q.add_bond(1, 2, 2)
    q.add_bond(2, 3, 1)
    q.add_bond(2, 4, 1)
    q.add_bond(3, 5, 1)
    q.add_bond(4, 6, 1)
    q.add_bond(5, 6, (1, 2))

    atom_fix = {2: (-1, False), 3: (1, False)}  # atom: (charge diff, radical)
    bonds_fix = ((1, 2, 8), (2, 3, 2))
    rules.append((q, atom_fix, bonds_fix))

    # Ferrocene covalent charge-free
    #  C5H5-(5)Fe(5)-C5H5
    #
    atoms = ({'atom': ListElement(['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Zr', 'Nb', 'Mo', 'Ru', 'Hf', 'W', 'Re', 'Ir']),
                     'neighbors': 10}, {'atom': 'C'}, {'atom': 'C'}, {'atom': 'C'}, {'atom': 'C'}, {'atom': 'C'},
             {'atom': 'C'}, {'atom': 'C'}, {'atom': 'C'}, {'atom': 'C'}, {'atom': 'C'})
    bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 5, 1), (1, 6, 1), (1, 7, 1), (1, 8, 1), (1, 9, 1), (1, 10, 1), (1, 11, 1),
             (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 6, 1), (6, 2, 1), (7, 8, 1), (8, 9, 1), (9, 10, 1), (10, 11, 1),
             (11, 7, 1))
    atom_fix = {1: {'charge': 2}, 2: {'charge': -1}, 3: {'hybridization': 2}, 4: {'hybridization': 2},
                5: {'hybridization': 2}, 6: {'hybridization': 2}, 7: {'charge': -1}, 8: {'hybridization': 2},
                9: {'hybridization': 2}, 10: {'hybridization': 2}, 11: {'hybridization': 2}}
    bonds_fix = (
    (1, 2, 8), (1, 3, 8), (1, 4, 8), (1, 5, 8), (1, 6, 8), (1, 7, 8), (1, 8, 8), (1, 9, 8), (1, 10, 8), (1, 11, 8),
    (3, 4, 2), (5, 6, 2), (8, 9, 2), (10, 11, 2))

    # Ferrocene covalent explicit H
    #  C5H5-(5)Fe(5)-C5H5
    #
    atoms = ({'atom': ListElement(['Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Zr', 'Nb', 'Mo', 'Ru', 'Hf', 'W', 'Re', 'Ir']),
                     'charge': 2, 'neighbors': 10}, {'atom': 'C', 'charge': -1}, {'atom': 'C'}, {'atom': 'C'},
             {'atom': 'C'}, {'atom': 'C'}, {'atom': 'C', 'charge': -1}, {'atom': 'C'}, {'atom': 'C'}, {'atom': 'C'},
             {'atom': 'C'})
    bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 5, 1), (1, 6, 1), (1, 7, 1), (1, 8, 1), (1, 9, 1), (1, 10, 1), (1, 11, 1),
             (2, 3, 1), (3, 4, 2), (4, 5, 1), (5, 6, 2), (6, 2, 1), (7, 8, 1), (8, 9, 2), (9, 10, 1), (10, 11, 2),
             (11, 7, 1))
    atom_fix = {}
    bonds_fix = (
    (1, 2, 8), (1, 3, 8), (1, 4, 8), (1, 5, 8), (1, 6, 8), (1, 7, 8), (1, 8, 8), (1, 9, 8), (1, 10, 8), (1, 11, 8))
    return rules


rules = Proxy(_rules)


__all__ = ['rules']
