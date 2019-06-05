# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from CachedMethods import cached_property
from collections import defaultdict
from typing import Dict


class Isomorphism:
    __slots__ = ()

    def __lt__(self, other):
        if len(self) >= len(other):
            return False
        return self.is_substructure(other)

    def __le__(self, other):
        return self.is_substructure(other)

    def __gt__(self, other):
        if len(self) <= len(other):
            return False
        return other.is_substructure(self)

    def __ge__(self, other):
        return other.is_substructure(self)

    def is_substructure(self, other) -> bool:
        """
        test self is substructure of other
        """
        try:
            next(self.get_mapping(other))
        except StopIteration:
            return False
        return True

    def is_equal(self, other) -> bool:
        """
        test self is structure of other
        """
        if len(self) != len(other):
            return False
        try:
            next(self.get_mapping(other))
        except StopIteration:
            return False
        return True

    def get_mapping(self, other) -> Dict[int, int]:
        """
        get self to other substructure mapping generator
        """
        size = len(self._atoms) - 1
        order, closures = self._compiled_query
        order_depth = {v[0]: k for k, v in enumerate(order)}

        o_atoms = other._atoms
        o_charges = other._charges
        o_radicals = other._radicals
        o_bonds = other._bonds

        stack = []
        path = []
        mapping = {}
        reversed_mapping = {}

        s_atom, s_charge, s_is_radical = order[0][2:-1]
        for n, o_atom in o_atoms.items():
            if s_atom == o_atom and s_charge == o_charges[n] and s_is_radical == o_radicals[n]:
                stack.append((n, 0))

        while stack:
            o_atom, depth = stack.pop()
            s_atom = order[depth][0]
            if depth == size:
                yield {s_atom: o_atom, **mapping}
            else:
                if len(path) != depth:
                    for x in path[depth:]:
                        del mapping[reversed_mapping[x]]
                    path = path[:depth]

                back = order[depth + 1][2]
                if back != s_atom:
                    fork = path[order_depth[back]]
                else:
                    fork = o_atom

                path.append(o_atom)
                mapping[s_atom] = o_atom
                reversed_mapping[o_atom] = s_atom

                lp = len(path)
                for o_n, o_bond in o_bonds[fork].items():
                    s_n, _, s_atom, s_charge, s_is_radical, s_bond = order[lp]
                    if o_n not in path and s_bond == o_bond and s_atom == o_atoms[o_n] and s_charge == o_charges[o_n] \
                            and s_is_radical == o_radicals[o_n] \
                            and all(bond == o_bonds[mapping[m]].get(o_n) for m, bond in closures[s_n]):
                        stack.append((o_n, lp))

    @cached_property
    def _compiled_query(self):
        atoms = self._atoms
        charges = self._charges
        radicals = self._radicals
        bonds = self._bonds

        start, atom = min(atoms.items())  # todo: optimize order
        stack = [(n, start, atoms[n], charges[n], radicals[n], bond) for n, bond in bonds[start].items()]
        seen = {start}
        order = [(start, atom, charges[start], radicals[start])]
        closures = defaultdict(list)
        while stack:
            front, back, *_ = atom = stack.pop()
            if front not in seen:
                order.append(atom)
                for n, bond in bonds[front].items():  # todo: optimize order
                    if n != back:
                        if n not in seen:
                            stack.append((n, front, atoms[n], charges[n], radicals[n], bond))
                        else:
                            closures[front].append((n, bond))
                seen.add(front)
        return order, closures


__all__ = ['Isomorphism']
