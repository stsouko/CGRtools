# -*- coding: utf-8 -*-
#
#  Copyright 2018-2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from abc import abstractmethod
from CachedMethods import cached_property
from collections import defaultdict
from itertools import permutations
from typing import Dict, Iterator, Any
from .._functions import lazy_product


frequency = {1: 10,  # H
             6: 9,  # C
             8: 8,  # O
             7: 7,  # N
             15: 6, 16: 6,  # P, S
             9: 5,  # F
             17: 4, 35: 4,  # Cl, Br
             53: 3,  # I
             5: 2, 14: 2,  # B, Si
             11: 1, 12: 1, 19: 1, 20: 1}  # Na, Mag,  K, Ca


def atom_frequency(x):
    return frequency.get(x.atomic_number, 0)


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
        Test self is substructure of other
        """
        try:
            next(self.get_mapping(other))
        except StopIteration:
            return False
        return True

    def is_equal(self, other) -> bool:
        """
        Test self is same structure as other
        """
        if len(self) != len(other):
            return False
        try:
            next(self.get_mapping(other))
        except StopIteration:
            return False
        return True

    @abstractmethod
    def get_mapping(self, other, *, automorphism_filter: bool = True,
                    optimize: bool = True, fallback: bool = False) -> Iterator[Dict[int, int]]:
        """
        Get self to other substructure mapping generator.

        :param automorphism_filter: Skip matches to same atoms.
        :param optimize: Morgan weights based automorphism preventing.
        :param fallback: Try without optimization then nothing matched.
        """
        if optimize:
            g = self.__components_mapping(other, other.atoms_order, automorphism_filter)
            m = next(g, None)
            if m is not None:
                yield m
                yield from g
                return
            elif not fallback:
                return
        yield from self.__components_mapping(other, {n: i for i, n in enumerate(other)}, automorphism_filter)

    def __components_mapping(self, other, o_order, automorphism_filter):
        components, closures = self._compiled_query
        o_atoms = other._atoms
        o_bonds = other._bonds

        seen = set()
        if len(components) == 1:
            for candidate in other.connected_components:
                for mapping in self._get_mapping(components[0], closures, o_atoms, o_bonds, set(candidate), o_order):
                    if automorphism_filter:
                        atoms = frozenset(mapping.values())
                        if atoms in seen:
                            continue
                        seen.add(atoms)
                    yield mapping
        else:
            for candidates in permutations((set(x) for x in other.connected_components), len(components)):
                mappers = [self._get_mapping(order, closures, o_atoms, o_bonds, component, o_order)
                           for order, component in zip(components, candidates)]
                for match in lazy_product(*mappers):
                    mapping = match[0].copy()
                    for m in match[1:]:
                        mapping.update(m)
                    if automorphism_filter:
                        atoms = frozenset(mapping.values())
                        if atoms in seen:
                            continue
                        seen.add(atoms)
                    yield mapping

    @staticmethod
    def _get_mapping(linear_query, query_closures, o_atoms, o_bonds, scope, groups):
        size = len(linear_query) - 1
        order_depth = {v[0]: k for k, v in enumerate(linear_query)}
        equal_cache = defaultdict(dict)

        stack = []
        path = []
        mapping = {}
        reversed_mapping = {}

        s_n, s_atom = linear_query[0]
        eqs = equal_cache[s_n]
        for n, o_atom in o_atoms.items():
            if n in scope:
                if s_atom == o_atom:
                    eqs[n] = True
                    stack.append((n, 0))
                else:
                    eqs[n] = False

        while stack:
            n, depth = stack.pop()
            current = linear_query[depth][0]
            if depth == size:
                yield {current: n, **mapping}
            else:
                if len(path) != depth:
                    for x in path[depth:]:
                        del mapping[reversed_mapping.pop(x)]
                    path = path[:depth]

                path.append(n)
                mapping[current] = n
                reversed_mapping[n] = current

                depth += 1
                s_n, back, s_atom, s_bond = linear_query[depth]
                if back != current:
                    n = path[order_depth[back]]

                eqs = equal_cache[s_n]
                uniq = set()
                for o_n, o_bond in o_bonds[n].items():
                    if o_n in scope and o_n not in reversed_mapping and s_bond == o_bond and groups[o_n] not in uniq:
                        uniq.add(groups[o_n])
                        if o_n in eqs:
                            if eqs[o_n]:
                                if all(bond == o_bonds[mapping[m]].get(o_n) for m, bond in query_closures[s_n]):
                                    stack.append((o_n, depth))
                        elif s_atom == o_atoms[o_n]:
                            eqs[o_n] = True
                            if all(bond == o_bonds[mapping[m]].get(o_n) for m, bond in query_closures[s_n]):
                                stack.append((o_n, depth))
                        else:
                            eqs[o_n] = False

    @cached_property
    def _compiled_query(self):
        return self.__compile_query(self._atoms, self._bonds, {n: atom_frequency(a) for n, a in self._atoms.items()})

    @staticmethod
    def __compile_query(atoms, bonds, atoms_frequencies):
        closures = defaultdict(list)
        components = []
        seen = set()
        while len(seen) < len(atoms):
            start = min(atoms.keys() - seen, key=lambda x: atoms_frequencies[x])
            seen.add(start)
            stack = [(n, start, atoms[n], bond) for n, bond in sorted(bonds[start].items(), reverse=True,
                                                                      key=lambda x: atoms_frequencies[x[0]])]
            order = [(start, atoms[start])]
            components.append(order)

            while stack:
                front, back, *_ = atom = stack.pop()
                if front not in seen:
                    order.append(atom)
                    for n, bond in sorted(bonds[front].items(), reverse=True, key=lambda x: atoms_frequencies[x[0]]):
                        if n != back:
                            if n not in seen:
                                stack.append((n, front, atoms[n], bond))
                            else:
                                closures[front].append((n, bond))
                    seen.add(front)
        return components, closures

    def is_automorphic(self):
        """
        Test for automorphism symmetry of graph.
        """
        try:
            next(self.get_automorphism_mapping())
        except StopIteration:
            return False
        return True

    def get_automorphism_mapping(self) -> Iterator[Dict[int, int]]:
        """
        Iterator of all possible automorphism mappings.
        """
        return self._get_automorphism_mapping(self.atoms_order, self._bonds,
                                              {n: atom_frequency(a) for n, a in self._atoms.items()})

    @classmethod
    def _get_automorphism_mapping(cls, atoms: Dict[int, int], bonds: Dict[int, Dict[int, Any]],
                                  atoms_frequencies: Dict[int, int]) -> Iterator[Dict[int, int]]:
        if len(atoms) == len(set(atoms.values())):
            return  # all atoms unique

        components, closures = cls.__compile_query(atoms, bonds, atoms_frequencies)
        groups = {x: n for n, x in enumerate(atoms)}
        mappers = [cls._get_mapping(order, closures, atoms, bonds, {x for x, *_ in order}, groups)
                   for order in components]
        if len(mappers) == 1:
            for mapping in mappers[0]:
                if any(k != v for k, v in mapping.items()):
                    yield mapping
        for match in lazy_product(*mappers):
            mapping = match[0]
            for m in match[1:]:
                mapping.update(m)
            if any(k != v for k, v in mapping.items()):
                yield mapping


__all__ = ['Isomorphism']
