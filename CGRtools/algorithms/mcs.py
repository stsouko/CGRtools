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
from abc import abstractmethod
from collections import defaultdict
from itertools import product, combinations, islice, permutations
from typing import Dict, Set, Optional, Iterator, Tuple, Type


class MCS:
    __slots__ = ()

    @abstractmethod
    def get_mcs_mapping(self, other, *, automorphism_filter: bool = True) -> Dict[int, int]:
        product_graph = self.__get_product(other)
        cliques = self.__clique(product_graph)

        mapping = [x for x in islice(cliques, 42)]
        max42 = max(len(x) for x in mapping)

        yield from (dict(x) for x in mapping if len(x) == max42)

        for x in cliques:
            if len(x) > max42:
                max42 = len(x)
                yield dict(x)
            elif len(x) == max42:
                yield dict(x)

    @staticmethod
    def __clique(graph: Dict[Tuple[int, int], Set[Tuple[int, int]]]) -> Iterator[Set[Tuple[int, int]]]:
        """
        clique search

        adopted from networkx algorithms.clique.find_cliques
        """
        subgraph = {x for x, y in graph.items() if y}  # skip isolated nodes
        if not subgraph:
            return  # empty or fully disconnected
        elif len(subgraph) == 2:  # dimer
            yield set(subgraph)
            return

        stack = []
        clique_atoms = [None]
        candidates = subgraph.copy()
        roots = candidates - graph[max(subgraph, key=lambda x: len(graph[x]))]

        while True:
            if roots:
                root = roots.pop()
                candidates.remove(root)
                clique_atoms[-1] = root
                neighbors = graph[root]
                neighbors_subgraph = subgraph & neighbors
                if not neighbors_subgraph:
                    yield set(clique_atoms)
                else:
                    neighbors_candidates = candidates & neighbors
                    if neighbors_candidates:
                        stack.append((subgraph, candidates, roots))
                        clique_atoms.append(None)
                        subgraph = neighbors_subgraph
                        candidates = neighbors_candidates
                        roots = candidates - graph[max(subgraph, key=lambda x: len(candidates & graph[x]))]
            elif not stack:
                return
            else:
                clique_atoms.pop()
                subgraph, candidates, roots = stack.pop()

    def __get_product(self, other) -> Dict[Tuple[int, int], Set[Tuple[int, int]]]:
        bonds = self._bonds
        o_bonds = other._bonds

        s_equal = defaultdict(list)  # equal self atoms
        for n, atom in self._atoms.items():
            s_equal[atom].append(n)
        p_equal = defaultdict(list)  # equal other atoms
        for n, atom in other._atoms.items():
            p_equal[atom].append(n)

        half_product = {}
        full_product = {}
        equal_atoms = {}
        for atom, ns in s_equal.items():
            ms = p_equal[atom]
            if ms:
                for nm in product(ns, ms):
                    full_product[nm] = set()
                    half_product[nm] = defaultdict(list)
                for n in ns:
                    equal_atoms[n] = ms  # memory save

        seen = set()
        for n, o_ns in equal_atoms.items():
            seen.add(n)
            for m, b in bonds[n].items():
                if m in equal_atoms and m not in seen:
                    o_ms = equal_atoms[m]
                    for o_n in o_ns:
                        ms = []
                        node1 = (n, o_n)
                        fms = full_product[node1]
                        for o_m, o_b in o_bonds[o_n].items():
                            if o_m in o_ms and b == o_b:
                                node2 = (m, o_m)
                                ms.append(node2)
                                full_product[node2].add(node1)
                                fms.add(node2)
                        if ms:
                            half_product[node1][m].extend(ms)

        seen = set()
        for n, moms in half_product.items():
            if moms:
                n1 = n[0]
                seen.add(n1)
                for m1, mom in moms.items():
                    if m1 not in seen:
                        for m in mom:
                            half_product[m][n1].append(n)
                for (m11, ms1), (m21, ms2) in combinations(moms.items(), 2):
                    if m11 not in seen:
                        if m21 in seen:
                            for m1, m2 in product(ms1, ms2):
                                if m2 not in full_product[m1] and m1[1] != m2[1]:
                                    half_product[m1][m21].append(m2)
                                    full_product[m1].add(m2)
                                    full_product[m2].add(m1)
                        else:
                            for m1, m2 in product(ms1, ms2):
                                if m1[1] != m2[1]:
                                    half_product[m1][m21].append(m2)
                                    full_product[m1].add(m2)
                                    full_product[m2].add(m1)
        return full_product


__all__ = ['MCS']
