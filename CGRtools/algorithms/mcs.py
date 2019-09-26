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
        yield from (dict(c) for c in cliques)
        #mapping = [(c, calc_bonds(c)) for c in islice(cliques, 42)]
        #max_bonds = max(x for _, x in mapping)
        #max_atoms = max(len(x) for x, y in mapping if y == max_bonds)
        #yield from (dict(x) for x, y in mapping if y == max_bonds and len(x) == max_atoms)

        #br = 0
        #for clique in cliques:
        #    if len(clique) < max_atoms:
        #        br += 1
        #        if br > 10:
        #            break
        #        continue
        #    elif br:
        #        br = 0
        #    bond_sum = calc_bonds(clique)
        #    if bond_sum > max_bonds:
        #        max_bonds = bond_sum
        #        max_atoms = len(clique)
        #        yield dict(clique)
        #    elif bond_sum == max_bonds:
        #        yield dict(clique)

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

        s_equal = defaultdict(set)  # equal self atoms
        for n, atom in self._atoms.items():
            s_equal[atom].add(n)
        p_equal = defaultdict(set)  # equal other atoms
        for n, atom in other._atoms.items():
            p_equal[atom].add(n)

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
                for mc in combinations(moms.values(), 2):
                    for m1, m2 in product(*mc):
                        if m1[1] != m2[1]:
                            full_product[m1].add(m2)
                            full_product[m2].add(m1)

        return full_product


__all__ = ['MCS']
