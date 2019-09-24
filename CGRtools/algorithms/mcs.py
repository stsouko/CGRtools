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
from itertools import product, combinations
from typing import Dict, List, Set, Hashable, Iterator, Tuple


class MCS:
    __slots__ = ()

    @abstractmethod
    def get_mcs_mapping(self, other) -> Dict[int, int]:
        bonds = self._bonds
        o_bonds = self._bonds
        product_graph = self.__get_product(other)
        mapping = {}
        for clique in self.__clique(product_graph):
            pass  # todo: bond filter
        return mapping

    @staticmethod
    def __clique(graph: Dict[Tuple[int, int], Set[Tuple[int, int]]]) -> Iterator[Dict[int, int]]:
        """
        clique search

        adopted from networkx algorithms.clique.find_cliques
        """
        subgraph = {x for x, y in graph.items() if y}  # skip isolated nodes
        if not subgraph:
            return  # empty or fully disconnected
        elif len(subgraph) == 2:  # dimer
            yield dict(subgraph)
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
                    yield dict(clique_atoms)
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

        s_equal = defaultdict(set)  # equal self atoms. from n - 1 to (n^2 - n)/2 complexity
        for n, atom in self._atoms.items():
            s_equal[atom].add(n)

        product_graph = {}
        equal_atoms = defaultdict(set)
        for n, atom in other._atoms.items():
            for m in s_equal.get(atom, ()):
                product_graph[(m, n)] = set()
                equal_atoms[m].add(n)

        for (n, ns_other), (m, ms_other) in combinations(equal_atoms.items(), 2):
            bond = bonds[n].get(m)
            if not bond:
                for o_n, o_m in product(ns_other, ms_other):
                    if o_n != o_m:
                        node1 = (n, o_n)
                        node2 = (m, o_m)
                        product_graph[node1].add(node2)
                        product_graph[node2].add(node1)
            else:
                for o_n, o_m in product(ns_other, ms_other):
                    if o_n != o_m:
                        o_bond = o_bonds[o_n].get(o_m)
                        if not o_bond or bond == o_bond:
                            node1 = (n, o_n)
                            node2 = (m, o_m)
                            product_graph[node1].add(node2)
                            product_graph[node2].add(node1)
        return product_graph


__all__ = ['MCS']
