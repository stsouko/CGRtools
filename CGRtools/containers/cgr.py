# -*- coding: utf-8 -*-
#
#  Copyright 2017-2019 Ramil Nugmanov <stsouko@live.ru>
#  Copyright 2018 Ravil Mukhametgaleev <sonic-mc@mail.ru>
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
from networkx.algorithms.isomorphism import GraphMatcher
from networkx.classes.function import frozen
from typing import List
from .common import BaseContainer
from ..algorithms import Morgan, SmilesCGR, CGRCompose, DepictCGR
from ..attributes import DynAtom, DynBond
from ..cache import cached_property


class CGRContainer(CGRCompose, Morgan, SmilesCGR, BaseContainer, DepictCGR):
    """
    storage for CGRs. has similar to molecules behavior
    """
    node_attr_dict_factory = DynAtom
    edge_attr_dict_factory = DynBond

    @cached_property
    def centers_list(self):
        """ get a list of lists of atoms of reaction centers
        """
        center = set()
        adj = defaultdict(set)
        for n, atom in self.atoms():
            if atom._reactant != atom._product:
                center.add(n)

        for n, m, bond in self.bonds():
            if bond._reactant != bond._product:
                adj[n].add(m)
                adj[m].add(n)
                center.add(n)
                center.add(m)

        out = []
        while center:
            n = center.pop()
            if n in adj:
                c = set(self.__plain_bfs(adj, n))
                out.append(list(c))
                center.difference_update(c)
            else:
                out.append([n])

        return out

    @cached_property
    def center_atoms(self):
        """ get list of atoms of reaction center (atoms with dynamic: bonds, charges, radicals).
        """
        nodes = set()
        for n, atom in self.atoms():
            if atom._reactant != atom._product:
                nodes.add(n)

        for n, m, bond in self.bonds():
            if bond._reactant != bond._product:
                nodes.add(n)
                nodes.add(m)

        return list(nodes)

    @cached_property
    def center_bonds(self):
        """ get list of bonds of reaction center (bonds with dynamic orders).
        """
        return [(n, m) for n, m, bond in self.bonds() if bond._reactant != bond._product]

    def reset_query_marks(self):
        """
        set or reset hyb and neighbors marks to atoms.
        """
        for i, atom in self.atoms():
            neighbors = 0
            hybridization = 1
            p_neighbors = 0
            p_hybridization = 1
            # hyb 1- sp3; 2- sp2; 3- sp1; 4- aromatic
            for j, bond in self._adj[i].items():
                isnth = self._node[j].element != 'H'

                order = bond.order
                if order:
                    if isnth:
                        neighbors += 1
                    if hybridization not in (3, 4):
                        if order == 4:
                            hybridization = 4
                        elif order == 3:
                            hybridization = 3
                        elif order == 2:
                            if hybridization == 2:
                                hybridization = 3
                            else:
                                hybridization = 2
                order = bond.p_order
                if order:
                    if isnth:
                        p_neighbors += 1
                    if p_hybridization not in (3, 4):
                        if order == 4:
                            p_hybridization = 4
                        elif order == 3:
                            p_hybridization = 3
                        elif order == 2:
                            if p_hybridization == 2:
                                p_hybridization = 3
                            else:
                                p_hybridization = 2

            atom._reactant._neighbors = neighbors
            atom._reactant._hybridization = hybridization
            atom._product._neighbors = p_neighbors
            atom._product._hybridization = p_hybridization
            atom.__dict__.clear()  # flush cache
        self.flush_cache()

    def substructure(self, atoms, meta=False, as_view=True):
        """
        create substructure containing atoms from nbunch list

        :param atoms: list of atoms numbers of substructure
        :param meta: if True metadata will be copied to substructure
        :param as_view: If True, the returned graph-view provides a read-only view
            of the original structure scaffold without actually copying any data
        """
        s = super().substructure(atoms, meta, as_view)
        if as_view:
            s.reset_query_marks = frozen
        return s

    @cached_property
    def aromatic_rings(self) -> List[List[int]]:
        """
        existed or formed aromatic rings atoms numbers
        """
        adj = self._bonds
        return [ring for ring in self.sssr if len(ring) in (5, 6, 7) and (
                adj[ring[0]][ring[-1]].order == 4 and all(adj[n][m].order == 4 for n, m in zip(ring, ring[1:])) or
                adj[ring[0]][ring[-1]].p_order == 4 and all(adj[n][m].p_order == 4 for n, m in zip(ring, ring[1:])))]

    def _matcher(self, other):
        """
        CGRContainer < CGRContainer
        """
        if isinstance(other, CGRContainer):
            return GraphMatcher(other, self, lambda x, y: x == y, lambda x, y: x == y)
        raise TypeError('only cgr-cgr possible')

    @staticmethod
    def __plain_bfs(adj, source):
        """modified NX fast BFS node generator"""
        seen = set()
        nextlevel = {source}
        while nextlevel:
            thislevel = nextlevel
            nextlevel = set()
            for v in thislevel:
                if v not in seen:
                    yield v
                    seen.add(v)
                    nextlevel.update(adj[v])


__all__ = ['CGRContainer']
