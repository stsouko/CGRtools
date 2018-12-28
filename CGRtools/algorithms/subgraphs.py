# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
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


class SubGraphs:
    def __and__(self, other):
        """
        substructure of graph
        """
        return self.substructure(other)

    def __sub__(self, other):
        """
        other nodes excluded substructure of graph
        :return graph or None
        """
        n = self._node.keys() - set(other)
        if n:
            return self.substructure(n)

    def substructure(self, atoms, meta=False):
        """
        create substructure containing atoms from nbunch list

        :param atoms: list of atoms numbers of substructure
        :param meta: if True metadata will be copied to substructure
        :return: container with substructure
        """
        s = self.subgraph(atoms).copy()
        if not meta:
            s.graph.clear()
        return s

    def augmented_substructure(self, atoms, dante=False, deep=1, meta=False):
        """
        get subgraph with atoms and their neighbors

        :param atoms: list of core atoms in graph
        :param dante: if True return list of graphs containing atoms, atoms + first circle, atoms + 1st + 2nd,
        etc up to deep or while new nodes available.
        :param deep: number of bonds between atoms and neighbors.
        :param meta: copy metadata to each substructure
        """
        nodes = [set(atoms)]
        for i in range(deep):
            n = {y for x in nodes[-1] for y in self._adj[x]} | nodes[-1]
            if n in nodes:
                break
            nodes.append(n)

        return [self.substructure(a, meta) for a in nodes] if dante else self.substructure(nodes[-1], meta)


__all__ = ['SubGraphs']
