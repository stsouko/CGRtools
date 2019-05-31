# -*- coding: utf-8 -*-
#
#  Copyright 2019 Ramil Nugmanov <stsouko@live.ru>
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


class Components:
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
        n = self._atoms.keys() - set(other)
        if n:
            return self.substructure(n)

    def substructure(self, atoms, *, meta=False, as_view=True):
        """
        create substructure containing atoms from nbunch list

        :param atoms: list of atoms numbers of substructure
        :param meta: if True metadata will be copied to substructure
        :param as_view: If True, the returned graph-view provides a read-only view
            of the original structure scaffold without actually copying any data.
        """
        s = self.subgraph(atoms)
        if as_view:
            s.add_atom = s.add_bond = s.delete_atom = s.delete_bond = frozen  # more informative exception
            return s
        s = s.copy()
        if not meta:
            s.graph.clear()
        return s

    def augmented_substructure(self, atoms, dante=False, deep=1, meta=False, as_view=True):
        """
        create substructure containing atoms and their neighbors

        :param atoms: list of core atoms in graph
        :param dante: if True return list of graphs containing atoms, atoms + first circle, atoms + 1st + 2nd,
            etc up to deep or while new nodes available
        :param deep: number of bonds between atoms and neighbors
        :param meta: copy metadata to each substructure
        :param as_view: If True, the returned graph-view provides a read-only view
            of the original graph without actually copying any data
        """
        nodes = [set(atoms)]
        for i in range(deep):
            n = {y for x in nodes[-1] for y in self._bonds[x]} | nodes[-1]
            if n in nodes:
                break
            nodes.append(n)
        if dante:
            return [self.substructure(a, meta, as_view) for a in nodes]
        else:
            return self.substructure(nodes[-1], meta, as_view)

    def split(self, meta=False):
        """
        split disconnected structure to connected substructures

        :param meta: copy metadata to each substructure
        :return: list of substructures
        """
        return [self.substructure(c, meta, False) for c in self.connected_components]

    @property
    def connected_components(self):
        if not self._atoms:
            return []
        atoms = set(self._atoms)
        components = []
        while atoms:
            start = atoms.pop()
            component = list(self.__component(start))
            components.append(component)
            atoms.difference_update(component)
        return components

    @property
    def connected_components_count(self):
        return len(self.connected_components)

    def __component(self, start):
        bonds = self._bonds
        seen = {start}
        queue = [start]
        while queue:
            start = queue.pop(0)
            yield start
            for i in bonds[start].keys() - seen:
                queue.append(i)
                seen.add(i)
