# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
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
from abc import ABC
from networkx import connected_components, Graph, relabel_nodes
from networkx.classes.function import frozen
from ..algorithms import Isomorphism, SSSR, Union
from ..cache import cached_property, cached_args_method
from ..periodictable import elements_list


class BaseContainer(Graph, Isomorphism, SSSR, Union, ABC):
    __slots__ = ('graph', '_node', '_adj', '_bonds', '_atoms', '_meta')

    def __init__(self, *args, **kwargs):
        """
        Empty data object initialization or conversion from another object type
        """
        super().__init__(*args, **kwargs)
        self._bonds = self._adj  # migration ad-hoc
        self._atoms = self._node
        self._meta = self.graph

    def __dir__(self):
        return [] or super().__dir__()

    def __getstate__(self):
        return {'node': self._atoms, 'adj': self._bonds, 'meta': self._meta}  # migration ad-hoc

    def __setstate__(self, state):
        if 'graph' in state:  # 3.0.10 compatibility.
            state['meta'] = state['graph']
        if 'node' in state:  # 3.1 compatibility.
            state['atoms'] = state['node']
            state['bonds'] = state['adj']

        self._meta = state['meta']
        self._atoms = state['atoms']
        self._bonds = state['bonds']

        self._adj = self._bonds  # migration ad-hoc
        self._node = self._atoms
        self.graph = self._meta

    def atom(self, n):
        return self._atoms[n]

    def bond(self, n, m):
        return self._bonds[n][m]

    @property
    def meta(self):
        return self._meta

    @cached_property
    def atoms_numbers(self):
        return list(self._atoms)

    @cached_property
    def atoms_count(self):
        return len(self)

    @cached_property
    def bonds_count(self):
        return self.size()

    def add_atom(self, atom, _map=None):
        """
        new atom addition
        """
        if _map is None:
            _map = max(self, default=0) + 1
        elif _map in self._atoms:
            raise KeyError('atom with same number exists')

        attr_dict = self.node_attr_dict_factory()
        if isinstance(atom, str):
            attr_dict.element = atom
        elif isinstance(atom, int):
            attr_dict.element = elements_list[atom - 1]
        else:
            attr_dict.update(atom)

        self._bonds[_map] = self.adjlist_inner_dict_factory()
        self._atoms[_map] = attr_dict
        self.flush_cache()
        return _map

    def add_bond(self, atom1, atom2, bond):
        """
        implementation of bond addition
        """
        if atom1 == atom2:
            raise KeyError('atom loops impossible')
        if atom1 not in self._atoms or atom2 not in self._atoms:
            raise KeyError('atoms not found')
        if atom1 in self._bonds[atom2]:
            raise KeyError('atoms already bonded')

        attr_dict = self.edge_attr_dict_factory()
        if isinstance(bond, int):
            attr_dict.order = bond
        else:
            attr_dict.update(bond)

        self._bonds[atom1][atom2] = self._bonds[atom2][atom1] = attr_dict
        self.flush_cache()

    def delete_atom(self, n):
        """
        implementation of atom removing
        """
        self.remove_node(n)
        self.flush_cache()

    def delete_bond(self, n, m):
        """
        implementation of bond removing
        """
        self.remove_edge(n, m)
        self.flush_cache()

    @cached_args_method
    def environment(self, atom):
        """
        pairs of (bond, atom) connected to atom

        :param atom: number
        :return: list
        """
        return tuple((bond, self._atoms[n]) for n, bond in self._bonds[atom].items())

    def substructure(self, atoms, meta=False, as_view=True):
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
            s._meta = s.graph
            s._bonds = s._adj
            s._atoms = s._node
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

    @cached_property
    def connected_components(self):
        return [list(x) for x in connected_components(self)]

    def split(self, meta=False):
        """
        split disconnected structure to connected substructures

        :param meta: copy metadata to each substructure
        :return: list of substructures
        """
        return [self.substructure(c, meta, False) for c in connected_components(self)]

    def remap(self, mapping, copy=False):
        return relabel_nodes(self, mapping, copy)

    def flush_cache(self):
        self.__dict__.clear()

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

    def atoms(self):
        """
        iterate over all atoms
        """
        return iter(self._atoms.items())

    def bonds(self):
        """
        iterate other all bonds
        """
        seen = set()
        for n, m_bond in self._bonds.items():
            seen.add(n)
            for m, bond in m_bond.items():
                if m not in seen:
                    yield n, m, bond

    @staticmethod
    def _get_subclass(name):
        """
        need for cyclic import solving
        """
        return next(x for x in BaseContainer.__subclasses__() if x.__name__ == name)


__all__ = ['BaseContainer']
