# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
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
from cached_property import cached_property
from functools import lru_cache
from networkx import Graph, relabel_nodes
from ..algorithms import Components, Isomorphism, Morgan, SSSR, SubGraphs
from ..periodictable import elements_list


class BaseContainer(Graph, Components, Isomorphism, Morgan, SSSR, SubGraphs, ABC):
    def __dir__(self):
        return [] or super().__dir__()

    def __getstate__(self):
        return {'graph': self.graph, '_node': self._node, '_adj': self._adj}

    def atom(self, n):
        return self._node[n]

    def bond(self, n, m):
        return self._adj[n][m]

    @property
    def meta(self):
        return self.graph

    @cached_property
    def atoms_numbers(self):
        return list(self._node)

    @cached_property
    def atoms_count(self):
        return self.order()

    @cached_property
    def bonds_count(self):
        return self.size()

    def add_atom(self, atom, _map=None):
        """
        new atom addition
        """
        if _map is None:
            _map = max(self, default=0) + 1
        elif _map in self._node:
            raise KeyError('atom with same number exists')

        attr_dict = self.node_attr_dict_factory()
        if isinstance(atom, str):
            attr_dict.element = atom
        elif isinstance(atom, int):
            attr_dict.element = elements_list[atom - 1]
        else:
            attr_dict.update(atom)

        self._adj[_map] = self.adjlist_inner_dict_factory()
        self._node[_map] = attr_dict
        self.flush_cache()
        return _map

    def add_bond(self, atom1, atom2, bond):
        """
        implementation of bond addition
        """
        if atom1 == atom2:
            raise KeyError('atom loops impossible')
        if atom1 not in self._node or atom2 not in self._node:
            raise KeyError('atoms not found')
        if atom1 in self._adj[atom2]:
            raise KeyError('atoms already bonded')

        attr_dict = self.edge_attr_dict_factory()
        if isinstance(bond, int):
            attr_dict.order = bond
        else:
            attr_dict.update(bond)

        self._adj[atom1][atom2] = self._adj[atom2][atom1] = attr_dict
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

    @lru_cache(1000)
    def environment(self, atom):
        """
        pairs of (bond, atom) connected to atom

        :param atom: number
        :return: list
        """
        return [(bond, self._node[n]) for n, bond in self._adj[atom].items()]

    def remap(self, mapping, copy=False):
        return relabel_nodes(self, mapping, copy)

    def flush_cache(self):
        self.environment.cache_clear()
        del self.__dict__['atoms_numbers']
        del self.__dict__['atoms_count']
        del self.__dict__['bonds_count']
        super().flush_cache()

    def _bonds(self):
        seen = set()
        for n, m_bond in self._adj.items():
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
