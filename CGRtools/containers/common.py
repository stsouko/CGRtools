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
from abc import ABC, abstractmethod
from networkx import Graph, relabel_nodes
from networkx.readwrite.json_graph import node_link_graph, node_link_data
from typing import Callable
from ..algorithms import hash_cgr_string, get_morgan


class BaseContainer(Graph, ABC):
    def __dir__(self):
        if self.__visible is None:
            self.__visible = [self.pickle.__name__, self.unpickle.__name__, self.copy.__name__, self.remap.__name__,
                              self.flush_cache.__name__, self.substructure.__name__,  self.get_morgan.__name__,
                              self.get_signature.__name__, self.get_signature_hash.__name__,
                              self.get_environment.__name__, self.mark.__name__, self.atom.__name__, self.bond.__name__,
                              self.stereo.__name__, self.add_atom.__name__, self.add_bond.__name__,
                              self.add_stereo.__name__, self.get_stereo.__name__,
                              self.delete_atom.__name__, self.delete_bond.__name__,
                              'meta', 'bonds_count', 'atoms_count', 'atom_numbers']
        return self.__visible

    def __getstate__(self):
        return {'graph': self.graph, '_node': self._node, '_adj': self._adj}

    def mark(self, n):
        return self._node[n].mark

    def atom(self, n):
        return self._node[n]

    def stereo(self, n, m=None):
        return (self._node[n] if m is None else self[n][m]).stereo

    def bond(self, n, m):
        return self[n][m]

    @property
    def atom_numbers(self):
        return list(self._node)

    @property
    def atoms_count(self):
        return self.order()

    @property
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
        if isinstance(atom, dict):
            self.add_node(_map, **atom)
        else:
            self.add_node(_map, atom=atom)
        self.flush_cache()
        return _map

    @abstractmethod
    def add_bond(self, atom1, atom2, bond):
        """
        implementation of bond addition
        """
        if atom1 not in self._node or atom2 not in self._node:
            raise KeyError('atoms not found')
        if self.has_edge(atom1, atom2):
            raise KeyError('atoms already bonded')
        if isinstance(bond, dict):
            self.add_edge(atom1, atom2, **bond)
        else:
            self.add_edge(atom1, atom2, bond=bond)
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

    @abstractmethod
    def add_stereo(self, *args, **kwargs):
        """
        implementation of stereo addition
        """
        pass

    def get_stereo(self, atom1, atom2):
        if self.__stereo_cache is None:
            self.__stereo_cache = self._prepare_stereo()
        return self.__stereo_cache.get((atom1, atom2), None)

    @abstractmethod
    def _prepare_stereo(self):
        """
        :return: dict of stereo labels on bonds
        """
        pass

    @property
    def meta(self):
        return self.graph

    def pickle(self):
        """ return json serializable Container
        """
        data = node_link_data(self, attrs=self.__attrs)
        data['class'] = type(self).__name__
        del data['multigraph'], data['directed']
        return data

    @classmethod
    def unpickle(cls, data):
        """convert json serializable Container into Molecule or CGR or Query object instance"""
        if 's_charge' in data['nodes'][0]:  # 2.8 compatibility
            is_cgr = not data['s_only']

            if data.get('is_query'):
                _class = 'QueryContainer'
            elif is_cgr:
                _class = 'CGRContainer'
            else:
                _class = 'MoleculeContainer'

            nodes = []
            for x in data['nodes']:
                n = {'atom': x['atom'], 'element': x['element'], 'charge': x['s_charge'], 'mark': x['mark'],
                     'x': x['s_x'], 'y': x['s_y'], 'z': x['s_z']}
                if 'isotope' in x:
                    n['isotope'] = x['isotope']
                if 's_hyb' in x:
                    n['hybridization'] = x['s_hyb']
                if 's_neighbors' in x:
                    n['neighbors'] = x['s_neighbors']
                if 's_stereo' in x:
                    n['stereo'] = x['s_stereo']
                if 's_radical' in x:
                    n['multiplicity'] = x['s_radical']
                if is_cgr:
                    n.update(p_x=x['p_x'], p_y=x['p_y'], p_z=x['p_z'], p_charge=x['p_charge'])
                    if 'p_hyb' in x:
                        n['p_hybridization'] = x['p_hyb']
                    if 'p_neighbors' in x:
                        n['p_neighbors'] = x['p_neighbors']
                    if 'p_stereo' in x:
                        n['p_stereo'] = x['p_stereo']
                    if 'p_radical' in x:
                        n['p_multiplicity'] = x['p_radical']
                nodes.append(n)

            bonds = []
            for x in data['bonds']:
                b = {'atom1': x['atom1'], 'atom2': x['atom2'], 'order': x['s_bond']}
                if 's_stereo' in x:
                    b['stereo'] = x['s_stereo']
                if is_cgr:
                    if 'p_bond' in x:
                        b['p_order'] = x['p_bond']
                    if 'p_stereo' in x:
                        b['p_stereo'] = x['p_stereo']
                bonds.append(b)

            data = {'graph': data['meta'], 'nodes': nodes, 'bonds': bonds, 'class': _class}

        graph = node_link_graph(data, multigraph=False, attrs=cls.__attrs)
        return _search_subclass(data['class'])(graph)

    def substructure(self, nbunch, meta=False):
        """
        create substructure containing atoms from nbunch list

        :param nbunch: list of atoms numbers of substructure
        :param meta: if True metadata will be copied to substructure
        :return: container with substructure
        """
        s = self.subgraph(nbunch).copy()
        if not meta:
            s.graph.clean()
        return s

    def get_environment(self, atoms, dante=False, deep=0):
        """
        get subgraph with atoms and their neighbors

        :param atoms: list of core atoms in graph
        :param dante: if True return list of graphs containing atoms, atoms + first circle, atoms + 1st + 2nd,
        etc up to deep or while new nodes available.
        :param deep: number of bonds between atoms and neighbors.
        """
        nodes = [set(atoms)]
        for i in range(deep):
            n = {y for x in nodes[-1] for y in self._adj[x]} | nodes[-1]
            if n in nodes:
                break
            nodes.append(n)

        return [self.substructure(a) for a in nodes] if dante else self.substructure(nodes[-1])

    def remap(self, mapping, copy=False):
        return relabel_nodes(self, mapping, copy=copy)

    def get_signature_hash(self, *args, **kwargs):
        return hash_cgr_string(self.get_signature(*args, **kwargs))

    def get_signature(self, isotope=False, stereo=False, hybridization=False, neighbors=False, element=True,
                      flush_cache=False,  weights=None):
        """
        return string representation of structure

        :param weights: dict of atoms in keys and orders in values
        :param isotope: set isotope marks
        :param stereo: set stereo marks
        :param hybridization: set hybridization mark of atom
        :param neighbors: set neighbors count mark of atom
        :param element: set elements marks
        :param flush_cache: recalculate signature if True
        """
        if flush_cache or self.__signatures is None:
            self.__signatures = {}

        k = (isotope, element, stereo, hybridization, neighbors)
        out = self.__signatures.get(k)
        if not out:
            sg = self._signature_generator(element, isotope, stereo, hybridization, neighbors)
            if not weights:
                weights = self.get_morgan(isotope, element, stereo, hybridization, neighbors, flush_cache=flush_cache)
            self.__signatures[k] = out = sg(self, weights)
        return out

    def get_morgan(self, isotope=False, element=True, stereo=False, hybridization=False, neighbors=False, labels=None,
                   flush_cache=False):
        if flush_cache or self.__weights is None:
            self.__weights = {}
        k = (isotope, element, stereo, hybridization, neighbors, labels)
        return self.__weights.get(k) or self.__weights.setdefault(k, get_morgan(self, isotope, element, stereo,
                                                                                hybridization, neighbors, labels))

    @abstractmethod
    def _signature_generator(self, *args, **kwargs) -> Callable:
        """
        container specific signature generation class
        """
        pass

    def flush_cache(self):
        self.__weights = self.__signatures = self.__pickle = self.__hash = self.__stereo_cache = None

    def __str__(self):
        return self.get_signature(isotope=True, stereo=True, hybridization=True, neighbors=True)

    def __repr__(self):
        if self.__pickle is None:
            self.__pickle = '%s.unpickle(%s)' % (self.__class__.__name__, self.pickle())
        return self.__pickle

    def __hash__(self):
        if self.__hash is None:
            self.__hash = int.from_bytes(self.get_signature_hash(isotope=True, stereo=True, hybridization=True,
                                                                 neighbors=True), 'big')
        return self.__hash

    def __eq__(self, other):
        return str(self) == str(other)

    __visible = __weights = __signatures = __pickle = __stereo_cache = __hash = None
    __attrs = dict(source='atom1', target='atom2', name='atom', link='bonds')


def _search_subclass(name, start=BaseContainer):
    for x in start.__subclasses__():
        if x.__name__ == name:
            return x
        d = _search_subclass(name, x)
        if d is not None:
            return d
