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
from collections import defaultdict
from itertools import chain
from networkx import Graph, relabel_nodes
from networkx.readwrite.json_graph import node_link_graph, node_link_data
from typing import Callable, Iterable
from warnings import warn
from ..algorithms import hash_cgr_string, get_morgan
from ..exceptions import InvalidData, InvalidAtom


class BaseContainer(Graph, ABC):
    def __init__(self, data=None, meta=None):
        """
        base structure container class

        :param data: MoleculeContainer [CGRContainer] or NX Graph object or other supported by NX for initialization
        :param meta: dictionary of metadata. like DTYPE-DATUM in RDF
        """
        super().__init__(data)
        self.__init()

        if meta is not None:
            if isinstance(meta, dict):
                self.__meta.update(meta)
            else:
                raise InvalidData('metadata can be dictionary')

    def __init(self):
        self.__atom_cache = {}
        self.__meta = {}
        self.__bond_cache = defaultdict(dict)

    def __dir__(self):
        if self.__visible is None:
            self.__visible = [self.pickle.__name__, self.unpickle.__name__, self.copy.__name__, self.remap.__name__,
                              self.flush_cache.__name__, self.substructure.__name__,  self.get_morgan.__name__,
                              self.get_signature.__name__, self.get_signature_hash.__name__, self.fix_data.__name__,
                              self.get_environment.__name__, self.mark.__name__, self.atom.__name__, self.bond.__name__,
                              self.stereo.__name__, self.add_atom.__name__, self.add_bond.__name__,
                              self.add_stereo.__name__, self.get_stereo.__name__,
                              self.delete_atom.__name__, self.delete_bond.__name__,
                              'meta', 'bonds_count', 'atoms_count', 'atom_numbers']
        return self.__visible

    def __getstate__(self):
        return {k: v for k, v in super().__getstate__().items()
                if not k.startswith('_BaseContainer__') and k != 'root_graph' or k == '_BaseContainer__meta'}

    def __setstate__(self, state):
        self.__dict__.update(state)
        self.root_graph = self
        self.__init()

    def mark(self, n):
        try:
            tmp = self.nodes[n]
        except KeyError:
            raise InvalidAtom('atom not found')
        return tmp['mark']

    def atom(self, n):
        if n not in self.__atom_cache:
            try:
                tmp = self.nodes[n]
            except KeyError:
                raise InvalidAtom('atom not found')

            self.__atom_cache[n] = self._atom_container(tmp)
        return self.__atom_cache[n]

    def stereo(self, n, m=None):
        try:
            tmp = self.nodes[n] if m is None else self[n][m]
        except KeyError:
            raise InvalidAtom('atom[s] or bond not found')

        return self._stereo_container(tmp)

    def bond(self, n, m):
        if m not in self.__bond_cache[n]:
            try:
                tmp = self[n][m]
            except KeyError:
                raise InvalidAtom('atom[s] or bond not found')

            self.__bond_cache[n][m] = self.__bond_cache[m][n] = self._bond_container(tmp)
        return self.__bond_cache[n][m]

    @property
    def atom_numbers(self):
        return list(self.nodes)

    @property
    def atoms_count(self):
        return self.order()

    @property
    def bonds_count(self):
        return self.size()

    @abstractmethod
    def add_atom(self, *args, **kwargs):
        """
        implementation of atom addition
        """
        pass

    @abstractmethod
    def add_bond(self, *args, **kwargs):
        """
        implementation of bond addition
        """
        pass

    def delete_atom(self, n):
        """
        implementation of atom removing
        """
        self.remove_node(n)

    def delete_bond(self, n, m):
        """
        implementation of bond removing
        """
        self.remove_edge(n, m)

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
        :return: dict of stere lables on bonds
        """
        pass

    @property
    def meta(self):
        return self.__meta

    @abstractmethod
    def pickle(self):
        """ return json serializable Container
        """
        node_marks = self._node_save
        edge_marks = self._edge_save

        g = Graph()
        g.add_nodes_from((n, {k: v for k, v in a.items() if k in node_marks}) for n, a in self.nodes(data=True))
        g.add_edges_from((n, m, {k: v for k, v in a.items() if k in edge_marks}) for n, m, a in self.edges(data=True))

        data = node_link_data(g, attrs=self.__attrs)
        data['meta'] = self.__meta
        return data

    @classmethod
    @abstractmethod
    def unpickle(cls, data) -> object:
        """convert json serializable Container into Molecule or CGR or Query object instance"""
        graph = node_link_graph(data, attrs=cls.__attrs)
        meta = data['meta']
        return graph, meta

    def copy(self):
        copy = super().copy()
        copy.meta.update(self.__meta)
        return copy

    def substructure(self, nbunch, meta=False):
        """
        create substructure containing atoms from nbunch list

        Notes: for prevent of data corruption in original structure, create copy of substructure. actual for templates

        :param nbunch: list of atoms numbers of substructure
        :param meta: if True metadata will be copied to substructure
        :return: Molecule or CGR container
        """
        return type(self)(self.subgraph(nbunch), self.meta if meta else None)

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
            n = set(chain.from_iterable(self.edges(nodes[i])))
            if not n or n in nodes:
                break
            nodes.append(n)

        if dante:
            centers = [self.substructure(a) for a in nodes]
        else:
            centers = self.substructure(nodes[-1])

        return centers

    def remap(self, mapping, copy=False):
        g = relabel_nodes(self, mapping, copy=copy)
        if copy:
            g.meta.update(self.meta)
        return g

    def get_signature_hash(self, *args, **kwargs):
        return hash_cgr_string(self.get_signature(*args, **kwargs))

    def get_signature(self, isotope=False, stereo=False, hybridization=False, neighbors=False, element=True,
                      flush_cache=False,  weights=None, hyb=False):
        """
        return string representation of structure

        :param weights: dict of atoms in keys and orders in values
        :param isotope: set isotope marks
        :param stereo: set stereo marks
        :param hybridization: set hybridization mark of atom
        :param neighbors: set neighbors count mark of atom
        :param element: set elements marks
        :param flush_cache: recalculate signature if True
        :param hyb: deprecated. see hybridization arg
        """
        if hyb:
            warn('attr hyb is deprecated, use hybridization instead', DeprecationWarning)
            hybridization = hyb

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

    def fix_data(self, copy=False, nodes_bunch=None, edges_bunch=None):
        g = self.copy() if copy else self
        for a in ((a for _, a in g.nodes(data=True)) if nodes_bunch is None else
                  (g.nodes[x] for x in g.nbunch_iter(nodes_bunch))):
            sp_new = self._attr_renew(a, self._node_marks)
            cm_new = self.__node_attr_clear(a)
            a.clear()
            a.update(sp_new)
            a.update(cm_new)

        for *_, a in g.edges(nbunch=edges_bunch, data=True):
            new = self._attr_renew(a, self._edge_marks)
            a.clear()
            a.update(new)

        if copy:
            return g
        self.flush_cache()

    def flush_cache(self):
        self.__weights = self.__signatures = self.__pickle = self.__hash = self.__stereo_cache = None

    def fresh_copy(self):
        """return a fresh copy graph with the same data structure but without atoms, bonds and metadata.
        """
        return type(self)()

    def __node_attr_clear(self, attr):
        new_attr = {}
        for s in self._node_base:
            ls = attr.get(s)
            if ls is not None:
                new_attr[s] = ls
        return new_attr

    @abstractmethod
    def _attr_renew(self, *args, **kwargs):
        """marks cleaning and fixing"""
        pass

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

    __meta = __visible = __bond_cache = __weights = __signatures = __pickle = __stereo_cache = None
    __hash = None
    __attrs = dict(source='atom1', target='atom2', name='atom', link='bonds')

    @property
    @abstractmethod
    def _node_save(self):
        pass

    @property
    @abstractmethod
    def _edge_save(self):
        pass

    @property
    @abstractmethod
    def _node_base(self) -> Iterable:
        pass

    @property
    @abstractmethod
    def _edge_marks(self):
        pass

    @property
    @abstractmethod
    def _node_marks(self):
        pass

    @classmethod
    @abstractmethod
    def _atom_container(cls, *args, **kwargs):
        pass

    @classmethod
    @abstractmethod
    def _stereo_container(cls, *args, **kwargs):
        pass

    @classmethod
    @abstractmethod
    def _bond_container(cls, *args, **kwargs):
        pass
