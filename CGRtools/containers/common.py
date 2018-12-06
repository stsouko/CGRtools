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
from hashlib import md5, sha256
from networkx import Graph, relabel_nodes, connected_components
from networkx.readwrite.json_graph import node_link_graph, node_link_data
from ..algorithms import Morgan, SSSR, StringCommon
from ..periodictable import elements_list


class BaseContainer(Graph, StringCommon, Morgan, SSSR, ABC):
    def __dir__(self):
        return list(self._visible) or super().__dir__()

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
        return next(x for x in BaseContainer.__subclasses__() if x.__name__ == data['class'])(graph)

    def environment(self, atom):
        """
        pairs of (bond, atom) connected to atom

        :param atom: number
        :return: list
        """
        return [(bond, self._node[n]) for n, bond in self._adj[atom].items()]

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

    def connected_components(self):
        return connected_components(self)

    def split(self, meta=False):
        """
        split disconnected structure to connected substructures

        :param meta: copy metadata to each substructure
        :return: list of substructures
        """
        return [self.substructure(c, meta) for c in connected_components(self)]

    def union(self, other):
        if not isinstance(other, BaseContainer):
            raise TypeError('BaseContainer subclass expected')
        if self._node.keys() & set(other):
            raise KeyError('mapping of graphs is not disjoint')

        # dynamic container resolving
        qc = next(x for x in BaseContainer.__subclasses__() if x.__name__ == 'QueryCGRContainer')
        qq = next(x for x in BaseContainer.__subclasses__() if x.__name__ == 'QueryContainer')
        cc = next(x for x in BaseContainer.__subclasses__() if x.__name__ == 'CGRContainer')

        if isinstance(self, qc):
            u = type(self)()
        elif isinstance(other, qc):
            u = type(other)()
        elif isinstance(self, cc):
            if isinstance(other, qq):  # force QueryCGRContainer
                u = qc()
            else:
                u = type(self)()
        elif isinstance(other, cc):
            if isinstance(self, qq):
                u = qc()
            else:
                u = type(other)()
        elif isinstance(self, qq):  # self has precedence
            u = type(self)()
        elif isinstance(other, qq):
            u = type(other)()
        else:
            u = type(self)()

        for n, a in self._node.items():
            u.add_atom(a, n)
        for n, a in other._node.items():
            u.add_atom(a, n)

        seen = set()
        for n, m_b in self._adj.items():
            seen.add(n)
            for m, b in m_b.items():
                if m not in seen:
                    u.add_bond(n, m, b)
        for n, m_b in other._adj.items():
            seen.add(n)
            for m, b in m_b.items():
                if m not in seen:
                    u.add_bond(n, m, b)
        return u

    def remap(self, mapping, copy=False):
        return relabel_nodes(self, mapping, copy)

    def get_signature_hash(self, *args, **kwargs):
        """
        concatenated md5 and sha256 hashes of cgr string
        :return: 48 bytes length string
        """
        bs = self.get_signature(*args, **kwargs).encode()
        return md5(bs).digest() + sha256(bs).digest()

    def get_signature(self, start=None, depth_limit=None, atom=True, isotope=False, stereo=False, hybridization=False,
                      neighbors=False, flush_cache=False,  weights=None):
        """
        return string representation of structure

        :param start: root atom map. need for augmented atom signatures
        :param depth_limit: dept of augmented atom signature
        :param weights: dict of atoms in keys and orders in values
        :param isotope: set isotope marks
        :param stereo: set stereo marks
        :param hybridization: set hybridization mark of atom
        :param neighbors: set neighbors count mark of atom
        :param atom: set elements marks
        :param flush_cache: recalculate signature if True
        """
        if flush_cache or self.__signatures is None:
            self.__signatures = {}

        if weights is None:
            k = (start, depth_limit, atom, isotope, stereo, hybridization, neighbors)
            if k in self.__signatures:
                return self.__signatures[k]
            weights = self.get_morgan(atom, isotope, stereo, hybridization, neighbors, flush_cache)
            sg = self._stringify(start, depth_limit, weights, atom, isotope, stereo, hybridization, neighbors)
            self.__signatures[k] = sg
            return sg
        else:
            return self._stringify(start, depth_limit, weights, atom, isotope, stereo, hybridization, neighbors)

    def get_morgan(self, atom=True, isotope=False, stereo=False, hybridization=False, neighbors=False,
                   flush_cache=False):
        if flush_cache or self.__weights is None:
            self.__weights = {}
        k = (atom, isotope, stereo, hybridization, neighbors)
        if k in self.__weights:
            return self.__weights[k]
        return self.__weights.setdefault(k, self._morgan(atom, isotope, stereo, hybridization, neighbors))

    @abstractmethod
    def _stringify(self, start=None, depth_limit=None, weights=None, atom=True, isotope=True, stereo=False,
                   hybridization=False, neighbors=False) -> str:
        """
        container specific signature generation
        """
        pass

    def flush_cache(self):
        self.__weights = self.__signatures = self.__pickle = self.__hash = None

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

    def __or__(self, other):
        """
        G | H is union of graphs
        """
        return self.union(other)

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

    __weights = __signatures = __pickle = __hash = None
    __attrs = dict(source='atom1', target='atom2', name='atom', link='bonds')
    _visible = ()


__all__ = ['BaseContainer']
