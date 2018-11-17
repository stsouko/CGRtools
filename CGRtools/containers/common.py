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
from hashlib import md5, sha256
from itertools import cycle
from networkx import Graph, relabel_nodes, connected_components
from networkx.readwrite.json_graph import node_link_graph, node_link_data
from ..algorithms.morgan import Morgan
from ..algorithms.strings import StringCommon
from ..attributes import Bond
from ..periodictable import elements_list, radical_unmap


class BaseContainer(Graph, StringCommon, Morgan, ABC):
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
            s.graph.clean()
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

    def union(self, g):
        if not isinstance(g, BaseContainer):
            raise TypeError('BaseContainer subclass expected')
        if self._node.keys() & set(g):
            raise KeyError('mapping of graphs is not disjoint')

        # dynamic container resolving
        qc = _search_subclass('QueryContainer')
        cc = _search_subclass('CGRContainer')

        if isinstance(self, qc):  # self has precedence
            u = type(self)()
        elif isinstance(g, qc):
            u = type(g)()
        elif isinstance(self, cc):
            u = type(self)()
        elif isinstance(g, cc):
            u = type(g)()
        else:
            u = type(self)()

        for n, a in self._node.items():
            u.add_atom(a, n)
        for n, a in g._node.items():
            u.add_atom(a, n)

        for n, m_b in self._adj.items():
            for m, b in m_b.items():
                u.add_bond(n, m, b)
        for n, m_b in g._adj.items():
            for m, b in m_b.items():
                u.add_bond(n, m, b)

    def compose(self, g, balance=False):
        """
        compose 2 graphs to CGR

        :param g: Molecule or CGR Container
        :return: CGRContainer or QueryContainer
        """
        if not isinstance(g, BaseContainer):
            raise TypeError('BaseContainer subclass expected')

        # dynamic container resolving
        qc = _search_subclass('QueryContainer')
        cc = _search_subclass('CGRContainer')

        # try to use custom containers
        if isinstance(self, qc):  # self has precedence
            h = type(self)()
        elif isinstance(g, qc):
            h = type(g)()
        elif isinstance(self, cc):
            h = type(self)()
        elif isinstance(g, cc):
            h = type(g)()
        else:  # use common CGR
            h = cc()

        common = self._node.keys() & g
        unique_reagent = self._node.keys() - common
        unique_product = g._node.keys() - common

        for n in common:
            h.add_atom(h.node_attr_dict_factory(self._node[n], g._node[n]), n)
        for n in unique_reagent:
            h.add_atom(h.node_attr_dict_factory(self._node[n], self._node[n]), n)
        for n in unique_product:
            h.add_atom(h.node_attr_dict_factory(g._node[n], g._node[n]), n)

        skin_reagent = defaultdict(list)
        seen = set()
        for n, m_b in self._adj.items():
            seen.add(n)
            if n in common:
                for m, b in m_b.items():
                    if m in seen:
                        continue
                    if m in common:
                        h.add_bond(n, m, h.edge_attr_dict_factory(b, g._adj[n][m]))
                    else:
                        skin_reagent[n].append(m)
                        h.add_bond(n, m, h.edge_attr_dict_factory(b, Bond(None, allow_none=True)))
            else:
                for m, b in m_b.items():
                    if m in seen:
                        continue
                    if m in common:
                        skin_reagent[m].append(n)
                        h.add_bond(n, m, h.edge_attr_dict_factory(b, Bond(None, allow_none=True)))
                    else:
                        h.add_bond(n, m, h.edge_attr_dict_factory(b, b))

        skin_product = defaultdict(list)
        seen = set()
        for n, m_b in g._adj.items():
            seen.add(n)
            if n in common:
                for m, b in m_b.items():
                    if m in seen:
                        continue
                    if m in common:
                        h.add_bond(n, m, h.edge_attr_dict_factory(self._adj[n][m], b))
                    else:
                        skin_product[n].append(m)
                        h.add_bond(n, m, h.edge_attr_dict_factory(Bond(None, allow_none=True), b))
            else:
                for m, b in m_b.items():
                    if m in seen:
                        continue
                    if m in common:
                        skin_product[m].append(n)
                        h.add_bond(n, m, h.edge_attr_dict_factory(Bond(None, allow_none=True), b))
                    else:
                        h.add_bond(n, m, h.edge_attr_dict_factory(b, b))

        if balance:
            """ calc unbalanced charges and radicals for skin atoms
            """
            meta = h.meta
            for n in (skin_reagent.keys() | skin_product.keys()):
                lost = skin_reagent[n]
                cycle_lost = cycle(lost)
                new = skin_product[n]
                cycle_new = cycle(new)
                atom = h._node[n]
                dr = atom.p_radical - atom.radical
                # radical balancing
                if dr > 0:  # radical added or increased.
                    for _, m in zip(range(dr), cycle_lost):  # homolysis
                        s_atom = h._node[m]
                        s_atom.p_multiplicity = radical_unmap[s_atom.p_radical + 1]
                        meta.setdefault('rule #14. atom lost. common atom radical added or increased. '
                                        'lost atom radical added', []).append((m, n))
                    for m in lost[dr:]:
                        meta.setdefault('rule #15. atom lost. common atom radical added or increased. '
                                        'lost atom radical unchanged', []).append((m, n))
                elif dr < 0:  # radical removed or decreased.
                    if n in skin_product:
                        for m in lost:
                            meta.setdefault('rule #20. atom lost. common atom radical removed or decreased. '
                                            'lost atom radical unchanged', []).append((m, n))
                    else:
                        for _, m in zip(range(-dr), cycle_lost):  # radical elimination
                            s_atom = h._node[m]
                            s_atom.p_multiplicity = radical_unmap[s_atom.p_radical + 1]
                            meta.setdefault('rule #21. atom lost. common atom radical removed or decreased. '
                                            'lost atom radical added', []).append((m, n))
                        for m in lost[-dr:]:
                            meta.setdefault('rule #20. atom lost. common atom radical removed or decreased. '
                                            'lost atom radical unchanged', []).append((m, n))
                else:
                    env = h.environment(n)
                    sv = atom.get_valence([(b.reagent, a.reagent) for b, a in env if b.order])
                    pv = atom.p_get_valence([(b.product, a.product) for b, a in env if b.p_order])
                    sh, ph = h.atom_total_h(n)

                    dv = pv - sv
                    dh = ph - sh
                    dc = atom.p_charge - atom.charge

                    if not (dv or dh or dc):  # common atom unchanged. Substitution, Elimination
                        for m in skins:
                            meta.setdefault('rule #1. atom lost. common atom unchanged. '
                                            'substitution, elimination, addition', []).append((m, n))
                    elif dv == dh == dc < 0:  # explicit hydrogen removing
                        for m in skins:
                            h._node[m].p_charge = 1
                            meta.setdefault('rule #4. atom lost. common atom deprotonation', []).append((m, n))
                    else:
                        for m in skins:
                            meta.setdefault('rule #5. atom lost. common atom changed. '
                                            'convert to reduction or oxidation', []).append((m, n))

                        pth = ph + sum(h.atom_total_h(x)[1] for x in skins)
                        if n in skin_product:
                            sth = sh + sum(h.atom_total_h(x)[0] for x in skin_product[n])
                        else:
                            sth = sh
                        dth = pth - sth

            for n, skins in skin_product.items():
                cycle_skins = cycle(skins)
                atom = h._node[n]
                dr = atom.p_radical - atom.radical
                # radical balancing
                if dr > 0:  # radical added or increased.
                    if n in skin_reagent:
                        for m in skins:
                            meta.setdefault('rule #16. atom new. common atom radical added or increased. '
                                            'new atom radical unchanged', []).append((m, n))
                    else:
                        for _, m in zip(range(dr), cycle_skins):  # radical addition
                            s_atom = h._node[m]
                            s_atom.multiplicity = radical_unmap[s_atom.radical + 1]
                            meta.setdefault('rule #17. atom new. common atom radical added or increased. '
                                            'new atom radical added', []).append((m, n))
                        for m in skins[dr:]:
                            meta.setdefault('rule #16. atom new. common atom radical added or increased. '
                                            'new atom radical unchanged', []).append((m, n))
                elif dr < 0:  # radical removed or decreased.
                    for _, m in zip(range(-dr), cycle_skins):  # recombination
                        s_atom = h._node[m]
                        s_atom.multiplicity = radical_unmap[s_atom.radical + 1]
                        meta.setdefault('rule #18. atom new. common atom radical removed or decreased. '
                                        'new atom radical added', []).append((m, n))
                    for m in skins[-dr:]:
                        meta.setdefault('rule #19. atom new. common atom radical removed or decreased. '
                                        'new atom radical unchanged', []).append((m, n))
                else:
                    env = h.environment(n)
                    sv = atom.get_valence([(b.reagent, a.reagent) for b, a in env if b.order])
                    pv = atom.p_get_valence([(b.product, a.product) for b, a in env if b.p_order])
                    sh, ph = h.atom_total_h(n)

                    dv = pv - sv
                    dh = ph - sh
                    dc = atom.p_charge - atom.charge

                    if not (dv or dh or dc):  # common atom unchanged. Substitution, Addition
                        for m in skins:
                            meta.setdefault('rule #2. atom new. common atom unchanged. '
                                            'substitution, elimination, addition', []).append((m, n))
                    elif dv == dh == dc > 0:  # explicit hydrogen addition
                        for m in skins:
                            h._node[m].charge = 1
                            h.meta.setdefault('rule #3. atom new. common atom protonation', []).append((m, n))
                    else:
                        for m in skins:
                            meta.setdefault('rule #6. atom new. common atom changed. '
                                            'convert to reduction or oxidation', []).append((m, n))

                        sth = sh + sum(h.atom_total_h(x)[0] for x in skins)
                        if n in skin_reagent:
                            pth = ph + sum(h.atom_total_h(x)[1] for x in skin_reagent[n])
                        else:
                            pth = ph
                        dth = pth - sth

            for n, sp in reverse_ext.items():

                # charge neutralization
                if dc > 0:
                    for _ in range(dc):
                        h.meta.setdefault('rule #7. charge neutralization. hydroxide radical added',
                                          []).append(h.add_atom(O(multiplicity=2), O(charge=-1)))
                elif dc < 0:
                    for _ in range(-dc):
                        h.meta.setdefault('rule #8. charge neutralization. hydrogen radical added',
                                          []).append(h.add_atom(H(multiplicity=2), H(charge=1)))

                # hydrogen balancing
                if dth > 0:
                    red_e = 0
                    for m in sp['products']:
                        if h.nodes[m]['element'] == 'H':  # set reduction H if explicit H count increased
                            h.nodes[m]['s_radical'] = 2
                            red_e += 1
                            h.meta.setdefault('rule #11. protonation. new explicit hydrogen radical added',
                                              []).append(m)

                    red = []
                    for _ in range(dth - red_e):  # add reduction agents
                        m = h.add_atom(H(multiplicity=2), H())
                        red.append(m)
                        h.meta.setdefault('rule #10. protonation. hydrogen radical added', []).append(m)
                    red = iter(red)

                    dih = sub(*h.atom_implicit_h(n))
                    if dih < 0:  # attach reduction H to central atom if implicit H atoms count increased
                        for _ in range(-dih):
                            m = next(red)
                            h.add_bond(m, n, None)
                            h.meta.setdefault('rule #12. protonation. new implicit hydrogen radical added',
                                              []).append(m)

                    for m in sp['reagents']:  # attach reduction H if detached group implicit H count increased
                        dih = sub(*h.atom_implicit_h(m))
                        if dih < 0:
                            for _ in range(-dih):
                                o = next(red)
                                h.add_bond(o, m, None)
                elif dth < 0:
                    oxo = []
                    for _ in range(-dth):
                        m = h.add_atom(O(multiplicity=2), O())
                        oxo.append(m)
                        h.meta.setdefault('rule #9. deprotonation. hydroxide radical added', []).append(m)
                    oxo = iter(oxo)

                    for m in sp['reagents']:
                        if h.nodes[m]['element'] == 'H':
                            o = next(oxo)
                            h.add_bond(o, m, None)
                            h.meta.setdefault('rule #13. hydrogen accepting by hydroxide radical added',
                                              []).append(m)

        return h

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
        self.__weights = self.__signatures = self.__pickle = self.__hash = self.__stereo_cache = None

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

    def __xor__(self, other):
        """
        G ^ H is CGR generation
        """
        return self.compose(other)

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

    __visible = __weights = __signatures = __pickle = __stereo_cache = __hash = None
    __attrs = dict(source='atom1', target='atom2', name='atom', link='bonds')


def _search_subclass(name, start=BaseContainer):
    for x in start.__subclasses__():
        if x.__name__ == name:
            return x
        d = _search_subclass(name, x)
        if d is not None:
            return d
