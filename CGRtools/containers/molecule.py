# -*- coding: utf-8 -*-
#
#  Copyright 2017 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
from collections import namedtuple, defaultdict
from itertools import chain
from networkx import Graph, relabel_nodes
from networkx.readwrite.json_graph import node_link_graph, node_link_data
from warnings import warn
from .. import InvalidData, InvalidAtom, InvalidStereo
from ..algorithms import get_morgan, CGRstring, hash_cgr_string, Valence, ValenceError, pyramid_volume
from ..periodictable import elements


class MoleculeContainer(Graph, Valence):
    def __init__(self, data=None, meta=None):
        Graph.__init__(self, data)
        Valence.__init__(self)
        if isinstance(meta, dict):
            self.__meta = meta

    def __dir__(self):
        if self.__visible is None:
            self.__visible = [self.pickle.__name__, self.unpickle.__name__, self.copy.__name__,
                              self.substructure.__name__, self.remap.__name__, self.add_stereo.__name__,
                              self.get_morgan.__name__, self.get_fear.__name__, self.get_fear_hash.__name__,
                              self.get_environment.__name__, self.fix_data.__name__, self.reset_query_marks.__name__,
                              self.atom.__name__, self.bond.__name__, self.add_atom.__name__, self.add_bond.__name__,
                              self.explicify_hydrogens.__name__, self.implicify_hydrogens.__name__,
                              'meta', 'bonds_count', 'atoms_count']  # properties names inaccessible
        return self.__visible

    def pickle(self):
        """ return json serializable CGR or Molecule
        """
        node_marks = self._node_save
        edge_marks = self._edge_save

        g = Graph()
        g.add_nodes_from((n, {k: v for k, v in a.items() if k in node_marks}) for n, a in self.nodes(data=True))
        g.add_edges_from((n, m, {k: v for k, v in a.items() if k in edge_marks}) for n, m, a in self.edges(data=True))

        data = node_link_data(g, attrs=self._attrs)
        data.update(meta=self.meta, s_only=True)
        return data

    @classmethod
    def unpickle(cls, data):
        """ convert json serializable CGR into MoleculeContainer object instance
        """
        if not data['s_only']:
            raise InvalidData('pickled data is invalid molecule. try cgr.unpickle')

        g = MoleculeContainer(node_link_graph(data, attrs=cls._attrs), data['meta'])
        g.fix_data()
        return g

    def copy(self):
        return self.__class__(self, self.meta)

    def substructure(self, nbunch, meta=False):
        sub = self.__class__(super().subgraph(nbunch), self.meta if meta else None)
        sub._fix_stereo()
        return sub

    def remap(self, mapping, copy=False):
        g = relabel_nodes(self, mapping, copy=copy)
        if copy:
            g = self.__class__(g, self.meta)
        return g

    @property
    def meta(self):
        if self.__meta is None:
            self.__meta = {}
        return self.__meta

    def get_fear_hash(self, weights=None, isotope=False, stereo=False, hyb=False, element=True):
        return hash_cgr_string(self.get_fear(weights=weights, isotope=isotope, stereo=stereo, hyb=hyb, element=element))

    def get_fear(self, weights=None, isotope=False, stereo=False, hyb=False, element=True, flush_cache=False):
        """
        :param weights: dict of atoms in keys and orders in values
        :param isotope: set isotope marks
        :param stereo: set stereo marks
        :param hyb: set hybridization mark of atom
        :param element: set elements marks
        :param flush_cache: recalculate fear if True
        :return: string representation of CGR
        """
        if flush_cache or self._fears is None:
            self._fears = {}
        k = (isotope, element, stereo, hyb)
        return self._fears.get(k) or self._fears.setdefault(k, CGRstring(isotope, stereo, hyb, element)
            (self, weights or self.get_morgan(isotope, element, stereo, flush_cache)))

    def add_atom(self, element, charge, radical=None, _map=None, mark='0', x=0, y=0, z=0):
        if element not in elements:
            raise InvalidData('element %s - not exists' % element)
        if _map is None:
            _map = max(self, default=0) + 1
        elif _map in self:
            raise InvalidData('mapping exists')
        if radical not in (None, 1, 2, 3):
            raise InvalidData('only monovalent (2), bivalent (1 singlet, 3 triplet) or None accepted')
        if not self._check_charge_radical(element, charge, radical=self._radical_map[radical]):
            raise InvalidData('charge and/or radical values impossible for this element')

        self.add_node(_map, element=element, s_charge=charge, mark=mark, s_x=x, s_y=y, s_z=z, map=_map)
        if radical:
            self.nodes[_map]['s_radical'] = radical

        self._weights = self._fears = self._pickle = None
        return _map

    def add_bond(self, atom1, atom2, mark, *, ignore=False):
        if atom1 not in self or atom2 not in self:
            raise InvalidData('atoms not found')
        if self.has_edge(atom1, atom2):
            raise InvalidData('bond exists')
        if mark not in (1, 2, 3, 4, 9):
            raise InvalidData('invalid bond mark')

        if not ignore:
            self._check_bonding(atom1, atom2, mark)
        self.add_edge(atom1, atom2, s_bond=mark)
        self._weights = self._fears = self._pickle = None

    def add_stereo(self, atom1, atom2, mark):
        if not self.has_edge(atom1, atom2):
            raise InvalidAtom('atom or bond not found')

        if self.nodes[atom1].get('s_stereo'):
            raise InvalidStereo('atom has stereo. change impossible')

        tmp = [(x, y['s_bond']) for x, y in self[atom1].items()]
        neighbors = [x for x, _ in tmp]

        if self.nodes[atom1]['s_z'] or any(self.nodes[x]['s_z'] for x in neighbors):
            raise InvalidStereo('molecule have 3d coordinates. bond up/down stereo unusable')

        neighbors_e = [self.nodes[x]['element'] for x in neighbors]
        implicit = self.atom_implicit_h(atom1)
        if implicit > 1 or implicit == 1 and 'H' in neighbors_e or neighbors_e.count('H') > 1:
            raise InvalidStereo('stereo impossible. too many H atoms')

        bonds = [x for _, x in tmp]
        total = implicit + len(neighbors)
        if total == 4:  # tetrahedron
            if any(x != 1 for x in bonds):
                raise InvalidStereo('only single bonded tetrahedron acceptable')

            weights = self.get_morgan(stereo=True)
            if len(neighbors) != len(set(weights[x] for x in neighbors)):
                raise InvalidStereo('stereo impossible. neighbors equivalent')

            order = sorted(neighbors, key=weights.get)
            vol = pyramid_volume(*((y['s_x'], y['s_y'], 0 if x != atom2 else mark) for x, y in
                                   ((x, self.nodes[x]) for x in (chain((atom1,), order) if implicit else order))))
            if not vol:
                raise InvalidStereo('unknown')

            self.nodes[atom1]['s_stereo'] = vol > 0 and 1 or -1
            self._weights = self._fears = self._pickle = None

    def bond(self, atom1, atom2):
        if self.__bond_cache is None:
            self.__bond_cache = defaultdict(dict)
        if atom2 not in self.__bond_cache[atom1]:
            try:
                tmp = self[atom1][atom2]
                res = self._bond_container(**{x: tmp.get(y) for x, y in self._bond_marks.items()})
                self.__bond_cache[atom1][atom2] = self.__bond_cache[atom2][atom1] = res
            except KeyError:
                raise InvalidAtom('atom or bond not found')
        return self.__bond_cache[atom1][atom2]

    def atom(self, n):
        if self.__atom_cache is None:
            self.__atom_cache = {}
        if n not in self.__atom_cache:
            try:
                tmp = self.nodes[n]
                self.__atom_cache[n] = self._atom_container(**{x: tmp.get(y) for x, y in self._atom_marks.items()})
            except KeyError:
                raise InvalidAtom('atom not found')
        return self.__atom_cache[n]

    def get_stereo(self, atom1, atom2):
        if self.__stereo_cache is None:
            self.__stereo_cache = {}
            nodes = list(self.nodes(data=True))
            while True:
                failed = []
                for i, tmp in enumerate(nodes, start=1):
                    n, attr = tmp
                    s = attr.get('s_stereo')
                    if s:
                        neighbors = list(self.neighbors(n))
                        len_n = len(neighbors)
                        if len_n in (3, 4):  # tetrahedron
                            if attr['s_z'] or any(self.nodes[x]['s_z'] for x in neighbors):
                                continue  # 3d molecules ignored

                            weights = self.get_morgan(stereo=True)
                            order = sorted(neighbors, key=weights.get)
                            for _ in range(len_n):
                                if (order[0], n) in self.__stereo_cache:
                                    order.append(order.pop(0))
                                else:
                                    failed.append(tmp)
                                    break
                            else:
                                failed.insert(0, tmp)
                                failed.extend(nodes[i:])
                                self.__stereo_cache = None
                                nodes = failed
                                break

                            if len_n == 4:
                                zero = self.nodes[order[0]]
                                zero = (zero['s_x'], zero['s_y'], 1)
                                first = self.nodes[order[1]]
                                first = (first['s_x'], first['s_y'], 0)
                            else:
                                zero = (attr['s_x'], attr['s_y'], 0)
                                first = self.nodes[order[0]]
                                first = (first['s_x'], first['s_y'], 1)

                            second = self.nodes[order[-2]]
                            third = self.nodes[order[-1]]
                            vol = pyramid_volume(zero, first,
                                                 (second['s_x'], second['s_y'], 0), (third['s_x'], third['s_y'], 0))

                            self.__stereo_cache[(n, order[0])] = 1 if vol > 0 and s == 1 or vol < 0 and s == -1 else -1
                else:
                    break

        return self.__stereo_cache.get((atom1, atom2), None)

    @property
    def bonds_count(self):
        return self.size()

    @property
    def atoms_count(self):
        return self.order()

    def get_environment(self, atoms, dante=False, deep=0):
        """ get subgraph with atoms and their neighbors
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

    def get_morgan(self, isotope=False, element=True, stereo=False, flush_cache=False):
        if flush_cache or self._weights is None:
            self._weights = {}
        k = (isotope, element, stereo)
        return self._weights.get(k) or self._weights.setdefault(k, get_morgan(self, isotope=isotope, element=element,
                                                                              stereo=stereo))

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
        self._weights = self._fears = self._pickle = None

    def reset_query_marks(self, copy=False):
        """
        set or reset hyb and neighbors marks to atoms.
        :param copy: if True return copy of graph and keep existing as is
        :return: graph if copy True else None
        """
        g = self.copy() if copy else self
        b, h, n = 's_bond', 's_hyb', 's_neighbors'
        for i, attr in g.nodes(data=True):
            label = dict(s_hyb=1, p_hyb=1, s_neighbors=0, p_neighbors=0)
            # hyb 1- sp3; 2- sp2; 3- sp1; 4- aromatic
            for node, bond in g[i].items():
                b_type = bond.get(b)
                if b_type and g.nodes[node]['element'] != 'H':
                    label[n] += 1
                if b_type in (1, None) or label[h] in (3, 4):
                    continue
                elif b_type == 4:
                    label[h] = 4
                elif b_type == 3 or (b_type == 2 and label[h] == 2):  # Если есть 3-я или две 2-х связи, то sp1
                    label[h] = 3
                elif b_type == 2 and label[h] != 3:
                    # Если есть 2-я связь, но до этого не было найдено другой 2-й, 3-й, или аром.
                    label[h] = 2

            attr.update(label)
        self._fears = self._pickle = None
        return g if copy else None

    def implicify_hydrogens(self):
        """
        remove explicit hydrogent if possible
        :return: number of removed hydrogens
        """
        explicit = {}
        c = 0
        for n, attr in self.nodes(data=True):
            if attr['element'] == 'H':
                m = next(self.neighbors(n))
                if self.nodes[m]['element'] != 'H':
                    explicit.setdefault(m, []).append(n)

        for n, h in explicit.items():
            atom = self.nodes[n]
            implicit = self._get_implicit_h(atom['element'], atom['s_charge'],
                                            [y['s_bond'] for x, y in self[n].items() if x not in h],
                                            radical=atom.get('s_radical', 0))
            if implicit:
                for x in h:
                    self.remove_node(x)
                    c += 1
        self._weights = self._fears = self._pickle = None
        return c

    def explicify_hydrogens(self):
        """
        add explicit hydrogens to atoms
        :return: number of added atoms
        """
        tmp = []
        for n, attr in self.nodes(data=True):
            if attr['element'] != 'H':
                for _ in range(self.atom_implicit_h(n)):
                    tmp.append(n)
        for n in tmp:
            self.add_bond(n, self.add_atom('H', 0), 1)

        self._weights = self._fears = self._pickle = None
        return len(tmp)

    def split_ions(self):
        for n, attr in self.nodes(data=True):
            if attr['element'] in self._metals:
                pass
        self._weights = self._fears = self._pickle = None

    def atom_implicit_h(self, atom):
        attr = self.nodes[atom]
        return self._get_implicit_h(attr['element'], attr['s_charge'], [x['s_bond'] for x in self[atom].values()],
                                    radical=self._radical_map[attr.get('s_radical')])

    def flush_cache(self):
        self._weights = self._fears = self._pickle = None

    @staticmethod
    def _attr_renew(attr, marks):
        new_attr = {}
        for s in marks:
            ls = attr.get(s)
            if ls is not None:
                new_attr[s] = ls
        return new_attr

    def _fix_stereo(self):
        tetrahedrons = []
        for atom, attr in self.nodes(data=True):
            neighbors = [self.nodes[x]['element'] for x in self.neighbors(atom)]
            implicit = self.atom_implicit_h(atom)
            if implicit > 1 or implicit == 1 and 'H' in neighbors or neighbors.count('H') > 1:
                continue

            total = implicit + len(neighbors)
            # tetrahedron
            if total == 4:
                tetrahedrons.append((atom, neighbors, implicit))

        if tetrahedrons:
            weights = self.get_morgan(stereo=True)
            for atom, neighbors, implicit in tetrahedrons:
                w_n = [weights[x] for x in neighbors]
                if len(w_n) != len(set(w_n)):
                    continue
                if implicit:
                    pass
                else:
                    pass

            self._weights = self._fears = self._pickle = None

    def _check_bonding(self, atom1, atom2, mark, label='s'):
        for atom, reverse in ((atom1, atom2), (atom2, atom1)):
            a = self.nodes[atom]
            lb = '%s_bond' % label
            lc = '%s_charge' % label
            lr = '%s_radical' % label
            try:
                if not self._check_valence(a['element'], a[lc],
                                           [x[lb] for x in self[atom].values() if x.get(lb)] + [mark],
                                           radical=self._radical_map[a.get(lr)]):
                    raise InvalidData('valence error')
            except ValenceError:
                tmp = [(y[lb], self.nodes[x]['element'])
                       for x, y in self[atom].items() if y.get(lb)] + [(mark, self.nodes[reverse]['element'])]
                if not self._check_valence(a['element'], a[lc], [x for x, _ in tmp], self._radical_map[a.get(lr)],
                                           neighbors=[x for _, x in tmp]):
                    raise InvalidData('valence error')

    @classmethod
    def __node_attr_clear(cls, attr):
        new_attr = {}
        for s in cls._node_base:
            ls = attr.get(s)
            if ls is not None:
                new_attr[s] = ls
        return new_attr

    def __str__(self):
        return self.get_fear(isotope=True, stereo=True)

    def __repr__(self):
        if self._pickle is None:
            self._pickle = '%s.unpickle(%s)' % (self.__class__.__name__, self.pickle())
        return self._pickle

    _attrs = dict(source='atom1', target='atom2', name='atom', link='bonds')
    _atom_marks = dict(charge='s_charge', stereo='s_stereo', neighbors='s_neighbors', hyb='s_hyb',
                       element='element', isotope='isotope', mark='mark', radical='s_radical')
    _bond_marks = dict(bond='s_bond', stereo='s_stereo')
    _atom_container = namedtuple('Atom', ['element', 'isotope', 'mark', 'charge', 'stereo', 'radical',
                                          'neighbors', 'hyb'])
    _bond_container = namedtuple('Bond', ['bond', 'stereo'])
    _node_base = ('element', 'isotope', 'mark', 's_x', 's_y', 's_z')
    _node_marks = ('s_neighbors', 's_hyb', 's_charge', 's_stereo', 's_radical')
    _node_save = _node_marks + _node_base
    _edge_save = _edge_marks = ('s_bond', 's_stereo')
    _radical_map = {1: 2, 2: 1, 3: 2, None: 0}
    _bond_map = {1: 1, 2: 2, 3: 3, 4: 1.5, 9: 1}
    __meta = __visible = __atom_cache = __bond_cache = __stereo_cache = _weights = _fears = _pickle = None

    def subgraph(self, *args, **kwargs):
        warn('subgraph name is deprecated. use sustructure instead', DeprecationWarning)
        return self.substructure(*args, **kwargs)
