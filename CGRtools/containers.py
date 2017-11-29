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
from itertools import chain, repeat, zip_longest
from networkx import Graph, relabel_nodes
from networkx.readwrite.json_graph import node_link_graph, node_link_data
from warnings import warn
from . import InvalidData, InvalidAtom, InvalidStereo
from .algorithms import get_morgan, CGRstring, hash_cgr_string, Valence, ValenceError, pyramid_volume
from .periodictable import elements

CGRTemplate = namedtuple('CGRTemplate', ['reagents', 'products', 'meta'])
MatchContainer = namedtuple('MatchContainer', ['mapping', 'meta', 'patch'])


class MergedReaction(namedtuple('MergedReaction', ['reagents', 'products'])):
    def copy(self):
        return MergedReaction(self.reagents.copy(), self.products.copy())


class MoleculeContainer(Graph, Valence):
    def __init__(self, meta=None):
        Graph.__init__(self)
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
        s_only = not isinstance(self, CGRContainer)
        node_marks = self._node_save
        edge_marks = self._edge_save

        g = Graph()
        g.add_nodes_from((n, {k: v for k, v in a.items() if k in node_marks}) for n, a in self.nodes(data=True))
        g.add_edges_from((n, m, {k: v for k, v in a.items() if k in edge_marks}) for n, m, a in self.edges(data=True))

        data = node_link_data(g, attrs=self.__attrs)
        data.update(meta=self.meta, s_only=s_only)
        return data

    @classmethod
    def unpickle(cls, data):
        """ convert json serializable CGR into MoleculeContainer or CGRcontainer object instance
        """
        g = node_link_graph(data, attrs=cls.__attrs)
        g.__class__ = MoleculeContainer if data['s_only'] else CGRContainer
        g.fix_data()
        g.meta.update(data['meta'])
        return g

    def copy(self):
        copy = super().copy()
        copy.__class__ = self.__class__
        copy.meta.update(self.meta)
        return copy

    def substructure(self, nbunch, meta=False):
        sub = super().subgraph(nbunch).copy()
        sub.__class__ = self.__class__
        if meta:
            sub.meta.update(self.meta)
        sub._fix_stereo()
        return sub

    def remap(self, mapping, copy=False):
        g = relabel_nodes(self, mapping, copy=copy)
        if copy:
            g.__class__ = self.__class__
        return g

    @property
    def meta(self):
        if self.__meta is None:
            self.__meta = {}
        return self.__meta

    def get_fear_hash(self, weights=None, isotope=False, stereo=False, hyb=False, element=True):
        return hash_cgr_string(self.get_fear(weights=weights, isotope=isotope, stereo=stereo, hyb=hyb, element=element))

    def get_fear(self, weights=None, isotope=False, stereo=False, hyb=False, element=True):
        """
        :param weights: dict of atoms in keys and orders in values 
        :param isotope: set isotope marks
        :param stereo: set stereo marks
        :param hyb: set hybridization mark of atom
        :param element: set elements marks
        :return: string representation of CGR
        """
        return CGRstring(isotope, stereo, hyb, element)(self, weights or get_morgan(self, isotope=isotope,
                                                                                    element=element))

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
        self.__weights = None
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
        self.__weights = None

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
                                   ((x, self.nodes[x])for x in (chain((atom1,), order) if implicit else order))))
            if not vol:
                raise InvalidStereo('unknown')

            self.nodes[atom1]['s_stereo'] = vol > 0 and 1 or -1
            self.__weights = None

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

    def get_morgan(self, isotope=False, element=True, stereo=False):
        if self.__weights is None:
            self.__weights = {}
        k = (isotope, element, stereo)
        return self.__weights.get(k) or self.__weights.setdefault(k, get_morgan(self, isotope=isotope, element=element,
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
        self.__weights = None

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
        self.__weights = None
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

        self.__weights = None
        return len(tmp)

    def split_ions(self):
        for n, attr in self.nodes(data=True):
            if attr['element'] in self._metals:
                pass
        self.__weights = None

    def atom_implicit_h(self, atom):
        attr = self.nodes[atom]
        return self._get_implicit_h(attr['element'], attr['s_charge'], [x['s_bond'] for x in self[atom].values()],
                                    radical=self._radical_map[attr.get('s_radical')])

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

            self.__weights = None

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

    __attrs = dict(source='atom1', target='atom2', name='atom', link='bonds')
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
    __meta = __visible = __atom_cache = __bond_cache = __weights = __stereo_cache = None

    def subgraph(self, *args, **kwargs):
        warn('subgraph name is deprecated. use sustructure instead', DeprecationWarning)
        return self.substructure(*args, **kwargs)


class CGRContainer(MoleculeContainer):
    def __dir__(self):
        if self.__visible is None:
            self.__visible = tmp = super().__dir__()
            tmp.append(self.get_center_atoms.__name__)
        return self.__visible

    def get_fear(self, weights=None, isotope=False, stereo=False, hyb=False, element=True):
        return CGRstring(isotope, stereo, hyb, element, True)(self, weights or get_morgan(self, isotope=isotope,
                                                                                          element=element))

    def get_center_atoms(self, stereo=False):
        """ get atoms of reaction center (dynamic bonds, stereo or charges).
        """
        nodes = set()
        for n, attr in self.nodes(data=True):
            if attr['s_charge'] != attr['p_charge'] or attr.get('s_radical') != attr.get('p_radical') or \
                   stereo and attr.get('s_stereo') != attr.get('p_stereo'):
                nodes.add(n)

        for *n, attr in self.edges(data=True):
            if attr.get('s_bond') != attr.get('p_bond') or stereo and attr.get('s_stereo') != attr.get('p_stereo'):
                nodes.update(n)

        return list(nodes)

    def add_atom(self, element, s_charge, p_charge=None, s_radical=None, p_radical=None, _map=None, mark='0',
                 s_x=0, s_y=0, s_z=0, p_x=None, p_y=None, p_z=None):
        if element not in elements:
            raise InvalidData('element %s - not exists' % element)
        if _map is None:
            _map = max(self, default=0) + 1
        elif _map in self:
            raise InvalidData('mapping exists')
        if s_radical not in (None, 1, 2, 3) or p_radical not in (None, 1, 2, 3):
            raise InvalidData('only monovalent (2), bivalent (1 singlet, 3 triplet) or None accepted')

        if p_charge is None:
            p_charge = s_charge
        if p_x is None:
            p_x = s_x
        if p_y is None:
            p_y = s_y
        if p_z is None:
            p_z = s_z

        if not (self._check_charge_radical(element, s_charge, radical=self._radical_map[s_radical]) and
                self._check_charge_radical(element, p_charge, radical=self._radical_map[p_radical])):
            raise InvalidData('charge and/or radical values impossible for this element')

        self.add_node(_map, element=element, s_charge=s_charge, p_charge=p_charge, mark=mark,
                      s_x=s_x, s_y=s_y, s_z=s_z, p_x=p_x, p_y=p_y, p_z=p_z, map=_map)

        if s_radical:
            self.nodes[_map]['s_radical'] = s_radical
        if p_radical:
            self.nodes[_map]['p_radical'] = p_radical
        return _map

    def add_bond(self, atom1, atom2, s_mark, p_mark, *, ignore=False):
        if atom1 not in self or atom2 not in self:
            raise InvalidData('atoms not found')
        if self.has_edge(atom1, atom2):
            raise InvalidData('bond exists')
        if not (s_mark or p_mark):
            raise InvalidData('empty bonds not allowed')

        if s_mark not in (1, 2, 3, 4, 9, None):
            raise InvalidData('invalid bond mark')
        elif not ignore and s_mark:
            self._check_bonding(atom1, atom2, s_mark)

        if p_mark not in (1, 2, 3, 4, 9, None):
            raise InvalidData('invalid bond mark')
        elif not ignore and p_mark:
            self._check_bonding(atom1, atom2, p_mark, label='p')

        if s_mark:
            self.add_edge(atom1, atom2, s_bond=s_mark)
        if p_mark:
            self.add_edge(atom1, atom2, p_bond=p_mark)

    def add_stereo(self, atom1, atom2, s_mark, p_mark):
        if self.has_edge(atom1, atom2):
            # todo: stereo calc
            pass

    def reset_query_marks(self, copy=False):
        """
        set or reset hyb and neighbors marks to atoms.
        :param copy: if True return copy of graph and keep existing as is
        :return: graph if copy True else None
        """
        g = self.copy() if copy else self
        for i in g:
            label = dict(s_hyb=1, p_hyb=1, sp_hyb=1, s_neighbors=0, p_neighbors=0, sp_neighbors=0)
            # hyb 1- sp3; 2- sp2; 3- sp1; 4- aromatic
            for b, h, n in (('s_bond', 's_hyb', 's_neighbors'), ('p_bond', 'p_hyb', 'p_neighbors')):
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

            for n, m, h in (('s_hyb', 'p_hyb', 'sp_hyb'), ('s_neighbors', 'p_neighbors', 'sp_neighbors')):
                label[h] = (label[n], label[m]) if label[n] != label[m] else label[n]

            g.nodes[i].update(label)
        return g if copy else None

    def implicify_hydrogens(self):
        """
        remove explicit hydrogent if possible
        :return: number of removed hydrogens
        """
        explicit = {}
        hydrogens = set()
        for_remove = []
        c = 0
        for n, attr in self.nodes(data=True):
            if attr['element'] == 'H':
                for m in self.neighbors(n):
                    if self.nodes[m]['element'] != 'H':
                        explicit.setdefault(m, []).append(n)
                    else:
                        hydrogens.add(m)
                        hydrogens.add(n)

        for n, h in explicit.items():
            atom = self.nodes[n]
            s_bonds = [y['s_bond'] for x, y in self[n].items() if x not in h and y.get('s_bond')]
            p_bonds = [y['p_bond'] for x, y in self[n].items() if x not in h and y.get('p_bond')]
            s_implicit = self._get_implicit_h(atom['element'], atom['s_charge'], s_bonds, atom.get('s_radical', 0))
            p_implicit = self._get_implicit_h(atom['element'], atom['p_charge'], p_bonds, atom.get('p_radical', 0))

            if not s_implicit and any(self[n][x].get('s_bond') for x in h):
                hydrogens.update(h)
            elif not p_implicit and any(self[n][x].get('p_bond') for x in h):
                hydrogens.update(h)
            else:
                for x in h:
                    for_remove.append(x)

        for x in for_remove:
            if x not in hydrogens:
                self.remove_node(x)
                c += 1
        return c

    def explicify_hydrogens(self):
        """
        add explicit hydrogens to atoms
        :return: number of added atoms
        """
        tmp = []
        for n, attr in self.nodes(data=True):
            if attr['element'] != 'H':
                si, pi = self.atom_implicit_h(n)
                if si or pi:
                    for s_mark, p_mark in zip_longest(repeat(1, si), repeat(1, pi)):
                        tmp.append((n, s_mark, p_mark))

        for n, s_mark, p_mark in tmp:
            self.add_bond(n, self.add_atom('H', 0), s_mark, p_mark)
        return len(tmp)

    def atom_implicit_h(self, atom):
        attr = self.nodes[atom]
        si = self._get_implicit_h(attr['element'], attr['s_charge'],
                                  [x['s_bond'] for x in self[atom].values() if x.get('s_bond')],
                                  radical=self._radical_map[attr.get('s_radical')])
        pi = self._get_implicit_h(attr['element'], attr['p_charge'],
                                  [x['p_bond'] for x in self[atom].values() if x.get('p_bond')],
                                  radical=self._radical_map[attr.get('p_radical')])
        return self.__implicit_container(si, pi)

    def _fix_stereo(self):
        pass

    @staticmethod
    def _attr_renew(attr, marks):
        new_attr = {}
        for s, p, sp in marks:
            ls = attr.get(s)
            lp = attr.get(p)

            if isinstance(ls, list):
                if isinstance(lp, list):
                    if ls == lp:
                        new_attr[sp] = new_attr[s] = new_attr[p] = list(set(ls))
                        continue
                    new_attr[sp] = list(set((x, y) for x, y in zip(ls, lp) if x != y))
                else:
                    new_attr[sp] = list(set((x, lp) for x in ls if x != lp))

                new_attr[s] = [x for x, _ in new_attr[sp]]
                new_attr[p] = [x for _, x in new_attr[sp]]
            elif isinstance(lp, list):
                new_attr[sp] = list(set((ls, x) for x in lp if x != ls))
                new_attr[s] = [x for x, _ in new_attr[sp]]
                new_attr[p] = [x for _, x in new_attr[sp]]
            elif ls != lp:
                if ls is not None:
                    new_attr[s] = ls
                if lp is not None:
                    new_attr[p] = lp
                new_attr[sp] = (ls, lp)
            elif ls is not None:
                new_attr[sp] = new_attr[s] = new_attr[p] = ls
        return new_attr

    _node_marks = tuple(('s_%s' % mark, 'p_%s' % mark, 'sp_%s' % mark)
                        for mark in ('charge', 'stereo', 'radical', 'neighbors', 'hyb'))
    _edge_marks = tuple(('s_%s' % mark, 'p_%s' % mark, 'sp_%s' % mark) for mark in ('bond', 'stereo'))

    __tmp1 = ('element', 'isotope', 'mark')
    __tmp2 = tuple(y for x in _node_marks for y in x[:2])
    __tmp3 = __tmp1 + __tmp2

    _node_base = __tmp1 + ('s_x', 's_y', 's_z', 'p_x', 'p_y', 'p_z')
    _node_save = __tmp2 + _node_base
    _edge_save = tuple(y for x in _edge_marks for y in x[:2])
    _atom_marks = {x: x for x in __tmp3}
    _bond_marks = {x: x for x in _edge_save}
    _atom_container = namedtuple('Atom', __tmp3)
    _bond_container = namedtuple('Bond', _edge_save)
    __implicit_container = namedtuple('ImplicitH', ('s_implicit', 'p_implicit'))
    __visible = None


class ReactionContainer(object):
    __slots__ = ('__reagents', '__products', '__reactants', '__meta')

    def __init__(self, substrats=None, products=None, reactants=None, meta=None, reagents=None):
        if substrats:
            warn('deprecated key. use reagents instead', DeprecationWarning)

        self.__reagents = substrats or reagents or []
        self.__products = products or []
        self.__reactants = reactants or []
        self.__meta = meta or {}

    def __getitem__(self, item):
        if item == 'substrats':
            warn('deprecated key. use reagents instead', DeprecationWarning)
            return self.__reagents
        elif item == 'reagents':
            return self.__reagents
        elif item == 'products':
            return self.__products
        elif item == 'reactants':
            return self.__reactants
        elif item == 'meta':
            return self.__meta
        else:
            raise Exception('Invalid key')

    def pickle(self):
        """ return json serializable reaction
        """
        return dict(reagents=[x.pickle() for x in self.__reagents], meta=self.meta,
                    products=[x.pickle() for x in self.__products], reactants=[x.pickle() for x in self.__reactants])

    @staticmethod
    def unpickle(data):
        """ convert json serializable reaction into ReactionContainer object instance 
        """
        return ReactionContainer(reagents=[MoleculeContainer.unpickle(x) for x in
                                           (data['reagents'] if 'reagents' in data else data['substrats'])],
                                 products=[MoleculeContainer.unpickle(x) for x in data['products']], meta=data['meta'],
                                 reactants=[MoleculeContainer.unpickle(x) for x in data.get('reactants', [])])

    @property
    def substrats(self):  # reverse compatibility
        warn('deprecated key. use reagents instead', DeprecationWarning)
        return self.__reagents

    @property
    def reagents(self):
        return self.__reagents

    @property
    def products(self):
        return self.__products

    @property
    def reactants(self):
        return self.__reactants

    @property
    def meta(self):
        return self.__meta

    def copy(self):
        return ReactionContainer(reagents=[x.copy() for x in self.__reagents], meta=self.__meta.copy(),
                                 products=[x.copy() for x in self.__products],
                                 reactants=[x.copy() for x in self.__reactants])
