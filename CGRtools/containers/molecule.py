# -*- coding: utf-8 -*-
#
#  Copyright 2017, 2018 Ramil Nugmanov <stsouko@live.ru>
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
from collections import namedtuple
from itertools import chain
from .common import BaseContainer
from ..algorithms import CGRstring, Valence, pyramid_volume
from ..exceptions import InvalidData, InvalidAtom, InvalidStereo, ValenceError
from ..periodictable import elements


class MoleculeContainer(BaseContainer, Valence):
    def __init__(self, *args, **kwargs):
        """
        storage for Molecules

        new molecule object creation or copy of data object
        """
        BaseContainer.__init__(self, *args, **kwargs)
        Valence.__init__(self)

    def __dir__(self):
        if self.__visible is None:
            self.__visible = tmp = super().__dir__()
            tmp.extend([self.explicify_hydrogens.__name__, self.implicify_hydrogens.__name__,
                        self.atom_implicit_h.__name__, self.reset_query_marks.__name__])
        return self.__visible

    def __getstate__(self):
        return {k: v for k, v in super().__getstate__().items()
                if not k.startswith(('_Valence__', '_MoleculeContainer__'))}

    def __setstate__(self, state):
        Valence.__init__(self)
        super().__setstate__(state)

    def pickle(self):
        """return json serializable Molecule"""
        data = super().pickle()
        data['s_only'] = True
        return data

    @classmethod
    def unpickle(cls, data):
        """ convert json serializable CGR into MoleculeContainer object instance
        """
        if not data['s_only']:
            raise InvalidData('pickled data is invalid molecule. try CGRContainer.unpickle')

        graph, meta = super().unpickle(data)
        g = cls(graph, meta)
        g.fix_data()
        return g

    def substructure(self, *args, **kwargs):
        sub = super().substructure(*args, **kwargs)
        sub._fix_stereo()
        return sub

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

        self.flush_cache()
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
        self.flush_cache()

    def delete_atom(self, n):
        """
        implementation of atom removing
        """
        # todo: stereo fixing
        self.remove_node(n)
        self.reset_query_marks()

    def delete_bond(self, n, m):
        """
        implementation of bond removing
        """
        # todo: stereo fixing
        self.remove_edge(n, m)
        self.reset_query_marks()

    def add_stereo(self, atom1, atom2, mark):
        if mark not in (1, -1):
            raise InvalidData('stereo mark invalid')
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
            self._tetrahedron_parse(atom1, atom2, mark, neighbors, bonds, implicit)
        else:
            raise InvalidStereo('unsupported stereo or stereo impossible. tetrahedron only supported')

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

    def reset_query_marks(self, copy=False):
        """
        set or reset hyb and neighbors marks to atoms.

        :param copy: if True return copy of graph and keep existing as is
        :return: graph if copy True else None
        """
        g = self.copy() if copy else self
        b, h, n = 's_bond', 's_hyb', 's_neighbors'
        for i, attr in g.nodes(data=True):
            label = dict(s_hyb=1, s_neighbors=0)
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

        if copy:
            return g
        self.flush_cache()

    def implicify_hydrogens(self):
        """
        remove explicit hydrogen if possible

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
        self.flush_cache()
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

        self.flush_cache()
        return len(tmp)

    def atom_implicit_h(self, atom):
        attr = self.nodes[atom]
        return self._get_implicit_h(attr['element'], attr['s_charge'], [x['s_bond'] for x in self[atom].values()],
                                    radical=self._radical_map[attr.get('s_radical')])

    def _tetrahedron_parse(self, atom1, atom2, mark, neighbors, bonds, implicit, label='s'):
        if any(x != 1 for x in bonds):
            raise InvalidStereo('only single bonded tetrahedron acceptable')

        weights = self.get_morgan(stereo=True, labels=label)
        if len(neighbors) != len(set(weights[x] for x in neighbors)):
            raise InvalidStereo('stereo impossible. neighbors equivalent')

        l_x = '%s_x' % label
        l_y = '%s_y' % label
        l_stereo = '%s_stereo' % label

        order = sorted(neighbors, key=weights.get)
        vol = pyramid_volume(*((y[l_x], y[l_y], 0 if x != atom2 else mark) for x, y in
                               ((x, self.nodes[x]) for x in (chain((atom1,), order) if implicit else order))))
        if not vol:
            raise InvalidStereo('unknown')

        self.nodes[atom1][l_stereo] = vol > 0 and 1 or -1
        self.flush_cache()

    @staticmethod
    def _attr_renew(attr, marks):
        new_attr = {}
        for s in marks:
            ls = attr.get(s)
            if ls is not None:
                new_attr[s] = ls
        return new_attr

    def _signature_generator(self, *args, **kwargs):
        return CGRstring(*args, **kwargs)

    def _fix_stereo(self):
        pass
        """
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

            self._weights = self._signatures = self._pickle = None
        """

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
    __visible = __stereo_cache = None
