# -*- coding: utf-8 -*-
#
#  Copyright 2017, 2018 Ramil Nugmanov <stsouko@live.ru>
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
from collections import defaultdict
from itertools import chain
from .common import BaseContainer
from ..algorithms import CGRstring, pyramid_volume, aromatize
from ..exceptions import InvalidData, InvalidAtom, InvalidStereo
from ..periodictable import Element, elements, H


class MoleculeContainer(BaseContainer):
    """
    storage for Molecules

    new molecule object creation or copy of data object
    """
    def __dir__(self):
        if self.__visible is None:
            self.__visible = super().__dir__() + [self.explicify_hydrogens.__name__, self.implicify_hydrogens.__name__,
                                                  self.atom_implicit_h.__name__, self.reset_query_marks.__name__,
                                                  self.aromatize.__name__, self.check_valence.__name__]
        return self.__visible

    def __getstate__(self):
        return {k: v for k, v in super().__getstate__().items() if not k.startswith('_MoleculeContainer__')}

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
        sub._fix_stereo_stage_2(self._fix_stereo_stage_1())
        return sub

    def add_atom(self, atom, _map=None, mark='0', x=0, y=0, z=0):
        if isinstance(atom, type):
            if not issubclass(atom, Element):
                TypeError('invalid type of atom')
            atom = atom()
        elif isinstance(atom, str):
            try:
                atom = elements[atom]()
            except KeyError:
                raise TypeError('invalid symbol of atom')
        elif not isinstance(atom, Element):
            raise TypeError('invalid type of atom')
        if _map is None:
            _map = max(self, default=0) + 1
        elif _map in self:
            raise InvalidData('mapping exists')

        self.add_node(_map, element=atom.symbol, s_charge=atom.charge, mark=mark, s_x=x, s_y=y, s_z=z, map=_map)
        if atom.multiplicity:
            self.nodes[_map]['s_radical'] = atom.multiplicity

        self.flush_cache()
        return _map

    def add_bond(self, atom1, atom2, mark=1, *, ignore=False):
        if atom1 not in self or atom2 not in self:
            raise InvalidData('atoms not found')
        if self.has_edge(atom1, atom2):
            raise InvalidData('bond exists')
        if mark not in (1, 2, 3, 4, 9):
            raise InvalidData('invalid bond mark')

        if not ignore:
            self.__check_bonding(atom1, atom2, mark)

        stereo = self._fix_stereo_stage_1()
        self.add_edge(atom1, atom2, s_bond=mark)
        self._fix_stereo_stage_2(stereo)
        self.flush_cache()

    def delete_atom(self, n):
        """
        implementation of atom removing
        """
        stereo = self._fix_stereo_stage_1()
        self.remove_node(n)
        self._fix_stereo_stage_2(stereo)
        self.reset_query_marks()

    def delete_bond(self, n, m):
        """
        implementation of bond removing
        """
        stereo = self._fix_stereo_stage_1()
        self.remove_edge(n, m)
        self._fix_stereo_stage_2(stereo)
        self.reset_query_marks()

    def add_stereo(self, atom1, atom2, mark):
        if mark not in (1, -1):
            raise InvalidData('stereo mark invalid')
        if not self.has_edge(atom1, atom2):
            raise InvalidAtom('atom or bond not found')

        if self.nodes[atom1].get('s_stereo'):
            raise self._stereo_exception3

        tmp = [(x, y['s_bond']) for x, y in self[atom1].items()]
        neighbors = [x for x, _ in tmp]

        if self.nodes[atom1]['s_z'] or any(self.nodes[x]['s_z'] for x in neighbors):
            raise self._stereo_exception1

        neighbors_e = [self.nodes[x]['element'] for x in neighbors]
        implicit = self.atom_implicit_h(atom1)
        if implicit > 1 or implicit == 1 and 'H' in neighbors_e or neighbors_e.count('H') > 1:
            raise self._stereo_exception4

        bonds = [x for _, x in tmp]
        total = implicit + len(neighbors)
        if total == 4:  # tetrahedron
            self._tetrahedron_parse(atom1, atom2, mark, neighbors, bonds, implicit)
        else:
            raise self._stereo_exception2

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
        explicit = defaultdict(list)
        c = 0
        for n, attr in self.nodes(data='element'):
            if attr == 'H':
                m = next(self.neighbors(n))
                if self.nodes[m]['element'] != 'H':
                    explicit[m].append(n)

        for n, h in explicit.items():
            atom = self.atom(n)
            len_h = len(h)
            for i in range(len_h, 0, -1):
                hi = h[:i]
                if atom.get_implicit_h([y['s_bond'] for x, y in self[n].items() if x not in hi]) == i:
                    for x in hi:
                        self.remove_node(x)
                        c += 1
                    break

        self.flush_cache()
        return c

    def explicify_hydrogens(self):
        """
        add explicit hydrogens to atoms

        :return: number of added atoms
        """
        tmp = []
        for n, attr in self.nodes(data='element'):
            if attr != 'H':
                for _ in range(self.atom_implicit_h(n)):
                    tmp.append(n)
        for n in tmp:
            self.add_bond(n, self.add_atom(H()), 1)

        self.flush_cache()
        return len(tmp)

    def aromatize(self) -> int:
        """
        convert structure to aromatic form

        :return: number of processed rings
        """
        res = aromatize(self)
        if res:
            self.flush_cache()
        return res

    def atom_implicit_h(self, atom):
        return self.atom(atom).get_implicit_h([x for *_, x in self.edges(atom, data='s_bond')])

    def atom_explicit_h(self, atom):
        return sum(self.nodes[x]['element'] == 'H' for x in self.neighbors(atom))

    def atom_total_h(self, atom):
        return self.atom_explicit_h(atom) + self.atom_implicit_h(atom)

    def check_valence(self):
        """
        check valences of all atoms

        :return: list of invalid atoms
        """
        errors = []
        for x in self.nodes():
            try:
                atom = self.atom(x)
            except InvalidAtom as e:
                errors.append('atom %d has error: %s' % (x, e))
            else:
                if atom.get_valence(*self._get_atom_environment(x)) is None:
                    errors.append('atom %d has invalid valence' % x)

        return errors

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

    def _prepare_stereo(self):
        _stereo_cache = {}
        nodes = list(self.nodes(data=True))
        while True:
            failed = []
            for i, tmp in enumerate(nodes, start=1):
                n, attr = tmp
                s = attr.get('s_stereo')
                if not s:
                    continue
                neighbors = list(self.neighbors(n))
                len_n = len(neighbors)
                if len_n in (3, 4):  # tetrahedron
                    if attr['s_z'] or any(self.nodes[x]['s_z'] for x in neighbors):
                        continue  # 3d molecules ignored

                    weights = self.get_morgan(stereo=True)
                    order = sorted(neighbors, key=weights.get)
                    for _ in range(len_n):
                        if (order[0], n) in _stereo_cache:
                            order.append(order.pop(0))
                        else:
                            failed.append(tmp)
                            break
                    else:
                        failed.insert(0, tmp)
                        failed.extend(nodes[i:])
                        _stereo_cache = None
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

                    _stereo_cache[(n, order[0])] = 1 if vol > 0 and s == 1 or vol < 0 and s == -1 else -1
            else:
                return _stereo_cache

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

    def _fix_stereo_stage_1(self):
        tmp = []
        for n, m in self.edges():
            s = self.get_stereo(n, m)
            if s:
                tmp.append((n, m, s))
            else:
                s = self.get_stereo(m, n)
                if s:
                    tmp.append((m, n, s))
        return tmp

    def _fix_stereo_stage_2(self, stereo):
        while True:
            failed_stereo = []
            for n, m, mark in stereo:
                try:
                    self.add_stereo(n, m, mark)
                except InvalidStereo:
                    failed_stereo.append((n, m, mark))
                except InvalidAtom:
                    continue
            if failed_stereo and len(stereo) > len(failed_stereo):
                stereo = failed_stereo
                continue
            break

    def __check_bonding(self, n, m, mark):
        for atom, reverse in ((n, m), (m, n)):
            bn, ng = self._get_atom_environment(atom)
            ng.append(self.nodes[reverse]['element'])
            bn.append(mark)

            if not self.atom(atom).get_valence(bn, ng):
                raise InvalidData('valence error')

    def _get_atom_environment(self, atom, label='s'):
        ng, bn = [], []
        label = '%s_bond' % label
        for a, attrs in self[atom].items():
            b = attrs.get(label)
            if b:
                bn.append(b)
                ng.append(self.nodes[a]['element'])
        return bn, ng

    @classmethod
    def _atom_container(cls, attrs):
        return elements[attrs['element']](charge=attrs['s_charge'], multiplicity=attrs.get('s_radical'),
                                          isotope=attrs.get('isotope'))

    @classmethod
    def _stereo_container(cls, attrs):
        return attrs.get('s_stereo')

    @classmethod
    def _bond_container(cls, attrs):
        return attrs['s_bond']

    _node_base = ('element', 'isotope', 'mark', 's_x', 's_y', 's_z')
    _node_marks = ('s_neighbors', 's_hyb', 's_charge', 's_stereo', 's_radical')
    _node_save = _node_marks + _node_base
    _edge_save = _edge_marks = ('s_bond', 's_stereo')

    __visible = None
    _stereo_exception1 = InvalidStereo('molecule have 3d coordinates. bond up/down stereo unusable')
    _stereo_exception2 = InvalidStereo('unsupported stereo or stereo impossible. tetrahedron only supported')
    _stereo_exception3 = InvalidStereo('atom has stereo. change impossible')
    _stereo_exception4 = InvalidStereo('stereo impossible. too many H atoms')
