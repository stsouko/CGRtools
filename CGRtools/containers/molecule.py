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
from collections import defaultdict, MutableMapping
from wrapt import ObjectProxy
from .common import BaseContainer
from ..algorithms import CGRstring, pyramid_volume, aromatize
from ..algorithms.morgan import initial_weights
from ..exceptions import InvalidData, InvalidAtom, InvalidStereo
from ..periodictable import Element, elements_classes, H


class AtomProxy(ObjectProxy):
    __slots__ = ('_self_nodes', '_self_node')

    def __init__(self, nodes, node, atom):
        super().__init__(atom)
        self._self_nodes = nodes
        self._self_node = node

    def copy(self):
        atom = self.__wrapped__
        return type(atom)(charge=atom.charge, multiplicity=atom.multiplicity, isotope=atom.isotope,
                          mapping=atom.mapping, mark=atom.mark, x=atom.x, y=atom.y, z=atom.z, stereo=atom.stereo)

    def update(self, value):
        nodes = self._self_nodes
        node = self._self_node
        if node not in nodes:
            raise KeyError('atom removed from molecule')
        atom = self.__wrapped__

        if isinstance(value, dict):
            if not value:  # ad-hoc for add_nodes_from method
                return
            if set(value) > self.__possible:
                raise ValueError('unknown atom attributes not allowed')

            if 'charge' in value or 'multiplicity' in value or 'isotope' in value or 'element' in value:
                old = dict(charge=atom.charge, multiplicity=atom.multiplicity, isotope=atom.isotope,
                           mapping=atom.mapping, mark=atom.mark, x=atom.x, y=atom.y, z=atom.z, stereo=atom.stereo)
                old.update(value)
                if 'element' in value:
                    del old['element']
                    if 'isotope' not in value:
                        del old['isotope']
                    value = elements_classes[value['element']](**old)
                else:
                    value = elements_classes[atom.symbol](**old)
                self.__wrapped__ = nodes[node] = value
            else:
                for k, v in value.items():
                    setattr(atom, k, v)
        elif isinstance(value, type):
            if not issubclass(value, Element):
                ValueError('only CGRtools.periodictable.Element subclasses allowed')
            if atom.number != value.number:
                self.__wrapped__ = nodes[node] = value(charge=atom.charge, multiplicity=atom.multiplicity,
                                                       mapping=atom.mapping, mark=atom.mark, stereo=atom.stereo,
                                                       x=atom.x, y=atom.y, z=atom.z)
        elif isinstance(value, Element):
            if atom.number != value.number or atom.charge != value.charge or \
                    atom.multiplicity != value.multiplicity or atom.isotope != value.isotope:
                value.neighbors = atom.neighbors
                value.hybridization = atom.hybridization
                self.__wrapped__ = nodes[node] = value
            else:
                atom.x = value.x
                atom.y = value.y
                atom.z = value.z
                atom.mapping = value.mapping
                atom.mark = value.mark
                atom.stereo = value.stereo
        else:
            raise ValueError('only CGRtools.periodictable.Element or dict of atom attributes allowed')

    __possible = {'element', 'charge', 'isotope', 'multiplicity', 'mapping', 'mark', 'x', 'y', 'z', 'stereo'}


class Atoms(MutableMapping):
    __slots__ = '_mapping'

    def __init__(self):
        self._mapping = {}

    def __delitem__(self, node):
        del self._mapping[node]

    def __len__(self):
        return len(self._mapping)

    def __iter__(self):
        return iter(self._mapping)

    def __contains__(self, node):
        return node in self._mapping

    def __setitem__(self, node, value):
        if not isinstance(node, int):  # ad-hoc for add_nodes_from method
            raise TypeError('invalid atom number')
        if isinstance(value, dict):
            if 'atom' not in value:
                if set(value) > self.__possible:
                    raise ValueError('unknown atom attributes not allowed')
                try:
                    value = elements_classes[value.get('element', 'A')](**{k: v for k, v in value.items()
                                                                           if k in self.__acceptable})
                except KeyError:
                    raise ValueError('invalid atom symbol')
                self._mapping[node] = value
                return
            value = value['atom']  # ad-hoc for add_node method

        if isinstance(value, type):
            if not issubclass(value, Element):
                ValueError('invalid type of atom')
            value = value()
        elif not isinstance(value, Element):
            raise ValueError('only CGRtools.periodictable.Element allowed')
        self._mapping[node] = value

    def __getitem__(self, node):
        try:
            atom = self._mapping[node]
        except KeyError:
            raise KeyError('atom not found')
        return AtomProxy(self._mapping, node, atom)

    def setdefault(self, node, value):
        if node not in self._mapping:
            self[node] = value
        return self[node]

    __acceptable = {'charge', 'isotope', 'multiplicity', 'mapping', 'mark', 'x', 'y', 'z', 'stereo'}
    __possible = __acceptable | {'element', 'hybridization', 'neighbors'}


class Bond:
    __slots__ = ('order', 'stereo')

    def __init__(self, order=1, stereo=None):
        self.order = order
        self.stereo = stereo

    def __setattr__(self, key, value):
        if key not in self.__possible:
            raise AttributeError('unknown bond attribute')
        acc = self.__acceptable[key]
        if value not in acc:
            raise ValueError(f'attribute {key} should be from acceptable list: {acc}')
        super().__setattr__(key, value)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __getitem__(self, key):
        return getattr(self, key)

    def __len__(self):
        return 2

    def __iter__(self):
        return iter(self.__possible)

    def __repr__(self):
        s = '' if self.stereo is None else f', stereo={self.stereo}'
        return f'{type(self).__name__}({self.order}{s})'

    def update(self, value):
        if isinstance(value, dict):
            if not value:
                return  # ad-hoc for add_edges_from method
            if 'bond' not in value:
                if set(value) > self.__possible:
                    raise ValueError('unknown bond attributes not allowed')
                for k, v in self.__acceptable.items():
                    if k in value and value[k] not in v:
                        raise ValueError(f'attribute {k} should be from acceptable list: {v}')
                for k, v in value.items():
                    super().__setattr__(k, v)
                return
            value = value['bond']
        if isinstance(value, Bond):
            self.order = value.order
            self.stereo = value.stereo
        else:
            raise TypeError('only Bond object or dict of bond attributes allowed')

    def copy(self):
        return type(self)(self.order, self.stereo)

    def items(self):
        """
        iterate other non-default bonds attrs-values pairs

        need for nx.readwrite.json_graph.node_link_data
        """
        def g():
            for k, d in self.__default.items():
                v = getattr(self, k)
                if v != d:
                    yield (k, v)
        return g()

    __possible = {'order', 'stereo'}
    __acceptable = {'order': {1, 2, 3, 4, 9}, 'stereo': {None, -1, 1, 0}}
    __default = {'order': 1, 'stereo': None}


class MoleculeContainer(BaseContainer):
    """
    storage for Molecules

    new molecule object creation or copy of data object
    """
    node_dict_factory = Atoms
    edge_attr_dict_factory = Bond

    def __dir__(self):
        if self.__visible is None:
            self.__visible = super().__dir__() + [self.explicify_hydrogens.__name__, self.implicify_hydrogens.__name__,
                                                  self.atom_implicit_h.__name__, self.reset_query_marks.__name__,
                                                  self.aromatize.__name__, self.check_valence.__name__]
        return self.__visible

    def substructure(self, *args, **kwargs):
        sub = super().substructure(*args, **kwargs)
        sub._fix_stereo_stage_2(self._fix_stereo_stage_1())
        return sub

    def add_bond(self, atom1, atom2, bond):
        stereo = self._fix_stereo_stage_1()
        super().add_bond(atom1, atom2, bond)
        self._fix_stereo_stage_2(stereo)

    def delete_atom(self, n):
        """
        implementation of atom removing
        """
        stereo = self._fix_stereo_stage_1()
        self.remove_node(n)
        self._fix_stereo_stage_2(stereo)

    def delete_bond(self, n, m):
        """
        implementation of bond removing
        """
        stereo = self._fix_stereo_stage_1()
        self.remove_edge(n, m)
        self._fix_stereo_stage_2(stereo)

    def add_stereo(self, atom1, atom2, mark):
        if mark not in (1, -1):
            raise InvalidData('stereo mark invalid')
        if not self.has_edge(atom1, atom2):
            raise InvalidAtom('atom or bond not found')

        if self.nodes[atom1].stereo:
            raise self._stereo_exception3

        if self._node[atom1].z or any(self._node[x].z for x in self._adj[atom1]):
            raise self._stereo_exception1

        implicit = self.atom_implicit_h(atom1)
        if implicit + self.atom_explicit_h(atom1) > 1:
            raise self._stereo_exception4

        total = implicit + len(self._adj[atom1])
        if total == 4 and all(x.order == 1 for x in self._adj[atom1].values()):  # tetrahedron
            self.__tetrahedron_parse(atom1, atom2, mark, implicit)
        else:
            raise self._stereo_exception2

    def reset_query_marks(self, copy=False):
        """
        set or reset hyb and neighbors marks to atoms.

        :param copy: if True return copy of graph and keep existing as is
        :return: graph if copy True else None
        """
        g = self.copy() if copy else self
        for i, atom in g._node.items():
            neighbors = 0
            hybridization = 1
            # hybridization 1- sp3; 2- sp2; 3- sp1; 4- aromatic
            for j, bond in g._adj[i].items():
                if g._node[j] != 'H':
                    neighbors += 1

                order = bond.order
                if order == 1 or hybridization in (3, 4):
                    continue
                elif order == 4:
                    hybridization = 4
                elif order == 3 or (order == 2 and hybridization == 2):  # Если есть 3-я или две 2-х связи, то sp1
                    hybridization = 3
                elif order == 2 and hybridization not in (3, 4):
                    # Если есть 2-я связь, но до этого не было найдено другой 2-й, 3-й, или аром.
                    hybridization = 2

            atom.neighbors = neighbors
            atom.hybridization = hybridization

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
        for n, atom in self._node.items():
            if atom == 'H':
                m = next(self.neighbors(n))
                if self._node[m] != 'H':
                    explicit[m].append(n)

        for n, h in explicit.items():
            atom = self._node[n]
            len_h = len(h)
            for i in range(len_h, 0, -1):
                hi = h[:i]
                if atom.get_implicit_h([y.order for x, y in self._adj[n].items() if x not in hi]) == i:
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
        for n, atom in self._node.items():
            if atom != 'H':
                for _ in range(atom.get_implicit_h([x.order for x in self._adj[n].values()])):
                    tmp.append(n)
        for n in tmp:
            self.add_bond(n, self.add_atom(H), Bond())

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
        return self._node[atom].get_implicit_h([x.order for x in self._adj[atom].values()])

    def atom_explicit_h(self, atom):
        return sum(self._node[x] == 'H' for x in self.neighbors(atom))

    def atom_total_h(self, atom):
        return self.atom_explicit_h(atom) + self.atom_implicit_h(atom)

    def check_valence(self):
        """
        check valences of all atoms

        :return: list of invalid atoms
        """
        return [f'atom {x} has invalid valence' for x, atom in self._node.items()
                if not atom.check_valence(self.environment(x))]

    def __tetrahedron_parse(self, atom1, atom2, mark, implicit):
        weights = self.get_morgan(stereo=True)

        neighbors = list(self._adj[atom1])
        if len(neighbors) != len(set(weights[x] for x in neighbors)):
            raise InvalidStereo('stereo impossible. neighbors equivalent')

        order = sorted(((x, self._node[x]) for x in neighbors), key=lambda x: weights[x[0]])
        if implicit:
            central = self._node[atom1]
            vol = pyramid_volume((central.x, central.y, 0), *((atom.x, atom.y, 0 if x != atom2 else mark)
                                                              for x, atom in order))
        else:
            vol = pyramid_volume(*((atom.x, atom.y, 0 if x != atom2 else mark) for x, atom in order))

        if not vol:
            raise InvalidStereo('unknown')

        self._node[atom1].stereo = vol > 0 and 1 or -1
        self.flush_cache()

    def _prepare_stereo(self):
        _stereo_cache = {}
        nodes = list(self._node.items())
        while True:
            failed = []
            for i, tmp in enumerate(nodes, start=1):
                n, atom = tmp
                s = atom.stereo
                if not s:
                    continue
                neighbors = list(self._adj[n])
                len_n = len(neighbors)
                if len_n in (3, 4):  # tetrahedron
                    if atom.z or any(self._node[x].z for x in neighbors):
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
                        zero = self._node[order[0]]
                        zero = (zero.x, zero.y, 1)
                        first = self._node[order[1]]
                        first = (first.x, first.y, 0)
                    else:
                        zero = (atom.x, atom.y, 0)
                        first = self._node[order[0]]
                        first = (first.x, first.y, 1)

                    second = self._node[order[-2]]
                    third = self._node[order[-1]]
                    vol = pyramid_volume(zero, first, (second.x, second.y, 0), (third.x, third.y, 0))

                    _stereo_cache[(n, order[0])] = 1 if vol > 0 and s == 1 or vol < 0 and s == -1 else -1
            else:
                return _stereo_cache

    def _signature_generator(self, *args, **kwargs):
        return CGRstring(*args, **kwargs)

    @property
    def _morgan_init(self):
        return initial_weights

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

    __visible = None
    _stereo_exception1 = InvalidStereo('molecule have 3d coordinates. bond up/down stereo unusable')
    _stereo_exception2 = InvalidStereo('unsupported stereo or stereo impossible. '
                                       'single bonded tetrahedron only supported')
    _stereo_exception3 = InvalidStereo('atom has stereo. change impossible')
    _stereo_exception4 = InvalidStereo('stereo impossible. too many H atoms')
