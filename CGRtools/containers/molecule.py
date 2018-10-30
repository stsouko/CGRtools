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
from collections.abc import MutableMapping
from .common import BaseContainer
from ..algorithms import CGRstring, pyramid_volume, aromatize
from ..algorithms.morgan import initial_weights
from ..exceptions import InvalidData, InvalidAtom, InvalidStereo
from ..periodictable import Element, elements_classes, H


class Atom(MutableMapping):
    __slots__ = '__atom'

    def __init__(self):
        super().__setattr__('_Atom__atom', None)

    def __getitem__(self, key):
        """
        dict like access to atom's attrs
        """
        if key == 'element':
            key = 'symbol'
        return getattr(self.__atom, key)

    def __setitem__(self, key, value):
        return setattr(self, key, value)

    def __getattr__(self, key):
        return getattr(self.__atom, key)

    def __setattr__(self, key, value):
        if key in self.__mutable:
            setattr(self.__atom, key, value)
        else:
            if key not in self.__possible:
                raise ValueError('unknown atom attributes not allowed')
            attrs = {k: getattr(self.__atom, k) for k in self.__acceptable}
            if key == 'element':
                if value == self.__atom.symbol:
                    return
                if value not in elements_classes or value == 'A':
                    raise ValueError('invalid atom symbol')

                del attrs['isotope']
                super().__setattr__('_Atom__atom', elements_classes[value](**attrs))
            else:
                attrs[key] = value
                super().__setattr__('_Atom__atom', elements_classes[self.__atom.symbol](**attrs))

    def __len__(self):
        return sum(1 for _ in self)

    def __iter__(self):
        """
        iterate other non-default atom's init attrs
        """
        if self.__atom is None:
            raise StopIteration
        yield 'element'
        if self.__atom.isotope != self.__atom.common_isotope:
            yield 'isotope'
        for k, d in self.__defaults.items():
            if d != getattr(self, k):
                yield k

    def __contains__(self, key):
        return key in self.__possible

    def __delitem__(self, key):
        raise NotImplemented('attribute deletion impossible')

    def __eq__(self, other):
        return self.__atom == other

    def __gt__(self, other):
        return self.__atom > other

    def __ge__(self, other):
        return self.__atom >= other

    def __lt__(self, other):
        return self.__atom < other

    def __le__(self, other):
        return self.__atom <= other

    def __repr__(self):
        return repr(self.__atom)

    def __hash__(self):
        return hash(self.__atom)

    def copy(self):
        atom = self.__atom
        return type(atom)(charge=atom.charge, multiplicity=atom.multiplicity, isotope=atom.isotope,
                          mapping=atom.mapping, mark=atom.mark, x=atom.x, y=atom.y, z=atom.z, stereo=atom.stereo)

    def update(self, *args, **kwargs):
        """
        update atom

        :param args: tuple with 1 or 0 elements. element can be dict of atom attrs or atom object or atom class or
        dict with 'atom' key with value equal to dict of atom attrs or atom object or atom class.
        atom keyed dict also may contain atom attrs. this attrs has precedence other 'atom' key value attrs.
        :param kwargs: atom attrs. has precedence other args[0]
        """
        if args:
            if len(args) > 1:
                raise TypeError('update expected at most 1 arguments')

            value = args[0]
            if isinstance(value, dict) and 'atom' in value:
                value.update(kwargs)
                kwargs = value
                value = value.pop('atom')
        elif 'atom' in kwargs:
            value = kwargs.pop('atom')
        else:
            value = kwargs
            kwargs = ()

        if isinstance(value, type):
            if not issubclass(value, Element):
                ValueError('only CGRtools.periodictable.Element subclasses allowed')
            if kwargs.keys() - self.__mutable.keys():
                raise KeyError('unmutable atom attributes not allowed')

            if self.__atom is None:
                super().__setattr__('_Atom__atom', value(**kwargs))
            else:
                attrs = {k: getattr(self.__atom, k) for k in self.__mutable}
                attrs.update(kwargs)
                super().__setattr__('_Atom__atom', value(**attrs))

        elif isinstance(value, Element):
            try:
                if not all(self.__mutable[k](v) for k, v in kwargs.items()):
                    raise ValueError('invalid value of mutable attribute')
            except KeyError:
                raise KeyError('unmutable atom attributes not allowed')

            attrs = {k: getattr(value, k) for k in self.__mutable}
            attrs.update(kwargs)
            if self.__atom is None or self.__atom != value:
                super().__setattr__('_Atom__atom', type(value)(charge=value.charge, isotope=value.isotope,
                                                               multiplicity=value.multiplicity, **attrs))
            else:
                for k, v in attrs.items():
                    setattr(self.__atom, k, v)
        else:
            if not isinstance(value, dict):
                try:
                    value = dict(value)
                except (TypeError, ValueError):
                    raise TypeError('invalid attrs sequence')

            value.update(kwargs)
            if not value:  # ad-hoc for add_nodes_from method
                return
            if value.keys() - self.__possible:
                raise ValueError('unknown atom attributes not allowed')

            if self.__atom is None:
                if 'element' not in value or value['element'] not in elements_classes or value['element'] == 'A':
                    raise ValueError('invalid atom symbol')

                super().__setattr__('_Atom__atom', elements_classes[value['element']](**{k: v for k, v in value.items()
                                                                                         if k != 'element'}))
            elif value.keys() - self.__mutable.keys():
                attrs = {k: getattr(self.__atom, k) for k in self.__acceptable}
                attrs.update(value)

                if 'element' in value:
                    if value['element'] not in elements_classes or value['element'] == 'A':
                        raise ValueError('invalid atom symbol')

                    del attrs['element']
                    if 'isotope' not in value:
                        del attrs['isotope']
                    super().__setattr__('_Atom__atom', elements_classes[value['element']](**attrs))
                else:
                    super().__setattr__('_Atom__atom', elements_classes[self.__atom.symbol](**attrs))
            else:
                if not all(self.__mutable[k](v) for k, v in value.items()):
                    raise ValueError('invalid value of mutable attribute')
                for k, v in value.items():
                    setattr(self.__atom, k, v)

    @staticmethod
    def __int_float(x):
        return isinstance(x, (float, int))

    __defaults = {'mapping': None, 'mark': '0', 'x': 0, 'y': 0, 'z': 0, 'stereo': None, 'charge': 0,
                  'multiplicity': None}
    __mutable = {'mapping': lambda x: x is None or isinstance(x, int), 'mark': lambda x: isinstance(x, str),
                 'x': __int_float, 'y': __int_float, 'z': __int_float, 'stereo': lambda x: x in {None, -1, 1, 0}}
    __acceptable = {'isotope', *__mutable}
    __possible = {'element', *__acceptable}


class Bond(MutableMapping):
    __slots__ = ('order', 'stereo')

    def __init__(self, order=1, stereo=None):
        self.order = order
        self.stereo = stereo

    def __getitem__(self, key):
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __setattr__(self, key, value):
        if key not in self.__acceptable:
            raise AttributeError('unknown bond attribute')
        if value not in self.__acceptable[key]:
            raise ValueError(f"attribute '{key}' value should be from acceptable list: {self.__acceptable[key]}")
        super().__setattr__(key, value)

    def __len__(self):
        return sum(1 for _ in self)

    def __iter__(self):
        """
        iterate other non-default bonds attrs

        need for nx.readwrite.json_graph.node_link_data
        """
        if self.order != self.__defaults['order']:
            yield 'order'
        if self.stereo != self.__defaults['stereo']:
            yield 'stereo'

    def __contains__(self, key):
        return key in self.__acceptable

    def __delitem__(self, key):
        raise NotImplemented('attribute deletion impossible')

    def __repr__(self):
        s = '' if self.stereo is None else f', stereo={self.stereo}'
        return f'{type(self).__name__}({self.order}{s})'

    def update(self, *args, **kwargs):
        """
        update bond
        :param args: tuple with 1 or 0 elements. element can be Bond object or dict of bond attrs or dict with 'bond'
        key with value equal to Bond object or dict of bond attrs.
        bond keyed dict also may contain bond attrs. this attrs has precedence other 'bond' key value attrs.
        :param kwargs: bond attrs. has precedence other args[0]
        """
        if args:
            if len(args) > 1:
                raise TypeError('update expected at most 1 arguments')

            value = args[0]
            if isinstance(value, dict) and 'bond' in value:
                value.update(kwargs)
                kwargs = value
                value = value.pop('bond')
        elif 'bond' in kwargs:
            value = kwargs.pop('bond')
        else:
            value = kwargs
            kwargs = ()

        if isinstance(value, Bond):
            try:
                if not all(v in self.__acceptable[k] for k, v in kwargs.items()):
                    raise ValueError('invalid attribute value')
            except KeyError:
                raise KeyError('unknown bond attributes not allowed')

            super().__setattr__('order', value.order)
            super().__setattr__('stereo', value.stereo)
            for k, v in kwargs.items():
                super().__setattr__(k, v)

        else:
            if not isinstance(value, dict):
                try:
                    value = dict(value)
                except (TypeError, ValueError):
                    raise TypeError('invalid attrs sequence')

            value.update(kwargs)
            if not value:
                return  # ad-hoc for add_edges_from method

            try:
                if not all(v in self.__acceptable[k] for k, v in kwargs.items()):
                    raise ValueError('invalid attribute value')
            except KeyError:
                raise KeyError('unknown bond attributes not allowed')

            for k, v in value.items():
                super().__setattr__(k, v)

    def copy(self):
        return type(self)(self.order, self.stereo)

    __acceptable = {'order': {1, 2, 3, 4, 9}, 'stereo': {None, -1, 1, 0}}
    __defaults = {'order': 1, 'stereo': None}


class MoleculeContainer(BaseContainer):
    """
    storage for Molecules

    new molecule object creation or copy of data object
    """
    node_attr_dict_factory = Atom
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
