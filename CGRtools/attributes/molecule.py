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
from collections.abc import MutableMapping
from ..periodictable import Element, elements_classes


class Atom(MutableMapping):
    __slots__ = ('_atom', '__hybridization', '__neighbors', '__color', '__stereo', '__mapping', '__mark',
                 '__x', '__y', '__z')

    def __init__(self, atom=None):
        if isinstance(atom, Atom):
            super().__setattr__('_atom', atom.atom)
            super().__setattr__('_Atom__x', atom.x)
            super().__setattr__('_Atom__y', atom.y)
            super().__setattr__('_Atom__z', atom.z)
            super().__setattr__('_Atom__mapping', atom.mapping)
            super().__setattr__('_Atom__stereo', atom.stereo)
            super().__setattr__('_Atom__mark', atom.mark)
            super().__setattr__('_Atom__color', atom.color)
            super().__setattr__('_Atom__hybridization', None)
            super().__setattr__('_Atom__neighbors', None)
        else:
            if atom is None:
                super().__setattr__('_atom', None)
            elif isinstance(atom, Element):
                super().__setattr__('_atom', atom)
            else:
                raise TypeError('invalid atom passed')

            super().__setattr__('_Atom__x', 0)
            super().__setattr__('_Atom__y', 0)
            super().__setattr__('_Atom__z', 0)
            super().__setattr__('_Atom__mapping', None)
            super().__setattr__('_Atom__stereo', None)
            super().__setattr__('_Atom__mark', '0')
            super().__setattr__('_Atom__color', None)
            super().__setattr__('_Atom__hybridization', None)
            super().__setattr__('_Atom__neighbors', None)

    @property
    def x(self):
        return self.__x

    @x.setter
    def x(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError('coordinates can be float')
        super().__setattr__('_Atom__x', float(value))

    @property
    def y(self):
        return self.__y

    @y.setter
    def y(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError('coordinates can be float')
        super().__setattr__('_Atom__y', float(value))

    @property
    def z(self):
        return self.__z

    @z.setter
    def z(self, value):
        if not isinstance(value, (float, int)):
            raise TypeError('coordinates can be float')
        super().__setattr__('_Atom__z', float(value))

    @property
    def mapping(self):
        return self.__mapping

    @mapping.setter
    def mapping(self, value):
        if value is None:
            super().__setattr__('_Atom__mapping', value)
        elif not isinstance(value, int):
            raise TypeError('mapping can be int')
        elif value >= 1000 or value < 0:
            raise ValueError('mapping can be in range 0-999')
        super().__setattr__('_Atom__mapping', value)

    @property
    def stereo(self):
        return self.__stereo

    @stereo.setter
    def stereo(self, value):
        if value is None:
            super().__setattr__('_Atom__stereo', None)
        elif not isinstance(value, int):
            raise TypeError('stereo can be int')
        elif value not in (-1, 1):
            raise ValueError('stereo can be: None, 1 or -1')
        super().__setattr__('_Atom__stereo', value)

    @property
    def neighbors(self):
        return self.__neighbors

    @neighbors.setter
    def neighbors(self, value):
        if value is None:
            super().__setattr__('_Atom__neighbors', None)
        elif not isinstance(value, int):
            raise TypeError('neighbors can be int')
        elif value < 0 or value > 998:
            raise ValueError('neighbors can be: None or in range 0-998')
        super().__setattr__('_Atom__neighbors', value)

    @property
    def hybridization(self):
        return self.__hybridization

    @hybridization.setter
    def hybridization(self, value):
        if value is None:
            super().__setattr__('_Atom__hybridization', None)
        elif not isinstance(value, int):
            raise TypeError('hybridization can be int')
        elif value not in (1, 2, 3, 4):
            raise ValueError('hybridization can be: None, 1, 2, 3 or 4')
        super().__setattr__('_Atom__hybridization', value)

    @property
    def mark(self):
        return self.__mark

    @mark.setter
    def mark(self, value):
        if not isinstance(value, str):
            raise TypeError('mark can be str')
        if not value or len(value) > 3:
            raise ValueError('mark can be 1-3 symbol length')
        super().__setattr__('_Atom__mark', value)

    @property
    def color(self):
        if self.__color:
            return self.__color.copy()

    @color.setter
    def color(self, value):
        if not value:
            super().__setattr__('_Atom__color', None)
        elif not isinstance(value, dict):
            raise TypeError('color can be dict')
        elif not all(isinstance(x, int) for x in value):
            raise TypeError('color keys can be int')
        elif not all(isinstance(x, str) for x in value.values()):
            raise TypeError('color values can be str')
        super().__setattr__('_Atom__color', value.copy())

    def __getitem__(self, key):
        """
        dict like access to atom's attrs
        """
        if key == 'element':
            return self._atom.symbol
        elif key in ('symbol', 'isotope', 'charge', 'multiplicity'):
            return getattr(self._atom, key)
        elif key in ('color', 'stereo', 'mapping', 'mark', 'x', 'y', 'z'):
            return getattr(self, key)
        raise KeyError('unknown atom attribute')

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __getattr__(self, key):
        if key == 'element':
            return self._atom.symbol
        return getattr(self._atom, key)

    def __setattr__(self, key, value):
        if key in ('color', 'stereo', 'mapping', 'mark', 'x', 'y', 'z', 'neighbors', 'hybridization'):
            super().__setattr__(key, value)
        elif key in ('element', 'symbol'):
            if value not in elements_classes or value == 'A':
                raise ValueError('invalid atom symbol')
            if self._atom is None:
                super().__setattr__('_atom', elements_classes[value]())
            elif value != self._atom.symbol:
                super().__setattr__('_atom', elements_classes[value](self._atom.charge, self._atom.multiplicity,
                                                                     self._atom.isotope))
        elif self._atom is None:
            raise TypeError('any atom not allowed')
        elif key in ('charge', 'isotope', 'multiplicity'):
            attrs = {'multiplicity': self._atom.multiplicity, 'isotope': self._atom.isotope,
                     'charge': self._atom.charge, key: value}
            super().__setattr__('_atom', elements_classes[self._atom.symbol](**attrs))
        else:
            raise KeyError('unknown atom attributes not allowed')

    def __len__(self):
        return sum(1 for _ in self)

    def __iter__(self):
        """
        iterate other non-default atom's init attrs
        """
        satom = self._atom
        if satom is None:
            raise StopIteration
        yield 'element'
        if satom.isotope != satom.common_isotope:
            yield 'isotope'
        if satom.charge:
            yield 'charge'
        if satom.multiplicity:
            yield 'multiplicity'
        for k, d in self.__defaults.items():
            if d != getattr(self, k):
                yield k

    def __delitem__(self, key):
        raise TypeError('attribute deletion impossible')

    def __eq__(self, other):
        if isinstance(other, Atom):
            return self._atom == other._atom
        return False

    def __str__(self):
        return self.format(isotope=True, stereo=True)

    def format(self, atom=True, isotope=False, stereo=False, hybridization=False, neighbors=False):
        smi = []

        if stereo and self.__stereo:
            smi.append(self._stereo_str[self.__stereo])
        if hybridization and self.__hybridization:
            smi.append(self._hybridization_str[self.__hybridization])
        if neighbors and self.__neighbors is not None:
            smi.append(str(self.__neighbors))
        if smi:
            smi.append(';')
            smi.insert(0, ';')

        satom = self._atom
        if atom:
            smi.insert(0, satom.symbol)
            if satom.charge:
                smi.append(self._charge_str[satom.charge])
            if satom.multiplicity:
                smi.append(self._multiplicity_str[satom.multiplicity])
        else:
            smi.insert(0, '*')

        if isotope and satom.isotope != satom.common_isotope:
            smi.insert(0, str(satom.isotope))

        if len(smi) != 1 or not atom:
            smi.insert(0, '[')
            smi.append(']')

        return ''.join(smi)

    @property
    def atom(self):
        """
        wrapped atom
        """
        return self._atom

    def copy(self):
        return type(self)(self)

    def update(self, *args, **kwargs):
        """
        update atom

        :param args: tuple with 1 or 0 elements. element can be dict of atom attrs or atom object or atom class.
        :param kwargs: atom attrs. has precedence other args[0]
        """
        if len(args) > 1:
            raise TypeError('update expected at most 1 arguments')
        elif args:
            value = args[0]
        else:
            value = kwargs
            kwargs = ()

        if isinstance(value, type):
            if not issubclass(value, Element):
                ValueError('only CGRtools.periodictable.Element subclasses allowed')
            if not all(self.__mutable[k](v) for k, v in kwargs.items()):
                raise ValueError('invalid value of attribute')
            for k, v in kwargs.items():
                super().__setattr__(k, v)
            super().__setattr__('_atom', value())

        elif isinstance(value, Element):
            if not all(self.__mutable[k](v) for k, v in kwargs.items()):
                raise ValueError('invalid value of attribute')
            for k, v in kwargs.items():
                super().__setattr__(k, v)
            super().__setattr__('_atom', value)
        else:
            if not isinstance(value, dict):
                try:
                    value = dict(value)
                except (TypeError, ValueError):
                    raise TypeError('invalid attrs sequence')

            value.update(kwargs)
            if not value:  # ad-hoc for add_nodes_from method
                return

            if self._atom is None and 'element' not in value:
                raise TypeError('any atom not allowed')

            mutable_keys = value.keys() & self.__mutable.keys()
            if not all(self.__mutable[k](value[k]) for k in mutable_keys):
                raise ValueError('invalid value of attribute')

            atom_keys = value.keys() - self.__mutable.keys()
            if atom_keys:
                attrs = {'multiplicity': self._atom.multiplicity, 'isotope': self._atom.isotope,
                         'charge': self._atom.charge, **{k: value[k] for k in atom_keys}}

                if 'element' in value:
                    if value['element'] not in elements_classes or value['element'] == 'A':
                        raise ValueError('invalid atom symbol')

                    del attrs['element']
                    if 'isotope' not in value:
                        del attrs['isotope']
                    super().__setattr__('_atom', elements_classes[value['element']](**attrs))
                else:
                    super().__setattr__('_atom', elements_classes[self._atom.symbol](**attrs))

            for k in mutable_keys:
                super().__setattr__(k, value[k])

    _hybridization_str = {4: 'a', 3: 't', 2: 'd', 1: 's', None: 'n'}
    _stereo_str = {1: '@', -1: '@@'}
    _multiplicity_str = {1: '*', 2: '*2', 3: '*3', None: 'n'}
    _charge_str = {-3: '-3', -2: '-2', -1: '-', 0: '0', 1: '+', 2: '+2', 3: '+3'}

    __defaults = {'mapping': None, 'mark': '0', 'x': 0, 'y': 0, 'z': 0, 'stereo': None, 'color': None}
    __mutable = {'mapping': lambda x: x is None or isinstance(x, int), 'mark': lambda x: isinstance(x, str),
                 'x': lambda x: isinstance(x, (float, int)), 'y': lambda x: isinstance(x, (float, int)),
                 'z': lambda x: isinstance(x, (float, int)), 'stereo': lambda x: x in (None, -1, 1),
                 'color': lambda x: x is None or isinstance(x, dict) and all(isinstance(y, int) for y in x) and
                                    all(isinstance(y, str) for y in x.values())}


class Bond(MutableMapping):
    __slots__ = ('order', 'stereo', '__allow_none')

    def __init__(self, order=None, stereo=None, allow_none=False):
        self.order = self._defaults['order'] if order is None else order
        self.stereo = self._defaults['stereo'] if stereo is None else stereo
        super().__setattr__('_Bond__allow_none', allow_none)

    def __getitem__(self, key):
        if key not in self._acceptable:
            raise KeyError('unknown bond attribute')
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __setattr__(self, key, value):
        if key not in self._acceptable:
            raise AttributeError('unknown bond attribute')
        if not self._acceptable[key](value, self.__allow_none):
            raise ValueError(f"attribute '{key}' value invalid")
        super().__setattr__(key, value)

    def __len__(self):
        return sum(1 for _ in self)

    def __iter__(self):
        """
        iterate other non-default bonds attrs

        need for nx.readwrite.json_graph.node_link_data
        """
        if self.order != self._defaults['order']:
            yield 'order'
        if self.stereo != self._defaults['stereo']:
            yield 'stereo'

    def __delitem__(self, key):
        raise TypeError('attribute deletion impossible')

    def __eq__(self, other):
        if isinstance(other, Bond):
            return self.order == other.order
        return False

    def __repr__(self):
        s = '' if self.stereo == self._defaults['stereo'] else f', stereo={self.stereo}'
        return f'{type(self).__name__}({self.order}{s})'

    def __hash__(self):
        return hash((self.order, self.stereo))

    def update(self, *args, **kwargs):
        """
        update bond

        :param args: tuple with 1 or 0 elements. element can be Bond object or dict of bond attrs.
        :param kwargs: bond attrs. has precedence other args[0]
        """
        if len(args) > 1:
            raise TypeError('update expected at most 1 arguments')
        elif args:
            value = args[0]
        else:
            value = kwargs
            kwargs = ()

        if isinstance(value, Bond):
            try:
                if not all(self._acceptable[k](v, self.__allow_none) for k, v in kwargs.items()):
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
                if not all(self._acceptable[k](v) for k, v in value.items()):
                    raise ValueError('invalid attribute value')
            except KeyError:
                raise KeyError('unknown bond attributes not allowed')

            for k, v in value.items():
                super().__setattr__(k, v)

    def copy(self):
        return type(self)(self.order, self.stereo, allow_none=self.__allow_none)

    _acceptable = {'order': lambda x, y: x in (1, 2, 3, 4, 9) or y and x is None,
                   'stereo': lambda x, _: x in (None, -1, 1)}
    _defaults = {'order': 1, 'stereo': None}


__all__ = ['Atom', 'Bond']
