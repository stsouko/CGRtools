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
from ..periodictable import Element, elements_classes, C


class Atom(MutableMapping):
    __slots__ = ('_atom', '__hybridization', '__neighbors', '__color', '__stereo', '__mapping', '__mark',
                 '__x', '__y', '__z', '_skip_checks')

    def __init__(self, *, skip_checks=False):
        super().__setattr__('_skip_checks', skip_checks)
        super().__setattr__('_atom', C())
        super().__setattr__('_Atom__x', 0.)
        super().__setattr__('_Atom__y', 0.)
        super().__setattr__('_Atom__z', 0.)
        super().__setattr__('_Atom__mapping', None)
        super().__setattr__('_Atom__stereo', None)
        super().__setattr__('_Atom__mark', None)
        super().__setattr__('_Atom__color', None)
        super().__setattr__('_Atom__hybridization', None)
        super().__setattr__('_Atom__neighbors', None)

    def __init_copy__(self, parent):
        super().__setattr__('_skip_checks', parent._skip_checks)
        super().__setattr__('_atom', parent._atom)
        super().__setattr__('_Atom__x', parent.x)
        super().__setattr__('_Atom__y', parent.y)
        super().__setattr__('_Atom__z', parent.z)
        super().__setattr__('_Atom__mapping', parent.mapping)
        super().__setattr__('_Atom__stereo', parent.stereo)
        super().__setattr__('_Atom__mark', parent.mark)
        super().__setattr__('_Atom__color', parent.color)
        super().__setattr__('_Atom__hybridization', None)
        super().__setattr__('_Atom__neighbors', None)

    def __setattr__(self, key, value):
        if key == 'element':
            value = self._element_check(value)
            if value != self._atom:
                super().__setattr__('_atom', value(self.charge, self.multiplicity))
        elif key in ('charge', 'isotope', 'multiplicity'):
            value = getattr(self, f'_{key}_check')(value)
            attrs = {'multiplicity': self.multiplicity, 'isotope': self.isotope,
                     'charge': self.charge, key: value}
            super().__setattr__('_atom', type(self._atom)(**attrs))
        elif key in ('neighbors', 'hybridization'):
            raise AttributeError('neighbors and hybridization change not allowed')
        else:
            if key in ('_neighbors', '_hybridization'):
                key = key[1:]
            if self._skip_checks:
                super().__setattr__(f'_Atom__{key}', value)
            else:
                super().__setattr__(key, value)

    def update(self, *args, **kwargs):
        """
        update atom

        :param args: tuple with 1 or 0 elements. element can be dict of atom attrs or atom object or atom class.
        :param kwargs: atom attrs. has precedence other args[0] if it's dict
        """
        if len(args) > 1:
            raise TypeError('update expected at most 1 arguments')
        elif args:
            value = args[0]
        else:
            value = kwargs
            kwargs = ()
        self._update(value, kwargs)

    def _update(self, value, kwargs):
        if isinstance(value, Atom):
            kwargs = self._check_kwargs(kwargs)
            super().__setattr__('_atom', value._atom)
            super().__setattr__('_Atom__x', value.x)
            super().__setattr__('_Atom__y', value.y)
            super().__setattr__('_Atom__z', value.z)
            super().__setattr__('_Atom__mapping', value.mapping)
            super().__setattr__('_Atom__stereo', value.stereo)
            super().__setattr__('_Atom__mark', value.mark)
            super().__setattr__('_Atom__color', value.color)
            for k, v in kwargs.items():
                super().__setattr__(k, v)
        if isinstance(value, type):
            if not issubclass(value, Element):
                ValueError('only CGRtools.periodictable.Element subclasses allowed')
            if not {'neighbors', 'hybridization', 'element'}.isdisjoint(kwargs):
                raise KeyError('neighbors and hybridization change not allowed. element override not allowed')
            kwargs = self._check_kwargs(kwargs)
            super().__setattr__('_atom', value())
            for k, v in kwargs.items():
                super().__setattr__(k, v)
        elif isinstance(value, Element):
            if not {'neighbors', 'hybridization', 'element'}.isdisjoint(kwargs):
                raise KeyError('neighbors and hybridization change not allowed. element override not allowed')
            kwargs = self._check_kwargs(kwargs)
            super().__setattr__('_atom', value)
            for k, v in kwargs.items():
                super().__setattr__(k, v)
        else:
            if not isinstance(value, dict):
                try:
                    value = dict(value)
                except (TypeError, ValueError):
                    raise TypeError('invalid attrs sequence')

            value.update(kwargs)
            if not value:  # ad-hoc for add_nodes_from method
                return
            if not {'neighbors', 'hybridization'}.isdisjoint(kwargs):
                raise KeyError('neighbors and hybridization change not allowed')
            mutable_keys = self._check_kwargs({k: value[k] for k in value.keys() -
                                               {'element', 'isotope', 'charge', 'multiplicity'}})

            if 'element' in value:
                element = self._element_check(value['element'])
                attrs = {'multiplicity': value.get('multiplicity', self.multiplicity),
                         'charge': value.get('charge', self.charge)}
                if 'isotope' in value:
                    attrs['isotope'] = value['isotope']
                super().__setattr__('_atom', element(**attrs))
            else:
                atom_keys = value.keys() & {'isotope', 'charge', 'multiplicity'}
                if atom_keys:
                    attrs = {'multiplicity': self.multiplicity, 'isotope': self.isotope,
                             'charge': self.charge, **{k: value[k] for k in atom_keys}}
                    super().__setattr__('_atom', type(self._atom)(**attrs))

            for k, v in mutable_keys.items():
                super().__setattr__(k, v)

    def stringify(self, atom=True, isotope=True, stereo=False, hybridization=False, neighbors=False):
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

        if atom:
            if self.element not in ('C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B'):
                atom = False
            smi.insert(0, self.element)
            if self.charge:
                smi.append(self._charge_str[self.charge])
            if self.multiplicity:
                smi.append(self._multiplicity_str[self.multiplicity])
            if isotope and self.isotope != self.common_isotope:
                smi.insert(0, str(self.isotope))
        else:
            smi.insert(0, '*')

        if len(smi) != 1 or not atom:
            smi.insert(0, '[')
            smi.append(']')

        return ''.join(smi)

    def weight(self, atom=True, isotope=False, stereo=False, hybridization=False, neighbors=False):
        weight = []
        if atom:
            weight.append(self.number)
            if isotope:
                weight.append(self.isotope)
            weight.append(self.charge)
            weight.append(self.multiplicity or 0)
        if stereo:
            weight.append(self.__stereo or 0)
        if hybridization:
            weight.append(self.__hybridization or 0)
        if neighbors:
            weight.append(self.__neighbors or -1)
        return tuple(weight)

    def __eq__(self, other):
        """
        == equality checks without stereo
        """
        if isinstance(other, Atom):
            return self._atom == other._atom
        return False

    def __ne__(self, other):
        """
        != equality checks with stereo
        """
        return self == other and self.stereo == other.stereo

    @property
    def x(self):
        return self.__x

    @x.setter
    def x(self, value):
        super().__setattr__('_Atom__x', self._x_check(value))

    @property
    def y(self):
        return self.__y

    @y.setter
    def y(self, value):
        super().__setattr__('_Atom__y', self._y_check(value))

    @property
    def z(self):
        return self.__z

    @z.setter
    def z(self, value):
        super().__setattr__('_Atom__z', self._z_check(value))

    @staticmethod
    def _x_check(x):
        return float(x)

    _z_check = _y_check = _x_check

    @property
    def mapping(self):
        return self.__mapping

    @mapping.setter
    def mapping(self, value):
        super().__setattr__('_Atom__mapping', self._mapping_check(value))

    @staticmethod
    def _mapping_check(x):
        if x is None:
            return None
        x = int(x)
        if 0 <= x <= 999:
            return x
        raise ValueError('mapping can be in range 0-999')

    @property
    def stereo(self):
        return self.__stereo

    @stereo.setter
    def stereo(self, value):
        super().__setattr__('_Atom__stereo', self._stereo_check(value))

    @staticmethod
    def _stereo_check(x):
        if x is None:
            return None
        x = int(x)
        if x in (-1, 1):
            return x
        raise ValueError('stereo can be: None, 1 or -1')

    @property
    def neighbors(self):
        return self.__neighbors

    @neighbors.setter
    def neighbors(self, value):
        super().__setattr__('_Atom__neighbors', self._neighbors_check(value))

    @staticmethod
    def _neighbors_check(x):
        if x is None:
            return None
        x = int(x)
        if 0 <= x <= 998:
            return x
        raise ValueError('neighbors can be: None or in range 0-998')

    @property
    def hybridization(self):
        return self.__hybridization

    @hybridization.setter
    def hybridization(self, value):
        super().__setattr__('_Atom__hybridization', self._hybridization_check(value))

    @staticmethod
    def _hybridization_check(x):
        if x is None:
            return None
        x = int(x)
        if x in (1, 2, 3, 4):
            return x
        raise ValueError('hybridization can be: None, 1, 2, 3 or 4')

    @property
    def mark(self):
        return self.__mark

    @mark.setter
    def mark(self, value):
        super().__setattr__('_Atom__mark', self._mark_check(value))

    @staticmethod
    def _mark_check(x):
        if x is None:
            return None
        x = int(x)
        if 1 <= x <= 999:
            return x
        raise ValueError('mark can be: None or in range 1-999')

    @property
    def color(self):
        if self.__color:
            return self.__color

    @color.setter
    def color(self, value):
        super().__setattr__('_Atom__color', self._color_check(value))

    @staticmethod
    def _color_check(x):
        if x is None:
            return None
        if not isinstance(x, dict):
            raise TypeError('color can be dict')
        elif not x:
            raise ValueError('empty colors dict not allowed')
        elif not all(isinstance(x, int) for x in x):
            raise TypeError('color keys can be int')
        elif not all(isinstance(x, str) for x in x.values()):
            raise TypeError('color values can be str')
        return x.copy()

    def __getitem__(self, key):
        """
        dict like access to atom's attrs
        """
        if key == 'element':
            return self._atom.symbol
        elif key in ('isotope', 'charge', 'multiplicity'):
            return getattr(self._atom, key)
        elif key in ('color', 'stereo', 'mapping', 'mark', 'x', 'y', 'z'):
            return getattr(self, key)
        raise KeyError('unknown atom attribute')

    def __setitem__(self, key, value):
        try:
            setattr(self, key, value)
        except AttributeError as e:
            raise KeyError from e

    def __getattr__(self, key):
        if key == 'element':
            return self._atom.symbol
        return getattr(self._atom, key)

    def __len__(self):
        return sum(1 for _ in self)

    def __iter__(self):
        """
        iterate other non-default atom's init attrs
        """
        yield 'element'
        if self.isotope != self.common_isotope:
            yield 'isotope'
        if self.charge:
            yield 'charge'
        if self.multiplicity:
            yield 'multiplicity'
        for k, d in self.__defaults.items():
            if d != getattr(self, k):
                yield k

    def __delitem__(self, key):
        raise TypeError('attribute deletion impossible')

    def __str__(self):
        return self.stringify(stereo=True)

    def copy(self):
        cls = type(self)
        copy = cls.__new__(cls)
        copy.__init_copy__(self)
        return copy

    @staticmethod
    def _element_check(x):
        if x != 'A' and x in elements_classes:
            return elements_classes[x]
        raise ValueError('invalid atom symbol')

    def _check_kwargs(self, kwargs):
        if not self._skip_checks:
            return {f'_Atom__{k}': getattr(self, f'_{k}_check')(v) for k, v in kwargs.items()}
        return {f'_Atom__{k}': v for k, v in kwargs.items()}

    _hybridization_str = {4: 'a', 3: 't', 2: 'd', 1: 's', None: 'n'}
    _stereo_str = {1: '@', -1: '@@'}
    _multiplicity_str = {1: '*', 2: '*2', 3: '*3', None: 'n'}
    _charge_str = {-3: '-3', -2: '-2', -1: '-', 0: '0', 1: '+', 2: '+2', 3: '+3'}
    __defaults = {'mapping': None, 'mark': None, 'x': 0., 'y': 0., 'z': 0., 'stereo': None, 'color': None}


class Bond(MutableMapping):
    __slots__ = ('order', 'stereo', '_skip_checks')

    def __init__(self, *, skip_checks=False):
        super().__setattr__('_skip_checks', skip_checks)
        super().__setattr__('order', 1)
        super().__setattr__('stereo', None)

    def __init_copy__(self, parent):
        super().__setattr__('_skip_checks', parent._skip_checks)
        super().__setattr__('order', parent.order)
        super().__setattr__('stereo', parent.stereo)

    def __setattr__(self, key, value):
        if not self._skip_checks:
            value = getattr(self, f'_{key}_check')(value)
        super().__setattr__(key, value)

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
        self._update(value, kwargs)

    def _update(self, value, kwargs):
        if isinstance(value, Bond):
            kwargs = self._check_kwargs(kwargs)
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

            for k, v in self._check_kwargs(value).items():
                super().__setattr__(k, v)

    def stringify(self, stereo=False):
        if stereo and self.stereo:
            return self._order_str[self.order] + self._stereo_str[self.stereo]
        return self._order_str[self.order]

    def weight(self, stereo=False):
        if stereo:
            return self.order, self.stereo or 0
        return self.order

    def __eq__(self, other):
        """
        == equality checks without stereo
        """
        if isinstance(other, Bond):
            return self.order == other.order
        return False

    def __ne__(self, other):
        """
        != equality checks with stereo
        """
        return self == other and self.stereo == other.stereo

    def __getitem__(self, key):
        if key not in self._defaults:
            raise KeyError('unknown bond attribute')
        return getattr(self, key)

    def __setitem__(self, key, value):
        try:
            setattr(self, key, value)
        except AttributeError as e:
            raise KeyError from e

    def __len__(self):
        return sum(1 for _ in self)

    def __iter__(self):
        """
        iterate other non-default bonds attrs

        need for nx.readwrite.json_graph.node_link_data
        """
        yield 'order'
        if self.stereo != self._defaults['stereo']:
            yield 'stereo'

    def __delitem__(self, key):
        raise TypeError('attribute deletion impossible')

    def __str__(self):
        return self.stringify(stereo=True)

    def copy(self):
        cls = type(self)
        copy = cls.__new__(cls)
        copy.__init_copy__(self)
        return copy

    def _check_kwargs(self, kwargs):
        if not self._skip_checks:
            kwargs = {k: getattr(self, f'_{k}_check')(v) for k, v in kwargs.items()}
        return kwargs

    @staticmethod
    def _order_check(x):
        if isinstance(x, int) and x in (1, 2, 3, 4, 9):
            return x
        raise ValueError('invalid order')

    @staticmethod
    def _stereo_check(x):
        if x is None:
            return None
        x = int(x)
        if x in (-1, 1):
            return x
        raise ValueError('stereo can be: None, 1 or -1')

    _defaults = {'order': 1, 'stereo': None}
    _order_str = {1: '-', 2: '=', 3: '#', 4: ':', 9: '~'}
    _stereo_str = {1: '@', -1: '@@'}


__all__ = ['Atom', 'Bond']
