# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
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


class Attribute(MutableMapping):
    __slots__ = '_skip_checks'

    def __setitem__(self, key, value):
        try:
            setattr(self, key, value)
        except AttributeError as e:
            raise KeyError from e

    def update(self, *args, **kwargs):
        if len(args) > 1:
            raise TypeError('update expected at most 1 arguments')
        elif args:
            value = args[0]
        else:
            value = kwargs
            kwargs = ()
        self._update(value, kwargs)

    def __delitem__(self, key):
        raise TypeError('attribute deletion impossible')

    def __len__(self):
        return sum(1 for _ in self)

    def copy(self):
        cls = type(self)
        copy = cls.__new__(cls)
        copy.__init_copy__(self)
        return copy

    def _check_kwargs(self, kwargs):
        if not self._skip_checks:
            return {k: getattr(self, f'_{k}_check')(v) for k, v in kwargs.items()}
        return kwargs

    @staticmethod
    def _stereo_check(x):
        if x is None:
            return
        elif x in (-1, 1):
            return x
        raise ValueError('stereo can be: None, 1 or -1')


class AtomAttribute(Attribute):
    @staticmethod
    def _isotope_check(x):
        if x is None:
            return
        elif isinstance(x, int) and x > 0:
            return x
        raise ValueError('invalid isotope')

    @staticmethod
    def _charge_check(x):
        if isinstance(x, int) and -3 <= x <= 3:
            return x
        raise ValueError('invalid charge')

    @staticmethod
    def _multiplicity_check(x):
        if x is None:
            return
        elif isinstance(x, int) and 1 <= x <= 3:
            return x
        raise ValueError('invalid multiplicity')

    @staticmethod
    def _x_check(x):
        return float(x)

    _z_check = _y_check = _x_check


class Atom(AtomAttribute):
    __slots__ = ('_atom', '__hybridization', '__neighbors', '__stereo', '__mapping', '__x', '__y', '__z')

    def __init__(self, *, skip_checks=False):
        super().__setattr__('_skip_checks', skip_checks)
        super().__setattr__('_atom', C())
        super().__setattr__('_Atom__x', 0.)
        super().__setattr__('_Atom__y', 0.)
        super().__setattr__('_Atom__z', 0.)
        super().__setattr__('_Atom__mapping', None)
        super().__setattr__('_Atom__stereo', None)
        super().__setattr__('_Atom__hybridization', None)
        super().__setattr__('_Atom__neighbors', None)

    def __init_copy__(self, parent):
        super().__setattr__('_skip_checks', parent._skip_checks)
        super().__setattr__('_atom', parent._atom)
        super().__setattr__('_Atom__x', parent.x)
        super().__setattr__('_Atom__y', parent.y)
        super().__setattr__('_Atom__z', parent.z)
        super().__setattr__('_Atom__mapping', parent.parsed_mapping)
        super().__setattr__('_Atom__stereo', parent.stereo)
        super().__setattr__('_Atom__hybridization', None)
        super().__setattr__('_Atom__neighbors', None)

    def __setattr__(self, key, value):
        if key == 'element':
            if not self._skip_checks:
                value = self._element_check(value)
            if value != self._atom:
                super().__setattr__('_atom', value(self.charge, self.multiplicity))
        elif key in ('charge', 'isotope', 'multiplicity'):
            if not self._skip_checks:
                value = getattr(self, f'_{key}_check')(value)
            attrs = {'multiplicity': self.multiplicity, 'isotope': self.isotope,
                     'charge': self.charge, key: value}
            super().__setattr__('_atom', type(self._atom)(**attrs))
        elif key == '_neighbors':
            super().__setattr__('_Atom__neighbors', value)
        elif key == '_hybridization':
            super().__setattr__('_Atom__hybridization', value)
        elif key == '_parsed_mapping':
            super().__setattr__('_Atom__mapping', value)
        elif self._skip_checks:
            super().__setattr__(f'_Atom__{key}', value)
        else:
            super().__setattr__(key, value)
        self.__dict__.clear()

    def _update(self, value, kwargs):
        if isinstance(value, Atom):
            if not {'element', 'isotope', 'charge', 'multiplicity'}.isdisjoint(kwargs):
                raise KeyError('element override not allowed')
            kwargs = self._check_kwargs(kwargs)
            super().__setattr__('_atom', value._atom)
            super().__setattr__('_Atom__x', value.x)
            super().__setattr__('_Atom__y', value.y)
            super().__setattr__('_Atom__z', value.z)
            super().__setattr__('_Atom__mapping', value.parsed_mapping)
            super().__setattr__('_Atom__stereo', value.stereo)
            if 'element' in kwargs:
                super().__setattr__('_atom',
                                    kwargs['element'](charge=kwargs.get('charge', self.charge),
                                                      multiplicity=kwargs.get('multiplicity', self.multiplicity),
                                                      isotope=kwargs.get('isotope')))
            elif kwargs.keys() & {'isotope', 'charge', 'multiplicity'}:
                super().__setattr__('_atom',
                                    type(self._atom)(charge=kwargs.get('charge', self.charge),
                                                     multiplicity=kwargs.get('multiplicity', self.multiplicity),
                                                     isotope=kwargs.get('isotope', self.isotope)))
            for k in ('stereo', 'x', 'y', 'z'):
                if k in kwargs:
                    super().__setattr__(f'_Atom__{k}', kwargs[k])

        elif isinstance(value, type):
            if not issubclass(value, Element):
                raise ValueError('only CGRtools.periodictable.Element subclasses allowed')
            if not {'element', 'isotope', 'charge', 'multiplicity'}.isdisjoint(kwargs):
                raise KeyError('element override not allowed')
            kwargs = self._check_kwargs(kwargs)
            super().__setattr__('_atom', value())
            for k, v in kwargs.items():
                super().__setattr__(f'_Atom__{k}', v)
        elif isinstance(value, Element):
            if not {'element', 'isotope', 'charge', 'multiplicity'}.isdisjoint(kwargs):
                raise KeyError('element override not allowed')
            kwargs = self._check_kwargs(kwargs)
            super().__setattr__('_atom', value)
            for k, v in kwargs.items():
                super().__setattr__(f'_Atom__{k}', v)
        else:
            if not isinstance(value, dict):
                try:
                    value = dict(value)
                except (TypeError, ValueError):
                    raise TypeError('invalid attrs sequence')

            value.update(kwargs)
            if not value:  # ad-hoc for add_nodes_from method
                return
            if not {'neighbors', 'hybridization', 'parsed_mapping'}.isdisjoint(kwargs):
                raise KeyError('neighbors hybridization and parsed_mapping change not allowed')

            value = self._check_kwargs(value)
            if 'element' in value:
                super().__setattr__('_atom',
                                    value['element'](charge=value.get('charge', self.charge),
                                                     multiplicity=value.get('multiplicity', self.multiplicity),
                                                     isotope=value.get('isotope')))
            elif value.keys() & {'isotope', 'charge', 'multiplicity'}:
                super().__setattr__('_atom',
                                    type(self._atom)(charge=value.get('charge', self.charge),
                                                     multiplicity=value.get('multiplicity', self.multiplicity),
                                                     isotope=value.get('isotope', self.isotope)))

            for k in ('stereo', 'x', 'y', 'z'):
                if k in value:
                    super().__setattr__(f'_Atom__{k}', value[k])
        self.__dict__.clear()

    def __int__(self):
        # 21b: 7b atom 9b isotope 3b charge 2b mult
        return self.number << 14 | self.isotope << 5 | (self.charge + 3) << 2 | (self.multiplicity or 0)

    def __eq__(self, other):
        """
        == equality checks without stereo
        """
        if isinstance(other, Atom):
            return self._atom == other._atom
        return False

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

    @property
    def stereo(self):
        return self.__stereo

    @stereo.setter
    def stereo(self, value):
        super().__setattr__('_Atom__stereo', self._stereo_check(value))

    @property
    def neighbors(self):
        return self.__neighbors

    @property
    def hybridization(self):
        return self.__hybridization

    @property
    def parsed_mapping(self):
        return self.__mapping

    def __getitem__(self, key):
        """
        dict like access to atom's attrs
        """
        if key == 'element':
            return self._atom.symbol
        elif key in ('isotope', 'charge', 'multiplicity'):
            return getattr(self._atom, key)
        elif key in ('stereo', 'x', 'y', 'z'):
            return getattr(self, key)
        raise KeyError('unknown atom attribute')

    def __getattr__(self, key):
        if key == 'element':
            value = self.__dict__['element'] = self._atom.symbol
        else:
            value = self.__dict__[key] = getattr(self._atom, key)
        return value

    def __iter__(self):
        yield 'element'
        if self.isotope != self.common_isotope:
            yield 'isotope'
        if self.charge:
            yield 'charge'
        if self.multiplicity:
            yield 'multiplicity'
        if self.stereo is not None:
            yield 'stereo'
        for k in ('x', 'y', 'z'):
            if getattr(self, k):
                yield k

    def __getstate__(self):
        return {'checks': self._skip_checks, 'atom': self._atom, 'x': self.__x, 'y': self.__y, 'z': self.__z,
                'mapping': self.__mapping, 'stereo': self.__stereo,
                'hybridization': self.__hybridization, 'neighbors': self.__neighbors}

    def __setstate__(self, state):
        super().__setattr__('_skip_checks', state['checks'])
        super().__setattr__('_atom', state['atom'])
        super().__setattr__('_Atom__x', state['x'])
        super().__setattr__('_Atom__y', state['y'])
        super().__setattr__('_Atom__z', state['z'])
        super().__setattr__('_Atom__mapping', state['mapping'])
        super().__setattr__('_Atom__stereo', state['stereo'])
        super().__setattr__('_Atom__hybridization', state['hybridization'])
        super().__setattr__('_Atom__neighbors', state['neighbors'])

    @staticmethod
    def _element_check(x):
        if x != 'A' and x in elements_classes:
            return elements_classes[x]
        raise ValueError('invalid atom symbol')


class Bond(Attribute):
    __slots__ = ('order', 'stereo')

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

    def _update(self, value, kwargs):
        if isinstance(value, Bond):
            kwargs = self._check_kwargs(kwargs)
            super().__setattr__('order', value.order)
            super().__setattr__('stereo', value.stereo)
            for k, v in kwargs.items():
                super().__setattr__(k, v)
        else:
            self._update_kwargs(value, kwargs)

    def __getitem__(self, key):
        if key not in ('order', 'stereo'):
            raise KeyError('unknown bond attribute')
        return getattr(self, key)

    def __iter__(self):
        yield 'order'
        if self.stereo:
            yield 'stereo'

    def __int__(self):
        return self.order

    def __eq__(self, other):
        """
        equality checks without stereo
        """
        if isinstance(other, Bond):
            return self.order == other.order
        return False

    def __getstate__(self):
        return {'checks': self._skip_checks, 'order': self.order, 'stereo': self.stereo}

    def __setstate__(self, state):
        super().__setattr__('_skip_checks', state['checks'])
        super().__setattr__('order', state['order'])
        super().__setattr__('stereo', state['stereo'])

    @staticmethod
    def _order_check(x):
        if x in (1, 2, 3, 4, 5):
            return x
        raise ValueError('invalid order')

    def _update_kwargs(self, value, kwargs):
        if not isinstance(value, dict):
            try:
                value = dict(value)
            except (TypeError, ValueError):
                raise TypeError('invalid attrs sequence')

        value.update(kwargs)
        if value:
            for k, v in self._check_kwargs(value).items():
                super().__setattr__(k, v)


__all__ = ['Atom', 'Bond']
