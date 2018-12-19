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

    def __str__(self):
        return self.stringify(stereo=True)

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

    _stereo_str = {1: '@', -1: '@@'}


class AtomAttribute(Attribute):
    @staticmethod
    def _x_check(x):
        return float(x)

    _z_check = _y_check = _x_check
    _hybridization_str = {4: 'a', 3: 't', 2: 'd', 1: 's', None: 'n'}
    _multiplicity_str = {1: '*', 2: '*2', 3: '*3', None: 'n'}
    _charge_str = {-3: '-3', -2: '-2', -1: '-', 0: '0', 1: '+', 2: '+2', 3: '+3'}


class BondAttribute(Attribute):
    __slots__ = ('order', 'stereo')

    def __init_copy__(self, parent):
        super().__setattr__('_skip_checks', parent._skip_checks)
        super().__setattr__('order', parent.order)
        super().__setattr__('stereo', parent.stereo)

    def __setattr__(self, key, value):
        if not self._skip_checks:
            value = getattr(self, f'_{key}_check')(value)
        super().__setattr__(key, value)

    def __getitem__(self, key):
        if key not in ('order', 'stereo'):
            raise KeyError('unknown bond attribute')
        return getattr(self, key)

    def __iter__(self):
        yield 'order'
        if self.stereo:
            yield 'stereo'

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

    def __getstate__(self):
        return {'checks': self._skip_checks, 'order': self.order, 'stereo': self.stereo}

    def __setstate__(self, state):
        super().__setattr__('_skip_checks', state['checks'])
        super().__setattr__('order', state['order'])
        super().__setattr__('stereo', state['stereo'])

    _order_str = {1: '-', 2: '=', 3: '#', 4: ':', 9: '~', None: '.'}


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
        super().__setattr__('_Atom__mapping', parent.mapping)
        super().__setattr__('_Atom__stereo', parent.stereo)
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
        elif key == '_neighbors':
            super().__setattr__('_Atom__neighbors', value)
        elif key == '_hybridization':
            super().__setattr__('_Atom__hybridization', value)
        elif self._skip_checks:
            super().__setattr__(f'_Atom__{key}', value)
        else:
            super().__setattr__(key, value)

    def _update(self, value, kwargs):
        if isinstance(value, Atom):
            if not {'element', 'isotope', 'charge', 'multiplicity'}.isdisjoint(kwargs):
                raise KeyError('element override not allowed')
            kwargs = self._check_kwargs(kwargs)
            super().__setattr__('_atom', value._atom)
            super().__setattr__('_Atom__x', value.x)
            super().__setattr__('_Atom__y', value.y)
            super().__setattr__('_Atom__z', value.z)
            super().__setattr__('_Atom__mapping', value.mapping)
            super().__setattr__('_Atom__stereo', value.stereo)
            if '_Atom__element' in kwargs:
                super().__setattr__('_atom',
                                    kwargs['_Atom__element'](charge=kwargs.get('_Atom__charge', self.charge),
                                                             multiplicity=kwargs.get('_Atom__multiplicity',
                                                                                     self.multiplicity),
                                                             isotope=kwargs.get('_Atom__isotope')))
            elif kwargs.keys() & {'_Atom__isotope', '_Atom__charge', '_Atom__multiplicity'}:
                super().__setattr__('_atom',
                                    type(self._atom)(charge=kwargs.get('_Atom__charge', self.charge),
                                                     multiplicity=kwargs.get('_Atom__multiplicity', self.multiplicity),
                                                     isotope=kwargs.get('_Atom__isotope', self.isotope)))
            for k in ('_Atom__stereo', '_Atom__mapping', '_Atom__x', '_Atom__y', '_Atom__z'):
                if k in kwargs:
                    super().__setattr__(k, kwargs[k])

        elif isinstance(value, type):
            if not issubclass(value, Element):
                ValueError('only CGRtools.periodictable.Element subclasses allowed')
            if not {'element', 'isotope', 'charge', 'multiplicity'}.isdisjoint(kwargs):
                raise KeyError('element override not allowed')
            kwargs = self._check_kwargs(kwargs)
            super().__setattr__('_atom', value())
            for k, v in kwargs.items():
                super().__setattr__(k, v)
        elif isinstance(value, Element):
            if not {'element', 'isotope', 'charge', 'multiplicity'}.isdisjoint(kwargs):
                raise KeyError('element override not allowed')
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

            value = self._check_kwargs(value)
            if '_Atom__element' in value:
                super().__setattr__('_atom', 
                                    value['_Atom__element'](charge=value.get('_Atom__charge', self.charge),
                                                            multiplicity=value.get('_Atom__multiplicity',
                                                                                   self.multiplicity),
                                                            isotope=value.get('_Atom__isotope')))
            elif value.keys() & {'_Atom__isotope', '_Atom__charge', '_Atom__multiplicity'}:
                super().__setattr__('_atom',
                                    type(self._atom)(charge=value.get('_Atom__charge', self.charge),
                                                     multiplicity=value.get('_Atom__multiplicity', self.multiplicity),
                                                     isotope=value.get('_Atom__isotope', self.isotope)))

            for k in ('_Atom__stereo', '_Atom__mapping', '_Atom__x', '_Atom__y', '_Atom__z'):
                if k in value:
                    super().__setattr__(k, value[k])

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

    @property
    def mapping(self):
        return self.__mapping

    @mapping.setter
    def mapping(self, value):
        super().__setattr__('_Atom__mapping', self._mapping_check(value))

    @staticmethod
    def _mapping_check(x):
        if x is None:
            return
        elif isinstance(x, int) and 0 <= x <= 999:
            return x
        raise ValueError('mapping can be in range 0-999')

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

    def __getitem__(self, key):
        """
        dict like access to atom's attrs
        """
        if key == 'element':
            return self._atom.symbol
        elif key in ('isotope', 'charge', 'multiplicity'):
            return getattr(self._atom, key)
        elif key in ('stereo', 'mapping', 'x', 'y', 'z', 'neighbors', 'hybridization'):
            return getattr(self, key)
        raise KeyError('unknown atom attribute')

    def __getattr__(self, key):
        if key == '__dict__':
            raise AttributeError()
        elif key == 'element':
            return self._atom.symbol
        return getattr(self._atom, key)

    def __iter__(self):
        yield 'element'
        if self.isotope != self.common_isotope:
            yield 'isotope'
        if self.charge:
            yield 'charge'
        if self.multiplicity:
            yield 'multiplicity'
        if self.mapping is not None:
            yield 'mapping'
        if self.stereo is not None:
            yield 'stereo'
        for k in ('x', 'y', 'z'):
            if getattr(self, k):
                yield k

    @staticmethod
    def _element_check(x):
        if x != 'A' and x in elements_classes:
            return elements_classes[x]
        raise ValueError('invalid atom symbol')

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

    def _check_kwargs(self, kwargs):
        if not self._skip_checks:
            return {f'_Atom__{k}': getattr(self, f'_{k}_check')(v) for k, v in kwargs.items()}
        return {f'_Atom__{k}': v for k, v in kwargs.items()}

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


class Bond(BondAttribute):
    def __init__(self, *, skip_checks=False):
        super(BondAttribute, self).__setattr__('_skip_checks', skip_checks)
        super(BondAttribute, self).__setattr__('order', 1)
        super(BondAttribute, self).__setattr__('stereo', None)

    def _update(self, value, kwargs):
        if isinstance(value, Bond):
            kwargs = self._check_kwargs(kwargs)
            super(BondAttribute, self).__setattr__('order', value.order)
            super(BondAttribute, self).__setattr__('stereo', value.stereo)
            for k, v in kwargs.items():
                super(BondAttribute, self).__setattr__(k, v)
        else:
            self._update_kwargs(value, kwargs)

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

    @staticmethod
    def _order_check(x):
        if x in (1, 2, 3, 4, 9):
            return x
        raise ValueError('invalid order')


__all__ = ['Atom', 'Bond']
