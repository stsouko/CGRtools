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
from .molecule import Bond, Atom
from ..periodictable import Element, elements_classes, elements_numbers


class QueryAtom(MutableMapping):
    __slots__ = ('_atom', '_skip_checks')

    def __init__(self, *, skip_checks=False):
        super().__setattr__('_skip_checks', skip_checks)
        super().__setattr__('_atom', self.__defaults.copy())

    def __init_copy__(self, parent):
        super().__setattr__('_atom', parent._atom.copy())

    def __setattr__(self, key, value):
        if not self._skip_checks:
            value = getattr(self, f'_{key}_check')(value)
        self._atom[key] = value

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
        self._update(value, kwargs)

    def _update(self, value, kwargs):
        if isinstance(value, QueryAtom):
            kwargs = self._check_kwargs(kwargs)
            self._atom.update(value)
            self._atom.update(kwargs)
        elif isinstance(value, Atom):
            kwargs = self._check_kwargs(kwargs)
            self._atom.update(element=(value.element,), charge=(value.charge,), stereo=value.stereo,
                              multiplicity=(value.multiplicity,) if value.multiplicity else (),
                              isotope=(value.isotope,) if value.isotope != value.common_isotope else (),
                              neighbors=(value.neighbors,) if value.neighbors else (),
                              hybridization=(value.hybridization,) if value.hybridization else ())
            self._atom.update(kwargs)
        elif isinstance(value, type):
            if not issubclass(value, Element):
                ValueError('only CGRtools.periodictable.Element subclasses allowed')
            if not {'charge', 'multiplicity', 'isotope', 'element'}.isdisjoint(kwargs):
                raise KeyError('charge, multiplicity, isotope and element override not allowed')
            kwargs = self._check_kwargs(kwargs)
            self._atom.update(element=(value.symbol,), charge=(0,), multiplicity=(), isotope=())
            self._atom.update(kwargs)
        elif isinstance(value, Element):
            if not {'charge', 'multiplicity', 'isotope', 'element'}.isdisjoint(kwargs):
                raise KeyError('charge, multiplicity, isotope and element override not allowed')
            kwargs = self._check_kwargs(kwargs)
            self._atom.update(element=(value.symbol,), charge=(value.charge,),
                              multiplicity=(value.multiplicity,) if value.multiplicity else (),
                              isotope=(value.isotope,) if value.isotope != value.common_isotope else ())
            self._atom.update(kwargs)
        else:
            if not isinstance(value, dict):
                try:
                    value = dict(value)
                except (TypeError, ValueError):
                    raise TypeError('invalid attrs sequence')

            value.update(kwargs)
            if not value:  # ad-hoc for add_nodes_from method
                return

            value = self._check_kwargs(value)
            self._atom.update(value)

    def stringify(self, atom=True, isotope=True, stereo=True, hybridization=True, neighbors=True):
        smi = []
        if stereo and self.stereo:
            smi.append(self._stereo_str[self.stereo])
        if hybridization:
            if len(self.hybridization) > 1:
                smi.append('<%s>' % ''.join(self._hybridization_str[x] for x in sorted(self.hybridization)))
            elif self.hybridization:
                smi.append(self._hybridization_str[self.hybridization[0]])
        if neighbors:
            if len(self.neighbors) > 1:
                smi.append('<%s>' % ''.join(str(x) for x in sorted(self.neighbors)))
            elif self.neighbors:
                smi.append(str(self.neighbors[0]))
        if smi:
            smi.append(';')
            smi.insert(0, ';')

        if atom:
            if self.element == ('A',):
                atom = False
                smi.insert(0, '*')
            elif len(self.element) > 1:
                atom = False
                smi.insert(0, ','.join(sorted(self.element, key=elements_numbers.get)))
            else:
                if self.element[0] not in ('C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B'):
                    atom = False
                smi.insert(0, self.element[0])

            if len(self.charge) > 1:
                smi.append('<%s>' % ''.join(self._charge_str[x] for x in sorted(self.charge)))
            elif self.charge != (0,):
                smi.append(self._charge_str[self.charge[0]])

            if len(self.multiplicity) > 1:
                smi.append('<%s>' % ''.join(self._multiplicity_str[x] for x in sorted(self.multiplicity)))
            elif self.multiplicity:
                smi.append(self._multiplicity_str[self.multiplicity[0]])

            if isotope:
                if len(self.isotope) > 1:
                    smi.insert(0, '<%s>' % ''.join(str(x) for x in sorted(self.isotope)))
                elif self.isotope:
                    smi.insert(0, str(self.isotope[0]))
        else:
            smi.insert(0, '*')

        if len(smi) != 1 or not atom:
            smi.insert(0, '[')
            smi.append(']')

        return ''.join(smi)

    def weight(self, atom=True, isotope=False, stereo=False, hybridization=False, neighbors=False):
        weight = []
        if atom:
            weight.append(tuple(sorted(elements_numbers[x] for x in self.element)))
            if isotope:
                weight.append(tuple(sorted(self.isotope)))
            weight.append(tuple(sorted(self.charge)))
            weight.append(self.multiplicity and tuple(sorted(self.multiplicity)))
        if stereo:
            weight.append(self.stereo or 0)
        if hybridization:
            weight.append(self.hybridization and tuple(sorted(self.hybridization)))
        if neighbors:
            weight.append(self.neighbors and tuple(sorted(self.neighbors)))
        return tuple(weight)

    def __eq__(self, other):
        if isinstance(other, QueryAtom):
            return all(sorted(self[attr]) == sorted(other[attr])
                       for attr in ('element', 'isotope', 'charge', 'multiplicity', 'hybridization', 'neighbors'))
        elif isinstance(other, Atom):
            return (other.element in self.element or self.element == ('A',)) and \
                    all(other[attr] in self[attr] for attr in
                        ('isotope', 'charge', 'multiplicity', 'hybridization', 'neighbors') if self[attr])
        return False

    def __ne__(self, other):
        if self == other:
            if self.stereo:
                return self.stereo == other.stereo
            return True
        return False

    def __getitem__(self, key):
        """
        dict like access to atom's attrs
        """
        return self._atom[key]

    def __getattr__(self, key):
        if key == '__dict__':
            raise AttributeError()
        try:
            return self._atom[key]
        except KeyError as e:
            raise AttributeError from e

    def __setitem__(self, key, value):
        try:
            setattr(self, key, value)
        except AttributeError as e:
            raise KeyError from e

    def __len__(self):
        return sum(1 for _ in self)

    def __iter__(self):
        return iter(self._atom)

    def __delitem__(self, key):
        raise TypeError('attribute deletion impossible')

    def __str__(self):
        return self.stringify()

    def copy(self):
        cls = type(self)
        copy = cls.__new__(cls)
        copy.__init_copy__(self)
        return copy

    @staticmethod
    def _element_check(x):
        if x == ('A',):
            return x
        elif x == 'A':
            return 'A',
        elif not isinstance(x, tuple):
            x = tuple(x)
        if 0 < len(x) == len(set(x)) and all(x != 'A' and x in elements_classes for x in x):
            return x
        raise ValueError('invalid element')

    @staticmethod
    def _stereo_check(x):
        if x is None:
            return None
        x = int(x)
        if x in (-1, 1):
            return x
        raise ValueError('stereo can be: None, 1 or -1')

    @staticmethod
    def _isotope_check(x):
        if isinstance(x, int):
            if x > 0:
                return x,
            raise ValueError('invalid isotope')
        elif not isinstance(x, tuple):
            x = tuple(x)
        if len(x) == len(set(x)) and all(isinstance(x, int) and x > 0 for x in x):
            return x
        raise ValueError('invalid isotope')

    @staticmethod
    def _charge_check(x):
        if isinstance(x, int):
            if -3 <= x <= 3:
                return x,
            raise ValueError('invalid charge')
        if not isinstance(x, tuple):
            x = tuple(x)
        if 0 < len(x) == len(set(x)) and all(isinstance(x, int) and -3 <= x <= 3 for x in x):
            return x
        raise ValueError('invalid charge')

    @staticmethod
    def _multiplicity_check(x):
        if isinstance(x, int):
            if 1 <= x <= 3:
                return x,
            raise ValueError('invalid multiplicity')
        if not isinstance(x, tuple):
            x = tuple(x)
        if len(x) == len(set(x)) and all(isinstance(x, int) and 1 <= x <= 3 for x in x):
            return x
        raise ValueError('invalid multiplicity')

    @staticmethod
    def _neighbors_check(x):
        if isinstance(x, int):
            if 0 <= x <= 998:
                return x,
            raise ValueError('invalid neighbors')
        if not isinstance(x, tuple):
            x = tuple(x)
        if len(x) == len(set(x)) and all(isinstance(x, int) and 0 <= x <= 998 for x in x):
            return x
        raise ValueError('invalid neighbors')

    @staticmethod
    def _hybridization_check(x):
        if isinstance(x, int):
            if 1 <= x <= 4:
                return x,
            raise ValueError('invalid hybridization')
        if not isinstance(x, tuple):
            x = tuple(x)
        if len(x) == len(set(x)) and all(isinstance(x, int) and 1 <= x <= 4 for x in x):
            return x
        raise ValueError('invalid hybridization')

    def _check_kwargs(self, kwargs):
        if not self._skip_checks:
            kwargs = {k: getattr(self, f'_{k}_check')(v) for k, v in kwargs.items()}
        return kwargs

    _hybridization_str = Atom._hybridization_str
    _stereo_str = Atom._stereo_str
    _multiplicity_str = Atom._multiplicity_str
    _charge_str = Atom._charge_str
    __defaults = {'element': ('A',), 'isotope': (), 'charge': (0,), 'multiplicity': (),
                  'neighbors': (), 'hybridization': (), 'stereo': None}


class QueryBond(Bond):
    def _update(self, value, kwargs):
        if isinstance(value, Bond):
            kwargs = self._check_kwargs(kwargs)
            if isinstance(value, QueryBond):
                super().__setattr__('order', value.order)
            else:
                super().__setattr__('order', (value.order,))
            super().__setattr__('stereo', value.stereo)
            for k, v in kwargs.items():
                super().__setattr__(k, v)
        else:
            super()._update(value, kwargs)

    def stringify(self, stereo=True):
        order = '<%s>' % ''.join(sorted(self._order_str[x] for x in sorted(self.order)))
        if stereo and self.stereo:
            return order + self._stereo_str[self.stereo]
        return order

    def weight(self, stereo=False):
        if stereo:
            return tuple(sorted(self.order)), self.stereo or 0
        return self.order

    def __eq__(self, other):
        """
        == equality checks. if stereo mark is presented in query, stereo also will be compared
        """
        if isinstance(other, QueryBond):
            return sorted(self.order) == sorted(other.order)
        elif isinstance(other, Bond):
            return other.order in self.order
        return False

    def __ne__(self, other):
        """
        != equality checks with stereo
        """
        if self == other:
            if self.stereo:
                return self.stereo == other.stereo
            return True
        return False

    @staticmethod
    def _order_check(x):
        if isinstance(x, int):
            if x in (1, 2, 3, 4, 9):
                return x,
            raise ValueError('invalid order')
        elif not isinstance(x, tuple):
            x = tuple(x)
        if x and all(x in (1, 2, 3, 4, 9) for x in x):
            return x
        raise ValueError('invalid order')

    _defaults = {'order': (1,), 'stereo': None}


__all__ = ['QueryAtom', 'QueryBond']
