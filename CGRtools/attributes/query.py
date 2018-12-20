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
from .molecule import Bond, Atom, AtomAttribute, BondAttribute
from ..periodictable import Element, elements_classes, elements_numbers


class QueryAtom(AtomAttribute):
    __slots__ = '_atom'

    def __init__(self, *, skip_checks=False):
        super().__setattr__('_skip_checks', skip_checks)
        super().__setattr__('_atom', {'element': ('A',), 'isotope': (), 'charge': (0,), 'multiplicity': (),
                                      'neighbors': (), 'hybridization': (), 'stereo': None, 'x': 0., 'y': 0., 'z': 0.})

    def __init_copy__(self, parent):
        super().__setattr__('_skip_checks', parent._skip_checks)
        super().__setattr__('_atom', parent._atom.copy())

    def __setattr__(self, key, value):
        if not self._skip_checks:
            value = getattr(self, f'_{key}_check')(value)
        self._atom[key] = value

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
                              hybridization=(value.hybridization,) if value.hybridization else (),
                              x=value.x, y=value.y, z=value.z)
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
                smi.append('<%s>' % ''.join(self._hybridization_str[x] for x in self.hybridization))
            elif self.hybridization:
                smi.append(self._hybridization_str[self.hybridization[0]])
        if neighbors:
            if len(self.neighbors) > 1:
                smi.append('<%s>' % ''.join(str(x) for x in self.neighbors))
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
                smi.insert(0, ','.join(self.element))
            else:
                if self.element[0] not in ('C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B'):
                    atom = False
                smi.insert(0, self.element[0])

            if len(self.charge) > 1:
                smi.append('<%s>' % ''.join(self._charge_str[x] for x in self.charge))
            elif self.charge != (0,):
                smi.append(self._charge_str[self.charge[0]])

            if len(self.multiplicity) > 1:
                smi.append('<%s>' % ''.join(self._multiplicity_str[x] for x in self.multiplicity))
            elif self.multiplicity:
                smi.append(self._multiplicity_str[self.multiplicity[0]])

            if isotope:
                if len(self.isotope) > 1:
                    smi.insert(0, '<%s>' % ','.join(str(x) for x in self.isotope))
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
            weight.append(tuple(elements_numbers[x] for x in self.element))
            if isotope:
                weight.append(self.isotope)
            weight.append(self.charge)
            weight.append(self.multiplicity)
        if stereo:
            weight.append(self.stereo or 0)
        if hybridization:
            weight.append(self.hybridization)
        if neighbors:
            weight.append(self.neighbors)
        return tuple(weight)

    def __eq__(self, other):
        if isinstance(other, QueryAtom):
            return all(set(self[attr]).issuperset(other[attr])
                       for attr in ('isotope', 'charge', 'multiplicity', 'hybridization', 'neighbors') if self[attr]) \
                   and (self.element == ('A',) or set(self.element).issuperset(other.element))
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

    def __iter__(self):
        return iter(self._atom)

    @staticmethod
    def _element_check(x):
        if x == ('A',):
            return x
        elif isinstance(x, str):
            if x in elements_classes:
                return x,
        elif 0 < len(x) == len(set(x)) and all(x != 'A' and x in elements_classes for x in x):
            return tuple(sorted(x, key=elements_numbers.get))
        raise ValueError('invalid element')

    @staticmethod
    def _isotope_check(x):
        if x is None:
            return ()
        elif isinstance(x, int):
            if x > 0:
                return x,
        elif len(x) == len(set(x)) and all(isinstance(x, int) and x > 0 for x in x):
            return tuple(sorted(x))
        raise ValueError('invalid isotope')

    def _charge_check(self, x):
        if isinstance(x, int):
            if -3 <= x <= 3:
                return x,
        elif x and all(isinstance(x, int) and -3 <= x <= 3 for x in x):
            if not self._skip_checks:
                if len(x) != len(set(x)):
                    raise ValueError('duplicates found')
                return tuple(sorted(x))
            return tuple(x)
        raise ValueError('invalid charge')

    def _multiplicity_check(self, x):
        if x is None:
            return ()
        elif isinstance(x, int):
            if 1 <= x <= 3:
                return x,
        elif all(isinstance(x, int) and 1 <= x <= 3 for x in x):
            if not self._skip_checks:
                if len(x) != len(set(x)):
                    raise ValueError('duplicates found')
                return tuple(sorted(x))
            return tuple(x)
        raise ValueError('invalid multiplicity')

    def _neighbors_check(self, x):
        if isinstance(x, int):
            if 0 <= x <= 998:
                return x,
        elif all(isinstance(x, int) and 0 <= x <= 998 for x in x):
            if not self._skip_checks:
                if len(x) != len(set(x)):
                    raise ValueError('duplicates found')
                return tuple(sorted(x))
            return tuple(x)
        raise ValueError('invalid neighbors')

    def _hybridization_check(self, x):
        if x in (1, 2, 3, 4):
            return x,
        elif all(x in (1, 2, 3, 4) for x in x):
            if not self._skip_checks:
                if len(x) != len(set(x)):
                    raise ValueError('duplicates found')
                return tuple(sorted(x))
            return tuple(x)
        raise ValueError('invalid hybridization')

    def __getstate__(self):
        return {'checks': self._skip_checks, 'atom': self._atom}

    def __setstate__(self, state):
        super().__setattr__('_skip_checks', state['checks'])
        super().__setattr__('_atom', state['atom'])


class QueryBond(BondAttribute):
    def __init__(self, *, skip_checks=False):
        super(BondAttribute, self).__setattr__('_skip_checks', skip_checks)
        super(BondAttribute, self).__setattr__('order', (1,))
        super(BondAttribute, self).__setattr__('stereo', None)

    def _update(self, value, kwargs):
        if isinstance(value, QueryBond):
            kwargs = self._check_kwargs(kwargs)
            super(BondAttribute, self).__setattr__('order', value.order)
            super(BondAttribute, self).__setattr__('stereo', value.stereo)
            for k, v in kwargs.items():
                super(BondAttribute, self).__setattr__(k, v)
        elif isinstance(value, Bond):
            kwargs = self._check_kwargs(kwargs)
            super(BondAttribute, self).__setattr__('order', (value.order,))
            super(BondAttribute, self).__setattr__('stereo', value.stereo)
            for k, v in kwargs.items():
                super(BondAttribute, self).__setattr__(k, v)
        else:
            self._update_kwargs(value, kwargs)

    def stringify(self, stereo=True):
        order = '<%s>' % ''.join(self._order_str[x] for x in self.order) \
            if len(self.order) > 1 else self._order_str[self.order[0]]

        if stereo and self.stereo:
            return order + self._stereo_str[self.stereo]
        return order

    def weight(self, stereo=False):
        if stereo:
            return self.order, self.stereo or 0
        return self.order

    def __eq__(self, other):
        """
        == equality checks. if stereo mark is presented in query, stereo also will be compared
        """
        if isinstance(other, QueryBond):
            return set(self.order).issuperset(other.order)
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

    def _order_check(self, x):
        if x in (1, 2, 3, 4, 9):
            return x,
        elif x and all(x in (1, 2, 3, 4, 9) for x in x):
            if not self._skip_checks:
                if len(x) != len(set(x)):
                    raise ValueError('duplicates found')
                return tuple(sorted(x))
            return tuple(x)
        raise ValueError('invalid order')


__all__ = ['QueryAtom', 'QueryBond']
