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
from .molecule import Atom, AtomAttribute
from ..cache import cached_property
from ..periodictable import Element, elements_classes


class QueryAtom(AtomAttribute):
    __slots__ = '_atom'

    def __init__(self, *, skip_checks=False):
        super().__setattr__('_skip_checks', skip_checks)
        super().__setattr__('_atom', {'element': None, 'isotope': None, 'charge': 0, 'multiplicity': None,
                                      'neighbors': (), 'hybridization': (), 'stereo': None,
                                      'x': 0., 'y': 0., 'z': 0.})

    def __init_copy__(self, parent):
        super().__setattr__('_skip_checks', parent._skip_checks)
        super().__setattr__('_atom', parent._atom.copy())

    def __setattr__(self, key, value):
        if not self._skip_checks:
            value = getattr(self, f'_{key}_check')(value)
        self._atom[key] = value
        self.__dict__.clear()

    def _update(self, value, kwargs):
        if isinstance(value, QueryAtom):
            kwargs = self._check_kwargs(kwargs)
            self._atom.update(value)
            self._atom.update(kwargs)
        elif isinstance(value, Atom):
            kwargs = self._check_kwargs(kwargs)
            self._atom.update(element=(value.element,), charge=value.charge, stereo=value.stereo,
                              multiplicity=value.multiplicity,
                              isotope=value.isotope if value.isotope != value.common_isotope else None,
                              neighbors=(value.neighbors,) if value.neighbors is not None else (),
                              hybridization=(value.hybridization,) if value.hybridization else (),
                              x=value.x, y=value.y, z=value.z)
            self._atom.update(kwargs)
        elif isinstance(value, type):
            if not issubclass(value, Element):
                raise ValueError('only CGRtools.periodictable.Element subclasses allowed')
            if not {'charge', 'multiplicity', 'isotope', 'element'}.isdisjoint(kwargs):
                raise KeyError('charge, multiplicity, isotope and element override not allowed')
            kwargs = self._check_kwargs(kwargs)
            self._atom.update(element=(value.symbol,) if value.symbol != 'A' else None,
                              charge=0, multiplicity=None, isotope=None)
            self._atom.update(kwargs)
        elif isinstance(value, Element):
            if not {'charge', 'multiplicity', 'isotope', 'element'}.isdisjoint(kwargs):
                raise KeyError('charge, multiplicity, isotope and element override not allowed')
            kwargs = self._check_kwargs(kwargs)
            self._atom.update(element=(value.symbol,) if value.symbol != 'A' else None,
                              charge=value.charge, multiplicity=value.multiplicity,
                              isotope=value.isotope if value.isotope != value.common_isotope else None)
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
        self.__dict__.clear()

    def __eq__(self, other):
        if isinstance(other, Atom):
            return (self.charge == other.charge and
                    (other.element in self.element_set if self.element else True) and
                    (self.isotope == other.isotope if self.isotope else True) and
                    (self.multiplicity == other.multiplicity if self.multiplicity else True) and
                    (other.neighbors in self.neighbors if self.neighbors else True) and
                    (other.hybridization in self.hybridization if self.hybridization else True))
        elif isinstance(other, QueryAtom):
            if self.element:
                if not other.element:
                    return False
                elif not self.element_set.issuperset(other.element_set):
                    return False
            if self.neighbors:
                if not other.neighbors:
                    return False
                elif not set(self.neighbors).issuperset(other.neighbors):
                    return False
            if self.hybridization:
                if not other.hybridization:
                    return False
                elif not set(self.hybridization).issuperset(other.hybridization):
                    return False

            return (self.charge == other.charge and
                    (self.isotope == other.isotope if self.isotope else True) and
                    (self.multiplicity == other.multiplicity if self.multiplicity else True))
        return False

    def __getitem__(self, key):
        """
        dict like access to atom's attrs
        """
        return self._atom[key]

    def __getattr__(self, key):
        try:
            value = self.__dict__[key] = self._atom[key]
        except KeyError as e:
            raise AttributeError from e
        return value

    def __iter__(self):
        return iter(self._atom)

    @cached_property
    def element_set(self):
        return set(self.element)

    @staticmethod
    def _element_check(x):
        if x is None:
            return
        elif isinstance(x, str):
            if x == 'A':
                return
            elif x in elements_classes:
                return x,
        else:
            y = tuple(x)
            if y == ('A',):
                return
            elif y and all(x != 'A' and x in elements_classes for x in y):
                return y

        raise ValueError(f'invalid element: {x}')

    def _neighbors_check(self, x):
        if x is None or x == ():
            return ()
        elif isinstance(x, int):
            if 0 <= x <= 998:
                return x,
        elif all(isinstance(x, int) and 0 <= x <= 998 for x in x):
            if not self._skip_checks and len(x) != len(set(x)):
                raise ValueError('duplicates found')
            return tuple(x)
        raise ValueError('invalid neighbors')

    def _hybridization_check(self, x):
        if x is None or x == ():
            return ()
        elif x in (1, 2, 3, 4):
            return x,
        elif all(x in (1, 2, 3, 4) for x in x):
            if not self._skip_checks and len(x) != len(set(x)):
                raise ValueError('duplicates found')
            return tuple(x)
        raise ValueError('invalid hybridization')

    def __getstate__(self):
        return {'checks': self._skip_checks, 'atom': self._atom}

    def __setstate__(self, state):
        super().__setattr__('_skip_checks', state['checks'])
        super().__setattr__('_atom', state['atom'])


__all__ = ['QueryAtom']
