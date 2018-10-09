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
from abc import ABC, ABCMeta, abstractmethod
from itertools import chain
from operator import ge, le, gt, lt
from os.path import splitext
from .data import (groups, periods, table_e, types, atom_valences_exceptions, atom_valences, common_isotope,
                   atom_charge_radical, atom_implicit_h, valence_electrons, electron_configuration, orbital_names)
from ..exceptions import InvalidAtom


class Element(ABC):
    def __init__(self, charge=0, multiplicity=None, isotope=None):
        self.__isotope = isotope
        self.__charge = charge
        self.__multiplicity = multiplicity

    @property
    def isotope(self):
        return self.__isotope or self.common_isotope

    @property
    def charge(self):
        return self.__charge

    @property
    def radical(self):
        return _radical_map[self.__multiplicity]

    @property
    def multiplicity(self):
        return self.__multiplicity

    @property
    @abstractmethod
    def symbol(self):
        pass

    @property
    @abstractmethod
    def common_isotope(self):
        pass

    def __repr__(self):
        r = ', multiplicity=%d' % self.__multiplicity if self.__multiplicity else ''
        i = ', isotope=%d' % self.__isotope if self.__isotope else ''
        return '{}({}{}{})'.format(self.__class__.__name__, self.__charge or '', r, i)


def arab2roman(number):
    roma = []
    for arabic, roman in [(1000, "M"), (900, "CM"), (500, "D"), (400, "CD"), (100, "C"), (90, "XC"), (50, "L"),
                          (40, "XL"), (10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I")]:
        factor, number = divmod(number, arabic)
        roma.append(roman * factor)
    return ''.join(roma)


def get_period(number):
    roma = arab2roman(number)

    class PeriodType(type):
        def __new__(mcs, cls_name, *args, **kwargs):
            return super().__new__(mcs, 'Period{}'.format(roma), *args, **kwargs)

        def __repr__(mcs):
            return "<class '{}.{}'>".format(splitext(__name__)[0], mcs.__name__)

    class _Period(metaclass=PeriodType):
        @property
        def number(self):
            return number

        def __repr__(self):
            return '{}()'.format(self.__class__.__name__)

    return _Period


def get_group(number):
    roma = arab2roman(number)

    class GroupType(type):
        def __new__(mcs, cls_name, *args, **kwargs):
            return super().__new__(mcs, 'Group{}'.format(roma), *args, **kwargs)

        def __repr__(mcs):
            return "<class '{}.{}'>".format(splitext(__name__)[0], mcs.__name__)

    class _Group(metaclass=GroupType):
        @property
        def number(self):
            return roma

        def __repr__(self):
            return '{}()'.format(self.__class__.__name__)

    return _Group


def get_element(symbol, number, _type, group, period):
    isotope = common_isotope[symbol]
    _electrons = valence_electrons[symbol]
    _configuration = ' '.join('%d%s%d' % (n, orbital_names[l], e) for (n, l), e in
                              electron_configuration[symbol].items())

    class ElementType(ABCMeta, type(group), type(period)):
        def __new__(mcs, cls_name, *args, **kwargs):
            return type.__new__(mcs, symbol, *args, **kwargs)

    class _Element(Element, group, period, metaclass=ElementType):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            if (self.symbol, self.charge, self.radical) not in atom_charge_radical:
                raise InvalidAtom('Invalid charge or number of unpaired electrons')

        @property
        def number(self):
            return number

        @property
        def symbol(self):
            return symbol

        @property
        def common_isotope(self):
            return isotope

        @property
        def electron_configuration(self):
            return _configuration

        @property
        def electrons(self):
            return _electrons

        @property
        def type(self):
            return _type

        def get_valence(self, bonds, neighbors):
            """
            get valence of atom for given environment include implicit hydrogens.
            :param bonds: bonds of atom
            :param neighbors: symbols of bonded atoms in same order as bonds
            :return: valence number or None if valence impossible
            """
            res = atom_valences.get((symbol, self.charge, self.radical, self.__bonds_sum(bonds)))
            if res is None:
                key = atom_valences_exceptions.get((symbol, self.charge, self.radical, len(bonds)))
                if key:
                    return key.get(tuple(sorted(zip(neighbors, bonds))))
            return res

        def get_implicit_h(self, bonds):
            return atom_implicit_h.get((symbol, self.charge, self.radical, self.__bonds_sum(bonds)), 0)

        @staticmethod
        def __bonds_sum(bonds):
            return int(sum(_bonds[x] for x in bonds))

        def __eq__(self, other):
            if isinstance(other, str):
                assert other in table_e, 'invalid element'
                return other == symbol
            elif isinstance(other, Element):
                return other.symbol == symbol and other.isotope == self.isotope and other.charge == self.charge and \
                       other.radical == self.radical
            return False

        def __gt__(self, other):
            return self.__ops(other, gt)

        def __ge__(self, other):
            return self.__ops(other, ge)

        def __lt__(self, other):
            return self.__ops(other, lt)

        def __le__(self, other):
            return self.__ops(other, le)

        def __ops(self, other, op):
            if isinstance(other, str):
                assert other in table_e, 'invalid element'
                return op(symbol, other)
            elif isinstance(other, Element):
                return op((symbol, self.isotope, self.charge, self.radical),
                          (other.symbol, other.isotope, other.charge, other.radical))
            raise TypeError('unorderable types {} and {}'.format(self.__class__.__name__, type(other)))

        def __hash__(self):
            return hash(symbol) ^ self.charge ^ self.isotope ^ self.radical

    return _Element


_bonds = {1: 1, 2: 2, 3: 3, 4: 1.5, 9: 1, None: 0}
_group_cache = {}
_period_cache = {}
_radical_map = {1: 2, 2: 1, 3: 2, None: 0}
_radical_unmap = {None: None, 0: None, 1: 2, 2: 3}
_groups = {z: x for x, y in enumerate(groups, start=1) for z in y}
_periods = {z: x for x, y in enumerate(periods, start=1) for z in y}
_types = {z: x for x, y in types.items() for z in y}

elements = {}

for n, s in enumerate(table_e, start=1):
    _g = _groups[n]
    _p = _periods[n]
    g = _group_cache.get(_g) or _group_cache.setdefault(_g, get_group(_g))
    p = _period_cache.get(_p) or _period_cache.setdefault(_p, get_period(_p))
    locals()[s] = elements[s] = get_element(s, n, _types[n], g, p)


for x in _group_cache.values():
    locals()[x.__name__] = x

for x in _period_cache.values():
    locals()[x.__name__] = x


__all__ = list(table_e) + [x.__name__ for x in chain(_group_cache.values(), _period_cache.values())] + \
          [Element.__name__, 'elements']
