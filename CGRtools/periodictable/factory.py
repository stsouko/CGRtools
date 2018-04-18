# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
from itertools import chain
from os.path import splitext
from .data import (groups, periods, table_e, table_i, types, electrons, atom_valences_exceptions, atom_valences,
                   atom_charge_radical, atom_implicit_h)


class Element:
    def __init__(self, charge=0, radical=0, isotope=None):
        assert (self.symbol, charge, radical) in atom_charge_radical, 'Invalid charge or number of unpaired electrons'
        self.__isotope = isotope
        self.__charge = charge
        self.__radical = radical

    @property
    def isotope(self):
        return self.__isotope or self.common_isotope

    @property
    def charge(self):
        return self.__charge

    @property
    def radical(self):
        return self.__radical

    def __repr__(self):
        r = ', %d' % self.__radical if self.__radical else ''
        i = ', %d' % self.__isotope if self.__isotope else ''
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
        def __new__(mcls, cls_name, *args, **kwargs):
            return super().__new__(mcls, 'Period{}'.format(roma), *args, **kwargs)

        def __repr__(mcls):
            return "<class '{}.{}'>".format(splitext(__name__)[0], mcls.__name__)

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
        def __new__(mcls, cls_name, *args, **kwargs):
            return super().__new__(mcls, 'Group{}'.format(roma), *args, **kwargs)

        def __repr__(mcls):
            return "<class '{}.{}'>".format(splitext(__name__)[0], mcls.__name__)

    class _Group(metaclass=GroupType):
        @property
        def number(self):
            return roma

        def __repr__(self):
            return '{}()'.format(self.__class__.__name__)

    return _Group


def get_element(symbol, number, isotope, _electrons, _type, group, period):
    class ElementType(type(group), type(period)):
        def __new__(mcls, cls_name, *args, **kwargs):
            return type.__new__(mcls, symbol, *args, **kwargs)

    class _Element(Element, group, period, metaclass=ElementType):
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
        def electrons(self):
            return _electrons

        @property
        def type(self):
            return _type

        def check_valence(self, bonds, neighbors):
            scrl = (symbol, self.charge, self.radical, len(bonds))
            if scrl in atom_valences_exceptions:
                return sorted(zip(neighbors, bonds)) in atom_valences_exceptions[scrl]

            return (symbol, self.charge, self.radical, self.__bonds_sum(bonds)) in atom_valences

        def get_implicit_h(self, bonds):
            return atom_implicit_h.get((symbol, self.charge, self.radical, self.__bonds_sum(bonds)), 0)

        @staticmethod
        def __bonds_sum(bonds):
            return int(sum(_bonds[x] for x in bonds))

    return _Element


_bonds = {1: 1, 2: 2, 3: 4, 4: 1.5, 9: 1}
_group_cache = {}
_period_cache = {}
isotopes = dict(zip(table_e, table_i))

_groups = {z: x for x, y in enumerate(groups, start=1) for z in y}
_periods = {z: x for x, y in enumerate(periods, start=1) for z in y}
_types = {z: x for x, y in types.items() for z in y}

elements = {}

for n, s in enumerate(table_e, start=1):
    _g = _groups[n]
    _p = _periods[n]
    g = _group_cache.get(_g) or _group_cache.setdefault(_g, get_group(_g))
    p = _period_cache.get(_p) or _period_cache.setdefault(_p, get_period(_p))
    e = next(n - x for x in reversed(electrons) if n > x)
    locals()[s] = elements[s] = get_element(s, n, isotopes[s], e, _types[n], g, p)


for x in _group_cache.values():
    locals()[x.__name__] = x

for x in _period_cache.values():
    locals()[x.__name__] = x


__all__ = list(table_e) + [x.__name__ for x in chain(_group_cache.values(), _period_cache.values())] + \
          [Element.__name__, 'elements', 'isotopes']
