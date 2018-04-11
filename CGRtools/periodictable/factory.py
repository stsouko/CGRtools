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
from .data import groups, periods, table_e, table_i, types, electrons


class Element:
    def __init__(self, isotope=None):
        self.__isotope = isotope

    @property
    def isotope(self):
        return self.__isotope or self.common_isotope

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, self.isotope or '')


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
            return "<class '{}.{}'>".format(__name__, mcls.__name__)

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
            return "<class '{}.{}'>".format(__name__, mcls.__name__)

    class _Group(metaclass=GroupType):
        @property
        def number(self):
            return roma

        def __repr__(self):
            return '{}()'.format(self.__class__.__name__)

    return _Group


def get_element(symbol, number, isotope, electrons, _type, group, period):
    class ElementType(type(group), type(period)):
        def __new__(mcls, cls_name, *args, **kwargs):
            return type.__new__(mcls, symbol, *args, **kwargs)

    class _Element(Element, group, period, metaclass=ElementType):
        @property
        def number(self):
            return number

        @property
        def common_isotope(self):
            return isotope

        @property
        def electrons(self):
            return electrons

        @property
        def type(self):
            return _type

    return _Element


_group_cache = {}
_period_cache = {}
_isotopes = dict(zip(table_e, table_i))

_groups = {z: x for x, y in enumerate(groups, start=1) for z in y}
_periods = {z: x for x, y in enumerate(periods, start=1) for z in y}
_types = {z: x for x, y in types.items() for z in y}

for n, s in enumerate(table_e, start=1):
    _g = _groups[n]
    _p = _periods[n]
    g = _group_cache.get(_g) or _group_cache.setdefault(_g, get_group(_g))
    p = _period_cache.get(_p) or _period_cache.setdefault(_p, get_period(_p))
    e = next(n - x for x in reversed(electrons) if n > x)
    locals()[s] = get_element(s, n, _isotopes[s], e, _types[n], g, p)


for x in _group_cache.values():
    locals()[x.__name__] = x

for x in _period_cache.values():
    locals()[x.__name__] = x


__all__ = list(table_e) + [x.__name__ for x in chain(_group_cache.values(),
                                                     _period_cache.values())] + [Element.__name__]
