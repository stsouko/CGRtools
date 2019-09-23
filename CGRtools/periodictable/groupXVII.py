# -*- coding: utf-8 -*-
#
#  Copyright 2019 Ramil Nugmanov <stsouko@live.ru>
#  Copyright 2019 Alexander Nikanshin <17071996sasha@gmail.com>
#  Copyright 2019 Tagir Akhmetshin <tagirshin@gmail.com>
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
from CachedMethods import FrozenDict
from .element import Element
from .groups import GroupXVII
from .periods import PeriodII, PeriodIII, PeriodIV, PeriodV, PeriodVI, PeriodVII


class F(Element, PeriodII, GroupXVII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 9

    @property
    def isotopes_distribution(self):
        return FrozenDict({19: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({19: 18.998403})

    @property
    def _common_valences(self):
        return 1,

    @property
    def _valences_exceptions(self):
        return (-1, False, 0, ()),


class Cl(Element, PeriodIII, GroupXVII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 17

    @property
    def isotopes_distribution(self):
        return FrozenDict({35: 0.7578, 37: 0.2422})

    @property
    def isotopes_masses(self):
        return FrozenDict({35: 34.968853, 37: 36.965903})

    @property
    def _common_valences(self):
        return 1,

    @property
    def _valences_exceptions(self):
        return ((-1, False, 0, ()),
                (-1, False, 0, ((1, 'Cl'), (1, 'I'))),  # [I-Cl-Cl]-

                (0, False, 0, ((1, 'O'), (2, 'O'))),  # HClO2
                (0, False, 0, ((1, 'O'), (2, 'O'), (2, 'O'))),  # HClO3
                (0, False, 0, ((1, 'O'), (2, 'O'), (2, 'O'), (2, 'O'))),  # HClO4

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'))),  # ClF3

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),  # ClF5
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (2, 'O'))),  # ClOF3
                (0, False, 0, ((1, 'F'), (2, 'O'), (2, 'O'))))  # ClO2F


class Br(Element, PeriodIV, GroupXVII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 35

    @property
    def isotopes_distribution(self):
        return FrozenDict({79: 0.5069, 81: 0.4931})

    @property
    def isotopes_masses(self):
        return FrozenDict({79: 78.918338, 81: 80.916291})

    @property
    def _common_valences(self):
        return 1,

    @property
    def _valences_exceptions(self):
        return ((-1, False, 0, ()),
                (-1, False, 0, ((1, 'Br'), (1, 'I'))),  # [I-Br-Br]-

                (0, False, 0, ((1, 'O'), (2, 'O'))),  # HBrO2
                (0, False, 0, ((1, 'O'), (2, 'O'), (2, 'O'))),  # HBrO3
                (0, False, 0, ((1, 'O'), (2, 'O'), (2, 'O'), (2, 'O'))),  # HBrO4

                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'))),  # Br(OX)3
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'))),  # BrF3

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),  # BrF5
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (2, 'O'))),  # BrOF3
                (0, False, 0, ((1, 'F'), (2, 'O'), (2, 'O'))),  # BrO2F

                (0, False, 0, ((1, 'F'), (2, 'O'), (2, 'O'), (2, 'O'))))  # BrO3F


class I(Element, PeriodV, GroupXVII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 53

    @property
    def isotopes_distribution(self):
        return FrozenDict({127: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({127: 126.904468})

    @property
    def _common_valences(self):
        return 1,

    @property
    def _valences_exceptions(self):
        return ((-1, False, 0, ()),
                (-1, False, 0, ((1, 'I'), (1, 'I'))),  # [I-I-I]-

                (0, False, 0, ((1, 'O'), (2, 'O'))),  # HIO2
                (0, False, 0, ((1, 'O'), (2, 'O'), (2, 'O'))),  # HIO3
                (0, False, 0, ((1, 'O'), (2, 'O'), (2, 'O'), (2, 'O'))),  # HIO4
                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (2, 'O'), (2, 'O'))),  # H3IO5
                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (2, 'O'))),  # H5IO6

                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'))),  # I(OX)3
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'))),  # IHal3
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'))),

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),  # IF5
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (2, 'O'))),  # IOF3
                (0, False, 0, ((1, 'F'), (2, 'O'), (2, 'O'))),  # IO2F

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),  # IF7
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (2, 'O'))),  # IOF5
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (2, 'O'), (2, 'O'))),  # IO2F3
                (0, False, 0, ((1, 'F'), (2, 'O'), (2, 'O'), (2, 'O'))))  # IO3F


class At(Element, PeriodVI, GroupXVII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 85

    @property
    def isotopes_distribution(self):
        return FrozenDict({210: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({210: 209.987155})

    @property
    def _common_valences(self):
        return 0, 1

    @property
    def _valences_exceptions(self):
        return ((1, False, 0, ()), (-1, False, 0, ()),
                (0, False, 0, ((1, 'O'), (2, 'O'), (2, 'O'))))


class Ts(Element, PeriodVII, GroupXVII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 117

    @property
    def isotopes_distribution(self):
        return FrozenDict({293: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({293: 293.0})

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()


__all__ = ['F', 'Cl', 'Br', 'I', 'At', 'Ts']
