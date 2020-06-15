# -*- coding: utf-8 -*-
#
#  Copyright 2019, 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Dayana Bashirova <dayana.bashirova@yandex.ru>
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
from .groups import GroupXII
from .periods import PeriodIV, PeriodV, PeriodVI, PeriodVII


class Zn(Element, PeriodIV, GroupXII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 30

    @property
    def isotopes_distribution(self):
        return FrozenDict({64: 0.4863, 66: 0.279, 67: 0.041, 68: 0.1875, 70: 0.0062, 62: 0., 69: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({64: 63.929147, 66: 65.926037, 67: 66.927131, 68: 67.924848, 70: 69.925325, 62: 61.934330,
                           69: 68.926550})

    @property
    def _common_valences(self):
        return 0, 2

    @property
    def _valences_exceptions(self):
        return ((2, False, 0, ()), (1, False, 0, ((1, 'C'),)),
                (-2, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),  # Zn[(OH)4]2-
                (-2, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'))),  # Zn[(CN)4]2-
                (-2, False, 0, ((1, 'N'), (1, 'N'), (1, 'N'), (1, 'N'))),  # Zn[(NCS)4]2-
                (-2, False, 0, ((1, 'S'), (1, 'S'), (1, 'S'), (1, 'S'))))  # Zn[(SCN)4]2-

    @property
    def atomic_radius(self):
        return 1.42


class Cd(Element, PeriodV, GroupXII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 48

    @property
    def isotopes_distribution(self):
        return FrozenDict({106: 0.0125, 108: 0.0089, 110: 0.1249, 111: 0.128, 112: 0.2413, 113: 0.1222, 114: 0.2873,
                           116: 0.0749})

    @property
    def isotopes_masses(self):
        return FrozenDict({106: 105.906458, 108: 107.904183, 110: 109.903006, 111: 110.904182, 112: 111.902757,
                           113: 112.904401, 114: 113.903358, 116: 115.904755})

    @property
    def _common_valences(self):
        return 0, 2

    @property
    def _valences_exceptions(self):
        return ((2, False, 0, ()),
                (-2, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),  # Cd[(OH)4]2-
                (-2, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'))),  # Cd[(CN)4]2-
                (-2, False, 0, ((1, 'N'), (1, 'N'), (1, 'N'), (1, 'N'))),  # Cd[(NCS)4]2-
                (-2, False, 0, ((1, 'S'), (1, 'S'), (1, 'S'), (1, 'S'))))  # Cd[(SCN)4]2-

    @property
    def atomic_radius(self):
        return 1.61


class Hg(Element, PeriodVI, GroupXII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 80

    @property
    def isotopes_distribution(self):
        return FrozenDict({196: 0.0015, 198: 0.0997, 199: 0.1687, 200: 0.231, 201: 0.1318, 202: 0.2986, 204: 0.0687,
                           197: 0., 203: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({196: 195.965815, 198: 197.966752, 199: 198.968262, 200: 199.968309, 201: 200.970285,
                           202: 201.970626, 204: 203.973476, 197: 196.967213, 203: 202.972873})

    @property
    def _common_valences(self):
        return 0, 2

    @property
    def _valences_exceptions(self):
        return ((-2, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'))),  # Hg[(CN)4]2-
                (-2, False, 0, ((1, 'N'), (1, 'N'), (1, 'N'), (1, 'N'))),  # Hg[(NCS)4]2-
                (-2, False, 0, ((1, 'S'), (1, 'S'), (1, 'S'), (1, 'S'))))  # Hg[(SCN)4]2-

    @property
    def atomic_radius(self):
        return 1.71


class Cn(Element, PeriodVII, GroupXII):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 112

    @property
    def isotopes_distribution(self):
        return FrozenDict({285: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({285: 285.177444})

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 1.71  # unknown, taken radius of previous element in group


__all__ = ['Zn', 'Cd', 'Hg', 'Cn']
