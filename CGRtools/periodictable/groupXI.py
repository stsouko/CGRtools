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
from .groups import GroupXI
from .periods import PeriodIV, PeriodV, PeriodVI, PeriodVII


class Cu(Element, PeriodIV, GroupXI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 29

    @property
    def isotopes_distribution(self):
        return FrozenDict({63: 0.6917, 65: 0.3083})

    @property
    def isotopes_masses(self):
        return FrozenDict({63: 62.929601, 65: 64.927794})

    @property
    def _common_valences(self):
        return 0, 1, 2

    @property
    def _valences_exceptions(self):
        return ((2, False, 0, ()),
                (-1, False, 0, ((1, 'Cl'), (1, 'Cl'))),  # CuCl2^-
                (-3, False, 0, ((1, 'S'), (1, 'S'))))  # CuS2^3- - это характерный комплекс для одновалентной меди


class Ag(Element, PeriodV, GroupXI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 47

    @property
    def isotopes_distribution(self):
        return FrozenDict({107: 0.51839, 109: 0.48161})

    @property
    def isotopes_masses(self):
        return FrozenDict({107: 106.905093, 109: 108.904756})

    @property
    def _common_valences(self):
        return 0, 1

    @property
    def _valences_exceptions(self):
        return ((-1, False, 0, ((1, 'Cl'), (1, 'Cl'))),  # AgCl2^1-
                (-1, False, 0, ((1, 'O'), (1, 'O'))),  # Ag(OH)2^1-
                (-1, False, 0, ((1, 'S'), (1, 'S'))),  # AgS2^1-
                (-1, False, 0, ((1, 'C'), (1, 'C'))))  # Ag(CN)2^1-


class Au(Element, PeriodVI, GroupXI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 79

    @property
    def isotopes_distribution(self):
        return FrozenDict({197: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({197: 196.966552})

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (-1, False, 0, ((1, 'Cl'), (1, 'Cl'))),  # AuCl2^1-
                (-1, False, 0, ((1, 'S'), (1, 'S'))),  # AuS2^1-
                (-1, False, 0, ((1, 'C'), (1, 'C'))),  # Au(CN)^1-
                (-1, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),  # Au(OH)4^1-
                (-1, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'))),  # Au(CN)4^1-
                (-1, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))))  # AuCl4^1-


class Rg(Element, PeriodVII, GroupXI):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 111

    @property
    def isotopes_distribution(self):
        return FrozenDict({282: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({282: 282.169127})

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()


__all__ = ['Cu', 'Ag', 'Au', 'Rg']
