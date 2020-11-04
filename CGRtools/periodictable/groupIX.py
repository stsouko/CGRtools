# -*- coding: utf-8 -*-
#
#  Copyright 2019, 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Tagir Akhmetshin <tagirshin@gmail.com>
#  Copyright 2019 Tansu Nasyrova <tansu.nasurova@gmail.com>
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
from .groups import GroupIX
from .periods import PeriodIV, PeriodV, PeriodVI, PeriodVII


class Co(Element, PeriodIV, GroupIX):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 27

    @property
    def isotopes_distribution(self):
        return FrozenDict({59: 1.0, 55: 0., 57: 0., 58: 0., 60: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({59: 58.9332, 55: 54.941999, 57: 56.936291, 58: 57.935753, 60: 59.933817})

    @property
    def _common_valences(self):
        return 0, 2, 3

    @property
    def _valences_exceptions(self):
        return ((2, False, 0, ()),
                (3, False, 0, ()),
                (0, False, 0, ((1, 'H'),)),
                (-3, False, 0, ((2, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),  # [CoO4]3-
                (-2, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),  # [CoF6]2-
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'), (1, 'H'))),  # HCo(CO)4

                (-1, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'))),  # [CoF3]-
                (-1, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (-1, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'))),  # [Co(NO3)3]-

                (-2, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),  # [CoF4]2-
                (-2, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (-2, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'), (1, 'Br'))),
                (-2, False, 0, ((1, 'I'), (1, 'I'), (1, 'I'), (1, 'I'))),
                (-2, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),  # [Co(OH)4]2-

                (-1, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'))),  # [Co(CN4)]-

                (0, False, 0, ((1, 'N'), (1, 'N'), (1, 'N'), (1, 'N'), (1, 'N'), (1, 'C'))),  # B12
                (-3, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'))),  # [Co(CN)6]3-
                (-3, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),  # [Co(OH)6]3-

                (-3, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),  # [CoCl5]3-
                (-3, False, 0, ((1, 'S'), (1, 'S'), (1, 'S'), (1, 'S'), (1, 'S'))),  # [Co(NCS)5]3-

                (-4, False, 0, ((1, 'S'), (1, 'S'), (1, 'S'), (1, 'S'), (1, 'S'), (1, 'S'))),  # [Co(NCS)6]4-
                (-4, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))))  # [Co(OH)6]4-

    @property
    def atomic_radius(self):
        return 1.52


class Rh(Element, PeriodV, GroupIX):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 45

    @property
    def isotopes_distribution(self):
        return FrozenDict({103: 1.0, 105: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({103: 102.905504, 105: 104.905694})

    @property
    def _common_valences(self):
        return 0, 3, 4

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ((2, 'O'),)),  # RhO
                (0, False, 0, ((1, 'O'), (1, 'O'))),  # Rh(OH)2
                (0, False, 0, ((2, 'S'),)),
                (0, False, 0, ((1, 'S'), (1, 'S'))),
                (-1, False, 0, ((1, 'Br'), (1, 'Br'), (1, 'Br'), (1, 'Br'))),  # [RhBr4]-
                (-3, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),  # [RhCl6]3-
                (-3, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'), (1, 'O'))),  # [Rh(NO2)6]3-
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'), (1, 'H'))),  # HRh(CO)4
                (0, False, 0, ((1, 'P'), (1, 'P'), (1, 'P'), (1, 'C'), (1, 'H'))))  # HRh(CO)[P(Ph)3]3

    @property
    def atomic_radius(self):
        return 1.73


class Ir(Element, PeriodVI, GroupIX):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 77

    @property
    def isotopes_distribution(self):
        return FrozenDict({191: 0.373, 193: 0.627, 192: 0.})

    @property
    def isotopes_masses(self):
        return FrozenDict({191: 190.960591, 193: 192.962924, 192: 191.962605})

    @property
    def _common_valences(self):
        return 0, 3, 4

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ((1, 'F'),)),
                (0, False, 0, ((1, 'Cl'),)),
                (0, False, 0, ((1, 'Br'),)),
                (0, False, 0, ((1, 'I'),)),
                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'S'), (1, 'S'))),
                (0, False, 0, ((2, 'S'),)),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (-3, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))))

    @property
    def atomic_radius(self):
        return 1.8


class Mt(Element, PeriodVII, GroupIX):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 109

    @property
    def isotopes_distribution(self):
        return FrozenDict({278: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({278: 278.15481})

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()

    @property
    def atomic_radius(self):
        return 1.8  # unknown, taken radius of previous element in group


__all__ = ['Co', 'Rh', 'Ir', 'Mt']
