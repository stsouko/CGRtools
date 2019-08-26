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
from .element import Element, FrozenDict
from .groups import GroupXV
from .periods import PeriodII, PeriodIII, PeriodIV, PeriodV, PeriodVI, PeriodVII


class N(Element, PeriodII, GroupXV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 7

    @property
    def isotopes_distribution(self):
        return FrozenDict({14: 0.99632, 15: 0.00368})

    @property
    def isotopes_masses(self):
        return FrozenDict({14: 14.003074, 15: 15.000109})

    @property
    def _common_valences(self):
        return 3,

    @property
    def _valences_exceptions(self):
        return ((-1, False, 2, ()), (1, False, 4, ()),
                (0, False, 1, ((2, 'O'), (1, 'C'), (1, 'C'))),  # O=[NH](C)C
                (0, False, 1, ((2, 'O'), (1, 'C'), (1, 'O'))),  # O=[NH](C)O
                (0, False, 1, ((2, 'O'), (1, 'C'), (1, 'N'))),  # O=[NH](C)N
                (0, False, 1, ((2, 'O'), (2, 'C'))),  # O=[NH]=C
                (0, False, 0, ((2, 'O'), (2, 'O'))),  # NO2
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'O'))),  # O-NO2
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'N'))),  # N-NO2
                (0, False, 0, ((2, 'O'), (2, 'O'), (1, 'C'))),  # C-NO2
                (0, False, 0, ((2, 'O'), (2, 'N'), (1, 'C'))),  # C-N(=O)=N
                (0, False, 0, ((2, 'C'), (2, 'O'), (1, 'O'))),  # O-N(=O)=C
                (0, False, 0, ((2, 'C'), (2, 'O'), (1, 'C'))),  # C-N(=O)=C
                (0, True, 0, ((2, 'O'),)),  # *NO
                (0, False, 0, ((2, 'N'), (2, 'N'))),  # N=N=N
                (0, False, 0, ((1, 'N'), (3, 'N'))),  # N-N#N
                (0, False, 0, ((2, 'C'), (3, 'N'))),  # C=N#N
                (0, False, 0, ((2, 'O'), (3, 'C'))),  # C#N=O
                (0, False, 0, ((2, 'N'), (3, 'N'))))  # N=N#N


class P(Element, PeriodIII, GroupXV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 15

    @property
    def isotopes_distribution(self):
        return FrozenDict({31: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({31: 30.973762})

    @property
    def _common_valences(self):
        return 3, 5

    @property
    def _valences_exceptions(self):
        return ((-1, False, 2, ()), (1, False, 4, (),),
                (-1, False, 0, ((1, 'C'), (1, 'C'), (1, 'C'), (1, 'C'))))  # Phosphonium ylide


class As(Element, PeriodIV, GroupXV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 33

    @property
    def isotopes_distribution(self):
        return FrozenDict({75: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({75: 74.921596})

    @property
    def _common_valences(self):
        return 0, 3, 5

    @property
    def _valences_exceptions(self):
        return (1, False, 4, ()),


class Sb(Element, PeriodV, GroupXV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 51

    @property
    def isotopes_distribution(self):
        return FrozenDict({121: 0.5721, 123: 0.4279})

    @property
    def isotopes_masses(self):
        return FrozenDict({121: 120.903818, 123: 122.904216})

    @property
    def _common_valences(self):
        return 0, 3, 5

    @property
    def _valences_exceptions(self):
        return (1, False, 4, ()),


class Bi(Element, PeriodVI, GroupXV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 83

    @property
    def isotopes_distribution(self):
        return FrozenDict({209: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({209: 208.980383})

    @property
    def _common_valences(self):
        return 0, 3

    @property
    def _valences_exceptions(self):
        return ((0, False, 0, ((1, 'Cl'),)),
                (0, False, 0, ((1, 'Br'),)),

                (0, False, 0, ((1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((1, 'Br'), (1, 'Br'))),
                (0, False, 0, ((1, 'I'), (1, 'I'))),
                (0, False, 0, ((1, 'S'), (1, 'S'))),
                (0, False, 0, ((2, 'S'),)),
                (0, False, 0, ((1, 'Se'), (1, 'Se'))),
                (0, False, 0, ((2, 'Se'),)),

                (0, False, 0, ((1, 'Cl'), (1, 'Cl'), (1, 'Cl'), (1, 'Cl'))),
                (0, False, 0, ((2, 'O'), (2, 'O'))),

                (0, False, 0, ((1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'), (1, 'F'))),
                (0, False, 0, ((1, 'O'), (2, 'O'), (2, 'O'))),
                (0, False, 0, ((1, 'O'), (1, 'O'), (1, 'O'), (2, 'O'))))


class Mc(Element, PeriodVII, GroupXV):
    __slots__ = ()

    @property
    def atomic_number(self):
        return 115

    @property
    def isotopes_distribution(self):
        return FrozenDict({289: 1.0})

    @property
    def isotopes_masses(self):
        return FrozenDict({289: 289.0})

    @property
    def _common_valences(self):
        return 0,

    @property
    def _valences_exceptions(self):
        return ()


__all__ = ['N', 'P', 'As', 'Sb', 'Bi', 'Mc']
