# -*- coding: utf-8 -*-
#
#  Copyright 2019 Ramil Nugmanov <stsouko@live.ru>
#  Copyright 2019 Tagir Akhmetshin <tagirshin@gmail.com>
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
from abc import ABC, abstractmethod
from collections import Mapping
from itertools import permutations
from typing import Optional, Tuple, Dict


class propertied_exceptions:
    def __init__(self, func):
        self.__doc__ = getattr(func, "__doc__")
        self.func = func

    def __get__(self, obj, cls):
        if obj is None:
            return
        try:
            exceptions = cls._cached_exceptions[obj.atomic_number]
        except KeyError:
            exceptions = cls._cached_exceptions[obj.atomic_number] = self.func(obj)

        obj.__dict__[self.func.__name__] = exceptions
        return exceptions


class Element(ABC):
    __slots__ = ('__isotope', '__charge', '__is_radical')
    _cached_exceptions = {}

    def __init__(self, charge: int = 0, isotope: Optional[int] = None, is_radical: bool = False):
        """
        element object with specified charge, isotope and multiplicity

        :param charge: formal charge of atom
        :param isotope: isotope number of element
        :param is_radical: multiplicity of atom
        """
        if isinstance(charge, int) and isinstance(is_radical, bool):
            if not charge or any(charge == c and is_radical == r for c, r, _ in self._valences_exceptions):
                self.__charge = charge
                self.__is_radical = is_radical
            else:
                raise ValueError('invalid charge and radical combination')
        else:
            raise TypeError('integer charge value and bool radical value required')

        if isotope is None:
            self.__isotope = None
        elif isinstance(isotope, int):
            if isotope in self.isotopes_distribution:
                self.__isotope = isotope
            else:
                raise ValueError('isotope number impossible or not stable')
        else:
            raise TypeError('integer isotope number required')

    @propertied_exceptions
    def all_exceptions(self):
        elements_classes = {x.__name__: x for x in Element.__subclasses__()}
        elements_numbers = {}
        exceptions = set()
        for charge, multiplicity, env in self.valences_exceptions:
            for x in permutations(env):
                tmp = []
                for bond, symbol in x:
                    try:
                        number = elements_numbers[symbol]
                    except KeyError:
                        number = elements_numbers[symbol] = elements_classes[symbol]().atomic_number
                    tmp.append((bond, number))
                exceptions.add((charge, multiplicity, tuple(tmp)))
        return frozenset(exceptions)

    @property
    def charge(self):
        return self.__charge

    @property
    def isotope(self):
        return self.__isotope

    @property
    def is_radical(self):
        return self.__is_radical

    @property
    def atomic_mass(self):
        mass = self.isotopes_masses
        if self.__isotope is None:
            return sum(x * mass[i] for i, x in self.isotopes_distribution.items())
        return mass[self.__isotope]

    def __eq__(self, other):
        return isinstance(other, Element) and self.__charge == other.charge and self.__isotope == other.isotope \
               and self.__is_radical == other.__is_radical

    def __hash__(self):
        """
        7bit | 9bit | 3bit | 1bit
        """
        return self.atomic_number << 13 | (self.__isotope or 0) << 4 | self.__charge + 3 << 1 | self.__is_radical

    @property
    @abstractmethod
    def atomic_number(self) -> int:
        """
        element number
        """

    @property
    @abstractmethod
    def isotopes_distribution(self) -> Dict[int, float]:
        """
        isotopes distribution in earth
        """

    @property
    @abstractmethod
    def isotopes_masses(self) -> Dict[int, float]:
        """
        isotopes distribution in earth
        """

    @property
    @abstractmethod
    def _common_valences(self) -> Tuple[int, ...]:
        """
        common valences of element
        """

    @property
    @abstractmethod
    def _valences_exceptions(self) -> Tuple[Tuple[int, bool, int, Tuple[Tuple[int, str], ...]]]:
        """
        exceptions in charges, radical state, implicit H count, and non H neighbors of element
        examples:
        (-1, False, 1, ()) - anion, not radical, has 1 implicit hydrogen if explicit atom not exists: [OH]- or C[O-]
        (0, True, 1, ()) - neutral, radical, has 1 implicit hydrogen if explicit atom not exists: [OH]* or C[O*]
        number of free electrons calculated as diff of default valence and number of connected atoms (include implicit)
        (0, False, 1, ((1, 'C'),)) - neutral, not radical, has 1 implicit hydrogen and 1 single bonded carbon:
        CO - alcohol or e.g. COC - ether
        (0, False, 0, ((1, 'O'), (2, 'O'))) - can be anion/cation or neutral: HNO2 or [NO2]-.
        state of neighbors atoms don't take into account. order of neighbors atoms don't take into account.
        use both: (0, False, 0, ((2, 'O'),)) and (0, False, 0, ((1, 'O'), (1, 'O'))) if chains possible
        use charge transfer for carbonyles, cyanides etc:
        (-1, False, 0, ((1, 'C'),)) - [M-]-C#[O+]
        user both for cyanates and isocyanates etc complexes:
        (-1, False, 0, ((1, 'O'),)) and (-1, False, 0, ((1, 'N'),))
        """


class FrozenDict(Mapping):
    __slots__ = '__d'

    def __init__(self, *args, **kwargs):
        self.__d = dict(*args, **kwargs)

    def __iter__(self):
        return iter(self.__d)

    def __len__(self):
        return len(self.__d)

    def __getitem__(self, key):
        return self.__d[key]
