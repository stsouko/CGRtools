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
from collections import Mapping, defaultdict
from typing import Optional, Tuple, Dict, Set, List
from weakref import ref
from ..exceptions import IsConnectedAtom, IsNotConnectedAtom, ValenceError


class class_cached_property:
    def __init__(self, func):
        self.__doc__ = getattr(func, '__doc__')
        self.func = func
        self.name = func.__name__

    def __get__(self, obj, cls):
        if obj is None:
            return self
        try:
            class_cache = cls.__class_cache__[type(obj)]
        except KeyError:
            cached = self.func(obj)
            cls.__class_cache__[type(obj)] = {self.name: cached}
        else:
            try:
                cached = class_cache[self.name]
            except KeyError:
                cached = class_cache[self.name] = self.func(obj)

        return cached


class Core(ABC):
    __slots__ = ('__isotope', '_graph', '_map')

    def __init__(self, isotope: Optional[int] = None):
        """
        element object with specified charge, isotope and multiplicity

        :param isotope: isotope number of element
        """
        if isinstance(isotope, int):
            if isotope not in self.isotopes_distribution:
                raise ValueError('isotope number impossible or not stable')
        elif isotope is not None:
            raise TypeError('integer isotope number required')
        self.__isotope = isotope

    def __repr__(self):
        return f'{self.__class__.__name__}({self.__isotope})'

    @property
    @abstractmethod
    def atomic_number(self) -> int:
        """
        element number
        """

    @property
    def isotope(self):
        return self.__isotope

    @property
    def atomic_mass(self):
        mass = self.isotopes_masses
        if self.__isotope is None:
            return sum(x * mass[i] for i, x in self.isotopes_distribution.items())
        return mass[self.__isotope]

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
    def charge(self):
        try:
            return self._graph()._charges[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def is_radical(self):
        try:
            return self._graph()._radicals[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def x(self):
        try:
            return self._graph()._plane[self._map][0]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def y(self):
        try:
            return self._graph()._plane[self._map][1]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def neighbors(self):
        try:
            return self._graph()._neighbors[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def hybridization(self):
        try:
            return self._graph()._hybridizations[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    def copy(self):
        """
        detached from graph copy of element
        """
        copy = object.__new__(self.__class__)
        copy._Core__isotope = self.__isotope
        return copy

    def _attach_to_graph(self, graph, _map):
        try:
            self._graph
        except AttributeError:
            self._graph = ref(graph)
            self._map = _map
        else:
            raise IsConnectedAtom

    def _change_map(self, _map):
        try:
            self._graph
        except AttributeError:
            raise IsNotConnectedAtom
        else:
            self._map = _map


class Element(Core):
    __slots__ = ()
    __class_cache__ = {}

    @classmethod
    def from_symbol(cls, symbol):
        """
        get Element class by its symbol
        """
        try:
            element = next(x for x in Element.__subclasses__() if x.__name__ == symbol)
        except StopIteration:
            raise ValueError(f'Element with symbol "{symbol}" not found')
        return element

    @classmethod
    def from_atomic_number(cls, number):
        """
        get Element class by its number
        """
        try:
            element = next(x for x in Element.__subclasses__() if x.atomic_number.fget(None) == number)
        except StopIteration:
            raise ValueError(f'Element with number "{number}" not found')
        return element

    def __eq__(self, other):
        return isinstance(other, Element) and self.atomic_number == other.atomic_number and \
               self.isotope == other.isotope

    def __int__(self):
        """
        21bit = 9bit | 7bit | 4bit | 1bit
        """
        return (self.isotope or 0) << 12 | self.atomic_number << 5 | self.charge + 4 << 1 | self.is_radical

    def valence_rules(self, charge: int, is_radical: bool, valence: int) -> \
            List[Tuple[Set[Tuple[int, 'Element']], Dict[Tuple[int, 'Element'], int], int]]:
        """
        valence rules for element with specific charge/radical state
        """
        try:
            return self._compiled_valence_rules[(charge, is_radical, valence)]
        except KeyError:
            raise ValenceError

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

    @class_cached_property
    def _compiled_charge_radical(self) -> Set[Tuple[int, bool]]:
        """
        exceptions in charges, radical state
        examples:
        (-1, False) - anion, not radical
        (0, True) - neutral radical
        """
        return {(c, r) for c, r, *_ in self._valences_exceptions}

    @class_cached_property
    def _compiled_valence_rules(self) -> \
            Dict[Tuple[int, bool, int], List[Tuple[Set[Tuple[int, 'Element']], Dict[Tuple[int, 'Element'], int], int]]]:
        """
        dictionary with key = (charge, is_radical, sum_of_bonds) and
        value = list of possible neighbors and implicit H count
        """
        elements_classes = {x.__name__: x for x in Element.__subclasses__()}

        rules = defaultdict(list)
        if self._common_valences[0]:  # atom has implicit hydrogens by default
            for valence in self._common_valences:
                for h in range(valence + 1):
                    rules[(0, False, valence - h)].append((set(), {}, h))  # any atoms and bonds possible
        else:
            for valence in self._common_valences:
                rules[(0, False, valence)].append((set(), {}, 0))  # any atoms and bonds possible

        for charge, is_radical, implicit, environment in self._valences_exceptions:
            explicit = sum(x for x, _ in environment)
            explicit_dict = defaultdict(int)
            explicit_set = set()
            for b, e in environment:
                be = (b, elements_classes[e])
                explicit_set.add(be)
                explicit_dict[be] += 1
            explicit_dict = dict(explicit_dict)

            if implicit:
                valence = explicit + implicit

                for h in range(implicit + 1):
                    rules[(charge, is_radical, valence - h)].append((explicit_set, explicit_dict, h))
            else:
                rules[(charge, is_radical, explicit)].append((explicit_set, explicit_dict, 0))
        return dict(rules)


class DynamicElement(Core):
    __slots__ = ('__p_charge', '__p_is_radical')

    @classmethod
    def from_symbol(cls, symbol):
        """
        get DynamicElement class by its symbol
        """
        try:
            element = next(x for x in DynamicElement.__subclasses__() if x.__name__ == symbol)
        except StopIteration:
            raise ValueError(f'DynamicElement with symbol "{symbol}" not found')
        return element

    @classmethod
    def from_atomic_number(cls, number):
        """
        get DynamicElement class by its number
        """
        try:
            element = next(x for x in DynamicElement.__subclasses__() if x.atomic_number.fget(None) == number)
        except StopIteration:
            raise ValueError(f'DynamicElement with number "{number}" not found')
        return element

    @property
    def p_charge(self):
        try:
            return self._graph()._p_charges[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def p_is_radical(self):
        try:
            return self._graph()._p_radicals[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def p_neighbors(self):
        try:
            return self._graph()._p_neighbors[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def p_hybridization(self):
        try:
            return self._graph()._p_hybridizations[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    def __eq__(self, other):
        return isinstance(other, DynamicElement) and self.atomic_number == other.atomic_number and \
               self.isotope == other.isotope

    def __int__(self):
        """
        26bit = 9bit | 7bit | 4bit | 4bit | 1bit| 1bit
        """
        return (self.isotope or 0) << 17 | self.atomic_number << 10 | self.charge + 4 << 6 | \
            self.__p_charge + 4 << 2 | self.is_radical << 1 | self.__p_is_radical


class QueryElement(Core):
    __slots__ = ()

    @classmethod
    def from_symbol(cls, symbol):
        """
        get Element class by its symbol
        """
        try:
            element = next(x for x in QueryElement.__subclasses__() if x.__name__ == symbol)
        except StopIteration:
            raise ValueError(f'Element with symbol "{symbol}" not found')
        return element

    @classmethod
    def from_atomic_number(cls, number):
        """
        get Element class by its number
        """
        try:
            element = next(x for x in QueryElement.__subclasses__() if x.atomic_number.fget(None) == number)
        except StopIteration:
            raise ValueError(f'Element with number "{number}" not found')
        return element

    def __eq__(self, other):
        return isinstance(other, (Element, QueryElement)) and self.atomic_number == other.atomic_number and \
            self.isotope == other.isotope

    def __int__(self):
        """
        21bit = 9bit | 7bit | 4bit | 1bit | 4bit | 2bit
        """
        return (self.isotope or 0) << 18 | self.atomic_number << 11 | self.charge + 4 << 7 | self.is_radical << 6 | \
            self.neighbors << 2 | self.hybridization - 1


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

    def __repr__(self):
        return repr(self.__d)
