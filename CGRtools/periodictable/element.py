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
from CachedMethods import class_cached_property
from collections import defaultdict
from collections.abc import Mapping
from typing import Optional, Tuple, Dict, Set, List
from weakref import ref
from ..exceptions import IsConnectedAtom, IsNotConnectedAtom, ValenceError


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


class Dynamic:
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


class Element(Core):
    __slots__ = ()
    __class_cache__ = {}

    @property
    def atomic_symbol(self):
        return self.__class__.__name__

    @property
    def implicit_hydrogens(self):
        try:
            return self._graph()._hydrogens[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def explicit_hydrogens(self):
        try:
            return self._graph()._explicit_hydrogens(self._map)
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def total_hydrogens(self):
        try:
            return self._graph()._total_hydrogens(self._map)
        except AttributeError:
            raise IsNotConnectedAtom

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
        """
        compare attached to molecules elements
        """
        return isinstance(other, Element) and self.atomic_number == other.atomic_number and \
            self.isotope == other.isotope and self.charge == other.charge and self.is_radical == other.is_radical

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


class DynamicElement(Core, Dynamic):
    __slots__ = ('__p_charge', '__p_is_radical')

    @property
    def atomic_symbol(self):
        return self.__class__.__name__[7:]

    @classmethod
    def from_symbol(cls, symbol):
        """
        get DynamicElement class by its symbol
        """
        try:
            element = next(x for x in DynamicElement.__subclasses__() if x.__name__[7:] == symbol)
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

    def __eq__(self, other):
        """
        compare attached to molecules dynamic elements
        """
        return isinstance(other, DynamicElement) and self.atomic_number == other.atomic_number and \
            self.isotope == other.isotope and self.charge == other.charge and self.is_radical == other.is_radical and \
            self.p_charge == other.p_charge and self.p_is_radical == other.p_is_radical

    def __int__(self):
        """
        26bit = 9bit | 7bit | 4bit | 4bit | 1bit| 1bit
        """
        return (self.isotope or 0) << 17 | self.atomic_number << 10 | self.charge + 4 << 6 | \
            self.p_charge + 4 << 2 | self.is_radical << 1 | self.p_is_radical


class QueryElement(Core):
    __slots__ = ()

    @property
    def atomic_symbol(self):
        return self.__class__.__name__[5:]

    @classmethod
    def from_symbol(cls, symbol):
        """
        get Element class by its symbol
        """
        try:
            element = next(x for x in QueryElement.__subclasses__() if x.__name__[5:] == symbol)
        except StopIteration:
            raise ValueError(f'QueryElement with symbol "{symbol}" not found')
        return element

    @classmethod
    def from_atomic_number(cls, number):
        """
        get Element class by its number
        """
        try:
            element = next(x for x in QueryElement.__subclasses__() if x.atomic_number.fget(None) == number)
        except StopIteration:
            raise ValueError(f'QueryElement with number "{number}" not found')
        return element

    def __eq__(self, other):
        """
        compare attached to molecules elements and query elements
        """
        if isinstance(other, Element):
            if self.atomic_number == other.atomic_number and self.isotope == other.isotope and \
                    self.charge == other.charge and self.is_radical == other.is_radical:
                if self.neighbors:
                    if other.neighbors in self.neighbors:
                        if self.hybridization:
                            if other.hybridization in self.hybridization:
                                return True
                        else:
                            return True
                elif self.hybridization:
                    if other.hybridization in self.hybridization:
                        return True
                else:
                    return True
        elif isinstance(other, QueryElement) and self.atomic_number == other.atomic_number and \
                self.isotope == other.isotope and self.charge == other.charge and self.is_radical == other.is_radical \
                and self.neighbors == other.neighbors and self.hybridization and other.hybridization:
            # equal query element has equal query marks
            return True
        return False

    def __int__(self):
        """
        21bit = 9bit | 7bit | 4bit | 1bit | 4bit | 4bit
        """
        return (self.isotope or 0) << 20 | self.atomic_number << 13 | self.charge + 4 << 9 | self.is_radical << 8 | \
            self._neighbors_bitmap.get(self.neighbors, 15) << 4 | self._hybridization_bitmap[self.hybridization]

    _hybridization_bitmap = {(): 0, (1,): 1, (2,): 2, (3,): 3, (4,): 4, (1, 2): 5, (1, 3): 6, (1, 4): 7, (2, 3): 8,
                             (2, 4): 9, (3, 4): 10, (1, 2, 3): 11, (1, 2, 4): 12, (1, 3, 4): 13, (2, 3, 4): 14,
                             (1, 2, 3, 4): 15}  # 4 bit
    _neighbors_bitmap = {(): 0, (1,): 1, (2,): 2, (3,): 3, (4,): 4, (5,): 5, (6,): 6,
                         (1, 2): 7, (2, 3): 8, (3, 4): 9, (4, 5): 10, (5, 6): 11,
                         (1, 2, 3): 12, (2, 3, 4): 13, (1, 2, 3, 4): 14}  # 15 is any other combination


class DynamicQueryElement(Core, Dynamic):
    __slots__ = ()

    @property
    def atomic_symbol(self):
        return self.__class__.__name__[12:]

    @classmethod
    def from_symbol(cls, symbol):
        """
        get Element class by its symbol
        """
        try:
            element = next(x for x in DynamicQueryElement.__subclasses__() if x.__name__[12:] == symbol)
        except StopIteration:
            raise ValueError(f'DynamicQueryElement with symbol "{symbol}" not found')
        return element

    @classmethod
    def from_atomic_number(cls, number):
        """
        get Element class by its number
        """
        try:
            element = next(x for x in DynamicQueryElement.__subclasses__() if x.atomic_number.fget(None) == number)
        except StopIteration:
            raise ValueError(f'DynamicQueryElement with number "{number}" not found')
        return element

    def __eq__(self, other):
        if isinstance(other, DynamicElement):
            if self.atomic_number == other.atomic_number and self.isotope == other.isotope and \
                    self.charge == other.charge and self.p_charge == other.p_charge and \
                    self.is_radical == other.is_radical and self.p_is_radical == other.p_is_radical:
                if self.neighbors:  # neighbors and p_neighbors all times paired
                    if (other.neighbors, other.p_neighbors) in zip(self.neighbors, self.p_neighbors):
                        if self.hybridization:
                            if (other.hybridization, other.p_hybridization) in zip(self.hybridization,
                                                                                   self.p_hybridization):
                                return True
                        else:
                            return True
                elif self.hybridization:
                    if (other.hybridization, other.p_hybridization) in zip(self.hybridization, self.p_hybridization):
                        return True
                else:
                    return True
        elif isinstance(other, DynamicQueryElement) and self.atomic_number == other.atomic_number and \
                self.isotope == other.isotope and self.charge == other.charge and self.p_charge == other.p_charge and \
                self.is_radical == other.is_radical and self.p_is_radical == other.p_is_radical and \
                self.neighbors == other.neighbors and self.hybridization and other.hybridization and \
                self.p_neighbors == other.p_neighbors and self.p_hybridization and other.p_hybridization:
            # equal query element has equal query marks
            return True
        return False

    def __int__(self):
        """
        36bit = 9bit | 7bit | 4bit | 4bit | 1bit | 1bit | 5bit | 5bit
        """
        return (self.isotope or 0) << 27 | self.atomic_number << 20 | self.charge + 4 << 16 | self.p_charge << 12 | \
            self.is_radical << 11 | self.p_is_radical << 10 | \
            self._hybridization_bitmap.get(tuple(zip(self.hybridization, self.p_hybridization)), 31) << 5 | \
            self._neighbors_bitmap.get(tuple(zip(self.neighbors, self.p_neighbors)), 31)

    _hybridization_bitmap = {(): 0, ((1, 1),): 1, ((2, 2),): 2, ((3, 3),): 3, ((4, 4),): 4,
                             ((1, 2),): 5, ((1, 3),): 6, ((1, 4),): 7, ((2, 3),): 8, ((2, 4),): 9, ((3, 4),): 10,
                             ((2, 1),): 11, ((3, 1),): 12, ((4, 1),): 13, ((3, 2),): 14,
                             ((4, 2),): 15, ((4, 3),): 16, }  # 31 is any other combination
    _neighbors_bitmap = {(): 0, ((1, 1),): 1, ((2, 2),): 2, ((3, 3),): 3, ((4, 4),): 4, ((5, 5),): 5, ((6, 6),): 6,
                         ((1, 2),): 7, ((2, 3),): 8, ((3, 4),): 9, ((4, 5),): 10, ((5, 6),): 11, ((2, 1),): 12,
                         ((3, 2),): 13, ((4, 3),): 14, ((5, 4),): 15, ((6, 5),): 16}  # 31 is any other combination


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
