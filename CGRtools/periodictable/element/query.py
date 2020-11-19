# -*- coding: utf-8 -*-
#
#  Copyright 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import Tuple, Dict, Type, List
from .core import Core
from .element import Element
from ...exceptions import IsNotConnectedAtom


class Query(Core):
    __slots__ = ()

    @property
    def neighbors(self) -> Tuple[int, ...]:
        try:
            return self._graph()._neighbors[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @neighbors.setter
    def neighbors(self, neighbors):
        try:
            g = self._graph()
            g._neighbors[self._map] = g._validate_neighbors(neighbors)
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    @Core.hybridization.setter
    def hybridization(self, hybridization):
        try:
            g = self._graph()
            g._hybridizations[self._map] = g._validate_hybridization(hybridization)
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    _hybridization_bitmap = {(): 0, (1,): 1, (2,): 2, (3,): 3, (4,): 4, (1, 2): 5, (1, 3): 6, (1, 4): 7, (2, 3): 8,
                             (2, 4): 9, (3, 4): 10, (1, 2, 3): 11, (1, 2, 4): 12, (1, 3, 4): 13, (2, 3, 4): 14,
                             (1, 2, 3, 4): 15}  # 4 bit
    _neighbors_bitmap = {(): 0, (1,): 1, (2,): 2, (3,): 3, (4,): 4, (5,): 5, (6,): 6,
                         (1, 2): 7, (2, 3): 8, (3, 4): 9, (4, 5): 10, (5, 6): 11,
                         (1, 2, 3): 12, (2, 3, 4): 13, (1, 2, 3, 4): 14}  # 15 is any other combination


class QueryElement(Query):
    __slots__ = ()

    @property
    def atomic_symbol(self) -> str:
        return self.__class__.__name__[5:]

    @classmethod
    def from_symbol(cls, symbol: str) -> Type['QueryElement']:
        """
        get Element class by its symbol
        """
        if symbol == 'A':
            return AnyElement
        try:
            element = next(x for x in QueryElement.__subclasses__() if x.__name__ == f'Query{symbol}')
        except StopIteration:
            raise ValueError(f'QueryElement with symbol "{symbol}" not found')
        return element

    @classmethod
    def from_atomic_number(cls, number: int) -> Type['QueryElement']:
        """
        get Element class by its number
        """
        if number == 0:
            return AnyElement
        try:
            element = next(x for x in QueryElement.__subclasses__() if x.atomic_number.fget(None) == number)
        except StopIteration:
            raise ValueError(f'QueryElement with number "{number}" not found')
        return element

    @Core.charge.setter
    def charge(self, charge):
        try:
            g = self._graph()
            # remove stereo data
            if self._map in g._atoms_stereo:
                del g._atoms_stereo[self._map]
            if g._cis_trans_stereo:
                for nm, path in g._stereo_cis_trans_paths.items():
                    if self._map in path and nm in g._cis_trans_stereo:
                        del g._cis_trans_stereo[nm]
            if g._allenes_stereo:
                for c, path in g._stereo_allenes_paths.items():
                    if self._map in path and c in g._allenes_stereo:
                        del g._allenes_stereo[c]
            g._charges[self._map] = g._validate_charge(charge)
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    @Core.is_radical.setter
    def is_radical(self, is_radical):
        try:
            g = self._graph()
            # remove stereo data
            if self._map in g._atoms_stereo:
                del g._atoms_stereo[self._map]
            if g._cis_trans_stereo:
                for nm, path in g._stereo_cis_trans_paths.items():
                    if self._map in path and nm in g._cis_trans_stereo:
                        del g._cis_trans_stereo[nm]
            if g._allenes_stereo:
                for c, path in g._stereo_allenes_paths.items():
                    if self._map in path and c in g._allenes_stereo:
                        del g._allenes_stereo[c]
            g._radicals[self._map] = g._validate_radical(is_radical)
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    def __eq__(self, other):
        """
        compare attached to molecules elements and query elements
        """
        if isinstance(other, Element):
            if self.atomic_number == other.atomic_number and self.charge == other.charge and \
                    self.is_radical == other.is_radical:
                if self.isotope and self.isotope != other.isotope:
                    return False
                if self.neighbors and other.neighbors not in self.neighbors:
                    return False
                if self.hybridization and other.hybridization not in self.hybridization:
                    return False
                return True
        elif isinstance(other, QueryElement) and self.atomic_number == other.atomic_number and \
                self.isotope == other.isotope and self.charge == other.charge and self.is_radical == other.is_radical \
                and self.neighbors == other.neighbors and self.hybridization and other.hybridization:
            # equal query element has equal query marks
            return True
        return False

    def __hash__(self):
        """
        21bit = 9bit | 7bit | 4bit | 1bit | 4bit | 4bit
        """
        return (self.isotope or 0) << 20 | self.atomic_number << 13 | self.charge + 4 << 9 | self.is_radical << 8 | \
            self._neighbors_bitmap.get(self.neighbors, 15) << 4 | self._hybridization_bitmap[self.hybridization]


class AnyElement(Query):
    __slots__ = ()

    def __init__(self, *args, **kwargs):
        super().__init__()

    @property
    def atomic_symbol(self) -> str:
        return 'A'

    @property
    def atomic_number(self) -> int:
        return 0

    @property
    def isotopes_distribution(self) -> Dict[int, float]:
        return {}

    @property
    def isotopes_masses(self) -> Dict[int, float]:
        return {}

    @property
    def atomic_radius(self):
        return 0.5

    @Core.charge.setter
    def charge(self, charge):
        try:
            g = self._graph()
            g._charges[self._map] = g._validate_charge(charge)
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    @Core.is_radical.setter
    def is_radical(self, is_radical):
        try:
            g = self._graph()
            g._radicals[self._map] = g._validate_radical(is_radical)
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    def __eq__(self, other):
        """
        Compare attached to molecules elements and query elements
        """
        if isinstance(other, Element):
            if self.charge == other.charge and self.is_radical == other.is_radical:
                if self.neighbors and other.neighbors not in self.neighbors:
                    return False
                if self.hybridization and other.hybridization not in self.hybridization:
                    return False
                return True
        elif isinstance(other, Query) and self.charge == other.charge and self.is_radical == other.is_radical \
                and self.neighbors == other.neighbors and self.hybridization and other.hybridization:
            return True
        return False

    def __hash__(self):
        """
        13bit = 4bit | 1bit | 4bit | 4bit
        """
        return self.charge + 4 << 9 | self.is_radical << 8 | \
            self._neighbors_bitmap.get(self.neighbors, 15) << 4 | self._hybridization_bitmap[self.hybridization]


class ListElement(AnyElement):
    __slots__ = ('_elements', '_numbers')

    def __init__(self, elements: List[str], *args, **kwargs):
        """
        Elements list
        """
        super().__init__()
        self._elements = tuple(elements)
        self._numbers = tuple(x.atomic_number.fget(None) for x in Element.__subclasses__() if x.__name__ in elements)

    @property
    def atomic_symbol(self) -> str:
        return ','.join(self._elements)

    def __eq__(self, other):
        """
        Compare attached to molecules elements and query elements
        """
        if isinstance(other, Element):
            if other.atomic_number in self._numbers:
                if self.charge != other.charge or self.is_radical != other.is_radical:
                    return False
                if self.neighbors and other.neighbors not in self.neighbors:
                    return False
                if self.hybridization and other.hybridization not in self.hybridization:
                    return False
                return True
        elif isinstance(other, Query) and self.charge == other.charge and self.is_radical == other.is_radical \
                and self.neighbors == other.neighbors and self.hybridization and other.hybridization:
            if isinstance(other, ListElement):
                return self._numbers == other._numbers
            if isinstance(other, AnyElement):
                return True
            return other.atomic_number in self._numbers
        return False

    def __hash__(self):
        """
        13bit = 4bit | 1bit | 4bit | 4bit
        """
        return 127 << 13 | self.charge + 4 << 9 | self.is_radical << 8 | \
            self._neighbors_bitmap.get(self.neighbors, 15) << 4 | self._hybridization_bitmap[self.hybridization]

    def __getstate__(self):
        state = super().__getstate__()
        state['elements'] = self._elements
        return state

    def __setstate__(self, state):
        self._elements = state['elements']
        self._numbers = tuple(x.atomic_number.fget(None) for x in Element.__subclasses__()
                              if x.__name__ in state['elements'])
        super().__setstate__(state)

    def __repr__(self):
        return f'{self.__class__.__name__}([{",".join(self._elements)}])'


__all__ = ['QueryElement', 'AnyElement', 'ListElement']
