# -*- coding: utf-8 -*-
#
#  Copyright 2020, 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import Tuple, Dict, Type, Union
from .core import Core
from .dynamic import Dynamic, DynamicElement
from .element import Element
from .query import QueryElement, AnyElement
from ..._functions import tuple_hash
from ...exceptions import IsNotConnectedAtom


class DynamicQuery(Dynamic):
    __slots__ = ()

    @property
    def neighbors(self) -> Tuple[int, ...]:
        """
        Number of neighbors of atom in reactant state.
        """
        try:
            return self._graph()._neighbors[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @property
    def p_neighbors(self) -> Tuple[int, ...]:
        """
        Number of neighbors of atom in product state.
        """
        try:
            return self._graph()._p_neighbors[self._map]
        except AttributeError:
            raise IsNotConnectedAtom

    @neighbors.setter
    def neighbors(self, neighbors):
        try:
            g = self._graph()
            neighbors = g._validate_neighbors(neighbors)
            neighbors, p_neighbors = g._validate_neighbors_pairing(neighbors, g._p_neighbors[self._map])
            g._neighbors[self._map] = neighbors
            g._p_neighbors[self._map] = p_neighbors
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    @p_neighbors.setter
    def p_neighbors(self, p_neighbors):
        try:
            g = self._graph()
            p_neighbors = g._validate_neighbors(p_neighbors)
            neighbors, p_neighbors = g._validate_neighbors_pairing(g._neighbors[self._map], p_neighbors)
            g._neighbors[self._map] = neighbors
            g._p_neighbors[self._map] = p_neighbors
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    @Core.hybridization.setter
    def hybridization(self, hybridization):
        try:
            g = self._graph()
            hybridization = g._validate_hybridization(hybridization)
            hybridization, p_hybridization = g._validate_hybridization_pairing(hybridization,
                                                                               g._p_hybridizations[self._map])
            g._hybridizations[self._map] = hybridization
            g._p_hybridizations[self._map] = p_hybridization
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom

    @Dynamic.p_hybridization.setter
    def p_hybridization(self, p_hybridization):
        try:
            g = self._graph()
            p_hybridization = g._validate_hybridization(p_hybridization)
            hybridization, p_hybridization = g._validate_hybridization_pairing(g._hybridizations[self._map],
                                                                               p_hybridization)
            g._hybridizations[self._map] = hybridization
            g._p_hybridizations[self._map] = p_hybridization
            g.flush_cache()
        except AttributeError:
            raise IsNotConnectedAtom


class DynamicQueryElement(DynamicQuery):
    __slots__ = ()

    @property
    def atomic_symbol(self) -> str:
        return self.__class__.__name__[12:]

    @classmethod
    def from_symbol(cls, symbol: str) -> Type[Union['DynamicQueryElement', 'DynamicAnyElement']]:
        """
        get Element class by its symbol
        """
        if symbol == 'A':
            return DynamicAnyElement
        try:
            element = next(x for x in DynamicQueryElement.__subclasses__() if x.__name__ == f'DynamicQuery{symbol}')
        except StopIteration:
            raise ValueError(f'DynamicQueryElement with symbol "{symbol}" not found')
        return element

    @classmethod
    def from_atomic_number(cls, number: int) -> Type[Union['DynamicQueryElement', 'DynamicAnyElement']]:
        """
        get Element class by its number
        """
        if number == 0:
            return DynamicAnyElement
        try:
            element = next(x for x in DynamicQueryElement.__subclasses__() if x.atomic_number.fget(None) == number)
        except StopIteration:
            raise ValueError(f'DynamicQueryElement with number "{number}" not found')
        return element

    @classmethod
    def from_atom(cls, atom: Union['Element', 'DynamicElement', 'DynamicQueryElement', 'DynamicAnyElement',
                                   'QueryElement', 'AnyElement']) -> Union['DynamicQueryElement', 'DynamicAnyElement']:
        """
        get DynamicQueryElement or DynamicAnyElement object from Element or DynamicElement or QueryElement object or
        copy of DynamicQueryElement or DynamicAnyElement
        """
        if isinstance(atom, (Element, DynamicElement, QueryElement, AnyElement)):
            return cls.from_atomic_number(atom.atomic_number)(atom.isotope)
        elif not isinstance(atom, (DynamicQueryElement, DynamicAnyElement)):
            raise TypeError('Element, DynamicElement, DynamicQueryElement, DynamicAnyElement,'
                            ' QueryElement or AnyElement expected')
        return atom.copy()

    def __eq__(self, other):
        if isinstance(other, DynamicElement):
            if self.atomic_number == other.atomic_number and \
                    self.charge == other.charge and self.p_charge == other.p_charge and \
                    self.is_radical == other.is_radical and self.p_is_radical == other.p_is_radical:
                if self.isotope and self.isotope != other.isotope:
                    return False
                if self.neighbors and (other.neighbors, other.p_neighbors) not in zip(self.neighbors, self.p_neighbors):
                    return False
                if self.hybridization and (other.hybridization, other.p_hybridization) not in zip(self.hybridization,
                                                                                                  self.p_hybridization):
                    return False
                return True
        elif isinstance(other, DynamicQueryElement) and self.atomic_number == other.atomic_number and \
                self.isotope == other.isotope and self.charge == other.charge and self.p_charge == other.p_charge and \
                self.is_radical == other.is_radical and self.p_is_radical == other.p_is_radical and \
                self.neighbors == other.neighbors and self.hybridization and other.hybridization and \
                self.p_neighbors == other.p_neighbors and self.p_hybridization and other.p_hybridization:
            # equal query element has equal query marks
            return True
        return False

    def __hash__(self):
        return tuple_hash((self.isotope or 0, self.atomic_number, self.charge, self.p_charge,
                           self.is_radical, self.p_is_radical, self.hybridization, self.p_hybridization,
                           self.neighbors, self.p_neighbors))


class DynamicAnyElement(DynamicQuery):  # except Hydrogen!
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

    def __eq__(self, other):
        if isinstance(other, DynamicElement):
            if self.charge == other.charge and self.p_charge == other.p_charge and \
                    self.is_radical == other.is_radical and self.p_is_radical == other.p_is_radical:
                if self.neighbors and (other.neighbors, other.p_neighbors) not in zip(self.neighbors, self.p_neighbors):
                    return False
                if self.hybridization and (other.hybridization, other.p_hybridization) not in zip(self.hybridization,
                                                                                                  self.p_hybridization):
                    return False
                return True
        elif isinstance(other, DynamicQuery) and self.charge == other.charge and self.p_charge == other.p_charge and \
                self.is_radical == other.is_radical and self.p_is_radical == other.p_is_radical and \
                self.neighbors == other.neighbors and self.hybridization and other.hybridization and \
                self.p_neighbors == other.p_neighbors and self.p_hybridization and other.p_hybridization:
            return True
            # equal query element has equal query marks
        return False

    def __hash__(self):
        return tuple_hash((self.charge, self.p_charge, self.is_radical, self.p_is_radical,
                           self.hybridization, self.p_hybridization, self.neighbors, self.p_neighbors))


__all__ = ['DynamicQueryElement', 'DynamicAnyElement']
