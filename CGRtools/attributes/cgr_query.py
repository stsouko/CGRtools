# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
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
from .cgr import DynAtom, DynAtomAttribute
from .molecule import Atom
from .query import QueryAtom
from ..cache import cached_property
from ..periodictable import Element


class DynQueryAtom(DynAtomAttribute):
    def __setattr__(self, key, value):
        if key in self._p_static:
            raise AttributeError(f'{key} is invalid')
        elif key.startswith('p_'):
            key = key[2:]
            value = getattr(self._reactant, f'_{key}_check')(value)
            if key in ('neighbors', 'hybridization'):
                if len(value) != len(self[key]):
                    raise ValueError(f'{key} lists in reactant and product should be equal size')
                value = sorted(zip(self[key], value))
                if len(value) != len(set(value)):
                    raise ValueError('duplicates found')
                self._reactant[key], value = zip(*value)
            self._product[key] = value
        else:
            value = getattr(self._reactant, f'_{key}_check')(value)
            if key in self._static:
                self._product[key] = value
            elif key in ('neighbors', 'hybridization'):
                if len(value) != len(self._product[key]):
                    raise ValueError(f'{key} lists in reactant and product should be equal size')
                value = sorted(zip(value, self._product[key]))
                if len(value) != len(set(value)):
                    raise ValueError('duplicates found')
                value, self._product[key] = zip(*value)
            self._reactant[key] = value
        self.__dict__.clear()

    def _update(self, value, kwargs):
        if isinstance(value, (DynQueryAtom, DynAtom)):
            r, p = self._split_check_kwargs(kwargs)
            p_value = value._product
            value = value._reactant
        elif isinstance(value, (QueryAtom, Atom, Element)):
            r, p = self._split_check_kwargs(kwargs)
            p_value = value
        elif isinstance(value, type):
            if not issubclass(value, Element):
                raise ValueError('only CGRtools.periodictable.Element subclasses allowed')
            r, p = self._split_check_kwargs(kwargs)
            p_value = value
        else:
            if not isinstance(value, dict):
                try:
                    value = dict(value)
                except (TypeError, ValueError):
                    raise TypeError('invalid attrs sequence')

            value.update(kwargs)
            if not value:  # ad-hoc for add_nodes_from method
                return
            value, p_value = self._split_check_kwargs(value)
            r = p = {}

        self._reactant._update(value, r)
        self._product._update(p_value, p)
        self.__dict__.clear()

    def __eq__(self, other):
        if isinstance(other, DynAtom):
            if self.neighbors:
                if (other.neighbors, other.p_neighbors) not in self.zip_neighbors:
                    return False
            if self.hybridization:
                if (other.hybridization, other.p_hybridization) not in self.zip_hybridization:
                    return False
            return ((other.element in self.element_set if self.element else True) and
                    (self.isotope == other.isotope if self.isotope else True) and
                    self.charge == other.charge and self.p_charge == other.p_charge and
                    self.multiplicity == other.multiplicity and self.p_multiplicity == other.p_multiplicity)
        elif isinstance(other, DynQueryAtom):
            if self.element:
                if not other.element:
                    return False
                elif not self.element_set.issuperset(other.element_set):
                    return False
            if self.neighbors:
                if not other.neighbors:
                    return False
                elif not self.zip_neighbors.issuperset(other.zip_neighbors):
                    return False
            if self.hybridization:
                if not other.hybridization:
                    return False
                elif not self.zip_hybridization.issuperset(other.zip_hybridization):
                    return False

            return ((self.isotope == other.isotope if self.isotope else True) and
                    self.charge == other.charge and self.p_charge == other.p_charge and
                    self.multiplicity == other.multiplicity and self.p_multiplicity == other.p_multiplicity)
        elif isinstance(other, (QueryAtom, Atom)):
            return self._reactant == other and self._product == other
        return False

    @cached_property
    def zip_neighbors(self):
        return set(zip(self.neighbors, self.p_neighbors))

    @cached_property
    def zip_hybridization(self):
        return set(zip(self.hybridization, self.p_hybridization))

    def _split_check_kwargs(self, kwargs):
        r, p = super()._split_check_kwargs(kwargs)
        for k in ('neighbors', 'hybridization'):
            if k in r:
                if k not in p:
                    raise ValueError(f'{k} attribute should be presented with p_{k}')
                elif len(r[k]) != len(p[k]):
                    raise ValueError(f'{k} lists in reactant and product should be equal size')
                value = sorted(zip(r[k], p[k]))
                if value:
                    if len(value) != len(set(value)):
                        raise ValueError('duplicates found')
                    r[k], p[k] = zip(*value)
        return r, p

    _factory = QueryAtom
    _static = {'element', 'isotope'}
    _p_static = {'p_element', 'p_isotope', 'p_element_set'}


__all__ = ['DynQueryAtom']
