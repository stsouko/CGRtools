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
from abc import ABC, abstractmethod
from CachedMethods import cached_property, cached_args_method
from typing import Dict, Optional, Tuple, Iterable, Iterator, Union, List, Type
from .bonds import Bond, DynamicBond
from ..algorithms.components import GraphComponents
from ..algorithms.isomorphism import Isomorphism
from ..algorithms.mcs import MCS
from ..algorithms.morgan import Morgan
from ..algorithms.sssr import SSSR
from ..exceptions import AtomNotFound
from ..periodictable.element import Core


class Graph(Isomorphism, MCS, SSSR, Morgan, GraphComponents, ABC):
    __slots__ = ('_atoms', '_bonds', '_meta', '_plane', '__dict__', '__weakref__', '_parsed_mapping', '_charges',
                 '_radicals')

    def __init__(self):
        """
        Empty data object initialization or conversion from another object type
        """
        self._atoms: Dict[int, Core] = {}
        self._charges: Dict[int, int] = {}
        self._radicals: Dict[int, bool] = {}
        self._plane: Dict[int, Tuple[float, float]] = {}
        self._bonds: Dict[int, Dict[int, Union[Bond, DynamicBond]]] = {}
        self._meta = {}
        self._parsed_mapping: Dict[int, int] = {}

    def __getstate__(self):
        return {'atoms': self._atoms, 'bonds': self._bonds, 'meta': self._meta, 'plane': self._plane,
                'parsed_mapping': self._parsed_mapping, 'charges': self._charges, 'radicals': self._radicals}

    def __setstate__(self, state):
        self._atoms = state['atoms']
        for n, a in state['atoms'].items():
            a._attach_to_graph(self, n)
        self._charges = state['charges']
        self._radicals = state['radicals']
        self._plane = state['plane']
        self._bonds = state['bonds']
        self._meta = state['meta']
        self._parsed_mapping = state['parsed_mapping']

    def __len__(self):
        return len(self._atoms)

    def __iter__(self):
        return iter(self._atoms)

    def __contains__(self, n: int):
        return n in self._atoms

    def __bool__(self):
        return bool(self._atoms)

    def atom(self, n: int) -> Core:
        return self._atoms[n]

    def has_atom(self, n: int) -> bool:
        return n in self._atoms

    def atoms(self) -> Iterator[Tuple[int, Core]]:
        """
        iterate over all atoms
        """
        return iter(self._atoms.items())

    @cached_property
    def atoms_count(self) -> int:
        return len(self._atoms)

    @cached_property
    def atoms_numbers(self) -> Tuple[int, ...]:
        return tuple(self._atoms)

    @cached_args_method
    def environment(self, atom: int) -> Tuple[Tuple[Union[Bond, DynamicBond], Core], ...]:
        """
        pairs of (bond, atom) connected to atom

        :param atom: number
        """
        atoms = self._atoms
        return tuple((bond, atoms[n]) for n, bond in self._bonds[atom].items())

    def bond(self, n: int, m: int) -> Union[Bond, DynamicBond]:
        return self._bonds[n][m]

    def has_bond(self, n: int, m: int) -> bool:
        try:
            self._bonds[n]  # check if atom exists
            return n in self._bonds[m]
        except KeyError:
            raise AtomNotFound

    def bonds(self) -> Iterator[Tuple[int, int, Union[Bond, DynamicBond]]]:
        """
        iterate other all bonds
        """
        seen = set()
        for n, m_bond in self._bonds.items():
            seen.add(n)
            for m, bond in m_bond.items():
                if m not in seen:
                    yield n, m, bond

    @cached_property
    def bonds_count(self) -> int:
        return sum(len(x) for x in self._bonds.values()) // 2

    @property
    def meta(self) -> Dict:
        return self._meta

    @abstractmethod
    def add_atom(self, atom, _map: Optional[int] = None, *, charge: int = 0,
                 is_radical: bool = False, xy: Tuple[float, float] = (0., 0.)) -> int:
        """
        new atom addition
        """
        if _map is None:
            _map = max(self._atoms, default=0) + 1
        elif not isinstance(_map, int):
            raise TypeError('mapping should be integer')
        elif _map in self._atoms:
            raise ValueError('atom with same number exists')

        if not isinstance(xy, tuple) or len(xy) != 2 or not isinstance(xy[0], float) or not isinstance(xy[1], float):
            raise TypeError('XY should be tuple with 2 float')

        charge = self._validate_charge(charge)
        is_radical = self._validate_radical(is_radical)

        self._atoms[_map] = atom
        self._charges[_map] = charge
        self._radicals[_map] = is_radical
        self._plane[_map] = xy
        self._bonds[_map] = {}
        atom._attach_to_graph(self, _map)
        self.__dict__.clear()
        return _map

    @abstractmethod
    def add_bond(self, n: int, m: int, bond):
        """
        new bond addition
        """
        if n == m:
            raise ValueError('atom loops impossible')
        if n not in self._bonds or m not in self._bonds:
            raise AtomNotFound('atoms not found')
        if n in self._bonds[m]:
            raise ValueError('atoms already bonded')

        self._bonds[n][m] = self._bonds[m][n] = bond
        self.__dict__.clear()

    @abstractmethod
    def delete_atom(self, n: int):
        """
        implementation of atom removing
        """
        del self._atoms[n]
        del self._charges[n]
        del self._radicals[n]
        del self._plane[n]
        sb = self._bonds
        for m in sb.pop(n):
            del sb[m][n]
        try:
            del self._parsed_mapping[n]
        except KeyError:
            pass
        self.__dict__.clear()

    def delete_bond(self, n: int, m: int):
        """
        implementation of bond removing
        """
        del self._bonds[n][m]
        del self._bonds[m][n]
        self.__dict__.clear()

    @abstractmethod
    def remap(self, mapping: Dict[int, int], *, copy: bool = False):
        if len(mapping) != len(set(mapping.values())):
            raise ValueError('mapping overlap')

        mg = mapping.get
        sp = self._plane
        sc = self._charges
        sr = self._radicals

        if copy:
            h = self.__class__()
            h._meta.update(self._meta)
            hb = h._bonds
            ha = h._atoms
            hc = h._charges
            hr = h._radicals
            hp = h._plane
            hm = h._parsed_mapping
            for n, atom in self._atoms.items():
                m = mg(n, n)
                hc[m] = sc[n]
                hr[m] = sr[n]
                hp[m] = sp[n]
                ha[m] = atom = atom.copy()
                atom._attach_to_graph(h, m)
        else:
            hb = {}
            ha = {}
            hc = {}
            hr = {}
            hp = {}
            hm = {}
            for n, atom in self._atoms.items():
                m = mg(n, n)
                hc[m] = sc[n]
                hr[m] = sr[n]
                hp[m] = sp[n]
                ha[m] = atom
                atom._change_map(m)  # change mapping number
            self._atoms = ha
            self._charges = hc
            self._radicals = hr
            self._plane = hp

        for n, m_bond in self._bonds.items():
            hb[mg(n, n)] = {mg(m, m): b for m, b in m_bond.items()}

        for n, m in self._parsed_mapping.items():
            hm[mg(n, n)] = m

        if copy:
            return h

        self._bonds = hb
        self._parsed_mapping = hm
        return self

    @abstractmethod
    def copy(self, *, meta: bool = True):
        """
        copy of graph

        :param meta: include metadata
        """
        copy = object.__new__(self.__class__)
        copy._meta = self._meta.copy() if meta else {}

        copy._charges = self._charges.copy()
        copy._radicals = self._radicals.copy()
        copy._plane = self._plane.copy()
        copy._parsed_mapping = self._parsed_mapping.copy()

        copy._bonds = cb = {n: {} for n in self._bonds}
        seen = set()
        for n, m_bond in self._bonds.items():
            seen.add(n)
            for m, bond in m_bond.items():
                if m not in seen:
                    cb[n][m] = cb[m][n] = bond.copy()

        copy._atoms = ca = {}
        for n, atom in self._atoms.items():
            ca[n] = atom = atom.copy()
            atom._attach_to_graph(copy, n)
        return copy

    @abstractmethod
    def substructure(self, atoms: Iterable[int], sub: Type['Graph'], *, meta: bool = False):
        if not atoms:
            raise ValueError('empty atoms list not allowed')
        if set(atoms) - self._atoms.keys():
            raise ValueError('invalid atom numbers')
        atoms = tuple(n for n in self._atoms if n in atoms)  # save original order
        sub = object.__new__(sub)

        sc = self._charges
        sr = self._radicals
        sp = self._plane
        sb = self._bonds

        sub._meta = self._meta.copy() if meta else {}
        sub._charges = {n: sc[n] for n in atoms}
        sub._radicals = {n: sr[n] for n in atoms}
        sub._plane = {n: sp[n] for n in atoms}
        sub._parsed_mapping = {n: m for n, m in self._parsed_mapping.items() if n in atoms}

        sub._bonds = cb = {n: {} for n in atoms}
        seen = set()
        for n in atoms:
            seen.add(n)
            for m, bond in sb[n].items():
                if m not in seen and m in atoms:
                    cb[n][m] = cb[m][n] = bond.copy()

        return sub, atoms

    def __and__(self, other):
        """
        substructure of graph
        """
        return self.substructure(other)

    def __sub__(self, other):
        """
        other nodes excluded substructure of graph
        :return graph or None
        """
        atoms = set(other)
        if atoms - self._atoms.keys():
            raise ValueError('invalid atom numbers')
        atoms = self._atoms.keys() - atoms
        if atoms:
            return self.substructure(atoms)
        raise ValueError('full substitution not allowed')

    def augmented_substructure(self, atoms: Iterable[int], deep: int = 1, **kwargs) -> 'Graph':
        """
        create substructure containing atoms and their neighbors

        :param atoms: list of core atoms in graph
        :param deep: number of bonds between atoms and neighbors
        :param meta: copy metadata to each substructure
        :param as_query: return Query object based on graph substructure. for Molecule and CGR only
        """
        return self.substructure(self.__augmented_substructure(atoms, deep)[-1], **kwargs)

    def augmented_substructures(self, atoms: Iterable[int], deep: int = 1, **kwargs) -> List['Graph']:
        """
        create list of substructures containing atoms and their neighbors

        :param atoms: list of core atoms in graph
        :param deep: number of bonds between atoms and neighbors
        :param meta: copy metadata to each substructure
        :param as_query: return Query object based on graph substructure. for Molecule and CGR only
        :return: list of graphs containing atoms, atoms + first circle, atoms + 1st + 2nd,
            etc up to deep or while new nodes available
        """
        return [self.substructure(a, **kwargs) for a in self.__augmented_substructure(atoms, deep)]

    def union(self, other: 'Graph') -> 'Graph':
        if self._atoms.keys() & other._atoms.keys():
            raise ValueError('mapping of graphs is not disjoint')

        u = self.copy(meta=False)
        u._charges.update(other._charges)
        u._radicals.update(other._radicals)
        u._plane.update(other._plane)
        u._parsed_mapping.update(other._parsed_mapping)
        return u

    def __or__(self, other):
        """
        G | H is union of graphs
        """
        return self.union(other)

    def split(self, meta: bool = False) -> List['Graph']:
        """
        split disconnected structure to connected substructures

        :param meta: copy metadata to each substructure
        :return: list of substructures
        """
        return [self.substructure(c, meta=meta) for c in self.connected_components]

    def flush_cache(self):
        self.__dict__.clear()

    @staticmethod
    def _validate_charge(charge):
        if not isinstance(charge, int):
            raise TypeError('formal charge should be int in range [-4, 4]')
        if charge > 4 or charge < -4:
            raise ValueError('formal charge should be in range [-4, 4]')
        return charge

    @staticmethod
    def _validate_radical(is_radical):
        if not isinstance(is_radical, bool):
            raise TypeError('radical state should be bool')
        return is_radical

    def __augmented_substructure(self, atoms, deep):
        atoms = set(atoms)
        if atoms - self._atoms.keys():
            raise ValueError('invalid atom numbers')
        nodes = [atoms]
        for i in range(deep):
            n = {y for x in nodes[-1] for y in self._bonds[x]} | nodes[-1]
            if n in nodes:
                break
            nodes.append(n)
        return nodes


__all__ = ['Graph']
