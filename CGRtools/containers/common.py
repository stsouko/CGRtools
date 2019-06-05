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
from typing import Dict, Optional, Tuple
from ..algorithms.isomorphism import Isomorphism


class Graph(Isomorphism, ABC):
    __slots__ = ('_atoms', '_bonds', '_meta', '_plane', '__dict__', '__weakref__', '_parsed_mapping', '_charges',
                 '_radicals')

    def __init__(self):
        """
        Empty data object initialization or conversion from another object type
        """
        self._atoms = {}
        self._charges = {}
        self._radicals = {}
        self._plane = {}
        self._bonds = {}
        self._meta = {}
        self._parsed_mapping = {}
        self.__dict__ = {}

    def __getstate__(self):
        return {'atoms': self._atoms, 'bonds': self._bonds, 'meta': self._meta, 'plane': self._plane,
                'parsed_mapping': self._parsed_mapping, 'charges': self._charges, 'radicals': self._radicals}

    def __setstate__(self, state):
        self._atoms = state['atoms']
        self._charges = state['charges']
        self._radicals = state['radicals']
        self._plane = state['plane']
        self._bonds = state['bonds']
        self._meta = state['meta']
        self._parsed_mapping = state['parsed_mapping']
        self.__dict__ = {}

    def __len__(self):
        return len(self._atoms)

    def __iter__(self):
        return iter(self._atoms)

    def atom(self, n):
        return self._atoms[n]

    def atoms(self):
        """
        iterate over all atoms
        """
        return iter(self._atoms.items())

    @cached_property
    def atoms_count(self):
        return len(self._atoms)

    @cached_property
    def atoms_numbers(self):
        return list(self._atoms)

    @cached_args_method
    def environment(self, atom):
        """
        pairs of (bond, atom) connected to atom

        :param atom: number
        :return: list
        """
        atoms = self._atoms
        return tuple((bond, atoms[n]) for n, bond in self._bonds[atom].items())

    def bond(self, n, m):
        return self._bonds[n][m]

    def bonds(self):
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
    def bonds_count(self):
        return sum(len(x) for x in self._bonds.values()) // 2

    @property
    def meta(self):
        return self._meta

    @abstractmethod
    def add_atom(self, atom, _map: Optional[int] = None, *, charge: int = 0,
                 is_radical: bool = False, xy: Tuple[float, float] = None) -> int:
        """
        new atom addition
        """
        if _map is None:
            _map = max(self._atoms, default=0) + 1
        elif not isinstance(_map, int):
            raise TypeError('mapping should be integer')
        elif _map in self._atoms:
            raise ValueError('atom with same number exists')

        if xy is None:
            xy = (0., 0.)
        elif not isinstance(xy, tuple) or len(xy) != 2 or not isinstance(xy[0], float) or not isinstance(xy[1], float):
            raise TypeError('XY should be tuple with 2 float')
        if not isinstance(charge, int) or charge > 4 or charge < -4:
            raise TypeError('formal charge should be int in range [-4, 4]')
        if not isinstance(is_radical, bool):
            raise TypeError('radical state should be bool')

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
            raise KeyError('atom loops impossible')
        if n not in self._bonds or m not in self._bonds:
            raise KeyError('atoms not found')
        if n in self._bonds[m]:
            raise KeyError('atoms already bonded')

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

    @abstractmethod
    def delete_bond(self, n: int, m: int):
        """
        implementation of bond removing
        """
        del self._bonds[n][m]
        del self._bonds[m][n]
        self.__dict__.clear()

    def substructure(self, atoms, *, meta=False, as_view=True):
        """
        create substructure containing atoms from atoms list

        :param atoms: list of atoms numbers of substructure
        :param meta: if True metadata will be copied to substructure
        :param as_view: If True, the returned graph-view provides a read-only view
            of the original structure scaffold without actually copying any data.
        """

    @abstractmethod
    def remap(self, mapping: Dict[int, int], *, copy=False):
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

        for n1, n2_bond in self._bonds.items():
            hb[mg(n1, n1)] = {mg(n2, n2): b for n2, b in n2_bond.items()}

        for n, m in self._parsed_mapping.items():
            hm[mg(n, n)] = m

        if copy:
            return h

        self._bonds = hb
        self._parsed_mapping = hm
        return self

    def flush_cache(self):
        self.__dict__.clear()

    @staticmethod
    def _get_subclass(name):
        """
        need for cyclic import solving
        """
        return next(x for x in Graph.__subclasses__() if x.__name__ == name)


__all__ = ['Graph']
