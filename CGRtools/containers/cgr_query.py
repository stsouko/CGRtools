# -*- coding: utf-8 -*-
#
#  Copyright 2019 Ramil Nugmanov <stsouko@live.ru>
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
from typing import List, Union, Tuple
from . import cgr, molecule, query  # cyclic imports resolve
from .common import Graph
from ..algorithms.smiles import QueryCGRSmiles
from ..periodictable import Element, DynamicElement, QueryElement, DynamicQueryElement


class QueryCGRContainer(Graph, QueryCGRSmiles):
    __slots__ = ('_neighbors', '_hybridization', '_p_neighbors', '_p_hybridization')

    def __init__(self):
        self._p_charges = {}
        self._p_radicals = {}
        self._neighbors = {}
        self._hybridization = {}
        self._p_neighbors = {}
        self._p_hybridization = {}
        super().__init__()

    def add_atom(self, atom: Union[DynamicQueryElement, DynamicElement, QueryElement, Element, int, str], *args,
                 p_charge: int = 0, p_is_radical: bool = False,
                 neighbors: Union[int, List[int], Tuple[int], None] = None,
                 hybridization: Union[int, List[int], Tuple[int], None] = None,
                 p_neighbors: Union[int, List[int], Tuple[int], None] = None,
                 p_hybridization: Union[int, List[int], Tuple[int], None] = None, **kwargs):
        if neighbors is None:
            neighbors = ()
        elif isinstance(neighbors, int):
            if neighbors < 0 or neighbors > 14:
                raise ValueError('neighbors should be in range [0, 14]')
            neighbors = (neighbors,)
        elif isinstance(neighbors, (tuple, list)):
            if not all(isinstance(n, int) for n in neighbors):
                raise TypeError('neighbors should be list or tuple of ints')
            if any(n < 0 or n > 14 for n in neighbors):
                raise ValueError('neighbors should be in range [0, 14]')
            neighbors = tuple(neighbors)
        else:
            raise TypeError('neighbors should be int or list or tuple of ints')
        if p_neighbors is None:
            p_neighbors = ()
        elif isinstance(p_neighbors, int):
            if p_neighbors < 0 or p_neighbors > 14:
                raise ValueError('neighbors should be in range [0, 14]')
            p_neighbors = (p_neighbors,)
        elif isinstance(p_neighbors, (tuple, list)):
            if not all(isinstance(n, int) for n in neighbors):
                raise TypeError('neighbors should be list or tuple of ints')
            if any(n < 0 or n > 14 for n in p_neighbors):
                raise ValueError('neighbors should be in range [0, 14]')
            p_neighbors = tuple(p_neighbors)
        else:
            raise TypeError('neighbors should be int or list or tuple of ints')
        if len(neighbors) != len(p_neighbors):
            raise ValueError('neighbors and p_neighbors should be same length')
        if len(set(zip(neighbors, p_neighbors))) != len(neighbors):
            raise ValueError('paired neighbors and p_neighbors should be unique')

        if hybridization is None:
            hybridization = ()
        elif isinstance(hybridization, int):
            if hybridization < 1 or hybridization > 4:
                raise ValueError('hybridization should be in range [1, 4]')
            hybridization = (hybridization,)
        elif isinstance(hybridization, (tuple, list)):
            if not all(isinstance(h, int) for h in hybridization):
                raise TypeError('hybridizations should be list or tuple of ints')
            if any(h < 1 or h > 4 for h in hybridization):
                raise ValueError('hybridizations should be in range [1, 4]')
            hybridization = tuple(hybridization)
        else:
            raise TypeError('hybridization should be int or list or tuple of ints')
        if p_hybridization is None:
            p_hybridization = ()
        elif isinstance(p_hybridization, int):
            if p_hybridization < 1 or p_hybridization > 4:
                raise ValueError('hybridization should be in range [1, 4]')
            p_hybridization = (p_hybridization,)
        elif isinstance(p_hybridization, (tuple, list)):
            if not all(isinstance(h, int) for h in p_hybridization):
                raise TypeError('hybridizations should be list or tuple of ints')
            if any(h < 1 or h > 4 for h in p_hybridization):
                raise ValueError('hybridizations should be in range [1, 4]')
            p_hybridization = tuple(p_hybridization)
        else:
            raise TypeError('hybridization should be int or list or tuple of ints')
        if len(hybridization) != len(p_hybridization):
            raise ValueError('hybridization and p_hybridization should be same length')
        if len(set(zip(hybridization, p_hybridization))) != len(hybridization):
            raise ValueError('paired hybridization and p_hybridization should be unique')

        if not isinstance(p_charge, int):
            raise TypeError('formal charge should be int in range [-4, 4]')
        if p_charge > 4 or p_charge < -4:
            raise ValueError('formal charge should be in range [-4, 4]')
        if not isinstance(p_is_radical, bool):
            raise TypeError('radical state should be bool')

        if not isinstance(atom, DynamicQueryElement):
            if isinstance(atom, (Element, QueryElement, DynamicElement)):
                atom = DynamicQueryElement.from_atomic_number(atom.atomic_number)(atom.isotope)
            if isinstance(atom, str):
                atom = DynamicQueryElement.from_symbol(atom)()
            elif isinstance(atom, int):
                atom = DynamicQueryElement.from_atomic_number(atom)()
            else:
                raise TypeError('QueryElement object expected')

        _map = super().add_atom(atom, *args, **kwargs)
        if neighbors:
            neighbors, p_neighbors = zip(*sorted(zip(neighbors, p_neighbors)))
        if hybridization:
            hybridization, p_hybridization = zip(*sorted(zip(hybridization, p_hybridization)))

        self._p_charges[_map] = p_charge
        self._p_radicals[_map] = p_is_radical
        self._neighbors[_map] = neighbors
        self._hybridization[_map] = hybridization
        self._p_neighbors[_map] = p_neighbors
        self._p_hybridization[_map] = p_hybridization
        return _map

    def add_bond(self, n, m, bond: Union['cgr.DynamicBond', 'molecule.Bond', int]):
        if not isinstance(bond, cgr.DynamicBond):
            if isinstance(bond, molecule.Bond):
                bond = object.__new__(cgr.DynamicBond)
                bond._DynamicBond__order = bond._DynamicBond__p_order = bond.order
            else:
                bond = cgr.DynamicBond(bond)
        super().add_bond(n, m, bond)

    def delete_atom(self, n):
        super().delete_atom(n)
        del self._p_charges[n]
        del self._p_radicals[n]
        del self._neighbors[n]
        del self._hybridization[n]
        del self._p_neighbors[n]
        del self._p_hybridization[n]

    def remap(self, mapping, *, copy=False) -> 'QueryCGRContainer':
        h = super().remap(mapping, copy=copy)
        mg = mapping.get
        spr = self._p_radicals
        sn = self._neighbors
        sh = self._hybridization
        spn = self._p_neighbors
        sph = self._p_hybridization

        if copy:
            hpc = h._p_charges
            hpr = h._p_radicals
            hn = h._neighbors
            hh = h._hybridization
            hpn = h._p_neighbors
            hph = h._p_hybridization
        else:
            hpc = {}
            hpr = {}
            hn = {}
            hh = {}
            hpn = {}
            hph = {}

        for n, c in self._p_charges.items():
            m = mg(n, n)
            hpc[m] = c
            hpr[m] = spr[n]
            hn[m] = sn[n]
            hpn[m] = spn[n]
            hh[m] = sh[n]
            hph[m] = sph[n]

        if copy:
            return h

        self._p_charges = hpc
        self._p_radicals = hpr
        self._neighbors = hn
        self._hybridization = hh
        self._p_neighbors = hpn
        self._p_hybridization = hph
        return self

    def copy(self, *, meta=True) -> 'QueryCGRContainer':
        copy = super().copy(meta=meta)
        copy._neighbors = self._neighbors.copy()
        copy._hybridization = self._hybridization.copy()
        copy._p_neighbors = self._p_neighbors.copy()
        copy._p_hybridization = self._p_hybridization.copy()
        copy._p_radicals = self._p_radicals.copy()
        copy._p_charges = self._p_charges.copy()
        return copy

    def substructure(self, atoms, *, meta=False) -> 'QueryCGRContainer':
        """
       create substructure containing atoms from atoms list

       :param atoms: list of atoms numbers of substructure
       :param meta: if True metadata will be copied to substructure
       """
        sub, atoms = super().substructure(atoms, meta=meta, sub=self.__class__)
        sa = self._atoms
        spc = self._p_charges
        spr = self._p_radicals
        sn = self._neighbors
        sh = self._hybridization
        spn = self._neighbors
        sph = self._hybridization

        sub._p_charges = {n: spc[n] for n in atoms}
        sub._p_radicals = {n: spr[n] for n in atoms}
        sub._neighbors = {n: sn[n] for n in atoms}
        sub._hybridization = {n: sh[n] for n in atoms}
        sub._neighbors = {n: spn[n] for n in atoms}
        sub._hybridization = {n: sph[n] for n in atoms}

        sub._atoms = ca = {}
        for n in atoms:
            ca[n] = atom = sa[n].copy()
            atom._attach_to_graph(sub, n)
        return sub

    def union(self, other) -> 'QueryCGRContainer':
        if isinstance(other, (QueryCGRContainer, cgr.CGRContainer)):
            u = super().union(other)
            u._p_charges.update(other._p_charges)
            u._p_radicals.update(other._p_radicals)

            if isinstance(other, QueryCGRContainer):
                u._neighbors.update(other._neighbors)
                u._hybridization.update(other._hybridization)
                u._p_neighbors.update(other._p_neighbors)
                u._p_hybridization.update(other._p_hybridization)

                ua = u._atoms
                for n, atom in other._atoms.items():
                    ua[n] = atom = atom.copy()
                    atom._attach_to_graph(u, n)
            else:  # CGRContainer
                un = u._neighbors
                uh = u._hybridization
                upn = u._p_neighbors
                uph = u._p_hybridization
                oh = other._hybridization
                opn = other._p_neighbors
                oph = other._p_hybridization
                for n, m in other._neighbors.items():
                    un[n] = (m,)
                    uh[n] = (oh[n])
                    upn[n] = (opn[n],)
                    uph[n] = (oph[n],)

                ua = u._atoms
                for n, atom in other._atoms.items():
                    ua[n] = atom = DynamicQueryElement.from_atomic_number(atom.atomic_number)(atom.isotope)
                    atom._attach_to_graph(u, n)

            ub = u._bonds
            for n in other._bonds:
                ub[n] = {}
            seen = set()
            for n, m_bond in other._bonds.items():
                seen.add(n)
                for m, bond in m_bond.items():
                    if m not in seen:
                        ub[n][m] = ub[m][n] = bond.copy()
            return u
        elif isinstance(other, (query.QueryContainer, molecule.MoleculeContainer)):
            u = super().union(other)
            u._p_charges.update(other._charges)
            u._p_radicals.update(other._radicals)

            if isinstance(other, query.QueryContainer):
                u._neighbors.update(other._neighbors)
                u._hybridization.update(other._hybridization)
                u._p_neighbors.update(other._neighbors)
                u._p_hybridization.update(other._hybridization)
            else:  # MoleculeContainer
                un = u._neighbors
                uh = u._hybridization
                upn = u._p_neighbors
                uph = u._p_hybridization
                oh = other._hybridization
                for n, m in other._neighbors.items():
                    un[n] = upn[n] = (m,)
                    uh[n] = uph[n] = (oh[n])

            ua = u._atoms
            for n, atom in other._atoms.items():
                ua[n] = atom = DynamicQueryElement.from_atomic_number(atom.atomic_number)(atom.isotope)
                atom._attach_to_graph(u, n)

            ub = u._bonds
            for n in other._bonds:
                ub[n] = {}
            seen = set()
            for n, m_bond in other._bonds.items():
                seen.add(n)
                for m, bond in m_bond.items():
                    if m not in seen:
                        ub[n][m] = ub[m][n] = bc = object.__new__(cgr.DynamicBond)
                        bc._DynamicBond__order = bc._DynamicBond__p_order = bond.order
            return u
        else:
            raise TypeError('Graph expected')

    def get_mapping(self, other):
        if isinstance(other, (QueryCGRContainer, cgr.CGRContainer)):
            return super().get_mapping(other)
        raise TypeError('CGRContainer or QueryCGRContainer expected')

    def __getstate__(self):
        return {'p_charges': self._p_charges, 'p_radicals': self._p_radicals, 'neighbors': self._neighbors,
                'hybridization': self._hybridization, 'p_neighbors': self._p_neighbors,
                'p_hybridization': self._p_hybridization, **super().__getstate__()}

    def __setstate__(self, state):
        super().__setstate__(state)
        self._neighbors = state['neighbors']
        self._hybridization = state['hybridization']
        self._p_neighbors = state['p_neighbors']
        self._p_hybridization = state['p_hybridization']
        self._p_charges = state['p_charges']
        self._p_radicals = state['p_radicals']


__all__ = ['QueryCGRContainer']
