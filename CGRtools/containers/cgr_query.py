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
from typing import List, Union, Tuple, Dict
from . import cgr, molecule, query  # cyclic imports resolve
from .bonds import Bond, DynamicBond
from .common import Graph
from ..algorithms.smiles import QueryCGRSmiles
from ..periodictable import Element, DynamicElement, QueryElement, DynamicQueryElement, AnyElement, DynamicAnyElement


class QueryCGRContainer(Graph, QueryCGRSmiles):
    __slots__ = ('_p_charges', '_p_radicals', '_neighbors', '_hybridizations', '_p_neighbors', '_p_hybridizations')

    def __init__(self):
        self._p_charges: Dict[int, int] = {}
        self._p_radicals: Dict[int, bool] = {}
        self._neighbors: Dict[int, Tuple[int, ...]] = {}
        self._hybridizations: Dict[int, Tuple[int, ...]] = {}
        self._p_neighbors: Dict[int, Tuple[int, ...]] = {}
        self._p_hybridizations: Dict[int, Tuple[int, ...]] = {}
        super().__init__()

    def add_atom(self, atom: Union[DynamicQueryElement, DynamicElement, QueryElement, Element,
                                   AnyElement, DynamicAnyElement, int, str], *args,
                 p_charge: int = 0, p_is_radical: bool = False,
                 neighbors: Union[int, List[int], Tuple[int, ...], None] = None,
                 hybridization: Union[int, List[int], Tuple[int, ...], None] = None,
                 p_neighbors: Union[int, List[int], Tuple[int, ...], None] = None,
                 p_hybridization: Union[int, List[int], Tuple[int, ...], None] = None, **kwargs):
        neighbors = self._validate_neighbors(neighbors)
        p_neighbors = self._validate_neighbors(p_neighbors)
        hybridization = self._validate_hybridization(hybridization)
        p_hybridization = self._validate_hybridization(p_hybridization)
        neighbors, p_neighbors = self._validate_neighbors_pairing(neighbors, p_neighbors)
        hybridization, p_hybridization = self._validate_hybridization_pairing(hybridization, p_hybridization)

        p_charge = self._validate_charge(p_charge)
        p_is_radical = self._validate_radical(p_is_radical)

        if not isinstance(atom, (DynamicQueryElement, DynamicAnyElement)):
            if isinstance(atom, AnyElement):
                atom = DynamicAnyElement()
            elif isinstance(atom, (Element, QueryElement, DynamicElement)):
                atom = DynamicQueryElement.from_atomic_number(atom.atomic_number)(atom.isotope)
            elif isinstance(atom, str):
                atom = DynamicQueryElement.from_symbol(atom)()
            elif isinstance(atom, int):
                atom = DynamicQueryElement.from_atomic_number(atom)()
            else:
                raise TypeError('QueryElement object expected')

        _map = super().add_atom(atom, *args, **kwargs)
        self._p_charges[_map] = p_charge
        self._p_radicals[_map] = p_is_radical
        self._neighbors[_map] = neighbors
        self._hybridizations[_map] = hybridization
        self._p_neighbors[_map] = p_neighbors
        self._p_hybridizations[_map] = p_hybridization
        return _map

    def add_bond(self, n, m, bond: Union[DynamicBond, Bond, int]):
        if not isinstance(bond, DynamicBond):
            if isinstance(bond, Bond):
                bond = object.__new__(DynamicBond)
                bond._DynamicBond__order = bond._DynamicBond__p_order = bond.order
            else:
                bond = DynamicBond(bond)
        super().add_bond(n, m, bond)

    def delete_atom(self, n):
        super().delete_atom(n)
        del self._p_charges[n]
        del self._p_radicals[n]
        del self._neighbors[n]
        del self._hybridizations[n]
        del self._p_neighbors[n]
        del self._p_hybridizations[n]

    def remap(self, mapping, *, copy=False) -> 'QueryCGRContainer':
        h = super().remap(mapping, copy=copy)
        mg = mapping.get
        spr = self._p_radicals
        sn = self._neighbors
        sh = self._hybridizations
        spn = self._p_neighbors
        sph = self._p_hybridizations

        if copy:
            hpc = h._p_charges
            hpr = h._p_radicals
            hn = h._neighbors
            hh = h._hybridizations
            hpn = h._p_neighbors
            hph = h._p_hybridizations
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
        self._hybridizations = hh
        self._p_neighbors = hpn
        self._p_hybridizations = hph
        return self

    def copy(self, **kwargs) -> 'QueryCGRContainer':
        copy = super().copy(**kwargs)
        copy._neighbors = self._neighbors.copy()
        copy._hybridizations = self._hybridizations.copy()
        copy._p_neighbors = self._p_neighbors.copy()
        copy._p_hybridizations = self._p_hybridizations.copy()
        copy._p_radicals = self._p_radicals.copy()
        copy._p_charges = self._p_charges.copy()
        return copy

    def substructure(self, atoms, **kwargs) -> 'QueryCGRContainer':
        """
        create substructure containing atoms from atoms list

        :param atoms: list of atoms numbers of substructure
        :param meta: if True metadata will be copied to substructure
        """
        sub, atoms = super().substructure(atoms, self.__class__, **kwargs)
        sa = self._atoms
        spc = self._p_charges
        spr = self._p_radicals
        sn = self._neighbors
        sh = self._hybridizations
        spn = self._p_neighbors
        sph = self._p_hybridizations

        sub._p_charges = {n: spc[n] for n in atoms}
        sub._p_radicals = {n: spr[n] for n in atoms}
        sub._neighbors = {n: sn[n] for n in atoms}
        sub._hybridizations = {n: sh[n] for n in atoms}
        sub._p_neighbors = {n: spn[n] for n in atoms}
        sub._p_hybridizations = {n: sph[n] for n in atoms}

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
                u._hybridizations.update(other._hybridizations)
                u._p_neighbors.update(other._p_neighbors)
                u._p_hybridizations.update(other._p_hybridizations)

                ua = u._atoms
                for n, atom in other._atoms.items():
                    ua[n] = atom = atom.copy()
                    atom._attach_to_graph(u, n)
            else:  # CGRContainer
                un = u._neighbors
                uh = u._hybridizations
                upn = u._p_neighbors
                uph = u._p_hybridizations
                oh = other._hybridizations
                opn = other._p_neighbors
                oph = other._p_hybridizations
                for n, m in other._neighbors.items():
                    un[n] = (m,)
                    uh[n] = (oh[n],)
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
                u._hybridizations.update(other._hybridizations)
                u._p_neighbors.update(other._neighbors)
                u._p_hybridizations.update(other._hybridizations)
            else:  # MoleculeContainer
                un = u._neighbors
                uh = u._hybridizations
                upn = u._p_neighbors
                uph = u._p_hybridizations
                oh = other._hybridizations
                for n, m in other._neighbors.items():
                    un[n] = upn[n] = (m,)
                    uh[n] = uph[n] = (oh[n],)

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
                        ub[n][m] = ub[m][n] = bc = object.__new__(DynamicBond)
                        bc._DynamicBond__order = bc._DynamicBond__p_order = bond.order
            return u
        else:
            raise TypeError('Graph expected')

    def get_mapping(self, other: Union['QueryCGRContainer', 'cgr.CGRContainer'], **kwargs):
        if isinstance(other, (QueryCGRContainer, cgr.CGRContainer)):
            return super().get_mapping(other, **kwargs)
        raise TypeError('CGRContainer or QueryCGRContainer expected')

    def get_mcs_mapping(self, other: Union['QueryCGRContainer', 'cgr.CGRContainer'], **kwargs):
        if isinstance(other, (QueryCGRContainer, cgr.CGRContainer)):
            return super().get_mcs_mapping(other, **kwargs)
        raise TypeError('CGRContainer or QueryCGRContainer expected')

    @staticmethod
    def _validate_neighbors(neighbors):
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
        return neighbors

    @staticmethod
    def _validate_hybridization(hybridization):
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
        return hybridization

    @staticmethod
    def _validate_neighbors_pairing(neighbors, p_neighbors):
        if len(neighbors) != len(p_neighbors):
            raise ValueError('neighbors and p_neighbors should be same length')
        if neighbors:
            if len(set(zip(neighbors, p_neighbors))) != len(neighbors):
                raise ValueError('paired neighbors and p_neighbors should be unique')
            neighbors, p_neighbors = zip(*sorted(zip(neighbors, p_neighbors)))
        return neighbors, p_neighbors

    @staticmethod
    def _validate_hybridization_pairing(hybridization, p_hybridization):
        if len(hybridization) != len(p_hybridization):
            raise ValueError('hybridization and p_hybridization should be same length')
        if hybridization:
            if len(set(zip(hybridization, p_hybridization))) != len(hybridization):
                raise ValueError('paired hybridization and p_hybridization should be unique')
            hybridization, p_hybridization = zip(*sorted(zip(hybridization, p_hybridization)))
        return hybridization, p_hybridization

    def __getstate__(self):
        return {'p_charges': self._p_charges, 'p_radicals': self._p_radicals, 'neighbors': self._neighbors,
                'hybridizations': self._hybridizations, 'p_neighbors': self._p_neighbors,
                'p_hybridizations': self._p_hybridizations, **super().__getstate__()}

    def __setstate__(self, state):
        super().__setstate__(state)
        self._neighbors = state['neighbors']
        self._hybridizations = state['hybridizations']
        self._p_neighbors = state['p_neighbors']
        self._p_hybridizations = state['p_hybridizations']
        self._p_charges = state['p_charges']
        self._p_radicals = state['p_radicals']


__all__ = ['QueryCGRContainer']
