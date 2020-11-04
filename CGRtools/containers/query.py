# -*- coding: utf-8 -*-
#
#  Copyright 2018-2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import List, Tuple, Union, Dict
from . import cgr, molecule  # cyclic imports resolve
from .bonds import Bond
from .common import Graph
from ..algorithms.calculate2d import Calculate2DMolecule
from ..algorithms.components import StructureComponents
from ..algorithms.depict import DepictQuery
from ..algorithms.smiles import QuerySmiles
from ..algorithms.stereo import QueryStereo
from ..periodictable import Element, QueryElement, AnyElement


class QueryContainer(QueryStereo, Graph, QuerySmiles, StructureComponents, DepictQuery, Calculate2DMolecule):
    __slots__ = ('_neighbors', '_hybridizations', '_atoms_stereo', '_cis_trans_stereo', '_allenes_stereo')

    def __init__(self):
        self._neighbors: Dict[int, Tuple[int, ...]] = {}
        self._hybridizations: Dict[int, Tuple[int, ...]] = {}
        self._atoms_stereo: Dict[int, bool] = {}
        self._allenes_stereo: Dict[int, bool] = {}
        self._cis_trans_stereo: Dict[Tuple[int, int], bool] = {}

        super().__init__()

    def add_atom(self, atom: Union[QueryElement, AnyElement, Element, int, str], *args,
                 neighbors: Union[int, List[int], Tuple[int, ...], None] = None,
                 hybridization: Union[int, List[int], Tuple[int, ...], None] = None, **kwargs):
        neighbors = self._validate_neighbors(neighbors)
        hybridization = self._validate_hybridization(hybridization)

        if not isinstance(atom, (QueryElement, AnyElement)):
            if isinstance(atom, Element):
                atom = QueryElement.from_atomic_number(atom.atomic_number)(atom.isotope)
            elif isinstance(atom, str):
                atom = QueryElement.from_symbol(atom)()
            elif isinstance(atom, int):
                atom = QueryElement.from_atomic_number(atom)()
            else:
                raise TypeError('QueryElement object expected')

        _map = super().add_atom(atom, *args, **kwargs)
        self._neighbors[_map] = neighbors
        self._hybridizations[_map] = hybridization
        return _map

    def add_bond(self, n, m, bond: Union[Bond, int]):
        if not isinstance(bond, Bond):
            bond = Bond(bond)

        sct = self._stereo_cis_trans_paths  # save
        sa = self._stereo_allenes_paths

        super().add_bond(n, m, bond)
        # remove stereo marks on bonded atoms and all its bonds
        if n in self._atoms_stereo:
            del self._atoms_stereo[n]
        if m in self._atoms_stereo:
            del self._atoms_stereo[m]
        if self._cis_trans_stereo:
            for nm, path in sct.items():
                if (n in path or m in path) and nm in self._cis_trans_stereo:
                    del self._cis_trans_stereo[nm]
        if self._allenes_stereo:
            for c, path in sa.items():
                if (n in path or m in path) and c in self._allenes_stereo:
                    del self._allenes_stereo[c]

    def delete_atom(self, n):
        bonds = set(self._bonds[n])  # save
        sct = self._stereo_cis_trans_paths
        sa = self._stereo_allenes_paths

        super().delete_atom(n)

        del self._neighbors[n]
        del self._hybridizations[n]

        sas = self._atoms_stereo
        if n in sas:
            del sas[n]
        for m in bonds:
            if m in sas:
                del sas[m]
        if self._cis_trans_stereo:
            for nm, path in sct.items():
                if not bonds.isdisjoint(path) and nm in self._cis_trans_stereo:
                    del self._cis_trans_stereo[nm]
        if self._allenes_stereo:
            for c, path in sa.items():
                if not bonds.isdisjoint(path) and c in self._allenes_stereo:
                    del self._allenes_stereo[c]

    def delete_bond(self, n, m):
        sct = self._stereo_cis_trans_paths  # save
        sa = self._stereo_allenes_paths

        super().delete_bond(n, m)

        if n in self._atoms_stereo:
            del self._atoms_stereo[n]
        if m in self._atoms_stereo:
            del self._atoms_stereo[m]
        if self._cis_trans_stereo:
            for nm, path in sct.items():
                if (n in path or m in path) and nm in self._cis_trans_stereo:
                    del self._cis_trans_stereo[nm]
        if self._allenes_stereo:
            for c, path in sa.items():
                if (n in path or m in path) and c in self._allenes_stereo:
                    del self._allenes_stereo[c]

    def remap(self, mapping, *, copy=False) -> 'QueryContainer':
        h = super().remap(mapping, copy=copy)
        mg = mapping.get
        sn = self._neighbors

        if copy:
            hn = h._neighbors
            hh = h._hybridizations
            has = h._atoms_stereo
            hal = h._allenes_stereo
            hcs = h._cis_trans_stereo
        else:
            hn = {}
            hh = {}
            has = {}
            hal = {}
            hcs = {}

        for n, hyb in self._hybridizations.items():
            m = mg(n, n)
            hn[m] = sn[n]
            hh[m] = hyb

        for n, stereo in self._atoms_stereo.items():
            has[mg(n, n)] = stereo
        for n, stereo in self._allenes_stereo.items():
            hal[mg(n, n)] = stereo
        for (n, m), stereo in self._cis_trans_stereo.items():
            hcs[(mg(n, n), mg(m, m))] = stereo

        if copy:
            return h

        self._neighbors = hn
        self._hybridizations = hh
        self._atoms_stereo = has
        self._allenes_stereo = hal
        self._cis_trans_stereo = hcs
        return self

    def copy(self, **kwargs) -> 'QueryContainer':
        copy = super().copy(**kwargs)
        copy._neighbors = self._neighbors.copy()
        copy._hybridizations = self._hybridizations.copy()
        copy._atoms_stereo = self._atoms_stereo.copy()
        copy._allenes_stereo = self._allenes_stereo.copy()
        copy._cis_trans_stereo = self._cis_trans_stereo.copy()
        return copy

    def substructure(self, atoms, **kwargs) -> 'QueryContainer':
        """
        create substructure containing atoms from atoms list

        :param atoms: list of atoms numbers of substructure
        :param meta: if True metadata will be copied to substructure
        """
        sub, atoms = super().substructure(atoms, self.__class__, **kwargs)
        sa = self._atoms
        sb = self._bonds
        sn = self._neighbors
        sh = self._hybridizations

        sub._neighbors = {n: sn[n] for n in atoms}
        sub._hybridizations = {n: sh[n] for n in atoms}

        lost = {n for n, a in sa.items() if a.atomic_number != 1} - set(atoms)  # atoms not in substructure
        not_skin = {n for n in atoms if lost.isdisjoint(sb[n])}
        sub._atoms_stereo = {n: s for n, s in self._atoms_stereo.items() if n in not_skin}
        sub._allenes_stereo = {n: s for n, s in self._allenes_stereo.items() if
                               not_skin.issuperset(self._stereo_allenes_paths[n]) and not_skin.issuperset(
                                       x for x in self._stereo_allenes[n] if x)}
        sub._cis_trans_stereo = {nm: s for nm, s in self._cis_trans_stereo.items() if
                                 not_skin.issuperset(self._stereo_cis_trans_paths[nm]) and not_skin.issuperset(
                                         x for x in self._stereo_cis_trans[nm] if x)}

        sub._atoms = ca = {}
        for n in atoms:
            atom = sa[n].copy()
            ca[n] = atom
            atom._attach_to_graph(sub, n)
        return sub

    def union(self, other, **kwargs) -> 'QueryContainer':
        if isinstance(other, (QueryContainer, molecule.MoleculeContainer)):
            u, other = super().union(other, **kwargs)
            if isinstance(other, QueryContainer):
                u._neighbors.update(other._neighbors)
                u._hybridizations.update(other._hybridizations)

                ua = u._atoms
                for n, atom in other._atoms.items():
                    atom = atom.copy()
                    ua[n] = atom
                    atom._attach_to_graph(u, n)
            else:
                un = u._neighbors
                uh = u._hybridizations
                oh = other._hybridizations
                for n, m_bond in other._bonds.items():
                    un[n] = (len(m_bond),)
                    uh[n] = (oh[n],)

                ua = u._atoms
                for n, atom in other._atoms.items():
                    atom = QueryElement.from_atomic_number(atom.atomic_number)(atom.isotope)
                    ua[n] = atom
                    atom._attach_to_graph(u, n)

            ub = u._bonds
            for n, m_bond in other._bonds.items():
                ub[n] = ubn = {}
                for m, bond in m_bond.items():
                    if m in ub:  # bond partially exists. need back-connection.
                        ubn[m] = ub[m][n]
                    else:
                        ubn[m] = bond.copy()

            u._atoms_stereo.update(other._atoms_stereo)
            u._allenes_stereo.update(other._allenes_stereo)
            u._cis_trans_stereo.update(other._cis_trans_stereo)
            return u
        elif isinstance(other, cgr.CGRContainer):
            raise TypeError('QueryContainer and CGRContainer unite impossible')
        elif isinstance(other, Graph):
            return other.union(self, **kwargs)
        else:
            raise TypeError('Graph expected')

    def get_mapping(self, other: Union['QueryContainer', 'molecule.MoleculeContainer'], **kwargs):
        if isinstance(other, (QueryContainer, molecule.MoleculeContainer)):
            return super().get_mapping(other, **kwargs)
        raise TypeError('MoleculeContainer or QueryContainer expected')

    def get_mcs_mapping(self, other: Union['QueryContainer', 'molecule.MoleculeContainer'], **kwargs):
        if isinstance(other, (QueryContainer, molecule.MoleculeContainer)):
            return super().get_mcs_mapping(other, **kwargs)
        raise TypeError('MoleculeContainer or QueryContainer expected')

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
            if len(set(neighbors)) != len(neighbors):
                raise ValueError('neighbors should be unique')
            neighbors = tuple(sorted(neighbors))
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
            if len(set(hybridization)) != len(hybridization):
                raise ValueError('hybridizations should be unique')
            hybridization = tuple(sorted(hybridization))
        else:
            raise TypeError('hybridization should be int or list or tuple of ints')
        return hybridization

    def __getstate__(self):
        return {'atoms_stereo': self._atoms_stereo, 'allenes_stereo': self._allenes_stereo,
                'cis_trans_stereo': self._cis_trans_stereo, 'neighbors': self._neighbors,
                'hybridizations': self._hybridizations, **super().__getstate__()}

    def __setstate__(self, state):
        if 'node' in state:  # 3.1 compatibility.
            state['atoms'] = {n: a.atom for n, a in state['node'].items()}
            state['charges'] = {n: a.charge for n, a in state['node'].items()}
            state['radicals'] = {n: a.is_radical for n, a in state['node'].items()}
            state['neighbors'] = {n: a.neighbors for n, a in state['node'].items()}
            state['hybridizations'] = {n: a.hybridization for n, a in state['node'].items()}
            state['parsed_mapping'] = {}
            state['bonds'] = bonds = {}
            for n, m_bond in state['adj'].items():
                bonds[n] = bn = {}
                for m, bond in m_bond.items():
                    if m in bonds:
                        bn[m] = bonds[m][n]
                    else:
                        bn[m] = Bond(bond.order)

            state['plane'] = {n: a.xy for n, a in state['node'].items()}
            state['atoms_stereo'] = {}
            state['allenes_stereo'] = {}
            state['cis_trans_stereo'] = {}
        elif 'allenes_stereo' not in state:  # <4.0.22
            state['atoms_stereo'] = {}  # flush existing stereo if exists.
            state['allenes_stereo'] = {}
            state['cis_trans_stereo'] = {}

        super().__setstate__(state)
        self._atoms_stereo = state['atoms_stereo']
        self._allenes_stereo = state['allenes_stereo']
        self._cis_trans_stereo = state['cis_trans_stereo']
        self._neighbors = state['neighbors']
        self._hybridizations = state['hybridizations']


__all__ = ['QueryContainer']
