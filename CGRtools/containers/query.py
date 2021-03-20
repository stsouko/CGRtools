# -*- coding: utf-8 -*-
#
#  Copyright 2018-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from . import molecule  # cyclic imports resolve
from .bonds import Bond, QueryBond
from .common import Graph
from ..algorithms.calculate2d import Calculate2DQuery
from ..algorithms.components import StructureComponents
from ..algorithms.depict import DepictQuery
from ..algorithms.smiles import QuerySmiles
from ..algorithms.stereo import Stereo
from ..periodictable import Element, QueryElement, AnyElement


class QueryContainer(Stereo, Graph, QuerySmiles, StructureComponents, DepictQuery, Calculate2DQuery):
    __slots__ = ('_neighbors', '_hybridizations', '_atoms_stereo', '_cis_trans_stereo', '_allenes_stereo',
                 '_hydrogens', '_rings_sizes', '_heteroatoms')

    def __init__(self):
        self._neighbors: Dict[int, Tuple[int, ...]] = {}
        self._hybridizations: Dict[int, Tuple[int, ...]] = {}
        self._atoms_stereo: Dict[int, bool] = {}
        self._allenes_stereo: Dict[int, bool] = {}
        self._cis_trans_stereo: Dict[Tuple[int, int], bool] = {}
        self._hydrogens: Dict[int, Tuple[int, ...]] = {}
        self._rings_sizes: Dict[int, Tuple[int, ...]] = {}
        self._heteroatoms: Dict[int, Tuple[int, ...]] = {}

        super().__init__()

    def add_atom(self, atom: Union[QueryElement, AnyElement, Element, int, str], *args,
                 neighbors: Union[int, List[int], Tuple[int, ...], None] = None,
                 hybridization: Union[int, List[int], Tuple[int, ...], None] = None,
                 hydrogens: Union[int, List[int], Tuple[int, ...], None] = None,
                 rings_sizes: Union[int, List[int], Tuple[int, ...], None] = None,
                 heteroatoms: Union[int, List[int], Tuple[int, ...], None] = None,
                 **kwargs):
        neighbors = self._validate_neighbors(neighbors)
        hybridization = self._validate_hybridization(hybridization)
        hydrogens = self._validate_neighbors(hydrogens)
        rings_sizes = self._validate_rings(rings_sizes)
        heteroatoms = self._validate_neighbors(heteroatoms)

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
        self._hydrogens[_map] = hydrogens
        self._rings_sizes[_map] = rings_sizes
        self._heteroatoms[_map] = heteroatoms
        return _map

    def add_bond(self, n, m, bond: Union[QueryBond, Bond, int, Tuple[int, ...]]):
        if isinstance(bond, Bond):
            bond = QueryBond.from_bond(bond)
        elif not isinstance(bond, QueryBond):
            bond = QueryBond(bond)

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
        del self._hydrogens[n]
        del self._rings_sizes[n]
        del self._heteroatoms[n]

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
        shg = self._hydrogens
        srs = self._rings_sizes
        sha = self._heteroatoms

        if copy:
            hn = h._neighbors
            hh = h._hybridizations
            has = h._atoms_stereo
            hal = h._allenes_stereo
            hcs = h._cis_trans_stereo
            hhg = h._hydrogens
            hrs = h._rings_sizes
            hha = h._heteroatoms
        else:
            hn = {}
            hh = {}
            has = {}
            hal = {}
            hcs = {}
            hhg = {}
            hrs = {}
            hha = {}

        for n, hyb in self._hybridizations.items():
            m = mg(n, n)
            hn[m] = sn[n]
            hh[m] = hyb
            hhg[m] = shg[n]
            hrs[m] = srs[n]
            hha[m] = sha[n]

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
        self._hydrogens = hhg
        self._rings_sizes = hrs
        self._heteroatoms = hha
        self._atoms_stereo = has
        self._allenes_stereo = hal
        self._cis_trans_stereo = hcs
        return self

    def copy(self, **kwargs) -> 'QueryContainer':
        copy = super().copy(**kwargs)
        copy._neighbors = self._neighbors.copy()
        copy._hybridizations = self._hybridizations.copy()
        copy._hydrogens = self._hydrogens.copy()
        copy._heteroatoms = self._heteroatoms.copy()
        copy._rings_sizes = self._rings_sizes.copy()
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
        sub, atoms = super().substructure(atoms, graph_type=self.__class__, atom_type=QueryElement,
                                          bond_type=QueryBond, **kwargs)
        sa = self._atoms
        sb = self._bonds
        sn = self._neighbors
        sh = self._hybridizations
        shg = self._hydrogens
        srs = self._rings_sizes
        sha = self._heteroatoms

        sub._neighbors = {n: sn[n] for n in atoms}
        sub._hybridizations = {n: sh[n] for n in atoms}
        sub._hydrogens = {n: shg[n] for n in atoms}
        sub._rings_sizes = {n: srs[n] for n in atoms}
        sub._heteroatoms = {n: sha[n] for n in atoms}

        lost = sa.keys() - set(atoms)  # atoms not in substructure
        not_skin = {n for n in atoms if lost.isdisjoint(sb[n])}
        sub._atoms_stereo = {n: s for n, s in self._atoms_stereo.items() if n in not_skin}
        sub._allenes_stereo = {n: s for n, s in self._allenes_stereo.items() if
                               not_skin.issuperset(self._stereo_allenes_paths[n]) and not_skin.issuperset(
                                       x for x in self._stereo_allenes[n] if x)}
        sub._cis_trans_stereo = {nm: s for nm, s in self._cis_trans_stereo.items() if
                                 not_skin.issuperset(self._stereo_cis_trans_paths[nm]) and not_skin.issuperset(
                                         x for x in self._stereo_cis_trans[nm] if x)}
        return sub

    def union(self, other, **kwargs) -> 'QueryContainer':
        if isinstance(other, QueryContainer):
            u, other = super().union(other, atom_type=QueryElement, bond_type=QueryBond, **kwargs)
            u._neighbors.update(other._neighbors)
            u._hybridizations.update(other._hybridizations)
            u._hydrogens.update(other._hydrogens)
            u._rings_sizes.update(other._rings_sizes)
            u._atoms_stereo.update(other._atoms_stereo)
            u._allenes_stereo.update(other._allenes_stereo)
            u._cis_trans_stereo.update(other._cis_trans_stereo)
            u._heteroatoms.update(other._heteroatoms)
            return u
        else:
            raise TypeError('QueryContainer expected')

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
    def _validate_rings(rings):
        if rings is None:
            rings = ()
        elif isinstance(rings, int):
            if rings < 3 and rings != 0:
                raise ValueError('rings should be greater or equal 3. ring equal to zero is no ring atom mark')
            rings = (rings,)
        elif isinstance(rings, (tuple, list)):
            if not all(isinstance(n, int) for n in rings):
                raise TypeError('rings should be list or tuple of ints')
            if any(n < 3 for n in rings):
                raise ValueError('rings should be greater or equal 3')
            if len(set(rings)) != len(rings):
                raise ValueError('rings should be unique')
            rings = tuple(sorted(rings))
        else:
            raise TypeError('rings should be int or list or tuple of ints')
        return rings

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
                'hybridizations': self._hybridizations, 'hydrogens': self._hydrogens,
                'rings_sizes': self._rings_sizes, 'heteroatoms': self._heteroatoms, **super().__getstate__()}

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
                        bn[m] = QueryBond(bond.order)

            state['plane'] = {n: a.xy for n, a in state['node'].items()}
        if 'allenes_stereo' not in state:  # <4.0.22
            state['atoms_stereo'] = {}  # flush existing stereo if exists.
            state['allenes_stereo'] = {}
            state['cis_trans_stereo'] = {}
        if 'hydrogens' not in state:  # <4.1.1
            state['hydrogens'] = {n: () for n in state['atoms']}
            state['rings_sizes'] = {n: () for n in state['atoms']}
        if isinstance(next((x for x in state['bonds'].values() for x in x.values()), None), Bond):  # <4.1.3
            bonds = {}
            for n, m_bond in state['bonds'].items():
                bonds[n] = bn = {}
                for m, bond in m_bond.items():
                    if m in bonds:
                        bn[m] = bonds[m][n]
                    else:
                        b = object.__new__(QueryBond)
                        b._QueryBond__order = (bond._Bond__order,)
                        bn[m] = b
            state['bonds'] = bonds
        if 'heteroatoms' not in state:  # <4.1.4
            state['heteroatoms'] = {n: () for n in state['atoms']}

        super().__setstate__(state)
        self._atoms_stereo = state['atoms_stereo']
        self._allenes_stereo = state['allenes_stereo']
        self._cis_trans_stereo = state['cis_trans_stereo']
        self._neighbors = state['neighbors']
        self._hybridizations = state['hybridizations']
        self._hydrogens = state['hydrogens']
        self._rings_sizes = state['rings_sizes']
        self._heteroatoms = state['heteroatoms']


__all__ = ['QueryContainer']
