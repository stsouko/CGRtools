# -*- coding: utf-8 -*-
#
#  Copyright 2017-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from CachedMethods import cached_args_method, cached_property
from collections import defaultdict, Counter
from typing import List, Union, Tuple, Optional, Dict
from . import cgr, query  # cyclic imports resolve
from .bonds import Bond, DynamicBond, QueryBond
from .common import Graph
from ..algorithms.aromatics import Aromatize
from ..algorithms.calculate2d import Calculate2DMolecule
from ..algorithms.components import StructureComponents
from ..algorithms.depict import DepictMolecule
from ..algorithms.huckel import Huckel
from ..algorithms.smiles import MoleculeSmiles
from ..algorithms.standardize import Standardize
from ..algorithms.stereo import MoleculeStereo
from ..algorithms.tautomers import Tautomers
from ..algorithms.x3dom import X3domMolecule
from ..exceptions import ValenceError, MappingError
from ..periodictable import Element, QueryElement


class MoleculeContainer(MoleculeStereo, Graph, Aromatize, Standardize, MoleculeSmiles, StructureComponents,
                        DepictMolecule, Calculate2DMolecule, Tautomers, Huckel, X3domMolecule):
    __slots__ = ('_conformers', '_hybridizations', '_atoms_stereo', '_hydrogens', '_cis_trans_stereo',
                 '_allenes_stereo')

    def __init__(self):
        self._conformers: List[Dict[int, Tuple[float, float, float]]] = []
        self._hybridizations: Dict[int, int] = {}
        self._hydrogens: Dict[int, Optional[int]] = {}
        self._atoms_stereo: Dict[int, bool] = {}
        self._allenes_stereo: Dict[int, bool] = {}
        self._cis_trans_stereo: Dict[Tuple[int, int], bool] = {}

        super().__init__()

    def add_atom(self, atom: Union[Element, int, str], *args, charge=0, is_radical=False, **kwargs):
        """
        Add new atom.
        """
        if not isinstance(atom, Element):
            if isinstance(atom, str):
                atom = Element.from_symbol(atom)()
            elif isinstance(atom, int):
                atom = Element.from_atomic_number(atom)()
            else:
                raise TypeError('Element object expected')

        _map = super().add_atom(atom, *args, charge=charge, is_radical=is_radical, **kwargs)
        self._hybridizations[_map] = 1
        self._conformers.clear()  # clean conformers. need full recalculation for new system

        if atom.atomic_number != 1:
            try:
                rules = atom.valence_rules(charge, is_radical, 0)
            except ValenceError:
                self._hydrogens[_map] = None
            else:
                for s, d, h in rules:
                    if h and not s:
                        self._hydrogens[_map] = h
                        break
                else:
                    self._hydrogens[_map] = 0
        else:
            self._hydrogens[_map] = 0
        return _map

    def add_bond(self, n, m, bond: Union[Bond, int]):
        """
        Connect atoms with bonds.

        For Thiele forms of molecule causes invalidation of internal state.
        Implicit hydrogens marks will not be set if atoms in aromatic rings.
        Call `kekule()` and `thiele()` in sequence to fix marks.
        """
        if not isinstance(bond, Bond):
            bond = Bond(bond)

        super().add_bond(n, m, bond)
        self._conformers.clear()  # clean conformers. need full recalculation for new system

        self._calc_implicit(n)
        self._calc_implicit(m)

        if bond.order != 1:  # 1 is neutral. 8 is rare. skip.
            self._calc_hybridization(n)
            self._calc_hybridization(m)

        if self._atoms[n].atomic_number != 1 and self._atoms[m].atomic_number != 1:  # not hydrogen
            # fix stereo if formed not to hydrogen bond
            self._fix_stereo()

    def delete_atom(self, n):
        """
        Remove atom.

        For Thiele forms of molecule causes invalidation of internal state.
        Implicit hydrogens marks will not be set if atoms in aromatic rings.
        Call `kekule()` and `thiele()` in sequence to fix marks.
        """
        old_bonds = self._bonds[n]  # save bonds
        isnt_hydrogen = self._atoms[n].atomic_number != 1
        super().delete_atom(n)

        del self._hybridizations[n]
        del self._hydrogens[n]
        self._conformers.clear()  # clean conformers. need full recalculation for new system

        for m in old_bonds:
            self._calc_implicit(m)

        for m in old_bonds:
            self._calc_hybridization(m)
        if isnt_hydrogen:  # hydrogen atom not used for stereo coding
            self._fix_stereo()

    def delete_bond(self, n, m):
        """
        Disconnect atoms.

        For Thiele forms of molecule causes invalidation of internal state.
        Implicit hydrogens marks will not be set if atoms in aromatic rings.
        Call `kekule()` and `thiele()` in sequence to fix marks.
        """
        super().delete_bond(n, m)
        self._conformers.clear()  # clean conformers. need full recalculation for new system

        self._calc_implicit(n)
        self._calc_implicit(m)
        self._calc_hybridization(n)
        self._calc_hybridization(m)

        if self._atoms[n].atomic_number != 1 and self._atoms[m].atomic_number != 1:
            self._fix_stereo()

    @cached_args_method
    def neighbors(self, n: int) -> int:
        """number of neighbors atoms excluding any-bonded"""
        return sum(b.order != 8 for b in self._bonds[n].values())

    @cached_args_method
    def heteroatoms(self, n: int) -> int:
        """
        Number of neighbored heteroatoms (not carbon or hydrogen)
        """
        atoms = self._atoms
        return sum(atoms[m].atomic_number not in (1, 6) for m in self._bonds[n])

    def remap(self, mapping, *, copy=False) -> 'MoleculeContainer':
        h = super().remap(mapping, copy=copy)
        mg = mapping.get
        shg = self._hydrogens

        if copy:
            hh = h._hybridizations
            hhg = h._hydrogens
            hc = h._conformers
            has = h._atoms_stereo
            hal = h._allenes_stereo
            hcs = h._cis_trans_stereo
        else:
            hh = {}
            hhg = {}
            hc = []
            has = {}
            hal = {}
            hcs = {}

        for n, hyb in self._hybridizations.items():
            m = mg(n, n)
            hh[m] = hyb
            hhg[m] = shg[n]

        hc.extend({mg(n, n): x for n, x in c.items()} for c in self._conformers)

        for n, stereo in self._atoms_stereo.items():
            has[mg(n, n)] = stereo
        for n, stereo in self._allenes_stereo.items():
            hal[mg(n, n)] = stereo
        for (n, m), stereo in self._cis_trans_stereo.items():
            hcs[(mg(n, n), mg(m, m))] = stereo

        if copy:
            return h

        self._hybridizations = hh
        self._hydrogens = hhg
        self._conformers = hc
        self._atoms_stereo = has
        self._allenes_stereo = hal
        self._cis_trans_stereo = hcs
        return self

    def copy(self, **kwargs) -> 'MoleculeContainer':
        copy = super().copy(**kwargs)
        copy._hybridizations = self._hybridizations.copy()
        copy._hydrogens = self._hydrogens.copy()
        copy._conformers = [c.copy() for c in self._conformers]
        copy._atoms_stereo = self._atoms_stereo.copy()
        copy._allenes_stereo = self._allenes_stereo.copy()
        copy._cis_trans_stereo = self._cis_trans_stereo.copy()
        return copy

    def substructure(self, atoms, *, as_query: bool = False, **kwargs) -> Union['MoleculeContainer',
                                                                                'query.QueryContainer']:
        """
        Create substructure containing atoms from atoms list.

        For Thiele forms of molecule In Molecule substructure causes invalidation of internal state.
        Implicit hydrogens marks will not be set if atoms in aromatic rings.
        Call `kekule()` and `thiele()` in sequence to fix marks.

        :param atoms: list of atoms numbers of substructure
        :param meta: if True metadata will be copied to substructure
        :param as_query: return Query object based on graph substructure
        """
        sub, atoms = super().substructure(atoms, graph_type=query.QueryContainer if as_query else self.__class__,
                                          atom_type=QueryElement if as_query else Element,
                                          bond_type=QueryBond if as_query else Bond, **kwargs)
        if as_query:
            sa = self._atoms
            sb = self._bonds
            sh = self._hybridizations
            shg = self._hydrogens
            sn = self.neighbors
            rs = self.atoms_rings_sizes.copy()

            lost = {n for n, a in sa.items() if a.atomic_number != 1} - set(atoms)  # atoms not in substructure
            not_skin = {n for n in atoms if lost.isdisjoint(sb[n])}
            sub._atoms_stereo = {n: s for n, s in self._atoms_stereo.items() if n in not_skin}
            sub._allenes_stereo = {n: s for n, s in self._allenes_stereo.items()
                                   if not_skin.issuperset(self._stereo_allenes_paths[n]) and
                                      not_skin.issuperset(x for x in self._stereo_allenes[n] if x)}
            sub._cis_trans_stereo = {nm: s for nm, s in self._cis_trans_stereo.items()
                                     if not_skin.issuperset(self._stereo_cis_trans_paths[nm]) and
                                        not_skin.issuperset(x for x in self._stereo_cis_trans[nm] if x)}

            sub._neighbors = {n: (sn(n),) for n in atoms}
            sub._hybridizations = {n: (sh[n],) for n in atoms}
            sub._hydrogens = {n: () if shg[n] is None else (shg[n],) for n in atoms}
            sub._rings_sizes = {n: rs.get(n, ()) for n in atoms}
            sub._heteroatoms = {n: () for n in atoms}
        else:
            sub._conformers = [{n: c[n] for n in atoms} for c in self._conformers]

            # recalculate query marks
            sub._hybridizations = {}
            sub._hydrogens = {}
            for n in atoms:
                sub._calc_hybridization(n)
                sub._calc_implicit(n)
            # fix_stereo will repair data
            sub._atoms_stereo = self._atoms_stereo
            sub._allenes_stereo = self._allenes_stereo
            sub._cis_trans_stereo = self._cis_trans_stereo
            sub._fix_stereo()
        return sub

    def union(self, other, **kwargs):
        if isinstance(other, MoleculeContainer):
            u, other = super().union(other, atom_type=Element, bond_type=Bond, **kwargs)
            u._conformers.clear()
            u._hybridizations.update(other._hybridizations)
            u._hydrogens.update(other._hydrogens)
            u._atoms_stereo.update(other._atoms_stereo)
            u._allenes_stereo.update(other._allenes_stereo)
            u._cis_trans_stereo.update(other._cis_trans_stereo)
            return u
        elif isinstance(other, Graph):
            return other.union(self, **kwargs)
        else:
            raise TypeError('MoleculeContainer expected')

    def compose(self, other: Union['MoleculeContainer', 'cgr.CGRContainer']) -> 'cgr.CGRContainer':
        """
        Compose 2 graphs to CGR.
        """
        sa = self._atoms
        sc = self._charges
        sr = self._radicals
        sp = self._plane
        sb = self._bonds

        bonds = []
        adj = defaultdict(lambda: defaultdict(lambda: [None, None]))

        if isinstance(other, MoleculeContainer):
            oa = other._atoms
            oc = other._charges
            or_ = other._radicals
            op = other._plane
            ob = other._bonds
            common = sa.keys() & other
            h = cgr.CGRContainer()
            atoms = h._atoms

            for n in sa.keys() - common:  # cleavage atoms
                h.add_atom(sa[n], n, charge=sc[n], is_radical=sr[n], xy=sp[n], p_charge=sc[n], p_is_radical=sr[n])
                for m, bond in sb[n].items():
                    if m not in atoms:
                        if m in common:  # bond to common atoms is broken bond
                            order = bond.order
                            bond = object.__new__(DynamicBond)
                            bond._DynamicBond__order, bond._DynamicBond__p_order = order, None
                        bonds.append((n, m, bond))
            for n in other._atoms.keys() - common:  # coupling atoms
                h.add_atom(oa[n], n, charge=oc[n], is_radical=or_[n], xy=op[n], p_charge=oc[n], p_is_radical=or_[n])
                for m, bond in ob[n].items():
                    if m not in atoms:
                        if m in common:  # bond to common atoms is formed bond
                            order = bond.order
                            bond = object.__new__(DynamicBond)
                            bond._DynamicBond__order, bond._DynamicBond__p_order = None, order
                        bonds.append((n, m, bond))
            for n in common:
                an = adj[n]
                for m, bond in sb[n].items():
                    if m in common:
                        an[m][0] = bond.order
                for m, bond in ob[n].items():
                    if m in common:
                        an[m][1] = bond.order
            for n in common:
                san = sa[n]
                if san.atomic_number != oa[n].atomic_number or san.isotope != oa[n].isotope:
                    raise MappingError(f'atoms with number {{{n}}} not equal')
                h.add_atom(san, n, charge=sc[n], is_radical=sr[n], xy=sp[n], p_charge=oc[n], p_is_radical=or_[n])
                for m, (o1, o2) in adj[n].items():
                    if m not in atoms:
                        bond = object.__new__(DynamicBond)
                        bond._DynamicBond__order, bond._DynamicBond__p_order = o1, o2
                        bonds.append((n, m, bond))
        elif isinstance(other, cgr.CGRContainer):
            oa = other._atoms
            oc = other._charges
            or_ = other._radicals
            opc = other._p_charges
            opr = other._p_radicals
            op = other._plane
            ob = other._bonds
            common = sa.keys() & other
            h = other.__class__()  # subclasses support
            atoms = h._atoms

            for n in sa.keys() - common:  # cleavage atoms
                h.add_atom(sa[n], n, charge=sc[n], is_radical=sr[n], xy=sp[n], p_charge=sc[n], p_is_radical=sr[n])
                for m, bond in sb[n].items():
                    if m not in atoms:
                        if m in common:  # bond to common atoms is broken bond
                            order = bond.order
                            bond = object.__new__(DynamicBond)
                            bond._DynamicBond__order, bond._DynamicBond__p_order = order, None
                        bonds.append((n, m, bond))
            for n in other._atoms.keys() - common:  # coupling atoms
                h.add_atom(oa[n].copy(), n, charge=oc[n], is_radical=or_[n], xy=op[n], p_charge=opc[n],
                           p_is_radical=opr[n])
                for m, bond in ob[n].items():
                    if m not in atoms:
                        if m in common:  # bond to common atoms is formed bond
                            order = bond.p_order
                            if order:  # skip broken bond. X>None => None>None
                                bond = object.__new__(DynamicBond)
                                bond._DynamicBond__order, bond._DynamicBond__p_order = None, order
                                bonds.append((n, m, bond))
                        else:
                            bonds.append((n, m, bond))
            for n in common:
                an = adj[n]
                for m, bond in sb[n].items():
                    if m in common:
                        an[m][0] = bond.order
                for m, bond in ob[n].items():
                    if m in an or m in common and bond.p_order:
                        # self has nm bond or other bond not broken
                        an[m][1] = bond.p_order
            for n in common:
                san = sa[n]
                if san.atomic_number != oa[n].atomic_number or san.isotope != oa[n].isotope:
                    raise MappingError(f'atoms with number {{{n}}} not equal')
                h.add_atom(san, n, charge=sc[n], is_radical=sr[n], xy=sp[n], p_charge=opc[n], p_is_radical=opr[n])
                for m, (o1, o2) in adj[n].items():
                    if m not in atoms:
                        bond = object.__new__(DynamicBond)
                        bond._DynamicBond__order, bond._DynamicBond__p_order = o1, o2
                        bonds.append((n, m, bond))
        else:
            raise TypeError('MoleculeContainer or CGRContainer expected')

        for n, m, bond in bonds:
            h.add_bond(n, m, bond)
        return h

    def __xor__(self, other):
        """
        G ^ H is CGR generation
        """
        return self.compose(other)

    def get_fast_mapping(self, other: 'MoleculeContainer') -> Optional[Dict[int, int]]:
        """
        Get self to other fast (suboptimal) structure mapping.
        Only one possible atoms mapping returned.
        Effective only for big molecules.
        """
        if isinstance(other, MoleculeContainer):
            if len(self) != len(other):
                return
            so = self.smiles_atoms_order
            oo = other.smiles_atoms_order
            if self != other:
                return
            return dict(zip(so, oo))
        raise TypeError('MoleculeContainer expected')

    def get_mapping(self, other: 'MoleculeContainer', **kwargs):
        if isinstance(other, MoleculeContainer):
            return super().get_mapping(other, **kwargs)
        raise TypeError('MoleculeContainer expected')

    def get_mcs_mapping(self, other: 'MoleculeContainer', **kwargs):
        if isinstance(other, MoleculeContainer):
            return super().get_mcs_mapping(other, **kwargs)
        raise TypeError('MoleculeContainer expected')

    @cached_property
    def molecular_charge(self) -> int:
        """
        Total charge of molecule
        """
        return sum(self._charges.values())

    @cached_property
    def is_radical(self) -> bool:
        """
        True if at least one atom is radical
        """
        return any(self._radicals.values())

    def __int__(self):
        """
        Total charge of molecule
        """
        return self.molecular_charge

    @cached_property
    def molecular_mass(self):
        return sum(x.atomic_mass for x in self._atoms.values())

    def __float__(self):
        return self.molecular_mass

    @cached_property
    def brutto(self) -> Dict[str, int]:
        """Counted atoms dict"""
        return Counter(x.atomic_symbol for x in self._atoms.values())

    @cached_args_method
    def _explicit_hydrogens(self, n: int) -> int:
        """
        Number of explicit hydrogen atoms connected to atom.

        Take into account any type of bonds with hydrogen atoms.
        """
        atoms = self._atoms
        return sum(atoms[m].atomic_number == 1 for m in self._bonds[n])

    @cached_args_method
    def _total_hydrogens(self, n: int) -> int:
        return self._hydrogens[n] + self._explicit_hydrogens(n)

    def _calc_implicit(self, n: int):
        atoms = self._atoms
        atom = atoms[n]
        if atom.atomic_number != 1:
            charge = self._charges[n]
            is_radical = self._radicals[n]
            explicit_sum = 0
            explicit_dict = defaultdict(int)
            for m, bond in self._bonds[n].items():
                order = bond.order
                if order == 4:  # aromatic rings not supported
                    self._hydrogens[n] = None
                    return
                elif order != 8:  # any bond used for complexes
                    explicit_sum += order
                    explicit_dict[(order, atoms[m].atomic_number)] += 1
            try:
                rules = atom.valence_rules(charge, is_radical, explicit_sum)
            except ValenceError:
                self._hydrogens[n] = None
                return
            for s, d, h in rules:
                if h and s.issubset(explicit_dict) and all(explicit_dict[k] >= c for k, c in d.items()):
                    self._hydrogens[n] = h
                    return
        self._hydrogens[n] = 0

    def _calc_hybridization(self, n: int):
        hybridization = 1
        for bond in self._bonds[n].values():
            order = bond.order
            if order == 4:
                self._hybridizations[n] = 4
                return
            elif order == 3:
                if hybridization != 3:
                    hybridization = 3
            elif order == 2:
                if hybridization == 2:
                    hybridization = 3
                elif hybridization == 1:
                    hybridization = 2
        self._hybridizations[n] = hybridization

    def __getstate__(self):
        return {'conformers': self._conformers, 'hydrogens': self._hydrogens, 'atoms_stereo': self._atoms_stereo,
                'allenes_stereo': self._allenes_stereo, 'cis_trans_stereo': self._cis_trans_stereo,
                **super().__getstate__()}

    def __setstate__(self, state):
        if '_BaseContainer__meta' in state:  # 2.8 reverse compatibility
            state['atoms'] = {n: Element.from_symbol(a['element'])(a.get('isotope')) for n, a in state['_node'].items()}
            state['charges'] = {n: a['s_charge'] for n, a in state['_node'].items()}
            state['radicals'] = {n: bool(a['s_radical']) for n, a in state['_node'].items()}
            state['plane'] = {n: (a['s_x'], a['s_y']) for n, a in state['node'].items()}

            state['bonds'] = bonds = {}
            for n, m_bond in state['_adj'].items():
                bonds[n] = bn = {}
                for m, bond in m_bond.items():
                    if m in bonds:
                        bn[m] = bonds[m][n]
                    else:
                        bn[m] = Bond(bond['s_bond'] if bond['s_bond'] != 9 else 5)

            state['conformers'] = []
            state['atoms_stereo'] = {}
            state['allenes_stereo'] = {}
            state['cis_trans_stereo'] = {}
            state['meta'] = state['_BaseContainer__meta']
            state['parsed_mapping'] = {}
        elif 'node' in state:  # 3.1 compatibility.
            state['atoms'] = {n: a.atom.copy() for n, a in state['node'].items()}
            state['charges'] = {n: a.charge for n, a in state['node'].items()}
            state['radicals'] = {n: a.is_radical for n, a in state['node'].items()}
            state['plane'] = {n: a.xy for n, a in state['node'].items()}
            state['parsed_mapping'] = {n: a.parsed_mapping for n, a in state['node'].items() if a.parsed_mapping}

            state['bonds'] = bonds = {}
            for n, m_bond in state['adj'].items():
                bonds[n] = bn = {}
                for m, bond in m_bond.items():
                    if m in bonds:
                        bn[m] = bonds[m][n]
                    else:
                        bn[m] = Bond(bond.order)

            state['conformers'] = []
            state['atoms_stereo'] = {}
            state['allenes_stereo'] = {}
            state['cis_trans_stereo'] = {}
        elif 'allenes_stereo' not in state:  # <4.0.22
            state['atoms_stereo'] = {}  # flush existing stereo if exists.
            state['allenes_stereo'] = {}
            state['cis_trans_stereo'] = {}

        super().__setstate__(state)
        self._conformers = state['conformers']
        self._atoms_stereo = state['atoms_stereo']
        self._allenes_stereo = state['allenes_stereo']
        self._cis_trans_stereo = state['cis_trans_stereo']

        if 'hydrogens' not in state:  # <4.0.38
            self._hydrogens = {}
            for n in state['atoms']:
                self._calc_implicit(n)
        else:
            self._hydrogens = state['hydrogens']

        # restore query marks
        self._hybridizations = {}
        for n in state['atoms']:
            self._calc_hybridization(n)


__all__ = ['MoleculeContainer']
