# -*- coding: utf-8 -*-
#
#  Copyright 2017-2019 Ramil Nugmanov <stsouko@live.ru>
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
from collections import defaultdict
from typing import List, Union
from . import cgr, query  # cyclic imports resolve
from .common import Graph
from ..algorithms.aromatics import Aromatize
from ..algorithms.depict import DepictMolecule
from ..algorithms.smiles import MoleculeSmiles
from ..exceptions import ValenceError, MappingError
from ..periodictable import Element, QueryElement


class Bond:
    __slots__ = ('__order',)

    def __init__(self, order):
        if not isinstance(order, int):
            raise TypeError('invalid order value')
        if order not in (1, 4, 2, 3, 8):
            raise ValueError('order should be from [1, 2, 3, 4, 8]')
        self.__order = order

    def __eq__(self, other):
        if isinstance(other, Bond):
            return self.__order == other.order
        return False

    def __repr__(self):
        return f'{self.__class__.__name__}({self.__order})'

    def __int__(self):
        return self.__order

    @property
    def order(self):
        return self.__order

    def copy(self):
        copy = object.__new__(self.__class__)
        copy._Bond__order = self.__order
        return copy


class MoleculeContainer(Graph, Aromatize, MoleculeSmiles, DepictMolecule):
    __slots__ = ('_conformers', '_neighbors', '_hybridizations', '_atoms_stereo', '_bonds_stereo', '_hydrogens')

    def __init__(self):
        self._conformers = []
        self._neighbors = {}
        self._hybridizations = {}
        self._hydrogens = {}
        self._atoms_stereo = {}
        self._bonds_stereo = {}
        super().__init__()

    def add_atom(self, atom: Union[Element, int, str], *args, charge=0, is_radical=False, **kwargs):
        if not isinstance(atom, Element):
            if isinstance(atom, str):
                atom = Element.from_symbol(atom)()
            elif isinstance(atom, int):
                atom = Element.from_atomic_number(atom)()
            else:
                raise TypeError('Element object expected')

        _map = super().add_atom(atom, *args, charge=charge, is_radical=is_radical, **kwargs)
        self._neighbors[_map] = 0
        self._hybridizations[_map] = 1
        self._bonds_stereo[_map] = {}
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
        if not isinstance(bond, Bond):
            bond = Bond(bond)

        super().add_bond(n, m, bond)
        self._conformers.clear()  # clean conformers. need full recalculation for new system
        sbs = self._bonds_stereo

        self._hydrogens[n] = self._calc_implicit(n)
        self._hydrogens[m] = self._calc_implicit(m)

        # calc query marks dynamically.
        if self._atoms[n].atomic_number != 1:  # not hydrogen
            self._neighbors[m] += 1
            self._hybridizations[m] = self._calc_hybridization(m)

            try:  # remove stereo marks on bonded atoms and all its bonds
                del self._atoms_stereo[m]
            except KeyError:
                pass
            if sbs[m]:
                for x in sbs[m]:
                    del sbs[x][m]  # remove incoming
                sbs[m] = {}  # remove outgoing
        if self._atoms[m].atomic_number != 1:  # not hydrogen
            self._neighbors[n] += 1
            self._hybridizations[n] = self._calc_hybridization(n)

            try:  # remove atom stereo state
                del self._atoms_stereo[n]
            except KeyError:
                pass
            if sbs[n]:
                for x in sbs[n]:
                    del sbs[x][n]
                sbs[n] = {}

    def delete_atom(self, n):
        old_bonds = self._bonds[n]  # save bonds
        atoms = self._atoms
        isnt_hydrogen = atoms[n].atomic_number != 1
        super().delete_atom(n)

        sn = self._neighbors
        sh = self._hybridizations
        shg = self._hydrogens
        sas = self._atoms_stereo
        sbs = self._bonds_stereo

        del sn[n]
        del sh[n]
        del shg[n]
        self._conformers.clear()  # clean conformers. need full recalculation for new system

        for m in old_bonds:
            shg[m] = self._calc_implicit(m)

        if isnt_hydrogen:  # neighbors query marks fix. ignore removed hydrogen
            for m in old_bonds:
                sh[m] = self._calc_hybridization(m)
                sn[m] -= 1

            # remove stereo marks on deleted atoms and all its neighbors
            try:
                del sas[n]
            except KeyError:
                pass
            for m in sbs.pop(n):
                del sbs[m][n]
                try:
                    del sas[m]
                except KeyError:
                    pass

    def delete_bond(self, n, m):
        super().delete_bond(n, m)
        self._conformers.clear()  # clean conformers. need full recalculation for new system

        atoms = self._atoms
        sbs = self._bonds_stereo
        sh = self._hybridizations
        sn = self._neighbors

        self._hydrogens[n] = self._calc_implicit(n)
        self._hydrogens[m] = self._calc_implicit(m)

        # neighbors query marks fix. ignore removed hydrogen
        if atoms[n].atomic_number != 1:
            sh[m] = self._calc_hybridization(m)
            sn[m] -= 1
            # remove stereo marks on unbonded atoms and all its bonds
            try:
                del self._atoms_stereo[m]
            except KeyError:
                pass
            if sbs[m]:
                for x in sbs[m]:
                    del sbs[x][m]
                sbs[m] = {}
        if atoms[m].atomic_number != 1:
            sh[n] = self._calc_hybridization(n)
            sn[n] -= 1

            try:
                del self._atoms_stereo[n]
            except KeyError:
                pass
            if sbs[n]:
                for x in sbs[n]:
                    del sbs[x][n]
                sbs[n] = {}

    def remap(self, mapping, *, copy=False) -> 'MoleculeContainer':
        h = super().remap(mapping, copy=copy)
        mg = mapping.get
        sn = self._neighbors
        shg = self._hydrogens
        sbs = self._bonds_stereo

        if copy:
            hn = h._neighbors
            hh = h._hybridizations
            hhg = h._hydrogens
            hc = h._conformers
            has = h._atoms_stereo
            hbs = h._bonds_stereo
        else:
            hn = {}
            hh = {}
            hhg = {}
            hc = []
            has = {}
            hbs = {}

        for n, hyb in self._hybridizations.items():
            m = mg(n, n)
            hn[m] = sn[n]
            hh[m] = hyb
            hhg[m] = shg[n]
            hbs[m] = {mg(x, x): s for x, s in sbs[n].items()}

        hc.extend({mg(n, n): x for n, x in c.items()} for c in self._conformers)

        for n, stereo in self._atoms_stereo.items():
            has[mg(n, n)] = stereo

        if copy:
            return h

        self._neighbors = hn
        self._hybridizations = hh
        self._hydrogens = hhg
        self._conformers = hc
        self._atoms_stereo = has
        self._bonds_stereo = hbs
        return self

    def copy(self, *, meta=True) -> 'MoleculeContainer':
        copy = super().copy(meta=meta)
        copy._neighbors = self._neighbors.copy()
        copy._hybridizations = self._hybridizations.copy()
        copy._hydrogens = self._hydrogens.copy()
        copy._conformers = [c.copy() for c in self._conformers]
        copy._atoms_stereo = self._atoms_stereo.copy()
        copy._bonds_stereo = {n: s.copy() for n, s in self._bonds_stereo.items()}
        return copy

    def substructure(self, atoms, *, as_query: bool = False, **kwargs) -> Union['MoleculeContainer',
                                                                                'query.QueryContainer']:
        """
        create substructure containing atoms from atoms list

        :param atoms: list of atoms numbers of substructure
        :param meta: if True metadata will be copied to substructure
        :param as_query: return Query object based on graph substructure
        """
        sub, atoms = super().substructure(atoms, query.QueryContainer if as_query else self.__class__, **kwargs)
        sa = self._atoms
        sb = self._bonds

        lost = {n for n, a in sa.items() if a.atomic_number != 1} - set(atoms)  # atoms not in substructure
        not_skin = {n for n in atoms if lost.isdisjoint(sb[n])}
        sub._atoms_stereo = {n: s for n, s in self._atoms_stereo.items() if n in not_skin}
        sub._bonds_stereo = cbs = {n: {} for n in atoms}
        for n, m_stereo in self._bonds_stereo.items():
            if n in not_skin:
                ns = cbs[n]
                for m, s in m_stereo.items():
                    if m in not_skin:
                        ns[m] = s

        if as_query:
            sub._atoms = ca = {}
            for n in atoms:
                atom = sa[n]
                ca[n] = atom = QueryElement.from_atomic_number(atom.atomic_number)(atom.isotope)
                atom._attach_to_graph(sub, n)

            sn = self._neighbors
            sh = self._hybridizations
            sub._neighbors = {n: (sn[n],) for n in atoms}
            sub._hybridizations = {n: (sh[n],) for n in atoms}
        else:
            sub._conformers = []
            sub._atoms = ca = {}
            for n in atoms:
                ca[n] = atom = sa[n].copy()
                atom._attach_to_graph(sub, n)

            # recalculate query marks
            sub._neighbors = sn = {}
            sub._hybridizations = sh = {}
            sub._hydrogens = shg = {}
            atoms = sub._atoms
            for n, m_bonds in sub._bonds.items():
                sn[n] = sum(atoms[m].atomic_number != 1 for m in m_bonds)
                sh[n] = sub._calc_hybridization(n)
                shg[n] = sub._calc_implicit(n)
        return sub

    def union(self, other):
        if isinstance(other, MoleculeContainer):
            u = super().union(other)
            u._conformers.clear()

            u._neighbors.update(other._neighbors)
            u._hybridizations.update(other._hybridizations)
            u._hydrogens.update(other._hydrogens)
            u._atoms_stereo.update(other._atoms_stereo)

            ub = u._bonds
            for n in other._bonds:
                ub[n] = {}
            seen = set()
            for n, m_bond in other._bonds.items():
                seen.add(n)
                for m, bond in m_bond.items():
                    if m not in seen:
                        ub[n][m] = ub[m][n] = bond.copy()

            us = u._bonds_stereo
            for n, m_stereo in other._bonds_stereo.items():
                us[n] = m_stereo.copy()

            ua = u._atoms
            for n, atom in other._atoms.items():
                ua[n] = atom = atom.copy()
                atom._attach_to_graph(u, n)
            return u
        elif isinstance(other, Graph):
            return other.union(self)
        else:
            raise TypeError('Graph expected')

    def compose(self, other):
        """
        compose 2 graphs to CGR

        :param other: Molecule or CGR Container
        :return: CGRContainer
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
            common = self._atoms.keys() & other
            h = cgr.CGRContainer()
            atoms = h._atoms

            for n in self._atoms.keys() - common:  # cleavage atoms
                h.add_atom(sa[n], n, charge=sc[n], is_radical=sr[n], xy=sp[n], p_charge=sc[n], p_is_radical=sr[n])
                for m, bond in sb[n].items():
                    if m not in atoms:
                        if m in common:  # bond to common atoms is broken bond
                            order = bond.order
                            bond = object.__new__(cgr.DynamicBond)
                            bond._DynamicBond__order, bond._DynamicBond__p_order = order, None
                        bonds.append((n, m, bond))
            for n in other._atoms.keys() - common:  # coupling atoms
                h.add_atom(oa[n], n, charge=oc[n], is_radical=or_[n], xy=op[n], p_charge=oc[n], p_is_radical=or_[n])
                for m, bond in ob[n].items():
                    if m not in atoms:
                        if m in common:  # bond to common atoms is formed bond
                            order = bond.order
                            bond = object.__new__(cgr.DynamicBond)
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
                        bond = object.__new__(cgr.DynamicBond)
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
            common = self._atoms.keys() & other
            h = other.__class__()  # subclasses support
            atoms = h._atoms

            for n in self._atoms.keys() - common:  # cleavage atoms
                h.add_atom(sa[n], n, charge=sc[n], is_radical=sr[n], xy=sp[n], p_charge=sc[n], p_is_radical=sr[n])
                for m, bond in sb[n].items():
                    if m not in atoms:
                        if m in common:  # bond to common atoms is broken bond
                            order = bond.order
                            bond = object.__new__(cgr.DynamicBond)
                            bond._DynamicBond__order, bond._DynamicBond__p_order = order, None
                        bonds.append((n, m, bond))
            for n in other._atoms.keys() - common:  # coupling atoms
                h.add_atom(oa[n], n, charge=oc[n], is_radical=or_[n], xy=op[n], p_charge=opc[n], p_is_radical=opr[n])
                for m, bond in ob[n].items():
                    if m not in atoms:
                        if m in common:  # bond to common atoms is formed bond
                            order = bond.p_order
                            if order:  # skip broken bond. X>None => None>None
                                bond = object.__new__(cgr.DynamicBond)
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
                        bond = object.__new__(cgr.DynamicBond)
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

    def get_mapping(self, other):
        if isinstance(other, MoleculeContainer):
            return super().get_mapping(other)
        raise TypeError('MoleculeContainer expected')

    def implicify_hydrogens(self):
        """
        remove explicit hydrogen if possible

        :return: number of removed hydrogens
        """
        atoms = self._atoms
        charges = self._charges
        radicals = self._radicals
        bonds = self._bonds
        explicit = defaultdict(list)
        for n, atom in atoms.items():
            if atom.atomic_number == 1:
                for m in bonds[n]:
                    if atoms[m].atomic_number != 1:
                        explicit[m].append(n)

        to_remove = set()
        for n, hs in explicit.items():
            atom = atoms[n]
            charge = charges[n]
            is_radical = radicals[n]
            len_h = len(hs)
            for i in range(len_h, 0, -1):
                hi = hs[:i]
                explicit_sum = 0
                explicit_dict = defaultdict(int)
                for m, bond in bonds[n].items():
                    if m not in hi:
                        explicit_sum += bond.order
                        explicit_dict[(bond.order, atoms[m].__class__)] += 1

                if any(s.issubset(explicit_dict) and all(explicit_dict[k] >= c for k, c in d.items()) and h >= i
                       for s, d, h in atom.valence_rules(charge, is_radical, explicit_sum)):
                    to_remove.update(hi)
                    break
        for n in to_remove:
            self.delete_atom(n)
        return len(to_remove)

    def explicify_hydrogens(self):
        """
        add explicit hydrogens to atoms

        :return: number of added atoms
        """
        to_add = []
        for n, h in self._hydrogens.items():
            try:
                to_add.extend([n] * h)
            except TypeError:
                raise ValenceError(f'atom {{{n}}} has valence error')
        for n in to_add:
            self.add_bond(n, self.add_atom('H'), 1)
        return len(to_add)

    def check_valence(self):
        """
        check valences of all atoms

        works only on molecules with aromatic rings in Kekule form
        :return: list of invalid atoms
        """
        atoms = self._atoms
        charges = self._charges
        radicals = self._radicals
        bonds = self._bonds
        errors = set(atoms)
        for n, atom in atoms.items():
            charge = charges[n]
            is_radical = radicals[n]
            explicit_sum = 0
            explicit_dict = defaultdict(int)
            for m, bond in bonds[n].items():
                order = bond.order
                if order == 4:  # aromatic rings not supported
                    break
                elif order != 8:  # any bond used for complexes
                    explicit_sum += order
                    explicit_dict[(order, atoms[m].__class__)] += 1
            else:
                try:
                    rules = atom.valence_rules(charge, is_radical, explicit_sum)
                except ValenceError:
                    pass
                else:
                    for s, d, h in rules:
                        if s.issubset(explicit_dict) and all(explicit_dict[k] >= c for k, c in d.items()):
                            errors.discard(n)
                            break
        return list(errors)

    @cached_property
    def aromatic_rings(self) -> List[List[int]]:
        """
        aromatic rings atoms numbers
        """
        bonds = self._bonds
        return [ring for ring in self.sssr if bonds[ring[0]][ring[-1]].order == 4
                and all(bonds[n][m].order == 4 for n, m in zip(ring, ring[1:]))]

    @cached_property
    def cumulenes(self) -> List[List[int]]:
        """
        alkenes, allenes and cumulenes atoms numbers
        """
        atoms = self._atoms
        bonds = self._bonds
        adj = defaultdict(set)  # carbon double bonds adjacency matrix
        for n, atom in atoms.items():
            if atom.atomic_number == 6:
                adj_n = adj[n].add
                b_sum = 0
                for m, bond in bonds[n].items():
                    b_sum += bond.order
                    if bond.order == 2 and atoms[m].atomic_number == 6:
                        adj_n(m)
                if b_sum > 4:
                    raise ValenceError(f'carbon atom: {n} has invalid valence = {b_sum}')
        if not adj:
            return []

        terminals = {x for x, y in adj.items() if len(y) == 1}
        cumulenes = []
        while terminals:
            m = terminals.pop()
            path = [m]
            cumulenes.append(path)
            while m not in terminals:
                n, m = m, adj[m].pop()
                adj[m].discard(n)
                path.append(m)
            terminals.discard(m)
        return cumulenes

    @cached_property
    def tetrahedrons(self) -> List[int]:
        """
        carbon sp3 atoms numbers
        """
        atoms = self._atoms
        bonds = self._bonds
        tetra = []
        for n, atom in atoms.items():
            if atom.atomic_number == 6:
                env = bonds[n]
                b_sum = sum(x.order for x in env.values())
                if b_sum > 4:
                    raise ValenceError(f'carbon atom: {n} has invalid valence = {b_sum}')
                elif all(x.order == 1 for x in env.values()):
                    tetra.append(n)
        return tetra

    @cached_args_method
    def _explicit_hydrogens(self, n) -> int:
        """
        number of explicit hydrogen atoms connected to atom.

        take into account any type of bonds with hydrogen atoms.
        """
        atoms = self._atoms
        return sum(atoms[m].atomic_number == 1 for m in self._bonds[n])

    @cached_args_method
    def _total_hydrogens(self, n) -> int:
        return self._hydrogens[n] + self._explicit_hydrogens(n)

    def _calc_implicit(self, n):
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
                    return
                elif order != 8:  # any bond used for complexes
                    explicit_sum += order
                    explicit_dict[(order, atoms[m].__class__)] += 1
            try:
                rules = atom.valence_rules(charge, is_radical, explicit_sum)
            except ValenceError:
                return
            for s, d, h in rules:
                if h and s.issubset(explicit_dict) and all(explicit_dict[k] >= c for k, c in d.items()):
                    return h
        return 0

    def _calc_hybridization(self, n):
        atoms = self._atoms
        hybridization = 1
        for m, bond in self._bonds[n].items():
            if atoms[m].atomic_number == 1:  # ignore hydrogen
                continue
            order = bond.order
            if order == 4:
                return 4
            elif order == 3:
                if hybridization != 3:
                    hybridization = 3
            elif order == 2:
                if hybridization == 2:
                    hybridization = 3
                elif hybridization == 1:
                    hybridization = 2
        return hybridization

    def __getstate__(self):
        return {'conformers': self._conformers, 'atoms_stereo': self._atoms_stereo, 'bonds_stereo': self._bonds_stereo,
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
            state['bonds_stereo'] = {n: {} for n in state['_node']}
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
            state['bonds_stereo'] = {n: {} for n in state['node']}

        super().__setstate__(state)
        self._conformers = state['conformers']
        self._atoms_stereo = state['atoms_stereo']
        self._bonds_stereo = state['bonds_stereo']

        # restore query and hydrogen marks
        self._neighbors = sn = {}
        self._hybridizations = sh = {}
        self._hydrogens = shg = {}
        atoms = state['atoms']
        for n, m_bonds in state['bonds'].items():
            sn[n] = sum(atoms[m].atomic_number != 1 for m in m_bonds)
            sh[n] = self._calc_hybridization(n)
            shg[n] = self._calc_implicit(n)


__all__ = ['MoleculeContainer']
