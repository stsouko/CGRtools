# -*- coding: utf-8 -*-
#
#  Copyright 2017-2019 Ramil Nugmanov <stsouko@live.ru>
#  Copyright 2018 Ravil Mukhametgaleev <sonic-mc@mail.ru>
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
from CachedMethods import cached_property
from collections import defaultdict
from typing import List, Union, Tuple, Dict, Optional
from . import cgr_query as query, molecule  # cyclic imports resolve
from .bonds import Bond, DynamicBond
from .common import Graph
from ..algorithms.depict import DepictCGR
from ..algorithms.smiles import CGRSmiles
from ..exceptions import MappingError
from ..periodictable import DynamicElement, Element, DynamicQueryElement


class CGRContainer(Graph, CGRSmiles, DepictCGR):
    __slots__ = ('_p_charges', '_p_radicals', '_neighbors', '_hybridizations', '_p_neighbors', '_p_hybridizations')

    def __init__(self):
        self._p_charges: Dict[int, int] = {}
        self._p_radicals: Dict[int, bool] = {}
        self._neighbors: Dict[int, int] = {}
        self._hybridizations: Dict[int, int] = {}
        self._p_neighbors: Dict[int, int] = {}
        self._p_hybridizations: Dict[int, int] = {}
        super().__init__()

    def add_atom(self, atom: Union[DynamicElement, Element, int, str], *args, p_charge: int = 0,
                 p_is_radical: bool = False, **kwargs):
        p_charge = self._validate_charge(p_charge)
        p_is_radical = self._validate_radical(p_is_radical)

        if not isinstance(atom, DynamicElement):
            if isinstance(atom, Element):
                atom = DynamicElement.from_atomic_number(atom.atomic_number)(atom.isotope)
            elif isinstance(atom, str):
                atom = DynamicElement.from_symbol(atom)()
            elif isinstance(atom, int):
                atom = DynamicElement.from_atomic_number(atom)()
            else:
                raise TypeError('DynamicElement object expected')

        _map = super().add_atom(atom, *args, **kwargs)
        self._p_charges[_map] = p_charge
        self._p_radicals[_map] = p_is_radical
        self._neighbors[_map] = 0
        self._hybridizations[_map] = 1
        self._p_neighbors[_map] = 0
        self._p_hybridizations[_map] = 1
        return _map

    def add_bond(self, n, m, bond: Union[DynamicBond, Bond, int]):
        if isinstance(bond, DynamicBond):
            order = bond.order
            p_order = bond.p_order
        elif isinstance(bond, Bond):
            order = p_order = bond.order
            bond = object.__new__(DynamicBond)
            bond._DynamicBond__order = bond._DynamicBond__p_order = order
        else:
            order = p_order = bond
            bond = DynamicBond(order, order)

        super().add_bond(n, m, bond)

        sh = self._hybridizations
        sph = self._p_hybridizations

        # calc query marks dynamically.
        if self._atoms[n].atomic_number != 1:  # not hydrogen
            if order:
                self._neighbors[m] += 1
                if order == 4:
                    if sh[m] != 4:
                        sh[m] = 4
                elif order == 3:
                    if sh[m] not in (3, 4):
                        sh[m] = 3
                elif order == 2:
                    if sh[m] == 1:
                        sh[m] = 2
                    elif sh[m] == 2:
                        sh[m] = 3
            if p_order:
                self._p_neighbors[m] += 1
                if p_order == 4:
                    if sph[m] != 4:
                        sph[m] = 4
                elif p_order == 3:
                    if sph[m] not in (3, 4):
                        sph[m] = 3
                elif p_order == 2:
                    if sph[m] == 1:
                        sph[m] = 2
                    elif sph[m] == 2:
                        sph[m] = 3
        if self._atoms[m].atomic_number != 1:  # not hydrogen
            if order:
                self._neighbors[n] += 1
                if order == 4:
                    if sh[n] != 4:
                        sh[n] = 4
                elif order == 3:
                    if sh[n] not in (3, 4):
                        sh[n] = 3
                elif order == 2:
                    if sh[n] == 1:
                        sh[n] = 2
                    elif sh[n] == 2:
                        sh[n] = 3
            if p_order:
                self._p_neighbors[n] += 1
                if p_order == 4:
                    if sph[n] != 4:
                        sph[n] = 4
                elif p_order == 3:
                    if sph[n] not in (3, 4):
                        sph[n] = 3
                elif p_order == 2:
                    if sph[n] == 1:
                        sph[n] = 2
                    elif sph[n] == 2:
                        sph[n] = 3

    def delete_atom(self, n):
        del self._p_charges[n]
        del self._p_radicals[n]

        atoms = self._atoms
        old_bonds = self._bonds[n]  # save bonds
        sn = self._neighbors
        sh = self._hybridizations
        spn = self._p_neighbors
        sph = self._p_hybridizations
        sb = self._bonds
        isnt_hydrogen = atoms[n].atomic_number != 1
        super().delete_atom(n)

        del sn[n]
        del sh[n]
        del spn[n]
        del sph[n]

        if isnt_hydrogen:  # neighbors query marks fix. ignore removed hydrogen
            for m, old_bond in old_bonds.items():
                if old_bond.order:
                    sn[m] -= 1
                if old_bond.p_order:
                    spn[m] -= 1

                hybridization = p_hybridization = 1
                for x, bond in sb[m].items():
                    if atoms[x].atomic_number == 1:  # ignore hydrogen
                        continue
                    order = bond.order
                    p_order = bond.p_order
                    if order:
                        if hybridization != 4:
                            if order == 4:
                                hybridization = 4
                            elif order == 3:
                                if hybridization != 3:
                                    hybridization = 3
                            elif order == 2:
                                if hybridization == 2:
                                    hybridization = 3
                                elif hybridization == 1:
                                    hybridization = 2
                    if p_order:
                        if p_hybridization != 4:
                            if p_order == 4:
                                p_hybridization = 4
                            elif p_order == 3:
                                if p_hybridization != 3:
                                    p_hybridization = 3
                            elif p_order == 2:
                                if p_hybridization == 2:
                                    p_hybridization = 3
                                elif p_hybridization == 1:
                                    p_hybridization = 2
                sh[m] = hybridization
                sph[m] = p_hybridization

    def delete_bond(self, n, m):
        old_bond = self._bonds[n][m]  # save bond
        super().delete_bond(n, m)

        atoms = self._atoms
        sh = self._hybridizations
        sn = self._neighbors
        sph = self._p_hybridizations
        spn = self._p_neighbors

        # neighbors query marks fix. ignore removed hydrogen
        if atoms[n].atomic_number != 1:
            if old_bond.order:
                sn[m] -= 1
            if old_bond.p_order:
                spn[m] -= 1

            hybridization = p_hybridization = 1
            for x, bond in self._bonds[m].items():
                if atoms[x].atomic_number == 1:  # ignore hydrogen
                    continue
                order = bond.order
                p_order = bond.p_order
                if order:
                    if hybridization != 4:
                        if order == 4:
                            hybridization = 4
                        elif order == 3:
                            if hybridization != 3:
                                hybridization = 3
                        elif order == 2:
                            if hybridization == 2:
                                hybridization = 3
                            elif hybridization == 1:
                                hybridization = 2
                if p_order:
                    if p_hybridization != 4:
                        if p_order == 4:
                            p_hybridization = 4
                        elif p_order == 3:
                            if p_hybridization != 3:
                                p_hybridization = 3
                        elif p_order == 2:
                            if p_hybridization == 2:
                                p_hybridization = 3
                            elif p_hybridization == 1:
                                p_hybridization = 2
            sh[m] = hybridization
            sph[m] = p_hybridization
        if atoms[m].atomic_number != 1:
            if old_bond.order:
                sn[n] -= 1
            if old_bond.p_order:
                spn[n] -= 1

            hybridization = p_hybridization = 1
            for x, bond in self._bonds[n].items():
                if atoms[x].atomic_number == 1:
                    continue
                order = bond.order
                p_order = bond.p_order
                if order:
                    if hybridization != 4:
                        if order == 4:
                            hybridization = 4
                        elif order == 3:
                            if hybridization != 3:
                                hybridization = 3
                        elif order == 2:
                            if hybridization == 2:
                                hybridization = 3
                            elif hybridization == 1:
                                hybridization = 2
                if p_order:
                    if p_hybridization != 4:
                        if p_order == 4:
                            p_hybridization = 4
                        elif p_order == 3:
                            if p_hybridization != 3:
                                p_hybridization = 3
                        elif p_order == 2:
                            if p_hybridization == 2:
                                p_hybridization = 3
                            elif p_hybridization == 1:
                                p_hybridization = 2
            sh[n] = hybridization
            sph[n] = p_hybridization

    def remap(self, mapping, *, copy=False) -> 'CGRContainer':
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

    def copy(self, **kwargs) -> 'CGRContainer':
        copy = super().copy(**kwargs)
        copy._neighbors = self._neighbors.copy()
        copy._hybridizations = self._hybridizations.copy()
        copy._p_neighbors = self._p_neighbors.copy()
        copy._p_hybridizations = self._p_hybridizations.copy()
        copy._p_radicals = self._p_radicals.copy()
        copy._p_charges = self._p_charges.copy()
        return copy

    def substructure(self, atoms, *, as_query: bool = False, **kwargs) -> Union['CGRContainer',
                                                                                'query.QueryCGRContainer']:
        """
        create substructure containing atoms from atoms list

        :param atoms: list of atoms numbers of substructure
        :param meta: if True metadata will be copied to substructure
        :param as_query: return Query object based on graph substructure
        """
        sub, atoms = super().substructure(atoms, query.QueryCGRContainer if as_query else self.__class__, **kwargs)
        sa = self._atoms
        spc = self._p_charges
        spr = self._p_radicals

        sub._p_charges = {n: spc[n] for n in atoms}
        sub._p_radicals = {n: spr[n] for n in atoms}

        if as_query:
            sub._atoms = ca = {}
            for n in atoms:
                atom = sa[n]
                ca[n] = atom = DynamicQueryElement.from_atomic_number(atom.atomic_number)(atom.isotope)
                atom._attach_to_graph(sub, n)

            sn = self._neighbors
            sh = self._hybridizations
            spn = self._p_neighbors
            sph = self._p_hybridizations
            sub._neighbors = {n: (sn[n],) for n in atoms}
            sub._hybridizations = {n: (sh[n],) for n in atoms}
            sub._p_neighbors = {n: (spn[n],) for n in atoms}
            sub._p_hybridizations = {n: (sph[n],) for n in atoms}
        else:
            sub._atoms = ca = {}
            for n in atoms:
                ca[n] = atom = sa[n].copy()
                atom._attach_to_graph(sub, n)

            # recalculate query marks
            sub._neighbors = sn = {}
            sub._hybridizations = sh = {}
            sub._p_neighbors = spn = {}
            sub._p_hybridizations = sph = {}
            atoms = sub._atoms
            for n, m_bonds in sub._bonds.items():
                neighbors = p_neighbors = 0
                hybridization = p_hybridization = 1
                for m, bond in m_bonds.items():
                    if atoms[m].atomic_number == 1:  # ignore hydrogen
                        continue
                    order = bond.order
                    p_order = bond.p_order
                    if order:
                        neighbors += 1
                        if hybridization != 4:
                            if order == 4:
                                hybridization = 4
                            elif order == 3:
                                if hybridization != 3:
                                    hybridization = 3
                            elif order == 2:
                                if hybridization == 2:
                                    hybridization = 3
                                elif hybridization == 1:
                                    hybridization = 2
                    if p_order:
                        p_neighbors += 1
                        if p_hybridization != 4:
                            if p_order == 4:
                                p_hybridization = 4
                            elif p_order == 3:
                                if p_hybridization != 3:
                                    p_hybridization = 3
                            elif p_order == 2:
                                if p_hybridization == 2:
                                    p_hybridization = 3
                                elif p_hybridization == 1:
                                    p_hybridization = 2
                sn[n] = neighbors
                sh[n] = hybridization
                spn[n] = p_neighbors
                sph[n] = p_hybridization
        return sub

    def union(self, other):
        if isinstance(other, CGRContainer):
            u = super().union(other)
            u._p_charges.update(other._p_charges)
            u._p_radicals.update(other._p_radicals)
            u._neighbors.update(other._neighbors)
            u._hybridizations.update(other._hybridizations)
            u._p_neighbors.update(other._p_neighbors)
            u._p_hybridizations.update(other._p_hybridizations)

            ub = u._bonds
            for n in other._bonds:
                ub[n] = {}
            seen = set()
            for n, m_bond in other._bonds.items():
                seen.add(n)
                for m, bond in m_bond.items():
                    if m not in seen:
                        ub[n][m] = ub[m][n] = bond.copy()

            ua = u._atoms
            for n, atom in other._atoms.items():
                ua[n] = atom = atom.copy()
                atom._attach_to_graph(u, n)
            return u
        elif isinstance(other, molecule.MoleculeContainer):
            u = super().union(other)
            u._p_charges.update(other._charges)
            u._p_radicals.update(other._radicals)
            u._neighbors.update(other._neighbors)
            u._hybridizations.update(other._hybridizations)
            u._p_neighbors.update(other._neighbors)
            u._p_hybridizations.update(other._hybridizations)

            ub = u._bonds
            for n, m_bond in other._bonds.items():
                ub[n] = {m: DynamicBond(b.order, b.order) for m, b in m_bond.items()}

            ua = u._atoms
            for n, atom in other._atoms.items():
                ua[n] = atom = DynamicElement.from_atomic_number(atom.atomic_number)(atom.isotope)
                atom._attach_to_graph(u, n)
            return u
        elif isinstance(other, Graph):  # Query or CGRQuery
            return other.union(self)
        else:
            raise TypeError('Graph expected')

    def compose(self, other: Union['molecule.MoleculeContainer', 'CGRContainer']) -> 'CGRContainer':
        """
        compose 2 graphs to CGR

        :param other: Molecule or CGR Container
        :return: CGRContainer
        """
        sa = self._atoms
        sc = self._charges
        sr = self._radicals
        spc = self._p_charges
        spr = self._p_radicals
        sp = self._plane
        sb = self._bonds

        bonds = []
        adj: Dict[int, Dict[int, List[Optional[int]]]] = defaultdict(lambda: defaultdict(lambda: [None, None]))
        h = self.__class__()  # subclasses support
        atoms = h._atoms

        if isinstance(other, molecule.MoleculeContainer):
            oa = other._atoms
            oc = other._charges
            or_ = other._radicals
            op = other._plane
            ob = other._bonds
            common = sa.keys() & other

            for n in sa.keys() - common:  # cleavage atoms
                h.add_atom(sa[n].copy(), n, charge=sc[n], is_radical=sr[n], xy=sp[n], p_charge=spc[n],
                           p_is_radical=spr[n])
                for m, bond in sb[n].items():
                    if m not in atoms:
                        if m in common:  # bond to common atoms is broken bond
                            order = bond.order
                            if order:  # skip formed bond. None>X => None>None
                                bond = object.__new__(DynamicBond)
                                bond._DynamicBond__order, bond._DynamicBond__p_order = order, None
                                bonds.append((n, m, bond))
                        else:
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
                for m, bond in ob[n].items():
                    if m in common:
                        an[m][1] = bond.order
                for m, bond in sb[n].items():
                    if m in an or m in common and bond.order:
                        an[m][0] = bond.order
            for n in common:
                san = sa[n]
                if san.atomic_number != oa[n].atomic_number or san.isotope != oa[n].isotope:
                    raise MappingError(f'atoms with number {{{n}}} not equal')
                h.add_atom(san.copy(), n, charge=sc[n], is_radical=sr[n], xy=sp[n], p_charge=oc[n], p_is_radical=or_[n])
                for m, (o1, o2) in adj[n].items():
                    if m not in atoms:
                        bond = object.__new__(DynamicBond)
                        bond._DynamicBond__order, bond._DynamicBond__p_order = o1, o2
                        bonds.append((n, m, bond))
        elif isinstance(other, CGRContainer):
            oa = other._atoms
            oc = other._charges
            or_ = other._radicals
            opc = other._p_charges
            opr = other._p_radicals
            op = other._plane
            ob = other._bonds
            common = sa.keys() & other

            for n in sa.keys() - common:  # cleavage atoms
                h.add_atom(sa[n].copy(), n, charge=sc[n], is_radical=sr[n], xy=sp[n], p_charge=spc[n],
                           p_is_radical=spr[n])
                for m, bond in sb[n].items():
                    if m not in atoms:
                        if m in common:  # bond to common atoms is broken bond
                            order = bond.order
                            if order:  # skip formed bond. None>X => None>None
                                bond = object.__new__(DynamicBond)
                                bond._DynamicBond__order, bond._DynamicBond__p_order = order, None
                                bonds.append((n, m, bond))
                        else:
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
                    if m in common and bond.order:  # skip formed bonds
                        an[m][0] = bond.order
                for m, bond in ob[n].items():
                    if m in an or m in common and bond.p_order:
                        # self has nm bond or other bond not broken
                        an[m][1] = bond.p_order
            for n in common:
                san = sa[n]
                if san.atomic_number != oa[n].atomic_number or san.isotope != oa[n].isotope:
                    raise MappingError(f'atoms with number {{{n}}} not equal')
                h.add_atom(san.copy(), n, charge=sc[n], is_radical=sr[n], xy=sp[n], p_charge=opc[n],
                           p_is_radical=opr[n])
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

    def get_mapping(self, other: 'CGRContainer', **kwargs):
        if isinstance(other, CGRContainer):
            return super().get_mapping(other, **kwargs)
        raise TypeError('CGRContainer expected')

    def get_mcs_mapping(self, other: 'CGRContainer', **kwargs):
        if isinstance(other, CGRContainer):
            return super().get_mcs_mapping(other, **kwargs)
        raise TypeError('CGRContainer expected')

    @cached_property
    def centers_list(self) -> Tuple[Tuple[int, ...], ...]:
        """ get a list of lists of atoms of reaction centers
        """
        radicals = self._radicals
        p_charges = self._p_charges
        p_radicals = self._p_radicals
        center = set()
        adj = defaultdict(set)

        for n, c in self._charges.items():
            if c != p_charges[n] or radicals[n] != p_radicals[n]:
                center.add(n)

        for n, m_bond in self._bonds.items():
            for m, bond in m_bond.items():
                if bond.order != bond.p_order:
                    adj[n].add(m)
        center.update(adj)

        out = []
        while center:
            n = center.pop()
            if n in adj:
                c = set(plain_bfs(adj, n))
                out.append(tuple(c))
                center.difference_update(c)
            else:
                out.append((n,))
        return tuple(out)

    @cached_property
    def center_atoms(self) -> Tuple[int, ...]:
        """ get list of atoms of reaction center (atoms with dynamic: bonds, charges, radicals).
        """
        radicals = self._radicals
        p_charges = self._p_charges
        p_radicals = self._p_radicals

        center = set()
        for n, c in self._charges.items():
            if c != p_charges[n] or radicals[n] != p_radicals[n]:
                center.add(n)

        for n, m_bond in self._bonds.items():
            if any(bond.order != bond.p_order for bond in m_bond.values()):
                center.add(n)

        return tuple(center)

    @cached_property
    def center_bonds(self) -> Tuple[Tuple[int, int], ...]:
        """ get list of bonds of reaction center (bonds with dynamic orders).
        """
        return tuple((n, m) for n, m, bond in self.bonds() if bond.order != bond.p_order)

    @cached_property
    def aromatic_rings(self) -> Tuple[Tuple[int, ...], ...]:
        """
        existed or formed aromatic rings atoms numbers
        """
        adj = self._bonds
        return tuple(ring for ring in self.sssr if
                     adj[ring[0]][ring[-1]].order == 4 and all(adj[n][m].order == 4 for n, m in zip(ring, ring[1:])) or
                     adj[ring[0]][ring[-1]].p_order == 4 and all(adj[n][m].p_order == 4 for n, m in zip(ring, ring[1:]))
                     )

    def decompose(self) -> Tuple['molecule.MoleculeContainer', 'molecule.MoleculeContainer']:
        """
        decompose CGR to pair of Molecules, which represents reactants and products state of reaction

        :return: tuple of two molecules
        """
        charges = self._charges
        p_charges = self._p_charges
        radicals = self._radicals
        p_radicals = self._p_radicals
        plane = self._plane

        reactants = molecule.MoleculeContainer()
        products = molecule.MoleculeContainer()

        for n, atom in self._atoms.items():
            atom = Element.from_atomic_number(atom.atomic_number)(atom.isotope)
            reactants.add_atom(atom, n, charge=charges[n], is_radical=radicals[n], xy=plane[n])
            products.add_atom(atom.copy(), n, charge=p_charges[n], is_radical=p_radicals[n], xy=plane[n])

        for n, m, bond in self.bonds():
            if bond.order:
                reactants.add_bond(n, m, bond.order)
            if bond.p_order:
                products.add_bond(n, m, bond.p_order)
        return reactants, products

    def __invert__(self):
        """
        decompose CGR
        """
        return self.decompose()

    def __getstate__(self):
        return {'p_charges': self._p_charges, 'p_radicals': self._p_radicals, **super().__getstate__()}

    def __setstate__(self, state):
        self._p_charges = state['p_charges']
        self._p_radicals = state['p_radicals']
        super().__setstate__(state)

        # restore query marks
        self._neighbors = sn = {}
        self._hybridizations = sh = {}
        self._p_neighbors = spn = {}
        self._p_hybridizations = sph = {}
        atoms = state['atoms']
        for n, m_bonds in state['bonds'].items():
            neighbors = p_neighbors = 0
            hybridization = p_hybridization = 1
            for m, bond in m_bonds.items():
                if atoms[m].atomic_number == 1:  # ignore hydrogen
                    continue
                order = bond.order
                p_order = bond.p_order
                if order:
                    neighbors += 1
                    if hybridization != 4:
                        if order == 4:
                            hybridization = 4
                        elif order == 3:
                            if hybridization != 3:
                                hybridization = 3
                        elif order == 2:
                            if hybridization == 2:
                                hybridization = 3
                            elif hybridization == 1:
                                hybridization = 2
                if p_order:
                    p_neighbors += 1
                    if p_hybridization != 4:
                        if p_order == 4:
                            p_hybridization = 4
                        elif p_order == 3:
                            if p_hybridization != 3:
                                p_hybridization = 3
                        elif p_order == 2:
                            if p_hybridization == 2:
                                p_hybridization = 3
                            elif p_hybridization == 1:
                                p_hybridization = 2
            sn[n] = neighbors
            sh[n] = hybridization
            spn[n] = p_neighbors
            sph[n] = p_hybridization


def plain_bfs(adj, source):
    """modified NX fast BFS node generator"""
    seen = set()
    nextlevel = {source}
    while nextlevel:
        thislevel = nextlevel
        nextlevel = set()
        for v in thislevel:
            if v not in seen:
                yield v
                seen.add(v)
                nextlevel.update(adj[v])


__all__ = ['CGRContainer']
