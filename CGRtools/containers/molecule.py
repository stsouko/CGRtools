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
from typing import List
from .common import Graph
from ..exceptions import ValenceError
from ..periodictable import Element


class Bond:
    __slots__ = ('order',)

    def __init__(self, order):
        self.order = order


class MoleculeContainer(Graph):
    __slots__ = ('_conformers', '_neighbors', '_hybridization', '_atoms_stereo', '_bonds_stereo')

    def __init__(self):
        self._conformers = []
        self._neighbors = {}
        self._hybridization = {}
        self._atoms_stereo = {}
        self._bonds_stereo = {}
        super().__init__()

    def add_atom(self, atom, _map=None):
        if not isinstance(atom, Element):
            if isinstance(atom, str):
                atom = Element.from_symbol(atom)()
            elif isinstance(atom, int):
                atom = Element.from_atomic_number(atom)()
            else:
                raise TypeError('Element object expected')

        _map = super().add_atom(atom, _map)
        self._neighbors[_map] = 0
        self._hybridization[_map] = 1
        self._bonds_stereo[_map] = {}
        self._conformers.clear()  # clean conformers. need full recalculation for new system
        return _map

    def add_bond(self, n, m, bond):
        if isinstance(bond, int):
            order = bond
            bond = Bond(order)
        elif isinstance(bond, Bond):
            order = bond.order
        else:
            raise TypeError('Bond object required')

        super().add_bond(n, m, bond)
        self._conformers.clear()  # clean conformers. need full recalculation for new system

        sbs = self._bonds_stereo
        sh = self._hybridization

        # calc query marks dynamically.
        if self._atoms[n].atomic_number != 1:  # not hydrogen
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
        if self._atoms[m].atomic_number != 1:  # not hydrogen
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
                # SO3 and other 3 double-bonded has hyb = 3

        # remove stereo marks on bonded atoms and all its bonds
        for x in sbs[n]:
            if x != m:
                del sbs[x][n]
        for x in sbs[m]:
            if x != n:
                del sbs[x][m]
        if sbs[n]:
            sbs[n] = {}
        if sbs[m]:
            sbs[m] = {}

        try:
            del self._atoms_stereo[n]
        except KeyError:
            pass
        try:
            del self._atoms_stereo[m]
        except KeyError:
            pass

    def delete_atom(self, n):
        bonds = self._bonds[n]  # save bonds
        atoms = self._atoms
        isnt_hydrogen = atoms[n].atomic_number != 1
        super().delete_atom(n)
        self._conformers.clear()  # clean conformers. need full recalculation for new system

        sn = self._neighbors
        sh = self._hybridization
        sas = self._atoms_stereo
        sbs = self._bonds_stereo
        sb = self._bonds

        del sn[n]
        del sh[n]

        if isnt_hydrogen:  # neighbors query marks fix. ignore removed hydrogen
            for m in bonds:
                hybridization = 1
                for x, bond in sb[m].items():
                    if atoms[x].atomic_number == 1:  # ignore hydrogen
                        continue
                    if hybridization in (3, 4):
                        continue
                    order = bond.order
                    if order == 4:
                        hybridization = 4
                    elif order == 3:
                        hybridization = 3
                    elif order == 2:
                        if hybridization == 2:
                            hybridization = 3
                        else:
                            hybridization = 2
                sh[m] = hybridization
                sn[m] -= 1

        # remove stereo marks on deleted atoms and all its neighbors
        for m in sbs.pop(n):
            del sbs[m][n]
        try:
            del sas[n]
        except KeyError:
            pass
        for m in bonds:
            try:
                del sas[m]
            except KeyError:
                pass

    def delete_bond(self, n, m):
        super().delete_bond(n, m)
        self._conformers.clear()  # clean conformers. need full recalculation for new system

        atoms = self._atoms
        sbs = self._bonds_stereo
        sh = self._hybridization
        sn = self._neighbors

        # neighbors query marks fix. ignore removed hydrogen
        if atoms[n].atomic_number != 1:
            hybridization = 1
            for x, bond in self._bonds[n].items():
                if atoms[x].atomic_number == 1:  # ignore hydrogen
                    continue
                if hybridization in (3, 4):
                    continue
                order = bond.order
                if order == 4:
                    hybridization = 4
                elif order == 3:
                    hybridization = 3
                elif order == 2:
                    if hybridization == 2:
                        hybridization = 3
                    else:
                        hybridization = 2
            sh[m] = hybridization
            sn[m] -= 1
        if atoms[m].atomic_number != 1:
            hybridization = 1
            for x, bond in self._bonds[m].items():
                if atoms[x].atomic_number == 1:
                    continue
                if hybridization in (3, 4):
                    continue
                order = bond.order
                if order == 4:
                    hybridization = 4
                elif order == 3:
                    hybridization = 3
                elif order == 2:
                    if hybridization == 2:
                        hybridization = 3
                    else:
                        hybridization = 2
            sh[n] = hybridization
            sn[n] -= 1

        # remove stereo marks on unbonded atoms and all its bonds
        for x in sbs[n]:
            if x != m:
                del sbs[x][n]
        for x in sbs[m]:
            if x != n:
                del sbs[x][m]
        if sbs[n]:
            sbs[n] = {}
        if sbs[m]:
            sbs[m] = {}

        try:
            del self._atoms_stereo[n]
        except KeyError:
            pass
        try:
            del self._atoms_stereo[m]
        except KeyError:
            pass

    def remap(self, mapping, *, copy=False):
        h = super().remap(mapping, copy=copy)
        mg = mapping.get
        sn = self._neighbors
        sh = self._hybridization
        sbs = self._bonds_stereo

        if copy:
            hn = h._neighbors
            hh = h._hybridization
            hc = h._conformers
            has = h._atoms_stereo
            hbs = h._bonds_stereo
        else:
            hn = {}
            hh = {}
            hc = []
            has = {}
            hbs = {}

        for n, hyb in sh.items():
            m = mg(n, n)
            hn[m] = sn[n]
            hh[m] = hyb
            hbs[m] = {mg(x, x): s for x, s in sbs[n].items()}

        hc.extend({mg(n, n): x for n, x in c.items()} for c in self._conformers)

        for n, stereo in self._atoms_stereo.items():
            has[mg(n, n)] = stereo

        if copy:
            return h

        self._neighbors = hn
        self._hybridization = hh
        self._conformers = hc
        self._atoms_stereo = has
        self._bonds_stereo = hbs
        return self

    def implicify_hydrogens(self):
        """
        remove explicit hydrogen if possible

        :return: number of removed hydrogens
        """
        atoms = self._atoms
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
                       for s, d, h in atom.valence_rules(explicit_sum)):
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
        atoms = self._atoms
        bonds = self._bonds
        to_add = []
        for n, atom in atoms.items():
            if atom.atomic_number != 1:
                explicit_sum = 0
                explicit_dict = defaultdict(int)
                for m, bond in bonds[n].items():
                    explicit_sum += bond.order
                    explicit_dict[(bond.order, atoms[m].__class__)] += 1
                for s, d, h in atom.valence_rules(explicit_sum):
                    if h and s.issubset(explicit_dict) and all(explicit_dict[k] >= c for k, c in d.items()):
                        to_add.extend([n] * h)
                        break
        for n in to_add:
            self.add_bond(n, self.add_atom('H'), 1)
        return len(to_add)

    @cached_args_method
    def atom_implicit_h(self, n) -> int:
        atoms = self._atoms
        if atoms[n].atomic_number != 1:
            explicit_sum = 0
            explicit_dict = defaultdict(int)
            for m, bond in self._bonds[n].items():
                explicit_sum += bond.order
                explicit_dict[(bond.order, atoms[m].__class__)] += 1
            for s, d, h in atoms[n].valence_rules(explicit_sum):
                if h and s.issubset(explicit_dict) and all(explicit_dict[k] >= c for k, c in d.items()):
                    return h
        return 0

    @cached_args_method
    def atom_explicit_h(self, n) -> int:
        atoms = self._atoms
        return sum(atoms[m].atomic_number == 1 for m in self._bonds[n])

    @cached_args_method
    def atom_total_h(self, n) -> int:
        return self.atom_implicit_h(n) + self.atom_explicit_h(n)

    def check_valence(self):
        """
        check valences of all atoms

        :return: list of invalid atoms
        """
        atoms = self._atoms
        bonds = self._bonds
        errors = set(atoms)
        for n, atom in atoms.items():
            explicit_sum = 0
            explicit_dict = defaultdict(int)
            for m, bond in bonds[n].items():
                explicit_sum += bond.order
                explicit_dict[(bond.order, atoms[m].__class__)] += 1
            try:
                rules = atoms[n].valence_rules(explicit_sum)
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
        return [ring for ring in self.sssr if len(ring) in (5, 6, 7) and bonds[ring[0]][ring[-1]].order == 4
                and all(bonds[n][m].order == 4 for n, m in zip(ring, ring[1:]))]

    @cached_property
    def cumulenes(self) -> List[List[int]]:
        """
        alkenes, allenes and cumulenes atoms numbers
        """
        atoms = self._atoms
        adj = defaultdict(set)  # carbon double bonds adjacency matrix
        for n, m_bond in self._bonds.items():
            if atoms[n].element == 'C':
                adj_n = adj[n].add
                for m, bond in m_bond.items():
                    if bond.order == 2 and atoms[m].element == 'C':
                        adj_n(m)
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
        return [n for n, env in self._bonds.items()
                if atoms[n].element == 'C' and all(x.order == 1 for x in env.values())]

    def _matcher(self, other):
        """
        return VF2 GraphMatcher

        MoleculeContainer < MoleculeContainer
        MoleculeContainer < CGRContainer
        """
        if isinstance(other, (self._get_subclass('CGRContainer'), MoleculeContainer)):
            return GraphMatcher(other, self, lambda x, y: x == y, lambda x, y: x == y)
        raise TypeError('only cgr-cgr possible')

    def __getstate__(self):
        return {'conformers': self._conformers, 'atoms_stereo': self._atoms_stereo, 'bonds_stereo': self._bonds_stereo,
                **super().__getstate__()}

    def __setstate__(self, state):
        if '_BaseContainer__meta' in state:  # 2.8 reverse compatibility
            state['atoms'] = {n: Element.from_symbol(atom['element'])(atom['s_charge'], atom.get('isotope'),
                                                               bool(atom['s_radical']))
                              for n, atom in state['_node'].items()}
            state['plane'] = {n: (atom['s_x'], atom['s_y']) for n, atom in state['node'].items()}

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
            state['bonds_stereo'] = {}
            state['meta'] = state['_BaseContainer__meta']
        elif 'node' in state:  # 3.0-3.1 compatibility.
            if 'graph' in state:  # 3.0.10 compatibility.
                state['meta'] = state['graph']
            state['atoms'] = {n: Element.from_atomic_number(a.number)(a.charge,
                                                                      a.isotope if a.common_isotope != a.isotope else
                                                                      None,
                                                                      bool(a.multiplicity))
                              for n, a in state['node'].items()}

            state['bonds'] = bonds = {}
            for n, m_bond in state['adj'].items():
                bonds[n] = bn = {}
                for m, bond in m_bond.items():
                    if m in bonds:
                        bn[m] = bonds[m][n]
                    else:
                        bn[m] = Bond(bond.order)

            state['plane'] = {n: (a.x, a.y) for n, a in state['node'].items()}
            state['conformers'] = []
            state['atoms_stereo'] = {}
            state['bonds_stereo'] = {}

        self._conformers = state['conformers']
        self._atoms_stereo = state['atoms_stereo']
        self._bonds_stereo = state['bonds_stereo']

        # restore query marks
        atoms = state['atoms']
        self._neighbors = sn = {}
        self._hybridization = sh = {}
        for n, m_bonds in state['bonds'].items():
            neighbors = 0
            hybridization = 1
            # hybridization 1- sp3; 2- sp2; 3- sp1; 4- aromatic
            for m, bond in m_bonds.items():
                if atoms[m].atomic_number == 1:  # ignore hydrogen
                    continue
                neighbors += 1
                if hybridization in (3, 4):
                    continue
                order = bond.order
                if order == 4:
                    hybridization = 4
                elif order == 3:
                    hybridization = 3
                elif order == 2:
                    if hybridization == 2:
                        hybridization = 3
                    else:
                        hybridization = 2
            sn[n] = neighbors
            sh[n] = hybridization


__all__ = ['MoleculeContainer']
