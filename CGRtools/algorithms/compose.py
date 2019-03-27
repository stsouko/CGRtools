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
from collections import defaultdict
from ..attributes import DynAtom, DynBond


class Compose:
    def __xor__(self, other):
        """
        G ^ H is CGR generation
        """
        return self.compose(other)

    def compose(self, other):
        """
        compose 2 graphs to CGR

        :param other: Molecule or CGR Container
        :return: CGRContainer
        """
        if not isinstance(other, Compose):
            raise TypeError('CGRContainer or MoleculeContainer [sub]class expected')

        cgr = self._get_subclass('CGRContainer')
        common = self._node.keys() & other

        if not common:
            if not (isinstance(self, cgr) or isinstance(other, cgr)):
                return cgr() | self | other
            return self | other

        unique_reactant = self._node.keys() - common
        unique_product = other._node.keys() - common

        h = cgr()
        atoms = h._node
        bonds = []
        common_adj = {n: {} for n in common}
        common_bonds = []

        r_atoms = {}
        r_skin = defaultdict(list)
        if isinstance(self, cgr):
            for n in unique_reactant:
                h.add_atom(self._node[n], n)
                for m, bond in self._adj[n].items():
                    if m not in atoms:
                        if m in common:  # bond to common atoms is broken bond
                            r_bond = bond._reactant
                            if r_bond is None:  # skip None>None
                                continue
                            r_skin[n].append(m)
                            bond = DynBond.__new__(DynBond)
                            bond.__init_copy__(r_bond, None)
                        bonds.append((n, m, bond))
            for n in common:
                r_atoms[n] = self._node[n]._reactant
                for m, bond in self._adj[n].items():
                    if m not in r_atoms and m in common:
                        tmp = [bond._reactant, None]
                        common_adj[n][m] = common_adj[m][n] = tmp
                        common_bonds.append((n, m, tmp))
        else:
            for n in unique_reactant:
                atom = DynAtom.__new__(DynAtom)  # add unique atom into CGR
                atom.__init_copy__(self._node[n], self._node[n])
                h.add_atom(atom, n)
                for m, r_bond in self._adj[n].items():  # unique atom neighbors
                    if m not in atoms:  # bond not analyzed yet
                        bond = DynBond.__new__(DynBond)
                        if m in common:  # bond to common atoms
                            r_skin[n].append(m)
                            bond.__init_copy__(r_bond, None)
                        else:  # bond static
                            bond.__init_copy__(r_bond, r_bond)
                        bonds.append((n, m, bond))
            for n in common:
                r_atoms[n] = self._node[n]
                for m, bond in self._adj[n].items():
                    if m not in r_atoms and m in common:  # analyze only common atoms bonds
                        tmp = [bond, None]  # reactant state only
                        common_adj[n][m] = common_adj[m][n] = tmp
                        common_bonds.append((n, m, tmp))

        p_atoms = {}
        p_skin = defaultdict(list)
        if isinstance(other, cgr):
            for n in unique_product:
                h.add_atom(other._node[n], n)
                for m, bond in other._adj[n].items():
                    if m not in atoms:
                        if m in common:  # bond to common atoms is new bond
                            p_bond = bond._product
                            if p_bond is None:  # skip None>None
                                continue
                            p_skin[n].append(m)
                            bond = DynBond.__new__(DynBond)
                            bond.__init_copy__(None, p_bond)
                        bonds.append((n, m, bond))
            for n in common:
                p_atoms[n] = other._node[n]._product
                n_bonds = common_adj[n]
                for m, bond in other._adj[n].items():
                    if m in n_bonds:
                        n_bonds[m][1] = bond._product
                    elif m not in p_atoms and m in common:  # new bond of reaction
                        p_bond = bond._product
                        if p_bond is None:  # skip None>None
                            continue
                        bond = DynBond.__new__(DynBond)
                        bond.__init_copy__(None, p_bond)
                        bonds.append((n, m, bond))
        else:
            for n in unique_product:
                atom = DynAtom.__new__(DynAtom)
                atom.__init_copy__(other._node[n], other._node[n])
                h.add_atom(atom, n)
                for m, p_bond in other._adj[n].items():
                    if m not in atoms:
                        bond = DynBond.__new__(DynBond)
                        if m in common:
                            p_skin[n].append(m)
                            bond.__init_copy__(None, p_bond)
                        else:
                            bond.__init_copy__(p_bond, p_bond)
                        bonds.append((n, m, bond))
            for n in common:
                p_atoms[n] = other._node[n]
                n_bonds = common_adj[n]
                for m, p_bond in other._adj[n].items():
                    if m in n_bonds:  # set product state of changed bond
                        n_bonds[m][1] = p_bond
                    elif m not in p_atoms and m in common:  # new bond of reaction
                        bond = DynBond.__new__(DynBond)
                        bond.__init_copy__(None, p_bond)
                        bonds.append((n, m, bond))

        for n, r_atom in r_atoms.items():  # prepare common DynAtom's
            p_atom = p_atoms[n]
            if r_atom.element != p_atom.element or r_atom.isotope != p_atom.isotope:
                raise ValueError('atom-to-atom mapping invalid')
            atom = DynAtom.__new__(DynAtom)
            atom.__init_copy__(r_atom, p_atom)
            h.add_atom(atom, n)

        for n, m, (r_bond, p_bond) in common_bonds:
            if r_bond is p_bond is None:  # skip None>None
                continue
            bond = DynBond.__new__(DynBond)
            bond.__init_copy__(r_bond, p_bond)
            h.add_bond(n, m, bond)

        for n, m, bond in bonds:
            h.add_bond(n, m, bond)

        return h


class CGRCompose(Compose):
    def __invert__(self):
        """
        decompose CGR
        """
        return self.decompose()

    def decompose(self):
        """
        decompose CGR to pair of Molecules, which represents reactants and products state of reaction

        :return: tuple of two molecules
        """
        mc = self._get_subclass('MoleculeContainer')
        reactants = mc()
        products = mc()

        for n, atom in self.atoms():
            reactants.add_atom(atom._reactant, n)
            products.add_atom(atom._product, n)

        for n, m, bond in self.bonds():
            if bond._reactant is not None:
                reactants.add_bond(n, m, bond._reactant)
            if bond._product is not None:
                products.add_bond(n, m, bond._product)
        return reactants, products


__all__ = ['Compose', 'CGRCompose']
