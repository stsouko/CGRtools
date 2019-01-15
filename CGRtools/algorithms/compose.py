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
from ..attributes import DynAtom, DynBond, Bond


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

        unique_reagent = self._node.keys() - common
        unique_product = other._node.keys() - common

        h = cgr()
        atoms = h._node
        bonds = []
        common_adj = {n: {} for n in common}
        common_bonds = []

        r_atoms = {}
        r_skin = defaultdict(list)
        if isinstance(self, cgr):
            for n in unique_reagent:
                h.add_atom(self._node[n], n)
                for m, bond in self._adj[n].items():
                    if m not in atoms:
                        if m in common:
                            r_skin[n].append(m)
                            r_bond = bond._reagent
                            bond = DynBond.__new__(DynBond)
                            bond.__init_copy__(r_bond, self.__none_bond)
                        bonds.append((n, m, bond))
            for n in common:
                r_atoms[n] = self._node[n]._reagent
                for m, bond in self._adj[n].items():
                    if m not in r_atoms and m in common:
                        tmp = [bond._reagent, self.__none_bond]
                        common_adj[n][m] = common_adj[m][n] = tmp
                        common_bonds.append((n, m, tmp))
        else:
            for n in unique_reagent:
                atom = DynAtom.__new__(DynAtom)
                atom.__init_copy__(self._node[n], self._node[n])
                h.add_atom(atom, n)
                for m, r_bond in self._adj[n].items():
                    if m not in atoms:
                        bond = DynBond.__new__(DynBond)
                        if m in common:
                            r_skin[n].append(m)
                            bond.__init_copy__(r_bond, self.__none_bond)
                        else:
                            bond.__init_copy__(r_bond, r_bond)
                        bonds.append((n, m, bond))
            for n in common:
                r_atoms[n] = self._node[n]
                for m, bond in self._adj[n].items():
                    if m not in r_atoms and m in common:
                        tmp = [bond, self.__none_bond]
                        common_adj[n][m] = common_adj[m][n] = tmp
                        common_bonds.append((n, m, tmp))

        p_atoms = {}
        p_skin = defaultdict(list)
        if isinstance(other, cgr):
            for n in unique_product:
                h.add_atom(other._node[n], n)
                for m, bond in other._adj[n].items():
                    if m not in atoms:
                        if m in common:
                            p_skin[n].append(m)
                            p_bond = bond._product
                            bond = DynBond.__new__(DynBond)
                            bond.__init_copy__(self.__none_bond, p_bond)
                        bonds.append((n, m, bond))
            for n in common:
                p_atoms[n] = other._node[n]._product
                n_bonds = common_adj[n]
                for m, bond in other._adj[n].items():
                    if m in n_bonds:
                        n_bonds[m][1] = bond._product
                    elif m not in p_atoms and m in common:
                        tmp = [self.__none_bond, bond._product]
                        n_bonds[m] = common_adj[m][n] = tmp
                        common_bonds.append((n, m, tmp))
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
                            bond.__init_copy__(self.__none_bond, p_bond)
                        else:
                            bond.__init_copy__(p_bond, p_bond)
                        bonds.append((n, m, bond))
            for n in common:
                p_atoms[n] = other._node[n]
                n_bonds = common_adj[n]
                for m, bond in other._adj[n].items():
                    if m in n_bonds:
                        n_bonds[m][1] = bond
                    elif m not in p_atoms and m in common:
                        tmp = [self.__none_bond, bond]
                        n_bonds[m] = common_adj[m][n] = tmp
                        common_bonds.append((n, m, tmp))

        for n, r_atom in r_atoms.items():
            p_atom = p_atoms[n]
            if r_atom.element != p_atom.element or r_atom.isotope != p_atom.isotope:
                raise ValueError('atom-to-atom mapping invalid')
            atom = DynAtom.__new__(DynAtom)
            atom.__init_copy__(r_atom, p_atom)
            h.add_atom(atom, n)

        for n, m, (r_bond, p_bond) in common_bonds:
            if r_bond.order is p_bond.order is None:
                continue
            bond = DynBond.__new__(DynBond)
            bond.__init_copy__(r_bond, p_bond)
            h.add_bond(n, m, bond)

        for n, m, bond in bonds:
            h.add_bond(n, m, bond)

        return h

    __none_bond = Bond(skip_checks=True)
    __none_bond.order = None


class CGRCompose(Compose):
    def __invert__(self):
        """
        decompose CGR
        """
        return self.decompose()

    def decompose(self):
        """
        decompose CGR to pair of Molecules, which represents reagents and products state of reaction

        :return: tuple of two molecules
        """
        mc = self._get_subclass('MoleculeContainer')
        reagents = mc()
        products = mc()

        for n, atom in self.atoms():
            reagents.add_atom(atom._reagent, n)
            products.add_atom(atom._product, n)

        for n, m, bond in self.bonds():
            if bond.order:
                reagents.add_bond(n, m, bond._reagent)
            if bond.p_order:
                products.add_bond(n, m, bond._product)
        return reagents, products


__all__ = ['Compose', 'CGRCompose']
