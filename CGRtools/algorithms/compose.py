# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
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


class CGRCompose:
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
        if not isinstance(other, CGRCompose):
            raise TypeError('CGRContainer or MoleculeContainer [sub]class expected')

        common = self._node.keys() & other
        if not common:
            return self.union(other)

        unique_reagent = self._node.keys() - common
        unique_product = other._node.keys() - common

        cgr = next(x for x in CGRCompose.__subclasses__() if x.__name__ == 'CGRContainer')

        h = cgr()
        atoms = h._node
        bonds = []

        r_atoms = {}
        r_bonds = []
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
                        r_bonds.append((n, m, bond._reagent))
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
                        r_bonds.append((n, m, bond))

        p_atoms = {}
        p_bonds = defaultdict(dict)
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
                for m, bond in other._adj[n].items():
                    if m not in p_atoms and m in common:
                        p_bonds[n][m] = p_bonds[m][n] = bond._product
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
                for m, bond in other._adj[n].items():
                    if m not in p_atoms and m in common:
                        p_bonds[n][m] = p_bonds[m][n] = bond

        for n, r_atom in r_atoms.items():
            p_atom = p_atoms[n]
            if r_atom.element != p_atom.element or r_atom.isotope != p_atom.isotope:
                raise ValueError('atom-to-atom mapping invalid')
            atom = DynAtom.__new__(DynAtom)
            atom.__init_copy__(r_atom, p_atom)
            h.add_atom(atom, n)
        for n, m, r_bond in r_bonds:
            p_bond = p_bonds[n][m]
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


__all__ = ['CGRCompose']
