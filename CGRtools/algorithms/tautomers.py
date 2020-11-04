# -*- coding: utf-8 -*-
#
#  Copyright 2020 Nail Samikaev <samikaevn@yandex.ru>
#  Copyright 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import TYPE_CHECKING, List, Tuple, Iterator
from ..containers.bonds import Bond

if TYPE_CHECKING:
    from CGRtools import MoleculeContainer


class Tautomers:
    __slots__ = ()

    def tautomerize(self) -> bool:
        """
        Convert structure to canonic tautomer form. Return True if structure changed.
        """
        # todo: implement
        return False

    def enumerate_tautomers(self) -> Iterator['MoleculeContainer']:
        """
        Enumerate all possible tautomeric forms of molecule. Supported hydrogen migration through delocalized chain.

        O=C-C[H] <-> [H]O-C=C
        O=C-C=C-C[H] <-> [H]O-C=C-C=C
        O=C-C[H]-C=C <-> [H]O-C=C-C=C
        O=C-C[H]-C=C <X> O=C-C=C-C[H]  not directly possible
        """
        yield self

        # making aromatic molecule
        mol = self.copy()
        mol.thiele()

        # generation of all tautomers
        seen = {self, mol}
        queue = [mol]
        while queue:
            for mol in queue.pop(0)._create_mol():
                if mol not in seen:
                    yield mol
                    seen.add(mol)
                    queue.append(mol)

    def _create_mol(self):
        atoms = self._atoms
        charges = self._charges
        radicals = self._radicals
        plane = self._plane
        bonds = self._bonds
        hybridizations = self._hybridizations
        hydrogens = self._hydrogens

        for path in self.__classic_case():
            mol = self.__class__()
            m_atoms = mol._atoms
            m_charges = mol._charges
            m_radicals = mol._radicals
            m_plane = mol._plane
            m_bonds = mol._bonds
            m_hybridizations = mol._hybridizations
            m_hydrogens = mol._hydrogens

            adj = {}
            for n, a in atoms.items():
                a = a.copy()
                m_atoms[n] = a
                a._attach_to_graph(mol, n)
                m_charges[n] = charges[n]
                m_radicals[n] = radicals[n]
                m_plane[n] = plane[n]
                m_bonds[n] = {}
                adj[n] = set()

            seen = set()
            for n, m, bond in path:
                if bond is not None:
                    m_bonds[n][m] = m_bonds[m][n] = Bond(bond)
                adj[n].add(m)
                adj[m].add(n)
                seen.add(n)
                seen.add(m)

            for n, ms in bonds.items():
                adjn = adj[n]
                for m, bond in ms.items():
                    if m not in adjn:
                        m_bonds[n][m] = m_bonds[m][n] = bond.copy()

            for n, h in hybridizations.items():
                if n in seen:
                    mol._calc_hybridization(n)
                    mol._calc_implicit(n)
                else:
                    m_hybridizations[n] = h
                    m_hydrogens[n] = hydrogens[n]
            mol.thiele()
            yield mol

    def __classic_case(self):
        bonds = self._bonds
        hydrogens = self._hydrogens

        # for each atom in entries
        for entry, has_hydrogen in self.__entries():
            seen = [entry]
            new_bonds = []

            # stack is neighbors
            stack = [(entry, n, b.order - 1, 1, True) for n, b in  # ketone case
                     bonds[entry].items() if 2 <= b.order <= 3]
            stack.extend((entry, n, b.order + 1, 0, False) for n, b in  # enol case
                         bonds[entry].items() if b.order <= 2 and has_hydrogen)

            while stack:
                previous, current, bond, depth, flag = stack.pop()

                # branch processing
                if flag:
                    if len(seen) > depth:
                        seen = seen[:depth]
                        new_bonds = new_bonds[:depth - 1]
                else:
                    if len(seen) > depth + 1:
                        seen = seen[:depth + 1]
                        new_bonds = new_bonds[:depth]

                # adding new bonds
                new_bonds.append((previous, current, bond))
                seen.append(current)

                # adding neighbors
                depth += 1
                diff = -1 if depth % 2 else 1
                stack.extend((current, nbg, bond, depth, flag) for nbg, bond in
                             ((nbg, bond.order + diff) for nbg, bond in
                              bonds[current].items() if nbg not in seen) if 1 <= bond <= 3)

                # time to yield
                if len(seen) % 2:
                    if has_hydrogen:  # enol
                        yield new_bonds
                    elif hydrogens[current]:  # ketone
                        yield new_bonds

    def __chain2circle(self):
        pass

    def __aromatic_case(self):
        pass

    def __entries(self) -> List[Tuple[int, bool]]:
        # possible: not radicals and not charged
        # O S Se -H[enol] (only 1 neighbor)
        # N -H[enol] (1 or 2 neighbor)
        # P -H[enol] (1 or 2 non-hydrogen bonds with order sum 1 or 2)

        # reverse (has double bonded neighbor):
        # O S Se (only 1 neighbor)
        # N P(1 or 2 neighbors)

        atoms = self._atoms
        charges = self._charges
        hydrogens = self._hydrogens
        neighbors = self.neighbors
        radicals = self._radicals
        hybridizations = self._hybridizations

        entries = []

        for n, a in atoms.items():
            if radicals[n] or charges[n]:
                continue
            if a.atomic_number in {8, 16, 34}:
                if neighbors(n) == 1:
                    entries.append((n, bool(hydrogens[n])))
            elif a.atomic_number in {7, 15}:
                if 0 < neighbors(n) < 3:
                    if hybridizations[n] == 1:  # amine
                        if hydrogens[n]:
                            entries.append((n, True))
                    elif hybridizations[n] == 2:  # imin
                        entries.append((n, False))
                        if hydrogens[n]:
                            entries.append((n, True))
                    elif hybridizations[n] == 3:  # nitrile
                        entries.append((n, False))
        return entries


__all__ = ['Tautomers']
