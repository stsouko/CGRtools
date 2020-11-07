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
from collections import deque
from typing import TYPE_CHECKING, List, Tuple, Iterator
from ..containers.bonds import Bond

if TYPE_CHECKING:
    from CGRtools import MoleculeContainer


class Tautomers:
    __slots__ = ()

    def tautomerize(self) -> bool:
        """
        Convert structure to canonical tautomeric form. Return True if structure changed.
        """
        canon = min(self.enumerate_tautomers(), key=lambda x: x.huckel_pi_electrons_energy)
        if canon != self:  # attach state of canonic tautomer to self
            # atoms, radicals state, parsed_mapping and plane are unchanged
            self._bonds = canon._bonds
            self._charges = canon._charges  # for zwitter-ionic tautomers
            self._hybridizations = canon._hybridizations
            self._hydrogens = canon._hydrogens
            self._atoms_stereo = canon._atoms_stereo
            self._allenes_stereo = canon._allenes_stereo
            self._cis_trans_stereo = canon._cis_trans_stereo
            self._conformers.clear()  # flush 3d
            self.flush_cache()
            return True
        return False

    def enumerate_tautomers(self) -> Iterator['MoleculeContainer']:
        """
        Enumerate all possible tautomeric forms of molecule. Supported hydrogen migration through delocalized chain.

        O=C-C[H] <-> [H]O-C=C
        O=C-C=C-C[H] <-> [H]O-C=C-C=C
        O=C-C[H]-C=C <-> [H]O-C=C-C=C
        O=C-C[H]-C=C <X> O=C-C=C-C[H]  not directly possible
        """
        copy = self.copy().thiele()
        yield copy
        seen = {self, copy}
        queue = deque([copy])
        while queue:
            current = queue.popleft()
            for mol in current._enumerate_zwitter_tautomers():
                if mol not in seen:
                    yield mol
                    seen.add(mol)
                    queue.append(mol)
            for mol in current._enumerate_ring_chain_tautomers():
                if mol not in seen:
                    yield mol
                    seen.add(mol)
                    queue.append(mol)
            for mol in current._enumerate_chain_tautomers():
                mol.thiele()
                if mol not in seen:
                    yield mol
                    seen.add(mol)
                    queue.append(mol)

    def _enumerate_zwitter_tautomers(self):
        return ()

    def _enumerate_ring_chain_tautomers(self):
        return ()

    def _enumerate_chain_tautomers(self):
        atoms = self._atoms
        charges = self._charges
        radicals = self._radicals
        plane = self._plane
        bonds = self._bonds
        hybridizations = self._hybridizations
        hydrogens = self._hydrogens

        for path in self.__enumerate_bonds():
            mol = self.__class__()
            mol._charges.update(charges)
            mol._radicals.update(radicals)
            mol._plane.update(plane)

            m_atoms = mol._atoms
            m_bonds = mol._bonds
            m_hybridizations = mol._hybridizations
            m_hydrogens = mol._hydrogens

            adj = {}
            for n, a in atoms.items():
                a = a.copy()
                m_atoms[n] = a
                a._attach_to_graph(mol, n)
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
            yield mol

    def __enumerate_bonds(self):
        bonds = self._bonds
        hydrogens = self._hydrogens

        # for each atom in entries
        for atom, hydrogen in self.__entries():
            path = [atom]
            new_bonds = []

            # stack is neighbors
            if hydrogen:  # enol
                stack = [(i, n.order, 1) for i, n in bonds[atom].items() if n.order == 1]
            else:  # ketone
                stack = [(i, n.order, 1) for i, n in bonds[atom].items() if n.order == 2]

            while stack:
                current, bond, depth = stack.pop()

                # branch processing
                if len(path) > depth:
                    path = path[:depth]
                    new_bonds = new_bonds[:depth - 1]

                # adding new bonds
                if bond == 1:
                    new_bonds.append((path[-1], current, 2))
                else:
                    new_bonds.append((path[-1], current, 1))
                path.append(current)

                # adding neighbors
                depth += 1
                nbg = [(x, y.order, depth) for x, y in bonds[current].items() if (y.order < 3) and (y.order != bond)]
                stack.extend(nbg)

                # time to yield
                if len(path) % 2:
                    if hydrogen:  # enol
                        yield new_bonds
                    elif hydrogens[current]:  # ketone
                        yield new_bonds

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
            if a.atomic_number in (8, 16, 34):  # O S Se
                if neighbors(n) == 1:
                    entries.append((n, bool(hydrogens[n])))
            elif a.atomic_number in (7, 15):
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
