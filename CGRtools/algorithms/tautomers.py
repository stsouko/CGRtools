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
        molecule = self.copy()
        molecule.thiele()

        # generation of all tautomers
        seen = {self, molecule}
        queue = [molecule]
        while queue:
            for molecule in queue.pop(0)._create_molecule():
                if molecule not in seen and not molecule.check_valence():
                    yield molecule
                    seen.add(molecule)
                    queue.append(molecule)

    def _create_molecule(self):
        base_atoms = self._atoms
        base_charges = self._charges
        base_radicals = self._radicals
        base_plane = self._plane
        base_bonds = self._bonds
        base_hybridizations = self._hybridizations
        base_hydrogens = self._hydrogens

        for changed_bonds in self.__changed_bonds():
            new_molecule = self.__class__()
            new_atoms = new_molecule._atoms
            new_charges = new_molecule._charges
            new_radicals = new_molecule._radicals
            new_plane = new_molecule._plane
            new_bonds = new_molecule._bonds
            new_hybridizations = new_molecule._hybridizations
            new_hydrogens = new_molecule._hydrogens

            adj = {}
            for n, a in base_atoms.items():
                a = a.copy()
                new_atoms[n] = a
                a._attach_to_graph(new_molecule, n)
                new_charges[n] = base_charges[n]
                new_radicals[n] = base_radicals[n]
                new_plane[n] = base_plane[n]
                new_bonds[n] = {}
                adj[n] = set()

            seen = set()
            for n, m, bond in changed_bonds:
                if bond is not None:
                    new_bonds[n][m] = new_bonds[m][n] = Bond(bond)
                adj[n].add(m)
                adj[m].add(n)
                seen.add(n)
                seen.add(m)

            for n, ms in base_bonds.items():
                adjn = adj[n]
                for m, bond in ms.items():
                    if m not in adjn:
                        new_bonds[n][m] = new_bonds[m][n] = bond.copy()

            for n, h in base_hybridizations.items():
                if n in seen:
                    new_molecule._calc_hybridization(n)
                    new_molecule._calc_implicit(n)
                else:
                    new_hybridizations[n] = h
                    new_hydrogens[n] = base_hydrogens[n]
            new_molecule.thiele()
            yield new_molecule

    def __changed_bonds(self):
        bonds = self._bonds
        hydrogens = self._hydrogens

        # entries for aromatic case
        entries = []

        # for each atom in entries
        for entry, has_hydrogen in self.__entries():
            seen = [entry]
            classic_case_bonds = []
            chain2circle_bonds = []
            aromatic_case_bonds = []

            # stack is neighbors
            stack = [(entry, n, b.order - 1, 1, True) for n, b in  # ketone case
                     bonds[entry].items() if 2 <= b.order <= 3]
            stack.extend((entry, n, b.order + 1, 0, False) for n, b in  # enol case
                         bonds[entry].items() if b.order <= 2 and has_hydrogen)

            # entries for aromatic case
            entries.append(entry)

            # main loop
            while stack:
                previous, current, bond, depth, flag = stack.pop()

                # branch processing
                if flag:
                    if len(seen) > depth:
                        seen = seen[:depth]
                        new_bonds = classic_case_bonds[:depth - 1]
                else:
                    if len(seen) > depth + 1:
                        seen = seen[:depth + 1]
                        new_bonds = classic_case_bonds[:depth]

                # adding new bonds
                classic_case_bonds.append((previous, current, bond))
                seen.append(current)

                # chain2circle_case
                if has_hydrogen and len(seen) >= 7 and bond == 1:
                    chain2circle_bonds.append((entry, previous, 1))
                    chain2circle_bonds.append((previous, current, bond))
                    yield chain2circle_bonds
                    chain2circle_bonds = []

                # adding neighbors
                depth += 1
                diff = -1 if depth % 2 else 1
                stack.extend((current, nbg, bond, depth, flag) for nbg, bond in
                             ((nbg, bond.order + diff) for nbg, bond in
                              bonds[current].items() if nbg not in seen) if 1 <= bond <= 3)

                # time to yield
                if len(seen) % 2:
                    if has_hydrogen:  # enol
                        yield classic_case_bonds
                    elif hydrogens[current]:  # ketone
                        yield classic_case_bonds

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
