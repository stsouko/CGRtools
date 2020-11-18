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
from itertools import product
from typing import TYPE_CHECKING, Iterator
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

    def enumerate_tautomers(self, *, treat_aromatic_rings=True) -> Iterator['MoleculeContainer']:
        """
        Enumerate all possible tautomeric forms of molecule. Supported hydrogen migration through delocalized chain.

        O=C-C[H] <-> [H]O-C=C
        O=C-C=C-C[H] <-> [H]O-C=C-C=C
        O=C-C[H]-C=C <-> [H]O-C=C-C=C
        O=C-C[H]-C=C <X> O=C-C=C-C[H]  not directly possible
        """
        copy = self.copy()
        if treat_aromatic_rings:
            copy.kekule()
            copy.thiele()  # prevent
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
                mol.kekule()
                mol.thiele()
                if mol not in seen:
                    yield mol
                    seen.add(mol)
                    queue.append(mol)

    def _enumerate_zwitter_tautomers(self):
        donors, acceptors = self.__h_donors_acceptors()

        sssr = self.sssr
        rings_count = self.rings_count
        atoms_rings = self.atoms_rings

        for d, a in product(donors, acceptors):
            mol = self.copy()

            # store cached sssr in new molecules for speedup
            mol.__dict__['rings_count'] = rings_count
            mol.__dict__['atoms_rings'] = atoms_rings
            mol.__dict__['sssr'] = sssr

            charges = mol._charges
            hydrogens = mol._hydrogens
            hydrogens[d] -= 1
            hydrogens[a] += 1
            charges[d] -= 1
            charges[a] += 1
            yield mol

    def _enumerate_ring_chain_tautomers(self):
        """
        [C:1]-[O,N,S,Se:1]-[H].[C:2]=[O,N,S,Se:2] >> [C:1]-[1]-[C:2]-[2]-[H] - exo 3-7
        [C:1]-[O,N,S,Se:1]-[H].[C:2]#[N:2] >> [C:1]-[1]-[C:2]=[2]-[H] - exo 5-7
        [C:1]-[O,N,S,Se:1]-[H].[C:2]=[N:2] >> [C:1]-[1]-[2]-[C:2]-[H] - endo 5-7
        """
        return ()

    def _enumerate_chain_tautomers(self):
        atoms = self._atoms
        charges = self._charges
        radicals = self._radicals
        plane = self._plane
        bonds = self._bonds
        hybridizations = self._hybridizations
        hydrogens = self._hydrogens

        sssr = self.sssr
        rings_count = self.rings_count
        atoms_rings = self.atoms_rings

        for path in self.__enumerate_bonds():
            mol = self.__class__()
            mol._charges.update(charges)
            mol._radicals.update(radicals)
            mol._plane.update(plane)

            # store cached sssr in new molecules for speedup
            mol.__dict__['rings_count'] = rings_count
            mol.__dict__['atoms_rings'] = atoms_rings
            mol.__dict__['sssr'] = sssr

            m_atoms = mol._atoms
            m_bonds = mol._bonds
            m_hybridizations = mol._hybridizations
            m_hydrogens = mol._hydrogens

            adj = {}
            for n, a in atoms.items():
                a = a.copy()
                m_atoms[n] = a
                a._attach_to_graph(mol, n)
                adj[n] = {}

            seen = set()
            for n, m, bond in path:
                adj[n][m] = adj[m][n] = Bond(bond)
                seen.add(n)
                seen.add(m)

            for n, ms in bonds.items():
                adjn = adj[n]
                bn = m_bonds[n] = {}
                for m, bond in ms.items():
                    if m in m_bonds:  # bond partially exists. need back-connection.
                        bn[m] = m_bonds[m][n]
                    elif m in adjn:
                        bn[m] = adjn[m]
                    else:
                        bn[m] = bond.copy()

            for n, h in hybridizations.items():
                if n in seen:
                    mol._calc_hybridization(n)
                    mol._calc_implicit(n)
                else:
                    m_hybridizations[n] = h
                    m_hydrogens[n] = hydrogens[n]
            yield mol

    def __enumerate_bonds(self):
        atoms = self._atoms
        bonds = self._bonds
        hydrogens = self._hydrogens
        hyb = self._hybridizations

        # for each atom in entries
        entries, forbidden = self.__entries()
        for atom, hydrogen in entries:
            path = []
            seen = {atom}
            # stack is neighbors
            if hydrogen:  # enol
                stack = [(atom, i, n.order + 1, 0, True) for i, n in bonds[atom].items() if n.order < 3]
            else:  # ketone
                stack = [(atom, i, n.order - 1, 0, False) for i, n in bonds[atom].items() if 1 < n.order < 4]

            while stack:
                last, current, bond, depth, order_up = stack.pop()

                if len(path) > depth:
                    seen.difference_update(x for _, x, _ in path[depth:])
                    path = path[:depth]

                path.append((last, current, bond))

                # adding neighbors
                depth += 1
                seen.add(current)
                diff = -1 if order_up else 1
                stack.extend((current, n, b, depth, not order_up) for n, b in
                             ((n, b.order + diff) for n, b in bonds[current].items() if
                              n not in seen and hyb[n] < 4) if 1 <= b <= 3)

                # time to yield
                if len(path) % 2 == 0:
                    if hydrogen:
                        if current in forbidden and path[0][0] != forbidden[current]:
                            continue
                        yield path
                    elif hydrogens[current]:
                        # allene carbon formation in small rings forbidden
                        if self.rings_count and hyb[current] == 2 and atoms[current].atomic_number == 6:
                            n = next(x for x, b in bonds[current].items() if x != last and b.order != 8)
                            common = set(self.atoms_rings.get(last, ())).intersection(self.atoms_rings.get(n, ()))
                            if common and min(len(x) for x in common) < 9:
                                continue
                        yield path

    def __entries(self):
        # possible: not radicals and not charged
        # O S Se -H[enol] (only 1 neighbor)
        # N -H[enol] (1 or 2 neighbor)
        # P -H[enol] (1 or 2 non-hydrogen bonds with order sum 1 or 2)

        # reverse (has double bonded neighbor): excluded heminal heteroatoms [CC(=O)O]
        # O S Se (only 1 neighbor)
        # N P(1 or 2 neighbors)
        atoms = self._atoms
        bonds = self._bonds
        charges = self._charges
        hydrogens = self._hydrogens
        neighbors = self.neighbors
        radicals = self._radicals
        hybridizations = self._hybridizations

        entries = []
        forbidden = {}
        for n, a in atoms.items():
            if radicals[n] or charges[n]:
                continue
            if a.atomic_number in (8, 16, 34):  # O S Se
                if neighbors(n) == 1:
                    m = next(m for m, b in bonds[n].items() if b.order != 8)
                    if hydrogens[n]:
                        # skip R=N-[OH]
                        if atoms[m].atomic_number == 7 and hybridizations[m] == 2 and not (charges[m] or radicals[m]):
                            continue
                        entries.append((n, True))
                    else:
                        if atoms[m].atomic_number == 6 and hybridizations[m] == 2:  # not R=C=O
                            c = False
                            h = None
                            for x, b in bonds[m].items():
                                if x != n and b.order == 1:  # skip first atom and any (8) bonds
                                    if atoms[x].atomic_number == 6:
                                        c = True
                                    else:
                                        h = x
                            if c and h:  # skip C-C(X)=O, X != C
                                forbidden[n] = h
                                continue
                        entries.append((n, False))
            elif a.atomic_number in (7, 15):
                if 0 < neighbors(n) < 3:
                    if hybridizations[n] == 1:  # amine
                        if hydrogens[n]:
                            entries.append((n, True))
                    elif hybridizations[n] == 2:  # imine
                        if hydrogens[n]:
                            entries.append((n, True))
                        m = next(m for m, b in bonds[n].items() if b.order == 2)
                        if atoms[m].atomic_number == 6 and hybridizations[m] == 2:  # not R=C=N-R
                            c = False
                            h = None
                            for x, b in bonds[m].items():
                                if x != n and b.order != 8:  # skip first atom and any (8) bonds
                                    if atoms[x].atomic_number == 6:
                                        c = True
                                    else:
                                        h = x
                            if c and h:  # skip C-C(X)=N-R, X != C
                                forbidden[n] = h
                                continue
                        entries.append((n, False))
                    elif hybridizations[n] == 3:  # nitrile
                        entries.append((n, False))
        return entries, forbidden

    def __h_donors_acceptors(self):
        # Ar - Se S O [N+] [P+] H-donor
        # [O-], [S-], [Se-], P, N - H-acceptor
        atoms = self._atoms
        bonds = self._bonds
        charges = self._charges
        hydrogens = self._hydrogens
        radicals = self._radicals
        hybridizations = self._hybridizations

        donors = []
        acceptors = []
        for n, a in atoms.items():
            if radicals[n]:  # skip radicals
                continue
            an = a.atomic_number
            if an in (8, 16, 34):
                if not charges[n] and hydrogens[n] == 1:
                    m = next(m for m, b in bonds[n].items() if b.order != 8)
                    hm = hybridizations[m]
                    if hm == 4:  # Ar[S,Se,O][H]
                        donors.append(n)
                    else:
                        man = atoms[m].atomic_number
                        if man in (8, 15, 16, 17, 33, 34, 35, 52, 53):  # oxo-acids, except Nitro, Boronic
                            if hm != 1:
                                donors.append(n)
                        elif man == 6 and hm == 2:
                            x = next(x for x, b in bonds[m].items() if b.order == 2)
                            if atoms[x].atomic_number in (8, 16, 34):  # carboxyl and S, Se analogs
                                donors.append(n)
                elif charges[n] == -1 and not any(charges[m] > 0 for m, b in bonds[n].items()
                                                  if b.order != 8):  # R[O,S,Se-]
                    # exclude zwitterions: nitro etc
                    acceptors.append(n)
            elif an in (7, 15):
                if not charges[n]:
                    hn = hybridizations[n]
                    if hn == 1:
                        ar = 0
                        for m, b in bonds[n].items():
                            if b.order == 8:
                                continue
                            hm = hybridizations[m]
                            if hm in (2, 3):  # pyrole? enamine or amide
                                break
                            elif hm == 4:
                                ar += 1
                        else:
                            if ar < 2:  # not Ar-N-Ar or pyrole
                                acceptors.append(n)
                    elif hn == 4:  # pyridine
                        acceptors.append(n)
                    elif hn == 2:
                        # imidazole, pyrroline, guanidine, imine
                        # any other cases ignored.
                        # N=CO >> NC=O - declined
                        # [Het]=N - declined
                        # R-[O,S,Se]-C=N, R != H - accepted (1,3-oxazole, etc)
                        m = next(m for m, b in bonds[n].items() if b.order == 2)
                        if atoms[m].atomic_number == 6 and hybridizations[m] == 2:  # SP2 carbon
                            for x, b in bonds[m].items():
                                if b.order != 8 and atoms[x].atomic_number in (8, 16, 34) and hydrogens[x]:
                                    break
                            else:
                                acceptors.append(n)
                elif charges[n] == 1 and hydrogens[n]:  # R[NH+] or =[N+H2?] - ammonia or imine-H+
                    donors.append(n)
        return donors, acceptors

    def __paths(self, donors, acceptors, minimal=3, maximal=7):
        bonds = self._bonds

        for start in donors:
            path = [start]
            seen = {start}
            stack = [(n, 1) for n in bonds[start]]

            while stack:
                current, depth = stack.pop()

                if len(path) > depth:
                    seen.difference_update(path[depth:])
                    path = path[:depth]

                path.append(current)
                seen.add(current)
                depth += 1
                if depth >= minimal and current in acceptors:
                    yield path.copy()
                if depth < maximal:
                    stack.extend((n, depth) for n in bonds[current] if n not in seen)

    def __rings(self):
        atoms = self._atoms
        bonds = self._bonds
        charges = self._charges
        hydrogens = self._hydrogens
        radicals = self._radicals
        neighbors = self.neighbors
        hybridizations = self._hybridizations

        for r in self.sssr:
            if len(r) > 7 or any(bonds[n][m].order == 8 for n, m in zip(r, r[1:])):
                continue
            lr = 1 - len(r)

            for i, n in enumerate(r):
                if atoms[n].atomic_number in (7, 8, 16, 34) and neighbors(n) == 2:  # N O S Se
                    # search for [C,N;R]-[N,O,S,Se:R]-[C;R]-,=[O,N;H!R] >> [C,N]-[N,O,S,Se;H].[C]=,#[O,N]
                    # or [C,N;R]-[N,O,S,Se:R]-[N;R]-[C;H!R] >> [C,N]-[N,O,S,Se;H].[N]=[C]
                    b = i - 1
                    a = lr + i
                    if charges[b] or charges[a] or radicals[b] or radicals[a]:
                        break
                    if atoms[b].atomic_number == 7:
                        if hybridizations[b] == 2:  # C=N-[N,O,S,Se:R]
                            if atoms[a].atomic_number in (6, 7):
                                ...
                        elif hybridizations[b] == 1:  # C-N-[n:R]
                            ...
            else:  # found
                ...
        return


__all__ = ['Tautomers']
