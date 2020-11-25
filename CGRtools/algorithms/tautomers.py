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
from CachedMethods import class_cached_property, cached_property
from collections import deque
from functools import reduce
from itertools import product
from operator import and_
from typing import TYPE_CHECKING, Iterator
from ..containers import query  # cyclic imports resolve
from ..containers.bonds import Bond
from ..periodictable import ListElement


if TYPE_CHECKING:
    from CGRtools import MoleculeContainer


class Tautomers:
    __slots__ = ()

    def tautomerize(self, *, prepare_molecules=True) -> bool:
        """
        Convert structure to canonical tautomeric form. Return True if structure changed.

        :param prepare_molecules: Standardize structures before. Aromatization and implicit hydrogens required.
        """
        canon = min(self.enumerate_tautomers(prepare_molecules=prepare_molecules),
                    key=lambda x: x.huckel_pi_electrons_energy)
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

    def enumerate_tautomers(self, *, prepare_molecules=True) -> Iterator['MoleculeContainer']:
        """
        Enumerate all possible tautomeric forms of molecule. Supported hydrogen migration through delocalized chain.

        O=C-C[H] <-> [H]O-C=C
        O=C-C=C-C[H] <-> [H]O-C=C-C=C
        O=C-C[H]-C=C <-> [H]O-C=C-C=C
        O=C-C[H]-C=C <X> O=C-C=C-C[H]  not directly possible

        :param prepare_molecules: Standardize structures before. Aromatization and implicit hydrogens required.
        """
        yield self.copy()
        atoms_stereo = self._atoms_stereo
        allenes_stereo = self._allenes_stereo
        cis_trans_stereo = self._cis_trans_stereo
        has_stereo = bool(atoms_stereo or allenes_stereo or cis_trans_stereo)
        bonds = self._bonds

        copy = self.copy()
        copy.clean_stereo()
        if prepare_molecules:
            copy.kekule()
            copy.implicify_hydrogens()
            copy.thiele()  # prevent

        entries = {}
        for n, h in copy._Tautomers__keto_enols[0]:
            if n in entries:
                del entries[n]
            else:
                entries[n] = h

        seen = {copy}
        queue = deque([copy])
        while queue:
            current = queue.popleft()

            for mol in current._enumerate_zwitter_tautomers():
                if mol not in seen:
                    seen.add(mol)
                    queue.append(mol)
                    if has_stereo:
                        mol = mol.copy()
                        mol._atoms_stereo.update(atoms_stereo)
                        mol._allenes_stereo.update(allenes_stereo)
                        mol._cis_trans_stereo.update(cis_trans_stereo)
                        mol._fix_stereo()
                    yield mol

            for mol, da in current._enumerate_ring_chain_tautomers():
                if mol not in seen:
                    seen.add(mol)
                    queue.append(mol)
                    if has_stereo:
                        mol = mol.copy()
                        mb = mol._bonds
                        mol._atoms_stereo.update((n, s) for n, s in atoms_stereo.items() if mb[n] == bonds[n])
                        mol._allenes_stereo.update(allenes_stereo)
                        mol._cis_trans_stereo.update(cis_trans_stereo)
                        mol._fix_stereo()
                    if da:
                        d, a = da
                        if mol._hydrogens[d]:  # imine
                            del entries[d]
                        else:
                            entries[d] = False
                        entries[a] = True
                    yield mol

            cur_entries = {}
            for mol in current._enumerate_keto_enol_tautomers():
                if mol not in seen:
                    seen.add(mol)
                    # prevent carbonyl migration
                    m_entries = mol._Tautomers__keto_enols[0]
                    if current is not copy:
                        if not cur_entries:
                            for n, h in current._Tautomers__keto_enols[0]:
                                if n in cur_entries:
                                    del cur_entries[n]
                                else:
                                    cur_entries[n] = h
                        if sum(h or -1 for n, h in m_entries if n in cur_entries and cur_entries[n] != h):
                            changes = [h or -1 for n, h in m_entries if n in entries and entries[n] != h]
                            if changes and not sum(changes):
                                continue

                    queue.append(mol)
                    if has_stereo:
                        mol = mol.copy()
                        mol._atoms_stereo.update(atoms_stereo)
                        mol._allenes_stereo.update(allenes_stereo)
                        mol._cis_trans_stereo.update(cis_trans_stereo)
                        mol._fix_stereo()
                    yield mol

    def _enumerate_zwitter_tautomers(self):
        donors, acceptors = self.__h_donors_acceptors

        sssr = self.sssr
        rings_count = self.rings_count
        atoms_rings = self.atoms_rings

        for d, a in product(donors, acceptors):
            mol = self.copy()

            # store cached sssr in new molecules for speedup
            mol.__dict__['rings_count'] = rings_count
            mol.__dict__['atoms_rings'] = atoms_rings
            mol.__dict__['sssr'] = sssr
            mol.__dict__['__cached_args_method_neighbors'] = self.__dict__['__cached_args_method_neighbors'].copy()

            charges = mol._charges
            hydrogens = mol._hydrogens
            hydrogens[d] -= 1
            hydrogens[a] += 1
            charges[d] -= 1
            charges[a] += 1
            yield mol

    def _enumerate_ring_chain_tautomers(self):
        bonds = self._bonds
        rings_count = self.rings_count
        atoms_rings = self.atoms_rings
        hyb = self._hybridizations

        for n, dnr, acc in self.__rings:
            mol = self.copy()
            bonds = mol._bonds
            hydrogens = mol._hydrogens
            hybridizations = mol._hybridizations

            del bonds[n][acc], bonds[acc][n]
            b = bonds[n][dnr]
            b._Bond__order = b._Bond__order + 1
            hydrogens[dnr] -= 1
            hydrogens[acc] += 1
            hybridizations[dnr] += 1
            hybridizations[n] += 1
            yield mol, (dnr, acc)

        donors, acceptors = self.__chains
        for d, (a, c) in product(donors, acceptors):
            path = self.__path(d, a)
            if not path:
                continue
            if len(path) == 4:  # skip 4-rings with sp2 and sp atoms
                if bonds[path[1]][path[2]] != 1:
                    continue
            elif rings_count:
                # skip 3 atoms in one aryl ring
                aryl = [x for x in path[1:-1] if hyb[x] == 4]
                if len(aryl) > 2:
                    for n, m, k in zip(aryl, aryl[1:], aryl[2:]):
                        if set(atoms_rings[n]) & set(atoms_rings[m]) & set(atoms_rings[k]):
                            continue
                # skip small rings with double bond
                if any(min(len(r) for r in atoms_rings[x]) < 7 for x in path[1:-1] if hyb[x] == 2 and x in atoms_rings):
                    continue

            mol = self.copy()
            bonds = mol._bonds
            hydrogens = mol._hydrogens
            hybridizations = mol._hybridizations

            b = bonds[a][c]
            b._Bond__order = b._Bond__order - 1
            bonds[d][a] = bonds[a][d] = Bond(1)
            hydrogens[c] += 1
            hydrogens[d] -= 1
            hybridizations[a] -= 1
            hybridizations[c] -= 1
            yield mol, None

    def _enumerate_keto_enol_tautomers(self):
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

        for path, hydrogen in self.__enumerate_keto_enol_tautomers():
            mol = self.__class__()
            mol._charges.update(charges)
            mol._radicals.update(radicals)
            mol._plane.update(plane)

            m_atoms = mol._atoms
            m_bonds = mol._bonds

            m_hybridizations = mol._hybridizations
            m_hydrogens = mol._hydrogens
            m_hybridizations.update(hybridizations)
            m_hydrogens.update(hydrogens)

            n = path[0][0]
            m = path[-1][1]
            if hydrogen:
                m_hydrogens[n] -= 1
                m_hydrogens[m] += 1
                m_hybridizations[n] += 1
                m_hybridizations[m] -= 1
            else:
                m_hydrogens[n] += 1
                m_hydrogens[m] -= 1
                m_hybridizations[n] -= 1
                m_hybridizations[m] += 1

            adj = {}
            for n, a in atoms.items():
                a = a.copy()
                m_atoms[n] = a
                a._attach_to_graph(mol, n)
                adj[n] = {}

            for n, m, bond in path:
                adj[n][m] = adj[m][n] = Bond(bond)

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

            mol.kekule()
            mol.thiele()
            # store cached sssr in new molecules for speedup
            mol.__dict__['rings_count'] = rings_count
            mol.__dict__['atoms_rings'] = atoms_rings
            mol.__dict__['sssr'] = sssr
            mol.__dict__['__cached_args_method_neighbors'] = self.__dict__['__cached_args_method_neighbors'].copy()
            yield mol

    def __enumerate_keto_enol_tautomers(self):
        atoms = self._atoms
        bonds = self._bonds
        hydrogens = self._hydrogens
        hyb = self._hybridizations

        # for each atom in entries
        entries, forbidden = self.__keto_enols
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
                        yield path, True
                    elif hydrogens[current]:
                        # allene carbon formation in small rings forbidden
                        if self.rings_count and hyb[current] == 2 and atoms[current].atomic_number == 6:
                            n = next(x for x, b in bonds[current].items() if x != last and b.order != 8)
                            common = set(self.atoms_rings.get(last, ())).intersection(self.atoms_rings.get(n, ()))
                            if common and min(len(x) for x in common) < 9:
                                continue
                        yield path, False

    @cached_property
    def __keto_enols(self):
        atoms = self._atoms
        bonds = self._bonds
        atoms_order = self.atoms_order
        connected_components = [set(x) for x in self.connected_components]

        entries = []
        forbidden = {}
        seen = set()
        for q, bl, wl, dnr, acc in self.__keto_enol_rules:
            components, closures = q._compiled_query
            for candidate in connected_components:
                for mapping in q._get_mapping(components[0], closures, atoms, bonds, candidate - seen, atoms_order):
                    n = mapping[1]
                    if n in seen:
                        continue
                    seen.add(n)
                    if bl:
                        seen.update(mapping[x] for x in bl)
                    if wl:
                        forbidden[n] = mapping[wl]
                    if dnr:
                        entries.append((n, True))
                    if acc:
                        entries.append((n, False))
        return entries, forbidden

    @cached_property
    def __h_donors_acceptors(self):
        atoms = self._atoms
        bonds = self._bonds
        atoms_order = self.atoms_order
        connected_components = [set(x) for x in self.connected_components]

        donors = []
        acceptors = []

        seen = set()
        for q, dnr, acc in self.__h_donor_acceptor_rules:
            components, closures = q._compiled_query
            for candidate in connected_components:
                for mapping in q._get_mapping(components[0], closures, atoms, bonds, candidate - seen, atoms_order):
                    n = mapping[1]
                    if n in seen:
                        continue
                    seen.add(n)
                    if dnr:
                        donors.append(n)
                    if acc:
                        acceptors.append(n)
        return donors, acceptors

    @cached_property
    def __rings(self):
        atoms = self._atoms
        bonds = self._bonds
        atoms_order = self.atoms_order
        atoms_rings = self.atoms_rings
        connected_components = [set(x) for x in self.connected_components]
        acceptable = {r for r in self.sssr if all(bonds[n][m].order != 8 for n, m in zip(r, r[1:]))}

        entries = []
        seen = set()
        for q, in_ring, dnr, acc in self.__ring_rules:
            components, closures = q._compiled_query
            for candidate in connected_components:
                for mapping in q._get_mapping(components[0], closures, atoms, bonds, candidate - seen, atoms_order):
                    n = mapping[1]
                    if n in seen:
                        continue
                    seen.add(n)
                    if n not in atoms_rings:
                        continue
                    in_ring = [mapping[x] for x in in_ring]
                    if any(x not in atoms_rings for x in in_ring):
                        continue
                    ring = reduce(and_, (set(atoms_rings[x]) for x in in_ring), set(atoms_rings[n])) & acceptable
                    if not ring:
                        continue
                    seen.update(in_ring)
                    entries.append((n, mapping[dnr], mapping[acc]))

        return entries

    @cached_property
    def __chains(self):
        atoms = self._atoms
        bonds = self._bonds
        atoms_order = self.atoms_order
        connected_components = [set(x) for x in self.connected_components]

        donors = []
        acceptors = []
        seen = set()
        for q, dnr in self.__chain_rules:
            components, closures = q._compiled_query
            for candidate in connected_components:
                for mapping in q._get_mapping(components[0], closures, atoms, bonds, candidate - seen, atoms_order):
                    n = mapping[1]
                    if n in seen:
                        continue
                    seen.add(n)
                    if dnr:
                        donors.append(n)
                    else:
                        m = mapping[2]
                        if m in seen:
                            continue
                        seen.add(m)
                        acceptors.append((n, m))

        return donors, acceptors

    @class_cached_property
    def __keto_enol_rules(self):
        rules = []  # query, black list [except first], allowed, H-donor, H-acceptor

        # C(=[O,S,Se])[O,S,Se][H,R] skip
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'Se']), neighbors=1)
        q.add_atom('C')
        q.add_atom(ListElement(['O', 'S', 'Se']), neighbors=(1, 2), hybridization=1)
        q.add_bond(1, 2, 2)
        q.add_bond(2, 3, 1)
        rules.append((q, None, 3, False, False))

        # C(=[O,S,Se])[O,S,Se-] skip
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'Se']), neighbors=1)
        q.add_atom('C')
        q.add_atom(ListElement(['O', 'S', 'Se']), charge=-1, hybridization=1)
        q.add_bond(1, 2, 2)
        q.add_bond(2, 3, 1)
        rules.append((q, None, None, False, False))

        # C(=[O,S,Se])N[H,R] skip
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'Se']), neighbors=1)
        q.add_atom('C')
        q.add_atom('N', hybridization=1)
        q.add_bond(1, 2, 2)
        q.add_bond(2, 3, 1)
        rules.append((q, None, 3, False, False))

        # C(=N)N[H,R] skip
        q = query.QueryContainer()
        q.add_atom('N', hybridization=2)
        q.add_atom('C')
        q.add_atom('N', hybridization=1)
        q.add_bond(1, 2, 2)
        q.add_bond(2, 3, 1)
        rules.append((q, None, 3, False, False))

        # C(=N)[O,S,Se][H,R] skip
        q = query.QueryContainer()
        q.add_atom('N', hybridization=2)
        q.add_atom('C')
        q.add_atom(ListElement(['O', 'S', 'Se']), neighbors=(1, 2), hybridization=1)
        q.add_bond(1, 2, 2)
        q.add_bond(2, 3, 1)
        rules.append((q, None, 3, False, False))

        # A=N-OH skip
        q = query.QueryContainer()
        q.add_atom('O', neighbors=1)
        q.add_atom('N', neighbors=2, hybridization=2)
        q.add_bond(1, 2, 1)
        rules.append((q, (2,), None, False, False))

        # [O,S,Se,N]=[PH]= skip
        q = query.QueryContainer()
        q.add_atom('P', neighbors=2, hybridization=3)
        q.add_atom(ListElement(['O', 'S', 'Se', 'N']))
        q.add_bond(1, 2, 2)
        rules.append((q, None, None, False, False))

        # A=[O,S,Se]
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'Se']), neighbors=1, hybridization=2)
        rules.append((q, None, None, False, True))

        # A-[O,S,Se]H
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'Se']), hybridization=1, neighbors=1)
        rules.append((q, None, None, True, False))

        # A#N
        q = query.QueryContainer()
        q.add_atom('N', hybridization=3)
        rules.append((q, None, None, False, True))

        # A=NH
        q = query.QueryContainer()
        q.add_atom('N', hybridization=2, neighbors=1)
        rules.append((q, None, None, True, True))

        # A=N
        q = query.QueryContainer()
        q.add_atom('N', hybridization=2, neighbors=2)
        rules.append((q, None, None, False, True))

        # A-NH
        q = query.QueryContainer()
        q.add_atom('N', hybridization=1, neighbors=(1, 2))
        rules.append((q, None, None, True, False))

        return rules

    @class_cached_property
    def __h_donor_acceptor_rules(self):
        rules = []  # query, H-donor, H-acceptor

        # [O,S,Se-]-[N+] skip
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'Se']), charge=-1)
        q.add_atom('N', charge=1)
        q.add_bond(1, 2, 1)
        rules.append((q, False, False))

        # [O,S,N,Se]-N= skip
        q = query.QueryContainer()
        q.add_atom('N', hybridization=2)
        q.add_atom(ListElement(['N', 'O', 'S', 'Se']))
        q.add_bond(1, 2, 1)
        rules.append((q, False, False))

        # [NH]-A= skip
        q = query.QueryContainer()
        q.add_atom('N', hybridization=1, neighbors=(1, 2))
        q.add_atom('A', hybridization=(2, 3))
        q.add_bond(1, 2, 1)
        rules.append((q, False, False))

        # [H][O,S,Se]-Ar
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'Se']), neighbors=1)
        q.add_atom('A', hybridization=4)
        q.add_bond(1, 2, 1)
        rules.append((q, True, False))

        # [H][O,S,Se]-[C,N,P,S,Cl,Se,Br,I]=O
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'Se']), neighbors=1)
        q.add_atom(ListElement(['C', 'N', 'P', 'S', 'Cl', 'Se', 'Br', 'I']))
        q.add_atom('O')
        q.add_bond(1, 2, 1)
        q.add_bond(2, 3, 2)
        rules.append((q, True, False))

        # [H][N+]
        q = query.QueryContainer()
        q.add_atom('N', charge=1, neighbors=(1, 2, 3), hybridization=1)
        rules.append((q, True, False))

        # [H][n+]
        q = query.QueryContainer()
        q.add_atom('N', charge=1, neighbors=2, hybridization=4)
        rules.append((q, True, False))

        # [H][N+]=
        q = query.QueryContainer()
        q.add_atom('N', charge=1, neighbors=(1, 2), hybridization=2)
        rules.append((q, True, False))

        # [O,S,Se-]-
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'Se']), charge=-1)
        rules.append((q, False, True))

        # N=,# or NR1-3
        q = query.QueryContainer()
        q.add_atom('N')
        rules.append((q, False, True))

        return rules

    @class_cached_property
    def __ring_rules(self):
        rules = []  # query, in ring [except first], H-donor, H-acceptor

        # -CC(OH)[O,S]-
        q = query.QueryContainer()
        q.add_atom('C', neighbors=3)  # entry
        q.add_atom('C')
        q.add_atom(ListElement(['O', 'S']), hybridization=1, neighbors=2)
        q.add_atom('O', neighbors=1)
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        rules.append((q, (2, 3), 4, 3))

        # -CC(OH)N(R)-
        q = query.QueryContainer()
        q.add_atom('C', neighbors=3)  # entry
        q.add_atom('C')
        q.add_atom('N', neighbors=3)
        q.add_atom('O', neighbors=1)
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        rules.append((q, (2, 3), 4, 3))

        # -CC(C)(OH)[O,S]-
        q = query.QueryContainer()
        q.add_atom('C')  # entry
        q.add_atom('C')
        q.add_atom(ListElement(['O', 'S']), hybridization=1, neighbors=2)
        q.add_atom('O', neighbors=1)
        q.add_atom('C')
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        q.add_bond(1, 5, 1)
        rules.append((q, (2, 3), 4, 3))

        # -CC(C)(OH)N(R)-
        q = query.QueryContainer()
        q.add_atom('C')  # entry
        q.add_atom('C')
        q.add_atom('N', neighbors=3)
        q.add_atom('O', neighbors=1)
        q.add_atom('C')
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        q.add_bond(1, 5, 1)
        rules.append((q, (2, 3), 4, 3))

        # -CC(NH)[O,S]-
        q = query.QueryContainer()
        q.add_atom('C', neighbors=3)  # entry
        q.add_atom('C')
        q.add_atom(ListElement(['O', 'S']), hybridization=1, neighbors=2)
        q.add_atom('N', neighbors=(1, 2), hybridization=1)
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        rules.append((q, (2, 3), 4, 3))

        # -CC(NH)N(R)-
        q = query.QueryContainer()
        q.add_atom('C', neighbors=3)  # entry
        q.add_atom('C')
        q.add_atom('N', neighbors=3)
        q.add_atom('N', neighbors=(1, 2), hybridization=1)
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        rules.append((q, (2, 3), 4, 3))

        # -CC(C)(NH)[O,S]-
        q = query.QueryContainer()
        q.add_atom('C')  # entry
        q.add_atom('C')
        q.add_atom(ListElement(['O', 'S']), hybridization=1, neighbors=2)
        q.add_atom('N', neighbors=(1, 2), hybridization=1)
        q.add_atom('C')
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        q.add_bond(1, 5, 1)
        rules.append((q, (2, 3), 4, 3))

        # -CC(C)(NH)N(R)-
        q = query.QueryContainer()
        q.add_atom('C')  # entry
        q.add_atom('C')
        q.add_atom('N', neighbors=3)
        q.add_atom('N', neighbors=(1, 2), hybridization=1)
        q.add_atom('C')
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        q.add_bond(1, 5, 1)
        rules.append((q, (2, 3), 4, 3))

        # -CC(=NH)[O,S]-
        q = query.QueryContainer()
        q.add_atom('C')  # entry
        q.add_atom('C')
        q.add_atom(ListElement(['O', 'S']), hybridization=1, neighbors=2)
        q.add_atom('N', neighbors=1)
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 2)
        rules.append((q, (2, 3), 4, 3))

        # -CC(=NH)N-
        q = query.QueryContainer()
        q.add_atom('C')  # entry
        q.add_atom('C')
        q.add_atom('N', hybridization=1)
        q.add_atom('N', neighbors=1)
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 2)
        rules.append((q, (2, 3), 4, 3))

        return rules

    @class_cached_property
    def __chain_rules(self):
        rules = []  # query, is donor

        # N-OH  imines azoles etc
        q = query.QueryContainer()
        q.add_atom('O', neighbors=1)
        q.add_atom('N')
        q.add_bond(1, 2, 1)
        rules.append((q, True))

        # C[O,S]H
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S']), neighbors=1)
        q.add_atom('C', hybridization=(1, 2, 4))
        q.add_bond(1, 2, 1)
        rules.append((q, True))

        # C-[NH]-C  2nd amines
        q = query.QueryContainer()
        q.add_atom('N', neighbors=2)
        q.add_atom('C', hybridization=(1, 4))
        q.add_atom('C', hybridization=1)
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        rules.append((q, True))

        # C-C(C)=[O,S,NH]
        q = query.QueryContainer()
        q.add_atom('C')
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=1)
        q.add_atom('C')
        q.add_atom('C')
        q.add_bond(1, 2, 2)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        rules.append((q, False))

        # C-C=[O,S,NH]
        q = query.QueryContainer()
        q.add_atom('C', neighbors=2)
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=1)
        q.add_atom('C')
        q.add_bond(1, 2, 2)
        q.add_bond(1, 3, 1)
        rules.append((q, False))

        # C-C(C)=NR
        q = query.QueryContainer()
        q.add_atom('C')
        q.add_atom('N', neighbors=2)
        q.add_atom('C')
        q.add_atom('C')
        q.add_bond(1, 2, 2)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        rules.append((q, False))

        # C-C=NR
        q = query.QueryContainer()
        q.add_atom('C', neighbors=2)
        q.add_atom('N', neighbors=2)
        q.add_atom('C')
        q.add_bond(1, 2, 2)
        q.add_bond(1, 3, 1)
        rules.append((q, False))

        # C-C#N
        q = query.QueryContainer()
        q.add_atom('C')
        q.add_atom('N')
        q.add_atom('C')
        q.add_bond(1, 2, 3)
        q.add_bond(1, 3, 1)
        rules.append((q, False))

        return rules

    def __path(self, donor, acceptor, minimal=4, maximal=7):
        bonds = self._bonds

        path = [donor]
        seen = {donor}
        stack = [(n, 1) for n in bonds[donor]]

        while stack:
            current, depth = stack.pop()

            if len(path) > depth:
                seen.difference_update(path[depth:])
                path = path[:depth]

            path.append(current)
            seen.add(current)
            depth += 1
            if depth >= minimal and current == acceptor:
                return path
            if depth < maximal:
                stack.extend((n, depth) for n in bonds[current] if n not in seen)


__all__ = ['Tautomers']
