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
from itertools import product, chain, repeat
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

        seen = {copy: None}
        queue = deque([copy])
        while queue:
            current = queue.popleft()

            for mol in current._enumerate_zwitter_tautomers():
                if mol not in seen:
                    seen[mol] = current
                    queue.append(mol)
                    if has_stereo:
                        mol = mol.copy()
                        mol._atoms_stereo.update(atoms_stereo)
                        mol._allenes_stereo.update(allenes_stereo)
                        mol._cis_trans_stereo.update(cis_trans_stereo)
                        mol._fix_stereo()
                    yield mol

            for mol in current._enumerate_ring_chain_tautomers():
                if mol not in seen:
                    seen[mol] = current
                    queue.append(mol)
                    if has_stereo:
                        mol = mol.copy()
                        mb = mol._bonds
                        mol._atoms_stereo.update((n, s) for n, s in atoms_stereo.items() if mb[n] == bonds[n])
                        mol._allenes_stereo.update(allenes_stereo)
                        mol._cis_trans_stereo.update(cis_trans_stereo)
                        mol._fix_stereo()
                    yield mol

            for mol, ket in current._enumerate_keto_enol_tautomers():
                if mol not in seen:
                    seen[mol] = current
                    # prevent carbonyl migration
                    if current is not copy and not ket:  # enol to ket potentially migrate ketone.
                        # search alpha hydroxy ketone inversion
                        before = seen[current]._sugar_groups
                        if any((k, e) in before for e, k in mol._sugar_groups):
                            continue

                    if has_stereo:
                        ster = mol.copy()
                        ster._atoms_stereo.update(atoms_stereo)
                        ster._allenes_stereo.update(allenes_stereo)
                        ster._cis_trans_stereo.update(cis_trans_stereo)
                        ster._fix_stereo()
                        yield ster
                    else:
                        yield mol
                    if len(mol.aromatic_rings) > len(current.aromatic_rings):
                        # found new aromatic ring. flush queue and start from it.
                        queue.clear()
                        queue.append(mol)
                        break
                    queue.append(mol)

    def _enumerate_zwitter_tautomers(self):
        atoms = self._atoms
        bonds = self._bonds
        atoms_order = self.atoms_order
        connected_components = [set(x) for x in self.connected_components]
        sssr = self.sssr
        rings_count = self.rings_count
        atoms_rings = self.atoms_rings
        atoms_rings_sizes = self.atoms_rings_sizes

        donors = []
        acceptors = []
        for q, acid in chain(zip(self.__acid_rules, repeat(True)), zip(self.__base_rules, repeat(False))):
            components, closures = q._compiled_query
            for candidate in connected_components:
                for mapping in q._get_mapping(components[0], closures, atoms, bonds, candidate, atoms_order):
                    n = mapping[1]
                    if acid:
                        donors.append(n)
                    else:
                        acceptors.append(n)

        for d, a in product(donors, acceptors):
            mol = self.copy()
            m_charges = mol._charges
            m_hydrogens = mol._hydrogens
            m_hydrogens[d] -= 1
            m_hydrogens[a] += 1
            m_charges[d] -= 1
            m_charges[a] += 1

            # store cached sssr in new molecules for speedup
            mol.__dict__['rings_count'] = rings_count
            mol.__dict__['atoms_rings'] = atoms_rings
            mol.__dict__['sssr'] = sssr
            mol.__dict__['atoms_rings_sizes'] = atoms_rings_sizes
            mol.__dict__['__cached_args_method_neighbors'] = self.__dict__['__cached_args_method_neighbors'].copy()
            yield mol

    def _enumerate_ring_chain_tautomers(self):
        atoms = self._atoms
        bonds = self._bonds
        atoms_order = self.atoms_order
        connected_components = [set(x) for x in self.connected_components]

        for q, t in chain(zip(self.__ring_rules, repeat(True)), zip(self.__chain_rules, repeat(False))):
            components, closures = q._compiled_query
            for candidate in connected_components:
                for mapping in q._get_mapping(components[0], closures, atoms, bonds, candidate, atoms_order):
                    n, m, k = mapping[1], mapping[2], mapping[3]
                    mol = self.copy()
                    m_bonds = mol._bonds
                    m_hydrogens = mol._hydrogens
                    m_hybridizations = mol._hybridizations

                    if t:
                        del m_bonds[n][k], m_bonds[k][n]
                        m_bonds[n][m]._Bond__order = 2
                        m_hydrogens[m] -= 1
                        m_hydrogens[k] += 1
                        m_hybridizations[m] += 1
                        m_hybridizations[n] += 1
                    else:
                        m_bonds[n][m]._Bond__order = 1
                        m_bonds[n][k] = m_bonds[k][n] = Bond(1)
                        m_hydrogens[k] -= 1
                        m_hydrogens[m] += 1
                        m_hybridizations[n] -= 1
                        m_hybridizations[m] -= 1
                    yield mol

    def _enumerate_keto_enol_tautomers(self):
        atoms = self._atoms
        bonds = self._bonds
        atoms_order = self.atoms_order
        connected_components = [set(x) for x in self.connected_components]

        sssr = self.sssr
        rings_count = self.rings_count
        atoms_rings = self.atoms_rings
        atoms_rings_sizes = self.atoms_rings_sizes
        neighbors = self.__dict__['__cached_args_method_neighbors']

        seen = set()
        for q, ket, *fix in self.__keto_enol_rules:
            components, closures = q._compiled_query
            for candidate in connected_components:
                for mapping in q._get_mapping(components[0], closures, atoms, bonds, candidate - seen, atoms_order):
                    a = mapping[1]
                    d = mapping[2]
                    if not ket:  # prevent 1-3 migration in 1-5 etc
                        if d in seen:
                            continue
                        seen.add(d)

                    mol = self.copy()
                    m_bonds = mol._bonds
                    m_hydrogens = mol._hydrogens
                    m_hybridizations = mol._hybridizations

                    for n, m, b in fix:
                        m_bonds[mapping[n]][mapping[m]]._Bond__order = b

                    m_hydrogens[a] += 1
                    m_hydrogens[d] -= 1
                    m_hybridizations[a] -= 1
                    m_hybridizations[d] += 1

                    # store cached sssr in new molecules for speedup
                    mol.__dict__['sssr'] = sssr
                    mol.__dict__['__cached_args_method_neighbors'] = neighbors.copy()
                    mol.kekule()
                    mol.thiele()
                    # restore after kekule-thiele
                    mol.__dict__['sssr'] = sssr
                    mol.__dict__['__cached_args_method_neighbors'] = neighbors.copy()
                    mol.__dict__['rings_count'] = rings_count
                    mol.__dict__['atoms_rings'] = atoms_rings
                    mol.__dict__['atoms_rings_sizes'] = atoms_rings_sizes
                    yield mol, ket

    @cached_property
    def _sugar_groups(self):
        atoms = self._atoms
        bonds = self._bonds
        atoms_order = self.atoms_order
        connected_components = [set(x) for x in self.connected_components]

        ek = []
        for q in self.__sugar_group_rules:
            components, closures = q._compiled_query
            for candidate in connected_components:
                for mapping in q._get_mapping(components[0], closures, atoms, bonds, candidate, atoms_order):
                    e, k = mapping[1], mapping[2]
                    ek.append((e, k))
        return ek

    @class_cached_property
    def __keto_enol_rules(self):
        rules = []  # query, is ketone, bond fix

        # first atom is acceptor of proton
        # second is donor
        # [C;X4,H]-[CH]=[O,S,N]
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=1)  # acceptor
        q.add_atom('C', hybridization=1, hydrogens=(1, 2, 3))  # donor
        q.add_atom('C', neighbors=2)
        q.add_bond(1, 3, 2)
        q.add_bond(2, 3, 1)
        rules.append((q, True, (1, 3, 1), (2, 3, 2)))

        # [C;X4,H]-[CH]=N[C;a,X4]
        q = query.QueryContainer()
        q.add_atom('N')  # acceptor
        q.add_atom('C', hybridization=1, hydrogens=(1, 2, 3))  # donor
        q.add_atom('C', neighbors=2)
        q.add_atom('C', hybridization=(1, 4))
        q.add_bond(1, 3, 2)
        q.add_bond(2, 3, 1)
        q.add_bond(1, 4, 1)
        rules.append((q, True, (1, 3, 1), (2, 3, 2)))

        # [C;X4,H]-C(C)=[O,S,N;H]
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=1)  # acceptor
        q.add_atom('C', hybridization=1, hydrogens=(1, 2, 3))  # donor
        q.add_atom('C')
        q.add_atom('C')
        q.add_bond(1, 3, 2)
        q.add_bond(2, 3, 1)
        q.add_bond(3, 4, 1)
        rules.append((q, True, (1, 3, 1), (2, 3, 2)))

        # [C;X4,H]-C(C)=N[C;a,X4]
        q = query.QueryContainer()
        q.add_atom('N')  # acceptor
        q.add_atom('C', hybridization=1, hydrogens=(1, 2, 3))  # donor
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C', hybridization=(1, 4))
        q.add_bond(1, 3, 2)
        q.add_bond(2, 3, 1)
        q.add_bond(3, 4, 1)
        q.add_bond(1, 5, 1)
        rules.append((q, True, (1, 3, 1), (2, 3, 2)))

        # C=[CH]-[O,S,N;H]
        q = query.QueryContainer()
        q.add_atom('C')  # acceptor
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=(1, 2))  # donor
        q.add_atom('C', neighbors=2)
        q.add_bond(1, 3, 2)
        q.add_bond(2, 3, 1)
        rules.append((q, False, (1, 3, 1), (2, 3, 2)))

        # C=C(C)-[O,S,N;H]
        q = query.QueryContainer()
        q.add_atom('C')  # acceptor
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=(1, 2))  # donor
        q.add_atom('C')
        q.add_atom('C')
        q.add_bond(1, 3, 2)
        q.add_bond(2, 3, 1)
        q.add_bond(4, 3, 1)
        rules.append((q, False, (1, 3, 1), (2, 3, 2)))

        # azo-Ar
        # C=N-[NH]-R >> C-N=N-R
        q = query.QueryContainer()
        q.add_atom('N', hydrogens=(1, 2))
        q.add_atom('N')
        q.add_atom('C', hybridization=2)
        q.add_bond(1, 2, 1)
        q.add_bond(2, 3, 2)
        # rules.append(q)

        # C=C([O,S,N]H)-[O,S,N]
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=(1, 2))
        q.add_atom('C')
        q.add_atom(ListElement(['O', 'S', 'N']), hybridization=1)
        q.add_atom('C', hybridization=2)
        q.add_bond(1, 2, 1)
        q.add_bond(2, 3, 1)
        q.add_bond(2, 4, 2)
        # rules.append(q)

        # [NH:R6]-[C:R6](=[O,S,NH])[C:R6]. alpha-pyridone
        q = query.QueryContainer()
        q.add_atom('C', rings_sizes=6)
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=1)
        q.add_atom('N', rings_sizes=6, neighbors=2, hybridization=1)
        q.add_atom('C', rings_sizes=6, hybridization=(2, 4))
        q.add_bond(1, 2, 2)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        # rules.append(q)

        return rules

    @class_cached_property
    def __acid_rules(self):
        rules = []

        # [H][O,S,Se]-Ar
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'Se']), neighbors=1)
        q.add_atom('A', hybridization=4)
        q.add_bond(1, 2, 1)
        rules.append(q)

        # [H][O,S,Se]-[C,N,P,S,Cl,Se,Br,I]=O
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'Se']), neighbors=1)
        q.add_atom(ListElement(['C', 'N', 'P', 'S', 'Cl', 'Se', 'Br', 'I']))
        q.add_atom('O')
        q.add_bond(1, 2, 1)
        q.add_bond(2, 3, 2)
        rules.append(q)

        # [H][N+]=,-,:
        q = query.QueryContainer()
        q.add_atom('N', charge=1, hydrogens=(1, 2, 3))
        rules.append(q)

        return rules

    @class_cached_property
    def __base_rules(self):
        rules = []

        # [O,S,Se-]-[C,N,P,S,Cl,Se,Br,I]=O
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'Se']), charge=-1)
        q.add_atom(ListElement(['C', 'N', 'P', 'S', 'Cl', 'Se', 'Br', 'I']))
        q.add_atom('O')
        q.add_bond(1, 2, 1)
        q.add_bond(2, 3, 2)
        rules.append(q)

        # [O,S,Se-]-Ar
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'Se']), charge=-1)
        q.add_atom('A', hybridization=4)
        q.add_bond(1, 2, 1)
        rules.append(q)

        # HN=CC
        q = query.QueryContainer()
        q.add_atom('N', neighbors=1)
        q.add_atom('C', neighbors=2)
        q.add_atom('C')
        q.add_bond(1, 2, 2)
        q.add_bond(2, 3, 1)
        rules.append(q)

        # HN=C(C)C
        q = query.QueryContainer()
        q.add_atom('N', neighbors=1)
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C')
        q.add_bond(1, 2, 2)
        q.add_bond(2, 3, 1)
        q.add_bond(2, 4, 1)
        rules.append(q)

        # CN=CC
        q = query.QueryContainer()
        q.add_atom('N')
        q.add_atom('C', neighbors=2)
        q.add_atom('C')
        q.add_atom('C')
        q.add_bond(1, 2, 2)
        q.add_bond(1, 3, 1)
        q.add_bond(2, 4, 1)
        rules.append(q)

        # CN=C(C)C
        q = query.QueryContainer()
        q.add_atom('N')
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C')
        q.add_bond(1, 2, 2)
        q.add_bond(1, 3, 1)
        q.add_bond(2, 4, 1)
        q.add_bond(2, 5, 1)
        rules.append(q)

        # H2N-C
        q = query.QueryContainer()
        q.add_atom('N', neighbors=1)
        q.add_atom('C', hybridization=1)
        q.add_bond(1, 2, 1)
        rules.append(q)

        # HN(C)C
        q = query.QueryContainer()
        q.add_atom('N', neighbors=2)
        q.add_atom('C', hybridization=(1, 4))
        q.add_atom('C', hybridization=1)
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        rules.append(q)

        # CN(C)C
        q = query.QueryContainer()
        q.add_atom('N')
        q.add_atom('C', hybridization=(1, 4))
        q.add_atom('C', hybridization=1)
        q.add_atom('C', hybridization=1)
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        rules.append(q)

        # :N:
        q = query.QueryContainer()
        q.add_atom('N', hybridization=4, hydrogens=0)
        rules.append(q)

        return rules

    @class_cached_property
    def __ring_rules(self):
        rules = []  # connection, H-donor, H-acceptor

        # aldehydes

        # C1([O,S,N;H])[O,S,N;!H]-CC1
        q = query.QueryContainer()
        q.add_atom('C', neighbors=3, rings_sizes=4)  # entry
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=(1, 2))
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=(2, 3), hybridization=1, hydrogens=0)
        q.add_atom('C', hybridization=1)
        q.add_atom('C', hybridization=1)
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        q.add_bond(3, 5, 1)
        q.add_bond(4, 5, 1)
        rules.append(q)

        # C1([O,S,N;H1])[O,S,N;!H]-[C;X4,a]:,-C-,=C1
        q = query.QueryContainer()
        q.add_atom('C', neighbors=3, rings_sizes=5)  # entry
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=(1, 2))
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=(2, 3), hybridization=1, hydrogens=0)
        q.add_atom('C')
        q.add_atom('C', hybridization=(1, 4))
        q.add_atom('C')
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        q.add_bond(3, 5, 1)
        q.add_bond(4, 6, (1, 2))
        q.add_bond(5, 6, (1, 4))
        rules.append(q)

        # C1([O,S,N;H])[O,S,N;!H]-[C;X4,a]:,-C-,=[C;X3,X4]-,=C1
        q = query.QueryContainer()
        q.add_atom('C', neighbors=3, rings_sizes=6)  # entry
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=(1, 2))
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=(2, 3), hybridization=1, hydrogens=0)
        q.add_atom('C')
        q.add_atom('C', hybridization=(1, 4))
        q.add_atom('C', hybridization=(1, 2))
        q.add_atom('C')
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        q.add_bond(3, 5, 1)
        q.add_bond(4, 6, (1, 2))
        q.add_bond(5, 7, (1, 4))
        q.add_bond(6, 7, (1, 2))
        rules.append(q)

        # C1([O,S,N;H])[O,S,N;!H]-[C;X4,a]:,-C-[S,O,N]-C1
        q = query.QueryContainer()
        q.add_atom('C', neighbors=3, rings_sizes=6)  # entry
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=(1, 2))
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=(2, 3), hybridization=1, hydrogens=0)
        q.add_atom('C')
        q.add_atom('C', hybridization=(1, 4))
        q.add_atom(ListElement(['S', 'O', 'N']))
        q.add_atom('C')
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        q.add_bond(3, 5, 1)
        q.add_bond(4, 6, 1)
        q.add_bond(5, 7, (1, 4))
        q.add_bond(6, 7, 1)
        rules.append(q)

        # ketones

        # C1([O,S,N;H])[O,S,N;!H]-CC1
        q = query.QueryContainer()
        q.add_atom('C', rings_sizes=4)  # entry
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=(1, 2))
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=(2, 3), hybridization=1, hydrogens=0)
        q.add_atom('C', hybridization=1)
        q.add_atom('C', hybridization=1)
        q.add_atom('C')
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        q.add_bond(3, 5, 1)
        q.add_bond(4, 5, 1)
        q.add_bond(1, 6, 1)
        rules.append(q)

        # C1([O,S,N;H1])[O,S,N;!H]-[C;X4,a]:,-C-,=C1
        q = query.QueryContainer()
        q.add_atom('C', rings_sizes=5)  # entry
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=(1, 2))
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=(2, 3), hybridization=1, hydrogens=0)
        q.add_atom('C')
        q.add_atom('C', hybridization=(1, 4))
        q.add_atom('C')
        q.add_atom('C')
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        q.add_bond(3, 5, 1)
        q.add_bond(4, 6, (1, 2))
        q.add_bond(5, 6, (1, 4))
        q.add_bond(1, 7, 1)
        rules.append(q)

        # C1([O,S,N;H])[O,S,N;!H]-[C;X4,a]:,-C-,=[C;X3,X4]-,=C1
        q = query.QueryContainer()
        q.add_atom('C', rings_sizes=6)  # entry
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=(1, 2))
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=(2, 3), hybridization=1, hydrogens=0)
        q.add_atom('C')
        q.add_atom('C', hybridization=(1, 4))
        q.add_atom('C', hybridization=(1, 2))
        q.add_atom('C')
        q.add_atom('C')
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        q.add_bond(3, 5, 1)
        q.add_bond(4, 6, (1, 2))
        q.add_bond(5, 7, (1, 4))
        q.add_bond(6, 7, (1, 2))
        q.add_bond(1, 8, 1)
        rules.append(q)

        # C1([O,S,N;H])[O,S,N;!H]-[C;X4,a]:,-C-[S,O,N]-C1
        q = query.QueryContainer()
        q.add_atom('C', rings_sizes=6)  # entry
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=(1, 2))
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=(2, 3), hybridization=1, hydrogens=0)
        q.add_atom('C')
        q.add_atom('C', hybridization=(1, 4))
        q.add_atom(ListElement(['S', 'O', 'N']))
        q.add_atom('C')
        q.add_atom('C')
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        q.add_bond(3, 5, 1)
        q.add_bond(4, 6, 1)
        q.add_bond(5, 7, (1, 4))
        q.add_bond(6, 7, 1)
        q.add_bond(1, 8, 1)
        rules.append(q)

        return rules

    @class_cached_property
    def __chain_rules(self):
        rules = []  # connection, H-acceptor, H-donor

        # aldehydes

        # C1([O,S,N;H])[O,S,N;!H]-CC1
        q = query.QueryContainer()
        q.add_atom('C', neighbors=2)  # entry
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=(1, 2), hybridization=2)
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=1)
        q.add_atom('C', hybridization=1)
        q.add_atom('C', hybridization=1)
        q.add_bond(1, 2, 2)
        q.add_bond(1, 4, 1)
        q.add_bond(3, 5, 1)
        q.add_bond(4, 5, 1)
        rules.append(q)

        # C1([O,S,N;H1])[O,S,N;!H]-[C;X4,a]:,-C-,=C1
        q = query.QueryContainer()
        q.add_atom('C', neighbors=2)  # entry
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=(1, 2), hybridization=2)
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=1)
        q.add_atom('C')
        q.add_atom('C', hybridization=(1, 4))
        q.add_atom('C')
        q.add_bond(1, 2, 2)
        q.add_bond(1, 4, 1)
        q.add_bond(3, 5, 1)
        q.add_bond(4, 6, (1, 2))
        q.add_bond(5, 6, (1, 4))
        rules.append(q)

        # C1([O,S,N;H])[O,S,N;!H]-[C;X4,a]:,-C-,=[C;X3,X4]-,=C1
        q = query.QueryContainer()
        q.add_atom('C', neighbors=2)  # entry
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=(1, 2), hybridization=2)
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=1)
        q.add_atom('C')
        q.add_atom('C', hybridization=(1, 4))
        q.add_atom('C', hybridization=(1, 2))
        q.add_atom('C')
        q.add_bond(1, 2, 2)
        q.add_bond(1, 4, 1)
        q.add_bond(3, 5, 1)
        q.add_bond(4, 6, (1, 2))
        q.add_bond(5, 7, (1, 4))
        q.add_bond(6, 7, (1, 2))
        rules.append(q)

        # C1([O,S,N;H])[O,S,N;!H]-[C;X4,a]:,-C-[S,O,N]-C1
        q = query.QueryContainer()
        q.add_atom('C', neighbors=2)  # entry
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=(1, 2), hybridization=2)
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=1)
        q.add_atom('C')
        q.add_atom('C', hybridization=(1, 4))
        q.add_atom(ListElement(['S', 'O', 'N']))
        q.add_atom('C')
        q.add_bond(1, 2, 2)
        q.add_bond(1, 4, 1)
        q.add_bond(3, 5, 1)
        q.add_bond(4, 6, 1)
        q.add_bond(5, 7, (1, 4))
        q.add_bond(6, 7, 1)
        rules.append(q)

        # ketones

        # C1([O,S,N;H])[O,S,N;!H]-CC1
        q = query.QueryContainer()
        q.add_atom('C')  # entry
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=(1, 2), hybridization=2)
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=1)
        q.add_atom('C', hybridization=1)
        q.add_atom('C', hybridization=1)
        q.add_atom('C')
        q.add_bond(1, 2, 2)
        q.add_bond(1, 4, 1)
        q.add_bond(3, 5, 1)
        q.add_bond(4, 5, 1)
        q.add_bond(1, 6, 1)
        rules.append(q)

        # C1([O,S,N;H1])[O,S,N;!H]-[C;X4,a]:,-C-,=C1
        q = query.QueryContainer()
        q.add_atom('C')  # entry
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=(1, 2), hybridization=2)
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=1)
        q.add_atom('C')
        q.add_atom('C', hybridization=(1, 4))
        q.add_atom('C')
        q.add_atom('C')
        q.add_bond(1, 2, 2)
        q.add_bond(1, 4, 1)
        q.add_bond(3, 5, 1)
        q.add_bond(4, 6, (1, 2))
        q.add_bond(5, 6, (1, 4))
        q.add_bond(1, 7, 1)
        rules.append(q)

        # C1([O,S,N;H])[O,S,N;!H]-[C;X4,a]:,-C-,=[C;X3,X4]-,=C1
        q = query.QueryContainer()
        q.add_atom('C')  # entry
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=(1, 2), hybridization=2)
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=1)
        q.add_atom('C')
        q.add_atom('C', hybridization=(1, 4))
        q.add_atom('C', hybridization=(1, 2))
        q.add_atom('C')
        q.add_atom('C')
        q.add_bond(1, 2, 2)
        q.add_bond(1, 4, 1)
        q.add_bond(3, 5, 1)
        q.add_bond(4, 6, (1, 2))
        q.add_bond(5, 7, (1, 4))
        q.add_bond(6, 7, (1, 2))
        q.add_bond(1, 8, 1)
        rules.append(q)

        # C1([O,S,N;H])[O,S,N;!H]-[C;X4,a]:,-C-[S,O,N]-C1
        q = query.QueryContainer()
        q.add_atom('C')  # entry
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=(1, 2), hybridization=2)
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=1)
        q.add_atom('C')
        q.add_atom('C', hybridization=(1, 4))
        q.add_atom(ListElement(['S', 'O', 'N']))
        q.add_atom('C')
        q.add_atom('C')
        q.add_bond(1, 2, 2)
        q.add_bond(1, 4, 1)
        q.add_bond(3, 5, 1)
        q.add_bond(4, 6, 1)
        q.add_bond(5, 7, (1, 4))
        q.add_bond(6, 7, 1)
        q.add_bond(1, 8, 1)
        rules.append(q)

        return rules

    @class_cached_property
    def __sugar_group_rules(self):
        rules = []

        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'N']), hydrogens=(1, 2))  # enol
        q.add_atom(ListElement(['O', 'N']))  # ketone
        q.add_atom('C', hybridization=1)
        q.add_atom('C', hybridization=2)
        q.add_bond(1, 3, 1)
        q.add_bond(3, 4, 1)
        q.add_bond(2, 4, 2)
        rules.append(q)

        return rules


__all__ = ['Tautomers']
