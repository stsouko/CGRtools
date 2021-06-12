# -*- coding: utf-8 -*-
#
#  Copyright 2020 Nail Samikaev <samikaevn@yandex.ru>
#  Copyright 2020, 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import TYPE_CHECKING, Iterator, Union
from ..containers import query  # cyclic imports resolve
from ..periodictable import ListElement


if TYPE_CHECKING:
    from CGRtools import MoleculeContainer


class Tautomers:
    __slots__ = ()

    def tautomerize(self: 'MoleculeContainer', *, prepare_molecules=True,
                    zwitter=True, keto_enol=True, limit: int = 1000) -> bool:
        """
        Convert structure to canonical tautomeric form. Return True if structure changed.
        """
        def key(m):  # inspired from https://github.com/mcs07/MolVS/
            # more aromatic rings is better
            # less charged atoms is better
            # keto-enol rule-based
            # smiles alphabetically sorted (descent)
            return len(m.aromatic_rings) * 100 - sum(x != 0 for x in m._charges.values()) * 5, str(m)

        canon = max(self.enumerate_tautomers(prepare_molecules=prepare_molecules, full=False, zwitter=zwitter,
                                             keto_enol=keto_enol, limit=limit), key=key)
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

    def enumerate_tautomers(self: 'MoleculeContainer', *, prepare_molecules=True, full=True,
                            zwitter=True, keto_enol=True, limit: int = 1000) -> Iterator['MoleculeContainer']:
        """
        Enumerate all possible tautomeric forms of molecule.

        :param prepare_molecules: Standardize structures before. Aromatization and implicit hydrogens required.
        :param full: Do full enumeration.
        :param zwitter: Enable acid-base tautomerization
        :param keto_enol: Enable keto-enol tautomerization
        :param limit: Maximum amount of generated structures.
        """
        if limit < 2:
            raise ValueError('limit should be greater or equal 2')
        yield self.copy()
        atoms_stereo = self._atoms_stereo
        allenes_stereo = self._allenes_stereo
        cis_trans_stereo = self._cis_trans_stereo
        has_stereo = bool(atoms_stereo or allenes_stereo or cis_trans_stereo)

        copy = self.copy()
        copy.clean_stereo()
        if prepare_molecules:
            copy.kekule()
            copy.implicify_hydrogens()
            copy.thiele()  # prevent

        seen = {copy: None}
        queue = deque([copy])
        counter = 1
        while queue:
            current = queue.popleft()

            if zwitter:
                for mol in current._enumerate_zwitter_tautomers(full):
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
                        counter += 1
                        if counter == limit:
                            return

            if keto_enol:
                for mol, ket in current._enumerate_keto_enol_tautomers(full):
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
                        counter += 1
                        if counter == limit:
                            return
                        if len(mol.aromatic_rings) > len(current.aromatic_rings):
                            # found new aromatic ring. flush queue and start from it.
                            queue.clear()
                            queue.append(mol)
                            break
                        queue.append(mol)

    def _enumerate_zwitter_tautomers(self: Union['MoleculeContainer', 'Tautomers'], full=True):
        sssr = self.sssr
        rings_count = self.rings_count
        atoms_rings = self.atoms_rings
        atoms_rings_sizes = self.atoms_rings_sizes
        if '__cached_args_method_neighbors' in self.__dict__:
            neighbors = self.__dict__['__cached_args_method_neighbors']
        else:
            neighbors = self.__dict__['__cached_args_method_neighbors'] = {}
        if '__cached_args_method_heteroatoms' in self.__dict__:
            heteroatoms = self.__dict__['__cached_args_method_heteroatoms']
        else:
            heteroatoms = self.__dict__['__cached_args_method_heteroatoms'] = {}

        donors = set()
        acceptors = set()
        for q, acid in chain(zip(self.__acid_rules if full else self.__stripped_acid_rules, repeat(True)),
                             zip(self.__base_rules if full else self.__stripped_base_rules, repeat(False))):
            for mapping in q.get_mapping(self, automorphism_filter=False):
                n = mapping[1]
                if acid:
                    donors.add(n)
                else:
                    acceptors.add(n)

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
            mol.__dict__['__cached_args_method_neighbors'] = neighbors.copy()
            mol.__dict__['__cached_args_method_heteroatoms'] = heteroatoms.copy()
            yield mol

    def _enumerate_keto_enol_tautomers(self: Union['MoleculeContainer', 'Tautomers'], full=True):
        sssr = self.sssr
        rings_count = self.rings_count
        atoms_rings = self.atoms_rings
        atoms_rings_sizes = self.atoms_rings_sizes
        if '__cached_args_method_neighbors' in self.__dict__:
            neighbors = self.__dict__['__cached_args_method_neighbors']
        else:
            neighbors = self.__dict__['__cached_args_method_neighbors'] = {}
        if '__cached_args_method_heteroatoms' in self.__dict__:
            heteroatoms = self.__dict__['__cached_args_method_heteroatoms']
        else:
            heteroatoms = self.__dict__['__cached_args_method_heteroatoms'] = {}

        seen = set()
        for q, ket, *fix in (self.__keto_enol_rules if full else self.__stripped_keto_enol_rules):
            for mapping in q.get_mapping(self, automorphism_filter=False):
                a = mapping[1]
                d = mapping[2]
                if not ket:  # prevent 1-3 migration in 1-5 etc
                    dv = (d, mapping[3])
                    if dv in seen:
                        continue
                    seen.add(dv)

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
                mol.__dict__['__cached_args_method_heteroatoms'] = heteroatoms.copy()
                k = mol.kekule()
                t = mol.thiele()
                # restore after kekule-thiele
                if k or t:
                    mol.__dict__['sssr'] = sssr
                    mol.__dict__['__cached_args_method_neighbors'] = neighbors.copy()
                    mol.__dict__['__cached_args_method_heteroatoms'] = heteroatoms.copy()
                mol.__dict__['rings_count'] = rings_count
                mol.__dict__['atoms_rings'] = atoms_rings
                mol.__dict__['atoms_rings_sizes'] = atoms_rings_sizes
                yield mol, ket

    @cached_property
    def _sugar_groups(self: Union['MoleculeContainer', 'Tautomers']):
        ek = []
        for q in self.__sugar_group_rules:
            for mapping in q.get_mapping(self, automorphism_filter=False):
                e, k = mapping[1], mapping[2]
                ek.append((e, k))
        return ek

    @class_cached_property
    def __stripped_keto_enol_rules(self):
        rules = []  # query, is ketone, bond fix
        # enoles are order-dependent

        # 2-pyridone. [O,S:1]=[C:3]1[N;H:2][C:4]=,:[C,N:5]-,:[C,N:6]=,:[C:7]1
        q = query.QueryContainer()  # first 3 atoms have special role!
        q.add_atom(ListElement(['O', 'S']), neighbors=1)  # 1st atom only acceptor
        q.add_atom('N', hydrogens=1)  # 2nd atom only donor
        q.add_atom('C')  # 3rd atom connected to donor
        q.add_atom('C')
        q.add_atom(ListElement(['C', 'N']))
        q.add_atom(ListElement(['C', 'N']))
        q.add_atom('C')
        q.add_bond(1, 3, 2)
        q.add_bond(2, 3, 1)
        q.add_bond(2, 4, 1)
        q.add_bond(3, 7, 1)
        q.add_bond(4, 5, (2, 4))
        q.add_bond(5, 6, (1, 4))
        q.add_bond(6, 7, (2, 4))
        rules.append((q, False, (1, 3, 1), (2, 3, 2)))

        # [C;a:10][N;H:2][N:3]=[C:4]1[C:5]=,:[C:6][C:7](=[O:1])[C:8]=,:[C:9]1
        q = query.QueryContainer()
        q.add_atom('O')
        q.add_atom('N', hydrogens=1, heteroatoms=1)
        q.add_atom('N')
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C', hybridization=4)
        q.add_bond(2, 3, 1)
        q.add_bond(3, 4, 2)
        q.add_bond(4, 5, 1)
        q.add_bond(5, 6, (2, 4))
        q.add_bond(6, 7, 1)
        q.add_bond(1, 7, 2)
        q.add_bond(7, 8, 1)
        q.add_bond(8, 9, (2, 4))
        q.add_bond(4, 9, 1)
        q.add_bond(2, 10, 1)
        rules.append((q, False, (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 6, 1), (6, 7, 2), (1, 7, 1)))

        # [C;a:10][N;H:2][N:3]=[C:4]1[C:5](=[O:1])[C:6]=,:[C:7]-,:[C:8]=,:[C:9]1
        q = query.QueryContainer()
        q.add_atom('O')
        q.add_atom('N', hydrogens=1, heteroatoms=1)
        q.add_atom('N')
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C', hybridization=4)
        q.add_bond(2, 3, 1)
        q.add_bond(3, 4, 2)
        q.add_bond(4, 5, 1)
        q.add_bond(5, 6, 1)
        q.add_bond(6, 7, (2, 4))
        q.add_bond(7, 8, (1, 4))
        q.add_bond(8, 9, (2, 4))
        q.add_bond(4, 9, 1)
        q.add_bond(1, 5, 2)
        q.add_bond(2, 10, 1)
        rules.append((q, False, (2, 3, 2), (3, 4, 1), (4, 5, 2), (1, 5, 1)))

        # [O,S,NH:1]=[C:7]([C,H])[C:6]=[C:5][C:4]=[C:3]([C,H])-[O,S,N;H:2]
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=1)  # acceptor
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=(1, 2), heteroatoms=0)  # donor
        q.add_atom('C', heteroatoms=1)
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C', heteroatoms=1)
        q.add_bond(2, 3, 1)
        q.add_bond(3, 4, 2)
        q.add_bond(4, 5, 1)
        q.add_bond(5, 6, 2)
        q.add_bond(6, 7, 1)
        q.add_bond(1, 7, 2)
        rules.append((q, False, (1, 7, 1), (2, 3, 2), (3, 4, 1), (4, 5, 2), (5, 6, 1), (6, 7, 2)))

        # todo: fix rules!!!

        # [C;X4,a:8][N:1]=[C:3]([C,H])[C:4]=[C:5][C:6]=[C:7]([C,H])-[O,S,N;H:2]
        q = query.QueryContainer()
        q.add_atom('N')  # acceptor
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=(1, 2), heteroatoms=0)  # donor
        q.add_atom('C', heteroatoms=1)
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C', heteroatoms=1)
        q.add_atom('C', hybridization=(1, 4))
        q.add_bond(1, 3, 2)
        q.add_bond(2, 7, 1)
        q.add_bond(3, 4, 1)
        q.add_bond(4, 5, 2)
        q.add_bond(5, 6, 1)
        q.add_bond(6, 7, 2)
        q.add_bond(1, 8, 1)
        rules.append((q, False, (1, 3, 1), (2, 7, 2), (3, 4, 2), (4, 5, 1), (5, 6, 2), (6, 7, 1)))

        # [O,S,NH:1]=[C:3]([C,H])[C:4]=[C:5]([C,H])-[O,S,N;H:2]
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=1)  # acceptor
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=(1, 2), heteroatoms=0)  # donor
        q.add_atom('C', heteroatoms=1)
        q.add_atom('C')
        q.add_atom('C', heteroatoms=1)
        q.add_bond(1, 3, 2)
        q.add_bond(2, 5, 1)
        q.add_bond(3, 4, 1)
        q.add_bond(4, 5, 2)
        rules.append((q, False, (1, 3, 1), (2, 5, 2), (3, 4, 2), (4, 5, 1)))

        # [C;X4,a:6][N:1]=[C:3]([C,H])[C:4]=[C:5]([C,H])-[O,S,N;H:2]
        q = query.QueryContainer()
        q.add_atom('N')  # acceptor
        q.add_atom(ListElement(['O', 'S', 'N']), hydrogens=(1, 2), heteroatoms=0)  # donor
        q.add_atom('C', heteroatoms=1)
        q.add_atom('C')
        q.add_atom('C', heteroatoms=1)
        q.add_atom('C', hybridization=(1, 4))
        q.add_bond(1, 3, 2)
        q.add_bond(2, 5, 1)
        q.add_bond(3, 4, 1)
        q.add_bond(4, 5, 2)
        q.add_bond(1, 6, 1)
        rules.append((q, False, (1, 3, 1), (2, 5, 2), (3, 4, 2), (4, 5, 1)))

        # [C;a:4]-[N:1]=[N:3]-[C;H:2]-[a:5]
        q = query.QueryContainer()
        q.add_atom('N')  # acceptor
        q.add_atom('C', hydrogens=1, hybridization=1)
        q.add_atom('N')
        q.add_atom('C', hybridization=4)
        q.add_atom('A', hybridization=4)
        q.add_bond(1, 3, 2)
        q.add_bond(1, 4, 1)
        q.add_bond(2, 3, 1)
        q.add_bond(2, 5, 1)
        rules.append((q, False, (1, 3, 1), (2, 3, 2)))

        # C=C([C,H])-[O,S,N;H]
        q = query.QueryContainer()
        q.add_atom('C')  # acceptor
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=1)  # donor
        q.add_atom('C', heteroatoms=1)
        q.add_bond(1, 3, 2)
        q.add_bond(2, 3, 1)
        rules.append((q, False, (1, 3, 1), (2, 3, 2)))

        # C=C([C,H])-[N;H]-N([C,H])[C,H]
        q = query.QueryContainer()
        q.add_atom('C')  # acceptor
        q.add_atom('N', hydrogens=1)  # donor
        q.add_atom('C', heteroatoms=1)
        q.add_atom('N', heteroatoms=1, hybridization=1)
        q.add_bond(1, 3, 2)
        q.add_bond(2, 3, 1)
        q.add_bond(2, 4, 1)
        rules.append((q, False, (1, 3, 1), (2, 3, 2)))

        # [N:1]=[C:3]-[N;H:2]
        q = query.QueryContainer()
        q.add_atom('N')  # acceptor
        q.add_atom('N', hydrogens=(1, 2))  # donor
        q.add_atom('C')
        q.add_bond(1, 3, 2)
        q.add_bond(2, 3, 1)
        rules.append((q, False, (1, 3, 1), (2, 3, 2)))

        # 2H‐pyrrole. [N:1]=1[C:2][C:3]=,:[C:4][C:5]=1
        q = query.QueryContainer()
        q.add_atom('N')  # acceptor
        q.add_atom('C', hydrogens=(1, 2))  # donor
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C')
        q.add_bond(1, 2, 1)
        q.add_bond(2, 3, 1)
        q.add_bond(3, 4, (2, 4))
        q.add_bond(4, 5, 1)
        q.add_bond(1, 5, 2)
        rules.append((q, False, (1, 5, 1)))  # thiele ad-hoc used!

        # [C:1]=[C:3]1[N;H:2][C:4]=,:[C,N:5]-,:[C,N:6]=,:[C:7]1
        q = query.QueryContainer()
        q.add_atom('C', hybridization=2)
        q.add_atom('N', hydrogens=1)
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom(ListElement(['C', 'N']))
        q.add_atom(ListElement(['C', 'N']))
        q.add_atom('C')
        q.add_bond(1, 3, 2)
        q.add_bond(2, 3, 1)
        q.add_bond(2, 4, 1)
        q.add_bond(3, 7, 1)
        q.add_bond(4, 5, (2, 4))
        q.add_bond(5, 6, (1, 4))
        q.add_bond(6, 7, (2, 4))
        rules.append((q, False, (1, 3, 1), (2, 3, 2)))

        # [O:1]=[C:3]1[C:4]:[C:5][C:2][C:6]:[C:7]1
        q = query.QueryContainer()
        q.add_atom('O')
        q.add_atom('C', hydrogens=(1, 2))
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C')
        q.add_atom('C')
        q.add_bond(1, 3, 2)
        q.add_bond(3, 4, 1)
        q.add_bond(4, 5, 4)
        q.add_bond(5, 2, 1)
        q.add_bond(2, 6, 1)
        q.add_bond(6, 7, 4)
        q.add_bond(3, 7, 1)
        rules.append((q, False, (1, 3, 1)))  # trick based on thiele search of aromatic rings

        return rules

    @class_cached_property
    def __keto_enol_rules(self):
        rules = list(self.__stripped_keto_enol_rules)  # query, is ketone, bond fix

        # first atom is acceptor of proton
        # second is donor
        # 1,3: [C;X4,H:2]-[C:3]([C,H])=[O,S,N;X1:1]
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=1)  # acceptor
        q.add_atom('C', hybridization=1, hydrogens=(1, 2, 3))  # donor
        q.add_atom('C', heteroatoms=1)
        q.add_bond(1, 3, 2)
        q.add_bond(2, 3, 1)
        rules.append((q, True, (1, 3, 1), (2, 3, 2)))

        # 1,3: [C;X4,H]-C([C,H])=N-[C;X4]
        q = query.QueryContainer()
        q.add_atom('N')  # acceptor
        q.add_atom('C', hybridization=1, hydrogens=(1, 2, 3))  # donor
        q.add_atom('C', heteroatoms=1)
        q.add_atom('C', hybridization=1, heteroatoms=1)
        q.add_bond(1, 3, 2)
        q.add_bond(2, 3, 1)
        q.add_bond(1, 4, 1)
        rules.append((q, True, (1, 3, 1), (2, 3, 2)))

        # 1,3: [C;X4,H]-C([C,H])=N-[C,N;a]
        q = query.QueryContainer()
        q.add_atom('N')  # acceptor
        q.add_atom('C', hybridization=1, hydrogens=(1, 2, 3))  # donor
        q.add_atom('C', heteroatoms=1)
        q.add_atom(ListElement(['C', 'N']), hybridization=4)
        q.add_bond(1, 3, 2)
        q.add_bond(2, 3, 1)
        q.add_bond(1, 4, 1)
        rules.append((q, True, (1, 3, 1), (2, 3, 2)))

        # 1,5: [C;X4,H:2]-[C:3]([C,H])=[C:5]([C,H])-[C:4]([C,H])=[O,S,N;X1:1]
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'N']), neighbors=1)  # acceptor
        q.add_atom('C', hybridization=1, hydrogens=(1, 2, 3))  # donor
        q.add_atom('C', heteroatoms=0)
        q.add_atom('C', heteroatoms=1)
        q.add_atom('C', heteroatoms=0)
        q.add_bond(1, 4, 2)
        q.add_bond(2, 3, 1)
        q.add_bond(3, 5, 2)
        q.add_bond(4, 5, 1)
        rules.append((q, True, (1, 4, 1), (2, 3, 2), (3, 5, 1), (4, 5, 2)))

        return rules

    @class_cached_property
    def __stripped_acid_rules(self):
        rules = []

        # Ammonia, H-imine,guanidine,amidine. [H][N+]=,-,:
        q = query.QueryContainer()
        q.add_atom('N', charge=1, hydrogens=(1, 2, 3, 4))
        rules.append(q)

        return rules

    @class_cached_property
    def __acid_rules(self):
        rules = list(self.__stripped_acid_rules)

        # Phenoles, [H][O,S,Se]-Ar
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'Se']), neighbors=1)
        q.add_atom(ListElement(['C', 'N']), hybridization=4)
        q.add_bond(1, 2, 1)
        rules.append(q)

        # Oxo-acids. [H][O,S,Se]-[C,N,P,S,Cl,Se,Br,I]=O
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'Se']), neighbors=1)
        q.add_atom(ListElement(['C', 'N', 'P', 'S', 'Cl', 'Se', 'Br', 'I']))
        q.add_atom('O')
        q.add_bond(1, 2, 1)
        q.add_bond(2, 3, 2)
        rules.append(q)

        # Nitro acid. [H]O-[N+](=O)[O-]
        q = query.QueryContainer()
        q.add_atom('O', neighbors=1)
        q.add_atom('N', charge=1)
        q.add_atom('O', charge=-1)
        q.add_atom('O')
        q.add_bond(1, 2, 1)
        q.add_bond(2, 3, 1)
        q.add_bond(2, 4, 2)
        rules.append(q)

        # Halogen acids
        q = query.QueryContainer()
        q.add_atom(ListElement(['F', 'Cl', 'Br', 'I']), neighbors=0)
        rules.append(q)
        return rules

    @class_cached_property
    def __stripped_base_rules(self):
        rules = []

        # Oxo-acid salts. [O,S,Se-]-[C,N,P,S,Cl,Se,Br,I]=O
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'Se']), charge=-1)
        q.add_atom(ListElement(['C', 'N', 'P', 'S', 'Cl', 'Se', 'Br', 'I']))
        q.add_atom('O')
        q.add_bond(1, 2, 1)
        q.add_bond(2, 3, 2)
        rules.append(q)

        # Phenole salts. [O,S,Se-]-Ar
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'S', 'Se']), charge=-1)
        q.add_atom(ListElement(['C', 'N']), hybridization=4)
        q.add_bond(1, 2, 1)
        rules.append(q)

        # Nitrate. [O-]-[N+](=O)[O-]
        q = query.QueryContainer()
        q.add_atom('O', charge=-1)
        q.add_atom('N', charge=1)
        q.add_atom('O', charge=-1)
        q.add_atom('O')
        q.add_bond(1, 2, 1)
        q.add_bond(2, 3, 1)
        q.add_bond(2, 4, 2)
        rules.append(q)

        # ions
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'F', 'Cl', 'Br', 'I']), charge=-1, neighbors=0)
        rules.append(q)

        return rules

    @class_cached_property
    def __base_rules(self):
        rules = list(self.__stripped_base_rules)

        # Guanidine
        q = query.QueryContainer()
        q.add_atom('N', heteroatoms=0)
        q.add_atom('C')
        q.add_atom('N')
        q.add_atom('N')
        q.add_bond(1, 2, 2)
        q.add_bond(2, 3, 1)
        q.add_bond(2, 4, 1)
        rules.append(q)

        # Oxo-guanidine, Amino-guanidine
        q = query.QueryContainer()
        q.add_atom('N')
        q.add_atom('C')
        q.add_atom('N')
        q.add_atom('N')
        q.add_atom(ListElement(['O', 'N']))
        q.add_bond(1, 2, 2)
        q.add_bond(2, 3, 1)
        q.add_bond(2, 4, 1)
        q.add_bond(1, 5, 1)
        rules.append(q)

        # O-alkyl-isourea, S-alkyl-isothiaurea
        q = query.QueryContainer()
        q.add_atom('N', heteroatoms=0)
        q.add_atom('C')
        q.add_atom('N')
        q.add_atom(ListElement(['O', 'S']), neighbors=2, hybridization=1, heteroatoms=0)
        q.add_bond(1, 2, 2)
        q.add_bond(2, 3, 1)
        q.add_bond(2, 4, 1)
        rules.append(q)

        # Dialkyl imidocarbonate
        q = query.QueryContainer()
        q.add_atom('N', heteroatoms=0)
        q.add_atom('C')
        q.add_atom('O', neighbors=2, heteroatoms=0)
        q.add_atom('O', neighbors=2, heteroatoms=0)
        q.add_bond(1, 2, 2)
        q.add_bond(2, 3, 1)
        q.add_bond(2, 4, 1)
        rules.append(q)

        # Amidine
        q = query.QueryContainer()
        q.add_atom('N', heteroatoms=0)
        q.add_atom('C', heteroatoms=2)
        q.add_atom('N')
        q.add_bond(1, 2, 2)
        q.add_bond(2, 3, 1)
        rules.append(q)

        # O-alkyl-imidate (oxazoline)
        q = query.QueryContainer()
        q.add_atom('N', heteroatoms=0)
        q.add_atom('C', heteroatoms=2)
        q.add_atom('O', neighbors=2, heteroatoms=0)
        q.add_bond(1, 2, 2)
        q.add_bond(2, 3, 1)
        rules.append(q)

        # Amidoxime. O-N=C([C,H])N
        q = query.QueryContainer()
        q.add_atom('N')
        q.add_atom('C', heteroatoms=2)
        q.add_atom('N')
        q.add_atom('O')
        q.add_bond(1, 2, 2)
        q.add_bond(2, 3, 1)
        q.add_bond(1, 4, 1)
        rules.append(q)

        # Oxime, Hydrazone. [O,N]-N=C([C,H])[C,H]
        q = query.QueryContainer()
        q.add_atom('N')
        q.add_atom('C', heteroatoms=1, hybridization=2)
        q.add_atom(ListElement(['N', 'O']))
        q.add_bond(1, 2, 2)
        q.add_bond(1, 3, 1)
        rules.append(q)

        # Imine. [C,H]N=C([C,H])[C,H]
        q = query.QueryContainer()
        q.add_atom('N', heteroatoms=0)
        q.add_atom('C', heteroatoms=1, hybridization=2)
        q.add_bond(1, 2, 2)
        rules.append(q)

        # Alkyl amine, Hydroxylamine, Hydrazine
        q = query.QueryContainer()
        q.add_atom('N', neighbors=1)
        q.add_atom(ListElement(['C', 'O', 'N']), hybridization=1, heteroatoms=1)
        q.add_bond(1, 2, 1)
        rules.append(q)

        # Dialkyl amine, Alkyl hydroxylamine, Alkyl hydrazine
        q = query.QueryContainer()
        q.add_atom('N', neighbors=2)
        q.add_atom(ListElement(['C', 'O', 'N']), hybridization=1, heteroatoms=1)
        q.add_atom('C', hybridization=1, heteroatoms=1)
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        rules.append(q)

        # Trialkyl amine, Dialkyl-hydroxylamine, Dialkyl-hydrazine
        q = query.QueryContainer()
        q.add_atom('N')
        q.add_atom(ListElement(['C', 'O', 'N']), hybridization=1, heteroatoms=1)
        q.add_atom('C', hybridization=1, heteroatoms=1)
        q.add_atom('C', hybridization=1, heteroatoms=1)
        q.add_bond(1, 2, 1)
        q.add_bond(1, 3, 1)
        q.add_bond(1, 4, 1)
        rules.append(q)

        # Pyridine. Imidazole. Triazole. :N:
        q = query.QueryContainer()
        q.add_atom('N', hybridization=4, hydrogens=0, neighbors=2)
        rules.append(q)

        return rules

    @class_cached_property
    def __sugar_group_rules(self):
        rules = []

        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'N']), hydrogens=(1, 2))  # enol
        q.add_atom(ListElement(['O', 'N']))  # ketone
        q.add_atom('C', hybridization=1, heteroatoms=1)
        q.add_atom('C', hybridization=2, heteroatoms=1)
        q.add_bond(1, 3, 1)
        q.add_bond(3, 4, 1)
        q.add_bond(2, 4, 2)
        rules.append(q)

        return rules


__all__ = ['Tautomers']
