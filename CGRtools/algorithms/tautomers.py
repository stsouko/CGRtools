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
from collections import deque, defaultdict
from itertools import product, chain, repeat
from typing import TYPE_CHECKING, Iterator, Union
from ..containers import query  # cyclic imports resolve
from ..exceptions import InvalidAromaticRing
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
            # keto-enol rule-based. not implemented
            # smiles alphabetically sorted (descent). bad rule!
            return len(m.aromatic_rings), -sum(x != 0 for x in m._charges.values()), str(m)

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

    def enumerate_tautomers(self: 'MoleculeContainer', *, prepare_molecules=True, full=True, zwitter=True,
                            keto_enol=True, hetero=True, limit: int = 1000) -> Iterator['MoleculeContainer']:
        """
        Enumerate all possible tautomeric forms of molecule.

        :param prepare_molecules: Standardize structures for correct processing.
        :param full: Do full enumeration.
        :param zwitter: Enable acid-base tautomerization
        :param keto_enol: Enable keto-enol tautomerization
        :param hetero: Enable hetero arenes tautomerization.
        :param limit: Maximum amount of generated structures.
        """
        if limit < 2:
            raise ValueError('limit should be greater or equal 2')

        atoms_stereo = self._atoms_stereo
        allenes_stereo = self._allenes_stereo
        cis_trans_stereo = self._cis_trans_stereo
        has_stereo = bool(atoms_stereo or allenes_stereo or cis_trans_stereo)

        copy = self.copy()
        copy.clean_stereo()
        if prepare_molecules:
            copy.kekule()
            copy.implicify_hydrogens()

        ster = copy.copy()
        ster.thiele()
        if has_stereo:
            ster._atoms_stereo.update(atoms_stereo)
            ster._allenes_stereo.update(allenes_stereo)
            ster._cis_trans_stereo.update(cis_trans_stereo)
            ster._fix_stereo()
        yield ster

        seen = {copy: None}
        queue = deque([copy])
        counter = 1
        while queue:
            current = queue.popleft()
            if keto_enol:
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
                        counter += 1
                        if counter == limit:
                            return
                        if len(mol.aromatic_rings) > len(current.aromatic_rings):
                            # found new aromatic ring. flush queue and start from it.
                            queue.clear()
                            queue.append(mol)
                            break
                        queue.append(mol)

            if zwitter:
                for mol in current._enumerate_zwitter_tautomers(full):
                    if mol not in seen:
                        seen[mol] = current
                        queue.append(mol)
                        if has_stereo:
                            mol = mol.copy()  # prevent seen hashtable breaking
                            mol._atoms_stereo.update(atoms_stereo)
                            mol._allenes_stereo.update(allenes_stereo)
                            mol._cis_trans_stereo.update(cis_trans_stereo)
                            mol._fix_stereo()
                        yield mol
                        counter += 1
                        if counter == limit:
                            return

            if hetero:
                for mol in current._enumerate_hetero_arene_tautomers():
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

    def _enumerate_keto_enol_tautomers(self: Union['MoleculeContainer', 'Tautomers']):
        if '__cached_args_method_neighbors' in self.__dict__:
            neighbors = self.__dict__['__cached_args_method_neighbors']
        else:
            neighbors = self.__dict__['__cached_args_method_neighbors'] = {}
        if '__cached_args_method_heteroatoms' in self.__dict__:
            heteroatoms = self.__dict__['__cached_args_method_heteroatoms']
        else:
            heteroatoms = self.__dict__['__cached_args_method_heteroatoms'] = {}

        for fix, ket in self.__enumerate_bonds():
            if ket:
                a = fix[-1][1]
                d = fix[0][0]
            else:
                a = fix[0][0]
                d = fix[-1][1]

            mol = self.copy()
            m_bonds = mol._bonds
            for n, m, b in fix:
                m_bonds[n][m]._Bond__order = b

            mol._hydrogens[a] += 1
            mol._hydrogens[d] -= 1
            mol._hybridizations[a] = 1
            mol._hybridizations[d] = 2

            # store cache in new molecules for speedup
            mol.__dict__['__cached_args_method_neighbors'] = neighbors.copy()
            mol.__dict__['__cached_args_method_heteroatoms'] = heteroatoms.copy()
            if mol.thiele(fix_tautomers=False, fix_metal_organics=False):
                mol.__dict__['__cached_args_method_neighbors'] = neighbors.copy()
                mol.__dict__['__cached_args_method_heteroatoms'] = heteroatoms.copy()
            yield mol, ket

    def _enumerate_zwitter_tautomers(self: Union['MoleculeContainer', 'Tautomers'], full=True):
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
            mol._hydrogens[d] -= 1
            mol._hydrogens[a] += 1
            mol._charges[d] -= 1
            mol._charges[a] += 1

            # store cache in new molecules for speedup
            mol.__dict__['__cached_args_method_neighbors'] = neighbors.copy()
            mol.__dict__['__cached_args_method_heteroatoms'] = heteroatoms.copy()
            yield mol

    def _enumerate_hetero_arene_tautomers(self: Union['MoleculeContainer', 'Tautomers']):
        if '__cached_args_method_neighbors' in self.__dict__:
            neighbors = self.__dict__['__cached_args_method_neighbors']
        else:
            neighbors = self.__dict__['__cached_args_method_neighbors'] = {}
        if '__cached_args_method_heteroatoms' in self.__dict__:
            heteroatoms = self.__dict__['__cached_args_method_heteroatoms']
        else:
            heteroatoms = self.__dict__['__cached_args_method_heteroatoms'] = {}

        atoms = self._atoms
        bonds = self._bonds
        hydrogens = self._hydrogens
        charges = self._charges
        radicals = self._radicals

        rings = defaultdict(list)  # aromatic skeleton
        for n, m_bond in bonds.items():
            for m, bond in m_bond.items():
                if bond.order == 4:
                    rings[n].append(m)
        if not rings:
            return

        acceptors = set()
        donors = set()
        single_bonded = set()
        for n, ms in rings.items():
            if len(ms) == 2:
                if atoms[n].atomic_number in (5, 7, 15) and not charges[n] and not radicals[n]:
                    # only neutral B, N, P
                    if hydrogens[n]:  # pyrrole
                        donors.add(n)
                    elif len(bonds[n]) == 2:  # pyridine
                        acceptors.add(n)
                    else:
                        single_bonded.add(n)
                elif charges[n] == -1 and atoms[n].atomic_number == 6:  # ferrocene
                    single_bonded.add(n)
            elif len(ms) == 3 and atoms[n].atomic_number in (5, 7, 15) and not charges[n] and not radicals[n]:
                single_bonded.add(n)
        if not donors or not acceptors:
            return

        atoms = set(rings)
        components = []
        while atoms:
            start = atoms.pop()
            component = {start: rings[start]}
            queue = deque([start])
            while queue:
                current = queue.popleft()
                for n in rings[current]:
                    if n not in component:
                        queue.append(n)
                        component[n] = rings[n]

            atoms.difference_update(component)
            if donors.isdisjoint(component) or acceptors.isdisjoint(component):
                continue
            components.append(component)

        if not components:
            return
        for component in components:
            for d, a in product(component.keys() & donors, component.keys() & acceptors):
                sb = component.keys() & single_bonded
                sb.add(a)  # now pyrrole
                try:
                    next(self._kekule_component(component, sb, ()))
                except InvalidAromaticRing:
                    continue
                mol = self.copy()
                mol._hydrogens[d] = 0
                mol._hydrogens[a] = 1

                # store cache in new molecules for speedup
                mol.__dict__['__cached_args_method_neighbors'] = neighbors.copy()
                mol.__dict__['__cached_args_method_heteroatoms'] = heteroatoms.copy()
                yield mol

    @cached_property
    def _sugar_groups(self: Union['MoleculeContainer', 'Tautomers']):
        ek = []
        for mapping in self.__sugar_group_rule.get_mapping(self, automorphism_filter=False):
            e, k = mapping[1], mapping[2]
            ek.append((e, k))
        return ek

    def __enumerate_bonds(self: Union['MoleculeContainer', 'Tautomers']):
        atoms = self._atoms
        bonds = self._bonds
        charge = self._charges
        hydrogens = self._hydrogens
        hybridizations = self._hybridizations
        heteroatoms = self.heteroatoms
        neighbors = self.neighbors

        # search neutral oxygen and nitrogen
        donors = []
        acceptors = []
        for n, a in atoms.items():
            if charge[n] or not neighbors(n) or heteroatoms(n):
                continue
            if a.atomic_number in (7, 8):
                if hybridizations[n] == 2:  # imine!
                    acceptors.append(n)
                elif hydrogens[n]:  # amine!
                    donors.append(n)

        for atom, hydrogen in chain(zip(donors, repeat(True)), zip(acceptors, repeat(False))):
            path = []
            seen = {atom}
            # stack is neighbors
            if hydrogen:  # enol
                stack = [(atom, n, 2, 0) for n in bonds[atom] if hybridizations[n] == 2]
            else:  # ketone
                stack = [(atom, n, 1, 0) for n, b in bonds[atom].items() if hybridizations[n] == 2 and b.order == 2]

            while stack:
                last, current, bond, depth = stack.pop()

                if len(path) > depth:
                    if len(path) % 2 == 0:
                        yield path, not hydrogen
                    seen.difference_update(x for _, x, _ in path[depth:])
                    path = path[:depth]

                path.append((last, current, bond))

                # adding neighbors
                depth += 1
                seen.add(current)
                if bond == 2:
                    next_bond = 1
                else:
                    next_bond = 2

                for n, b in bonds[current].items():
                    if n == last:
                        continue
                    elif n in seen:  # aromatic ring destruction. pyridine double bonds shift
                        if not stack:
                            path = [None]
                    elif b.order == bond and atoms[n].atomic_number == 6:
                        hb = hybridizations[n]
                        if hb == 2:  # grow up
                            stack.append((current, n, next_bond, depth))
                        elif hydrogen:
                            if hb == 3:  # OC=CC=C=A case
                                cp = path.copy()
                                cp.append((current, n, 1))
                                yield cp, True  # ketone found
                        elif hb == 1 and hydrogens[n]:  # ketone >> enole
                            cp = path.copy()
                            cp.append((current, n, 2))
                            yield cp, False

            if len(path) % 2 == 0:  # time to yield
                yield path, not hydrogen

    @class_cached_property
    def __sugar_group_rule(self):
        q = query.QueryContainer()
        q.add_atom(ListElement(['O', 'N']), hydrogens=(1, 2))  # enol
        q.add_atom(ListElement(['O', 'N']))  # ketone
        q.add_atom('C', hybridization=1, heteroatoms=1)
        q.add_atom('C', hybridization=2, heteroatoms=1)
        q.add_bond(1, 3, 1)
        q.add_bond(3, 4, 1)
        q.add_bond(2, 4, 2)
        return q

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


__all__ = ['Tautomers']
