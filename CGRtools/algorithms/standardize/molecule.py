# -*- coding: utf-8 -*-
#
#  Copyright 2018-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2018 Tagir Akhmetshin <tagirshin@gmail.com>
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
from CachedMethods import class_cached_property
from collections import defaultdict
from typing import List, TYPE_CHECKING, Union, Tuple
from ...containers import query  # cyclic imports resolve
from ...containers.bonds import Bond
from ...exceptions import ValenceError
from ...periodictable import ListElement, H


if TYPE_CHECKING:
    from CGRtools import MoleculeContainer


class Standardize:
    __slots__ = ()

    def canonicalize(self: 'MoleculeContainer', *, logging=False) -> \
            Union[bool, List[Tuple[Tuple[int, ...], int, str]]]:
        """
        Convert molecule to canonical forms of functional groups and aromatic rings without explicit hydrogens.

        :param logging: return log.
        """
        k = self.kekule()
        s = self.standardize(fix_stereo=False, logging=logging)
        h = self.implicify_hydrogens(fix_stereo=False)
        t = self.thiele()
        if logging:
            if k:
                s.insert(0, ((), -1, 'kekulized'))
            if h:
                s.append(((), -1, 'implicified'))
            if t:
                s.append(((), -1, 'aromatized'))
            return s
        return k or s or h or t

    def standardize(self: Union['MoleculeContainer', 'Standardize'], *, fix_stereo=True, logging=False) -> \
            Union[bool, List[Tuple[Tuple[int, ...], int, str]]]:
        """
        Standardize functional groups. Return True if any non-canonical group found.

        :param logging: return list of fixed atoms with matched rules.
        """
        neutralized = self.neutralize(fix_stereo=False, logging=logging)
        hs, log = self.__standardize()
        if hs:
            if not neutralized:
                self.flush_cache()
            for n in hs:
                self._calc_implicit(n)
            # second round. need for intersected groups.
            hs, lg = self.__standardize()
            log.extend(lg)
            for n in hs:
                self._calc_implicit(n)
            if fix_stereo:
                self._fix_stereo()
            if logging:
                if neutralized:
                    log.append((tuple(neutralized), -1, 'neutralized'))
                return log
            return True
        if neutralized:
            if fix_stereo:
                self._fix_stereo()
            if logging:
                log.append((tuple(neutralized), -1, 'neutralized'))
                return log
            return True
        if logging:
            return log
        return False

    def neutralize(self: Union['MoleculeContainer', 'Standardize'], *, fix_stereo=True, logging=False) -> \
            Union[bool, List[int]]:
        """
        Transform biradical or dipole resonance structures into neutral form. Return True if structure form changed.

        :param logging: return list of changed atoms.
        """
        atoms = self._atoms
        charges = self._charges
        radicals = self._radicals
        bonds = self._bonds
        entries, exits, rads, constrains = self.__entries()
        hs = set()
        while len(rads) > 1:
            n = rads.pop()
            for path in self.__find_delocalize_path(n, rads, constrains):
                l, m, b = path[-1]
                if b == 1:  # required pi-bond
                    continue
                try:
                    atoms[m].valence_rules(charges[m], False, sum(int(y) for x, y in bonds[m].items() if x != l) + b)
                except ValenceError:
                    continue
                self.__patch_path(path)
                radicals[n] = radicals[m] = False
                rads.discard(m)
                hs.add(n)
                hs.update(x for _, x, _ in path)
                break  # path found
            # path not found. atom n keep as is
        while entries and exits:
            n = entries.pop()
            for path in self.__find_delocalize_path(n, exits, constrains):
                l, m, b = path[-1]
                try:
                    atoms[m].valence_rules(charges[m] - 1, radicals[m],
                                           sum(int(y) for x, y in bonds[m].items() if x != l) + b)
                except ValenceError:
                    continue
                self.__patch_path(path)
                charges[n] = charges[m] = 0
                exits.discard(m)
                hs.add(n)
                hs.update(x for _, x, _ in path)
                break  # path from negative atom to positive atom found.
            # path not found. keep negative atom n as is
        if hs:
            self.flush_cache()
            for n in hs:
                self._calc_implicit(n)
                self._calc_hybridization(n)
            if fix_stereo:
                self._fix_stereo()
            if logging:
                return list(hs)
            return True
        if logging:
            return []
        return False

    def remove_hydrogen_bonds(self: 'MoleculeContainer', *, keep_to_terminal=True, fix_stereo=True) -> int:
        """Remove hydrogen bonds marked with 8 (any) bond

        :param keep_to_terminal: Keep any bonds to terminal hydrogens
        :return: removed bonds count
        """
        bonds = self._bonds
        hg = defaultdict(set)
        for n, atom in self._atoms.items():
            if atom.atomic_number == 1:
                for m, b in bonds[n].items():
                    if b.order == 8:
                        hg[n].add(m)

        if keep_to_terminal:
            for n, ms in hg.items():
                if len(bonds[n]) == len(ms):  # H~A or A~H~A etc case
                    m = ms.pop()
                    if m in hg:  # H~H case
                        hg[m].discard(n)

        seen = set()
        c = 0
        for n, ms in hg.items():
            seen.add(n)
            for m in ms:
                if m in seen:
                    continue
                c += 1
                del bonds[n][m], bonds[m][n]
        if c:
            self.flush_cache()
            if fix_stereo:
                self._fix_stereo()
        return c

    def implicify_hydrogens(self: 'MoleculeContainer', *, fix_stereo=True) -> int:
        """
        Remove explicit hydrogen if possible.
        Works only with Kekule forms of aromatic structures.
        Keeps isotopes of hydrogen.

        :return: number of removed hydrogens
        """
        atoms = self._atoms
        charges = self._charges
        radicals = self._radicals
        bonds = self._bonds
        explicit = defaultdict(list)
        for n, atom in atoms.items():
            if atom.atomic_number == 1 and (atom.isotope is None or atom.isotope == 1):
                if len(bonds[n]) > 1:
                    raise ValenceError(f'Hydrogen atom {{{n}}} has invalid valence. Try to use remove_hydrogen_bonds()')
                for m, b in bonds[n].items():
                    if b.order == 1:
                        if atoms[m].atomic_number != 1:  # not H-H
                            explicit[m].append(n)
                    elif b.order != 8:
                        raise ValenceError(f'Hydrogen atom {{{n}}} has invalid valence {{{b.order}}}.')

        to_remove = set()
        for n, hs in explicit.items():
            atom = atoms[n]
            charge = charges[n]
            is_radical = radicals[n]
            len_h = len(hs)
            for i in range(len_h, 0, -1):
                hi = hs[:i]
                explicit_sum = 0
                explicit_dict = defaultdict(int)
                for m, bond in bonds[n].items():
                    if m not in hi and bond.order != 8:
                        explicit_sum += bond.order
                        explicit_dict[(bond.order, atoms[m].atomic_number)] += 1
                try:
                    rules = atom.valence_rules(charge, is_radical, explicit_sum)
                except ValenceError:
                    break
                if any(s.issubset(explicit_dict) and all(explicit_dict[k] >= c for k, c in d.items()) and h >= i
                       for s, d, h in rules):
                    to_remove.update(hi)
                    break
        for n in to_remove:
            self.delete_atom(n)
        if to_remove and fix_stereo:
            self._fix_stereo()
        return len(to_remove)

    def explicify_hydrogens(self: 'MoleculeContainer', *, fix_stereo=True, start_map=None, return_maps=False) -> int:
        """
        Add explicit hydrogens to atoms.

        :return: number of added atoms
        """
        hydrogens = self._hydrogens
        to_add = []
        for n, h in hydrogens.items():
            try:
                to_add.extend([n] * h)
            except TypeError:
                raise ValenceError(f'atom {{{n}}} has valence error')

        if return_maps:
            log = []
        if to_add:
            bonds = self._bonds
            m = start_map
            for n in to_add:
                m = self.add_atom(H(), _map=m)
                bonds[n][m] = bonds[m][n] = Bond(1)
                hydrogens[n] = 0
                if return_maps:
                    log.append((n, m))
                m += 1

            if fix_stereo:
                self._fix_stereo()
            if return_maps:
                return log
            return len(to_add)
        if return_maps:
            return log
        return 0

    def check_valence(self: 'MoleculeContainer') -> List[int]:
        """
        Check valences of all atoms.

        Works only on molecules with aromatic rings in Kekule form.
        :return: list of invalid atoms
        """
        atoms = self._atoms
        charges = self._charges
        radicals = self._radicals
        bonds = self._bonds
        hydrogens = self._hydrogens
        errors = set(atoms)
        for n, atom in atoms.items():
            charge = charges[n]
            is_radical = radicals[n]
            explicit_sum = 0
            explicit_dict = defaultdict(int)
            for m, bond in bonds[n].items():
                order = bond.order
                if order == 4:  # aromatic rings not supported
                    break
                elif order != 8:  # any bond used for complexes
                    explicit_sum += order
                    explicit_dict[(order, atoms[m].atomic_number)] += 1
            else:
                try:
                    rules = atom.valence_rules(charge, is_radical, explicit_sum)
                except ValenceError:
                    pass
                else:
                    hs = hydrogens[n]
                    for s, d, h in rules:
                        if h == hs and s.issubset(explicit_dict) and all(explicit_dict[k] >= c for k, c in d.items()):
                            errors.discard(n)
                            break
        return list(errors)

    def clean_isotopes(self: 'MoleculeContainer') -> bool:
        """
        Clean isotope marks from molecule.
        Return True if any isotope found.
        """
        atoms = self._atoms
        isotopes = [x for x in atoms.values() if x.isotope]
        if isotopes:
            for i in isotopes:
                i._Core__isotope = None
            self.flush_cache()
            self._fix_stereo()
            return True
        return False

    def __standardize(self: Union['MoleculeContainer', 'Standardize']):
        atom_map = {'charge': self._charges, 'is_radical': self._radicals, 'hybridization': self._hybridizations}
        bonds = self._bonds
        hs = set()
        log = []
        flush = False
        for r, (pattern, atom_fix, bonds_fix) in enumerate(self.__standardize_compiled_rules):
            seen = set()
            for mapping in pattern.get_mapping(self, automorphism_filter=False):
                match = set(mapping.values())
                if not match.isdisjoint(seen):  # skip intersected groups
                    continue
                seen.update(match)
                for n, fix in atom_fix.items():
                    n = mapping[n]
                    for key, value in fix.items():
                        atom_map[key][n] = value
                for n, m, b in bonds_fix:
                    n = mapping[n]
                    m = mapping[m]
                    if m in bonds[n]:
                        bonds[n][m]._Bond__order = b
                        if b == 8:  # expected original molecule don't contain `any` bonds or these bonds not changed
                            flush = True
                    else:
                        flush = True
                        bonds[n][m] = bonds[m][n] = Bond(b)
                log.append((tuple(match), r, str(pattern)))
                # flush cache
                if flush:
                    try:
                        del self.__dict__['__cached_args_method_neighbors']
                    except KeyError:  # already flushed before
                        pass
                    flush = False
            hs.update(seen)
        return hs, log

    def __patch_path(self: 'MoleculeContainer', path):
        bonds = self._bonds
        for n, m, b in path:
            bonds[n][m]._Bond__order = b

    def __find_delocalize_path(self: 'MoleculeContainer', start, finish, constrains):
        bonds = self._bonds
        stack = [(start, n, 0, b.order + 1) for n, b in bonds[start].items() if n in constrains and b.order < 3]
        path = []
        seen = {start}
        while stack:
            last, current, depth, order = stack.pop()
            if len(path) > depth:
                seen.difference_update(x for _, x, _ in path[depth:])
                path = path[:depth]

            path.append((last, current, order))

            if current in finish:
                if depth:  # one bonded ignored. we search double bond transfer! A=A-A >> A-A=A.
                    yield path
                continue  # stop grow

            depth += 1
            seen.add(current)
            diff = -1 if depth % 2 else 1
            stack.extend((current, n, depth, b) for n, b in  # I want 3.8
                         ((n, b.order + diff) for n, b in bonds[current].items() if n not in seen and n in constrains)
                         if 1 <= b <= 3)

    def __entries(self: 'MoleculeContainer'):
        charges = self._charges
        radicals = self._radicals
        atoms = self._atoms
        hybs = self._hybridizations
        bonds = self._bonds

        transfer = set()
        entries = set()
        exits = set()
        rads = set()
        for n, a in atoms.items():
            if a.atomic_number not in {5, 6, 7, 8, 14, 15, 16, 33, 34, 52} or hybs[n] == 4:
                # filter non-organic set, halogens and aromatics
                continue
            if charges[n] == -1:
                if len(bonds[n]) == 4 and a.atomic_number == 5:  # skip boron
                    continue
                entries.add(n)
            elif charges[n] == 1:
                if len(bonds[n]) == 4 and a.atomic_number == 7:  # skip boron
                    continue
                exits.add(n)
            elif radicals[n]:
                rads.add(n)
            transfer.add(n)
        return entries, exits, rads, transfer

    @class_cached_property
    def __standardize_compiled_rules(self):
        rules = []
        for atoms, bonds, atom_fix, bonds_fix in self.__standardize_rules():
            q = query.QueryContainer()
            for a in atoms:
                q.add_atom(**a)
            for n, m, b in bonds:
                q.add_bond(n, m, b)
            rules.append((q, atom_fix, bonds_fix))
        return rules

    @staticmethod
    def __standardize_rules():
        rules = []

        # Triazine-like
        #
        #   B - N        [B-] - [N+]
        #  /     \       /        \
        # N       B >> [N+]       [B-]
        #  \     /       \        /
        #   B - N        [B-] - [N+]
        #
        atoms = ({'atom': 'B', 'neighbors': 4, 'hybridization': 1}, {'atom': 'N', 'neighbors': 4, 'hybridization': 1},
                 {'atom': 'B', 'neighbors': 4, 'hybridization': 1}, {'atom': 'N', 'neighbors': 4, 'hybridization': 1},
                 {'atom': 'B', 'neighbors': 4, 'hybridization': 1}, {'atom': 'N', 'neighbors': 4, 'hybridization': 1})
        bonds = ((1, 2, 1), (2, 3, 1), (3, 4, 1), (4, 5, 1), (5, 6, 1), (1, 6, 1))
        atom_fix = {1: {'charge': -1}, 2: {'charge': 1}, 3: {'charge': -1}, 4: {'charge': 1}, 5: {'charge': -1},
                    6: {'charge': 1}}
        bonds_fix = ()
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Boron
        #
        #     A   A             A      A
        #     |   |             |      |
        # A - B = N - A >> A - [B-] - [N+] - A
        #     |   |             |      |
        #     A   A             A      A
        #
        atoms = ({'atom': 'B', 'neighbors': 4, 'hybridization': 2}, {'atom': 'N', 'neighbors': 4, 'hybridization': 2})
        bonds = ((1, 2, 2),)
        atom_fix = {1: {'charge': -1, 'hybridization': 1}, 2: {'charge': 1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Boron
        #
        #      С            С
        #     //           //
        # B - N   >> B- - [N+]
        #      \            \
        #       A            A
        #
        atoms = ({'atom': 'B'}, {'atom': 'N', 'neighbors': 3, 'hybridization': 2}, {'atom': 'C'})
        bonds = ((1, 2, 1), (2, 3, 2))
        atom_fix = {1: {'charge': -1}, 2: {'charge': 1}}
        bonds_fix = ()
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Boron
        #
        #     A   A             A      A
        #     |   |             |      |
        # A - B - N - A >> A - [B-] - [N+] - A
        #     |   |             |      |
        #     A   A             A      A
        #
        atoms = ({'atom': 'B', 'neighbors': 4, 'hybridization': 1}, {'atom': 'N', 'neighbors': 4, 'hybridization': 1})
        bonds = ((1, 2, 1),)
        atom_fix = {1: {'charge': -1}, 2: {'charge': 1}}
        bonds_fix = ()
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Boron
        #
        #     A    A        A       A
        #     |   //        |      //
        # A - B - N >> A - [B-] - [N+]
        #     |    \        |       \
        #     A     A       A        A
        #
        atoms = ({'atom': 'B', 'neighbors': 4, 'hybridization': 1}, {'atom': 'N', 'neighbors': 3, 'hybridization': 2})
        bonds = ((1, 2, 1),)
        atom_fix = {1: {'charge': -1}, 2: {'charge': 1}}
        bonds_fix = ()
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Boron-Nitrogen
        #
        #   A      A       A   A
        #   |      |       |   |
        #  [B-] = [N+] >>  B - N
        #   |      |       |   |
        #   A      A       A   A
        #
        atoms = ({'atom': 'B', 'charge': -1, 'neighbors': 3, 'hybridization': 2},
                 {'atom': 'N', 'charge': 1, 'neighbors': 3, 'hybridization': 2})
        bonds = ((1, 2, 2),)
        atom_fix = {1: {'charge': 0, 'hybridization': 1}, 2: {'charge': 0, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        #
        # A   H   A     A     H   A
        #  \ / \ /       \  .. \ /
        #   B   B    >>   B     B
        #  / \ / \       / \  .. \
        # A   H   A     A   H     A
        #
        atoms = ({'atom': 'B', 'neighbors': 4}, {'atom': 'B', 'neighbors': 4},
                 {'atom': 'H', 'neighbors': 2}, {'atom': 'H', 'neighbors': 2})
        bonds = ((1, 3, 1), (1, 4, 1), (2, 3, 1), (2, 4, 1))
        atom_fix = {}
        bonds_fix = ((1, 3, 8), (2, 4, 8))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        #
        #        [A-]                 A
        #         |                   |
        # [A-] - [B+3] - [A-] >> A - [B-] - A
        #         |                   |
        #        [A-]                 A
        #
        atoms = ({'atom': 'B', 'charge': 3, 'neighbors': 4}, {'atom': 'A', 'charge': -1}, {'atom': 'A', 'charge': -1},
                 {'atom': 'A', 'charge': -1}, {'atom': 'A', 'charge': -1})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 5, 1))
        atom_fix = {1: {'charge': -1}, 2: {'charge': 0}, 3: {'charge': 0}, 4: {'charge': 0}, 5: {'charge': 0}}
        bonds_fix = ()
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        #
        #        A             A
        #        |             |
        # [A-] - B - A >> A - [B-] - A
        #        |             |
        #        A             A
        #
        atoms = ({'atom': 'B', 'neighbors': 4}, {'atom': 'A', 'charge': -1},
                 {'atom': 'A'}, {'atom': 'A'}, {'atom': 'A'})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 5, 1))
        atom_fix = {1: {'charge': -1}, 2: {'charge': 0}}
        bonds_fix = ()
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        #
        #      A             A
        #      |             |
        # A -  B - A >> A - [B-] - A
        #      |             |
        #      A             A
        #
        atoms = ({'atom': 'B', 'neighbors': 4, 'hybridization': 1},)
        bonds = ()
        atom_fix = {1: {'charge': -1}}
        bonds_fix = ()
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Ammonium
        #
        #     A         A
        #     |         |
        #     N    >>  [N+]
        #   / | \     / | \
        #  A  A  A   A  A  A
        #
        atoms = ({'atom': 'N', 'neighbors': 4, 'hybridization': 1},)
        bonds = ()
        atom_fix = {1: {'charge': 1}}
        bonds_fix = ()
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitrone
        #
        #      O          [O-]
        #     //          /
        # C = N  >> C = [N+]
        #      \          \
        #       A          A
        #
        atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'O', 'neighbors': 1}, {'atom': 'C'}, {'atom': 'A'})
        bonds = ((1, 2, 2), (1, 3, 2), (1, 4, 1))
        atom_fix = {1: {'charge': 1, 'hybridization': 2}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitrone like
        #
        #      N          [N-]
        #     //          /
        # C = N  >> C = [N+]
        #      \         \
        #       A         A
        #
        atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'N', 'neighbors': (1, 2), 'hybridization': 2},
                 {'atom': 'C'}, {'atom': 'A'})
        bonds = ((1, 2, 2), (1, 3, 2), (1, 4, 1))
        atom_fix = {1: {'charge': 1, 'hybridization': 2}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitrone
        #
        # N = N = O >> N = [N+] - [O-]
        #     |             |
        #     A             A
        #
        atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'O', 'neighbors': 1}, {'atom': 'N'}, {'atom': 'A'})
        bonds = ((1, 2, 2), (1, 3, 2), (1, 4, 1))
        atom_fix = {1: {'charge': 1, 'hybridization': 2}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitrone
        #
        # [N-] - [N+] = O >> N = [N+] - [O-]
        #         |               |
        #         A               A
        #
        atoms = ({'atom': 'N', 'charge': 1, 'neighbors': 3}, {'atom': 'O', 'neighbors': 1},
                 {'atom': 'N', 'charge': -1, 'hybridization': 1, 'neighbors': (1, 2)}, {'atom': 'A'})
        bonds = ((1, 2, 2), (1, 3, 1), (1, 4, 1))
        atom_fix = {2: {'charge': -1, 'hybridization': 1}, 3: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 1), (1, 3, 2))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitrone
        #
        # [N+] = N = N - ? >> [N+] = [N+] - [N-] - ?
        #        |                    |
        #        A                    A
        #
        atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'N', 'neighbors': (1, 2), 'hybridization': 2},
                 {'atom': 'N', 'charge': 1, 'hybridization': 2, 'neighbors': 3}, {'atom': 'A'})
        bonds = ((1, 2, 2), (1, 3, 2), (1, 4, 1))
        atom_fix = {1: {'charge': 1, 'hybridization': 2}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitrone [CN(=O)=N(=O)C]
        #
        # [N+] = N = O >> [N+] = [N+] - O-
        #        |                |
        #        A                A
        #
        atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'O', 'neighbors': 1},
                 {'atom': 'N', 'charge': 1, 'hybridization': 2, 'neighbors': 3}, {'atom': 'A'})
        bonds = ((1, 2, 2), (1, 3, 2), (1, 4, 1))
        atom_fix = {1: {'charge': 1, 'hybridization': 2}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitro
        #
        #      O          O
        #     //         //
        # C = N  >> C - [N+]
        #      \         \
        #       OH       [O-]
        #
        atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'O', 'neighbors': 1}, {'atom': 'O', 'neighbors': 1},
                 {'atom': 'C', 'hybridization': 2})
        bonds = ((1, 2, 1), (1, 3, 2), (1, 4, 2))
        atom_fix = {1: {'charge': 1, 'hybridization': 2}, 2: {'charge': -1}, 4: {'hybridization': 1}}
        bonds_fix = ((1, 4, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitro
        #
        #       [O-]          [O-]
        #       /             /
        # C = [N+ ] >>  C - [N+]
        #      \             \\
        #       OH            O
        #
        atoms = ({'atom': 'N', 'charge': 1, 'neighbors': 3}, {'atom': 'O', 'charge': -1, 'neighbors': 1},
                 {'atom': 'O', 'neighbors': 1}, {'atom': 'C', 'hybridization': 2})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 2))
        atom_fix = {3: {'hybridization': 2}, 4: {'hybridization': 1}}
        bonds_fix = ((1, 3, 2), (1, 4, 1))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitro
        #
        #      O           [O-]
        #     //           /
        # A - N   >> A - [N+]
        #     \\          \\
        #      O           O
        #
        atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'O', 'neighbors': 1},
                 {'atom': 'O', 'neighbors': 1}, {'atom': 'A'})
        bonds = ((1, 2, 2), (1, 3, 2), (1, 4, 1))
        atom_fix = {1: {'charge': 1, 'hybridization': 2}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitro
        #
        #         O              [O-]
        #        //              /
        # [A-] - N   >> [A-] - [N+]
        #        \\             \\
        #         O              O
        #
        atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'O', 'neighbors': 1},
                 {'atom': 'O', 'neighbors': 1}, {'atom': 'A', 'charge': -1})
        bonds = ((1, 2, 2), (1, 3, 2), (1, 4, 1))
        atom_fix = {1: {'charge': 1, 'hybridization': 2}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitro
        #
        # O : N : O      O = [N+] - [O-]
        #     |      >>       |
        #     A               A
        #
        atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'O', 'neighbors': 1},
                 {'atom': 'O', 'neighbors': 1}, {'atom': 'A'})
        bonds = ((1, 2, 4), (1, 3, 4), (1, 4, 1))
        atom_fix = {1: {'charge': 1, 'hybridization': 2}, 2: {'charge': -1, 'hybridization': 1},
                    3: {'hybridization': 1}}
        bonds_fix = ((1, 2, 1), (1, 3, 2))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # N-same
        #
        #      N           [N-]
        #     //           /
        # A - N   >> A - [N+]
        #     \\          \\
        #      N           N
        #
        atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'N', 'hybridization': 2, 'neighbors': (1, 2)},
                 {'atom': 'N', 'hybridization': 2}, {'atom': 'A'})
        bonds = ((1, 2, 2), (1, 3, 2), (1, 4, 1))
        atom_fix = {1: {'charge': 1, 'hybridization': 2}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitrite
        #
        #   O        [O-]
        #  //        /
        # [N-]  >>  N
        #  \\       \\
        #   O        O
        #
        atoms = ({'atom': 'N', 'neighbors': 2, 'charge': -1}, {'atom': 'O', 'neighbors': 1},
                 {'atom': 'O', 'neighbors': 1})
        bonds = ((1, 2, 2), (1, 3, 2))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Azide
        #
        #  A - [N-] - [N+] # N  >> A - N = [N+] = [N-]
        #
        atoms = ({'atom': 'N', 'charge': 1, 'neighbors': 2},
                 {'atom': 'N', 'charge': -1, 'neighbors': 2, 'hybridization': 1}, {'atom': 'N', 'neighbors': 1})
        bonds = ((1, 2, 1), (1, 3, 3))
        atom_fix = {2: {'charge': 0, 'hybridization': 2}, 3: {'charge': -1, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2), (1, 3, 2))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Azide
        #
        #  A - [N+] # N = [N-]  >> A - N = [N+] = [N-]
        #
        atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'N', 'charge': 1, 'neighbors': 2},
                 {'atom': 'N', 'charge': -1, 'neighbors': 1}, {'atom': 'A'})
        bonds = ((1, 2, 3), (1, 3, 2), (2, 4, 1))
        atom_fix = {1: {'charge': 1}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Azide
        #
        # A - N = N # N >> A - N = [N+] = [N-]
        #
        atoms = ({'atom': 'N', 'neighbors': 2},
                 {'atom': 'N', 'neighbors': 2, 'hybridization': 2}, {'atom': 'N', 'neighbors': 1})
        bonds = ((1, 2, 2), (1, 3, 3))
        atom_fix = {1: {'charge': 1}, 3: {'charge': -1, 'hybridization': 2}}
        bonds_fix = ((1, 3, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Azide
        #
        # A - N = N = N >> A - N = [N+] = [N-]
        #
        atoms = ({'atom': 'N', 'neighbors': 2},
                 {'atom': 'N', 'neighbors': 2, 'hybridization': 2}, {'atom': 'N', 'neighbors': 1})
        bonds = ((1, 2, 2), (1, 3, 2))
        atom_fix = {1: {'charge': 1}, 3: {'charge': -1}}
        bonds_fix = ()
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Azide
        #
        # A - NH - N # N >> A - N = [N+] = [N-]
        #
        atoms = ({'atom': 'N', 'neighbors': 2},
                 {'atom': 'N', 'neighbors': 2, 'hybridization': 1}, {'atom': 'N', 'neighbors': 1})
        bonds = ((1, 2, 1), (1, 3, 3))
        atom_fix = {1: {'charge': 1}, 2: {'hybridization': 2}, 3: {'charge': -1, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2), (1, 3, 2))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Azide
        #
        # [N-] # N = N - A >> [N-] == [N+] == N - A
        #
        atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'N', 'neighbors': 1, 'hybridization': 3, 'charge': -1},
                 {'atom': 'N', 'hybridization': 2, 'neighbors': 2})
        bonds = ((1, 2, 3), (1, 3, 2))
        atom_fix = {1: {'charge': 1}, 2: {'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Azide
        #
        # [N-] == N # N >> [N-] == [N+] == [N-]
        #
        atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'N', 'neighbors': 1, 'charge': -1},
                 {'atom': 'N', 'neighbors': 1})
        bonds = ((1, 2, 2), (1, 3, 3))
        atom_fix = {1: {'charge': 1}, 3: {'hybridization': 2, 'charge': -1}}
        bonds_fix = ((1, 3, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitrile oxide
        #
        # - C # N = O >> - C # [N+] - [O-]
        #
        atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'O', 'neighbors': 1}, {'atom': 'C'})
        bonds = ((1, 2, 2), (1, 3, 3))
        atom_fix = {1: {'charge': 1}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # N-oxide
        #
        #    |           |
        #  - N - ? >> - [N+] - ?
        #    \\          |
        #     O         [O-]
        #
        atoms = ({'atom': 'N', 'neighbors': (3, 4), 'hybridization': 2}, {'atom': 'O', 'neighbors': 1})
        bonds = ((1, 2, 2),)
        atom_fix = {1: {'charge': 1, 'hybridization': 1}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # N-oxide radical
        #
        #    |         |
        #  - N*  >>  - N
        #    \\        |
        #     O        O*
        #
        atoms = ({'atom': 'N', 'neighbors': 3, 'hybridization': 2, 'is_radical': True}, {'atom': 'O', 'neighbors': 1})
        bonds = ((1, 2, 2),)
        atom_fix = {1: {'hybridization': 1, 'is_radical': False}, 2: {'is_radical': True, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitrous oxide
        #
        # O = N # N >> O = [N+] = [N-]
        #
        atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'N', 'neighbors': 1}, {'atom': 'O', 'neighbors': 1})
        bonds = ((1, 2, 3), (1, 3, 2))
        atom_fix = {1: {'charge': 1}, 2: {'charge': -1, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitrous oxide
        #
        # [O-] - [N+] # N >> O = [N+] = [N-]
        #
        atoms = ({'atom': 'N', 'neighbors': 2, 'charge': 1}, {'atom': 'N', 'neighbors': 1},
                 {'atom': 'O', 'neighbors': 1, 'charge': -1, 'hybridization': 1})
        bonds = ((1, 2, 3), (1, 3, 1))
        atom_fix = {2: {'charge': -1, 'hybridization': 2}, 3: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2), (1, 3, 2))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Diazo
        #
        #  [C-] - [N+] # N  >>  C = [N+] = [N-]
        #
        atoms = ({'atom': 'N', 'charge': 1, 'neighbors': 2}, {'atom': 'N', 'neighbors': 1},
                 {'atom': 'C', 'charge': -1, 'hybridization': 1})
        bonds = ((1, 2, 3), (1, 3, 1))
        atom_fix = {2: {'charge': -1, 'hybridization': 2}, 3: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2), (1, 3, 2))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Diazo
        #
        # C = N # N >> C = [N+] = [N-]
        #
        atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'N', 'neighbors': 1}, {'atom': 'C', 'hybridization': 2})
        bonds = ((1, 2, 3), (1, 3, 2))
        atom_fix = {1: {'charge': 1}, 2: {'charge': -1, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Diazo
        #
        # C # N == N - ? >> [C-] = [N+] = N - ?
        #
        atoms = ({'atom': 'C', 'neighbors': 1, 'hybridization': 3}, {'atom': 'N', 'neighbors': 2, 'hybridization': 3},
                 {'atom': 'N', 'neighbors': (1, 2), 'hybridization': 2})
        bonds = ((1, 2, 3), (2, 3, 2))
        atom_fix = {1: {'charge': -1, 'hybridization': 2}, 2: {'charge': 1}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Diazo
        #
        # A - C # N = N - ? >> A - [C-] = [N+] = N - ?
        #
        atoms = ({'atom': 'C', 'neighbors': 2, 'hybridization': 3}, {'atom': 'N', 'neighbors': 2, 'hybridization': 3},
                 {'atom': 'N', 'neighbors': (1, 2), 'hybridization': 2}, {'atom': 'A'})
        bonds = ((1, 2, 3), (2, 3, 2), (1, 4, 1))
        atom_fix = {1: {'charge': -1, 'hybridization': 2}, 2: {'charge': 1}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Diazo
        #
        # C           C
        #  \           \
        #   N # N >>   [N+] = [N-]
        #  /           /
        # C           C
        #
        atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'N', 'neighbors': 1}, {'atom': 'C'}, {'atom': 'C'})
        bonds = ((1, 2, 3), (1, 3, 1), (1, 4, 1))
        atom_fix = {1: {'charge': 1, 'hybridization': 2}, 2: {'charge': -1, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Diazo
        #
        #    |           |
        #  - N - ? >> - [N+] - ?
        #    \\          |
        #     N         [N-]
        #
        atoms = ({'atom': 'N', 'neighbors': (3, 4), 'hybridization': 2}, {'atom': 'N', 'hybridization': 2})
        bonds = ((1, 2, 2),)
        atom_fix = {1: {'charge': 1, 'hybridization': 1}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Diazonium
        #
        #  C - N = [N+]  >>  C - [N+] # N
        #
        atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'N', 'charge': 1, 'neighbors': 1},
                 {'atom': 'C', 'hybridization': 1})
        bonds = ((1, 2, 2), (1, 3, 1))
        atom_fix = {1: {'charge': 1, 'hybridization': 3}, 2: {'charge': 0, 'hybridization': 3}}
        bonds_fix = ((1, 2, 3),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Isocyanate
        #
        #  [N+] - [C-] = O  >>  N = C = O
        #
        atoms = ({'atom': 'C', 'charge': -1, 'neighbors': 2},
                 {'atom': 'N', 'charge': 1, 'neighbors': (1, 2), 'hybridization': 1}, {'atom': 'O', 'neighbors': 1})
        bonds = ((1, 2, 1), (1, 3, 2))
        atom_fix = {1: {'charge': 0, 'hybridization': 3}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitroso
        #
        # - [N+] - [O-]  >>  - N = O
        #
        atoms = ({'atom': 'N', 'charge': 1, 'neighbors': 2, 'hybridization': 1},
                 {'atom': 'O', 'charge': -1, 'neighbors': 1})
        bonds = ((1, 2, 1),)
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Iminium
        #
        #  [C+] - N(R2)  >> C = [N+](R2)
        #
        atoms = ({'atom': 'N', 'neighbors': 3, 'hybridization': 1},
                 {'atom': 'C', 'charge': 1, 'hybridization': 1, 'neighbors': (1, 2, 3)})
        bonds = ((1, 2, 1),)
        atom_fix = {1: {'charge': 1, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitrilium
        #
        #  [C+] = N(R)  >>  C # [N+](R)
        #
        atoms = ({'atom': 'N', 'neighbors': (1, 2), 'hybridization': 2},
                 {'atom': 'C', 'charge': 1, 'neighbors': (1, 2), 'hybridization': 2})
        bonds = ((1, 2, 2),)
        atom_fix = {1: {'charge': 1, 'hybridization': 3}, 2: {'charge': 0, 'hybridization': 3}}
        bonds_fix = ((1, 2, 3),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Amide
        #
        #      [O,S]H    [O,S]
        #      /         //
        # N = C  >> NH - C
        #
        atoms = ({'atom': 'C', 'hybridization': 2}, {'atom': ListElement(['O', 'S']), 'neighbors': 1},
                 {'atom': 'N', 'hybridization': 2})
        bonds = ((1, 2, 1), (1, 3, 2))
        atom_fix = {2: {'hybridization': 2}, 3: {'hybridization': 1}}
        bonds_fix = ((1, 2, 2), (1, 3, 1))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Phosphonic
        #
        #       O                   O
        #       |                   |
        #  C - [P+] - [O-]  >>  C - P = O
        #       |                   |
        #       O                   O
        #
        atoms = ({'atom': 'P', 'charge': 1, 'neighbors': 4}, {'atom': 'O', 'charge': -1, 'neighbors': 1},
                 {'atom': 'O'}, {'atom': 'O'}, {'atom': 'C'})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 5, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Phosphonium ylide
        #
        #       C                   C
        #       |                   |
        #  C - [P-] - [C+]  >>  C - P = C
        #       |                   |
        #       C                   C
        #
        atoms = ({'atom': 'P', 'charge': -1, 'neighbors': 4}, {'atom': 'C', 'charge': 1, 'neighbors': (1, 2, 3)},
                 {'atom': 'C'}, {'atom': 'C'},
                 {'atom': 'C'})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 5, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Hexafluorophosphate
        #
        #   F   F         F   F
        #    \ /           \ /
        # F - P - F >> F - [P-] - F
        #    / \           / \
        #   F   F         F   F
        #
        atoms = ({'atom': 'P', 'neighbors': 6}, {'atom': 'F', 'neighbors': 1}, {'atom': 'F', 'neighbors': 1},
                 {'atom': 'F', 'neighbors': 1}, {'atom': 'F', 'neighbors': 1}, {'atom': 'F', 'neighbors': 1},
                 {'atom': 'F', 'neighbors': 1})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 5, 1), (1, 6, 1), (1, 7, 1))
        atom_fix = {1: {'charge': -1}}
        bonds_fix = ()
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Sulfodioxide
        #
        #       C                   C
        #       |                   |
        #  O = [S+] - [O-]  >>  O = S = O
        #       |                   |
        #       C                   C
        #
        atoms = ({'atom': 'S', 'charge': 1, 'neighbors': 4},
                 {'atom': 'O', 'charge': -1, 'neighbors': 1}, {'atom': 'O', 'neighbors': 1},
                 {'atom': 'C'}, {'atom': 'C'})
        bonds = ((1, 2, 1), (1, 3, 2), (1, 4, 1), (1, 5, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Sulfoxonium ylide
        #
        #       C                   C
        #       |                   |
        #  C = [S+] - [O-]  >>  C = S = O
        #       |                   |
        #       C                   C
        #
        atoms = ({'atom': 'S', 'charge': 1, 'neighbors': 4},
                 {'atom': 'O', 'charge': -1, 'neighbors': 1}, {'atom': 'C'}, {'atom': 'C'}, {'atom': 'C'})
        bonds = ((1, 2, 1), (1, 3, 2), (1, 4, 1), (1, 5, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Sulfoxide
        #
        #           C            C
        #          /            /
        # [O-] - [S+]  >>  O = S
        #          \            \
        #           C            C
        #
        atoms = ({'atom': 'S', 'charge': 1, 'neighbors': 3},
                 {'atom': 'O', 'charge': -1, 'neighbors': 1}, {'atom': 'C'}, {'atom': 'C'})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Sulfone
        #
        #       [O-]          O
        #       /            //
        # A = [S+2]  >>  A = S
        #       \            \\
        #       [O-]          O
        #
        atoms = ({'atom': 'S', 'charge': 2, 'neighbors': 3}, {'atom': 'O', 'charge': -1, 'neighbors': 1},
                 {'atom': 'O', 'charge': -1, 'neighbors': 1}, {'atom': 'A'})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 2))
        atom_fix = {1: {'charge': 0, 'hybridization': 3}, 2: {'charge': 0, 'hybridization': 2},
                    3: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2), (1, 3, 2))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Sulfone
        #
        #          A                  A
        #          |                  |
        # [O-] - [S+2] - [O-] >>  O = S = O
        #          |                  |
        #          A                  A
        #
        atoms = ({'atom': 'S', 'charge': 2, 'neighbors': 4}, {'atom': 'O', 'charge': -1, 'neighbors': 1},
                 {'atom': 'O', 'charge': -1, 'neighbors': 1}, {'atom': 'A'}, {'atom': 'A'})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 5, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 3}, 2: {'charge': 0, 'hybridization': 2},
                    3: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2), (1, 3, 2))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Sulfonium ylide
        #
        #  C - [S-] - [C+]  >>  C - S = C
        #       |                   |
        #       C                   C
        #
        atoms = ({'atom': 'S', 'charge': -1, 'neighbors': 3},
                 {'atom': 'C', 'charge': 1, 'hybridization': 1, 'neighbors': (1, 2, 3)}, {'atom': 'C'}, {'atom': 'C'})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Sulfite
        #
        #  O              O
        #  \\             \\
        #  [S-] - A  >>    S - A
        #  //             /
        #  O            [O-]
        #
        atoms = ({'atom': 'S', 'charge': -1, 'neighbors': 3}, {'atom': 'O', 'neighbors': 1},
                 {'atom': 'O', 'neighbors': 1}, {'atom': 'A'})
        bonds = ((1, 2, 2), (1, 3, 2), (1, 4, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Sulfite
        #
        #           A           A
        #          /           /
        # [O-] - [S+]  >> O = S
        #          \           \
        #           N           N
        #
        atoms = ({'atom': 'O', 'charge': -1, 'neighbors': 1}, {'atom': 'N', 'hybridization': 1},
                 {'atom': 'S', 'charge': 1, 'neighbors': 3}, {'atom': 'A'})
        bonds = ((1, 3, 1), (2, 3, 1), (2, 4, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 3: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 3, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Sulfine
        #
        #  [O-] - [S+] = C  >>  O = S = C
        #
        atoms = ({'atom': 'S', 'charge': 1, 'neighbors': 2},
                 {'atom': 'O', 'charge': -1, 'neighbors': 1}, {'atom': 'C'})
        bonds = ((1, 2, 1), (1, 3, 2))
        atom_fix = {1: {'charge': 0, 'hybridization': 3}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # amide S
        #
        #      A1,3                 A1,3
        #      |                    |
        #  N = S - [OH]  >>  [NH] - S = O
        #
        atoms = ({'atom': 'S', 'neighbors': (3, 5), 'hybridization': 2},
                 {'atom': 'O', 'neighbors': 1}, {'atom': 'N', 'neighbors': (1, 2), 'hybridization': 2})
        bonds = ((1, 2, 1), (1, 3, 2))
        atom_fix = {2: {'hybridization': 2}, 3: {'hybridization': 1}}
        bonds_fix = ((1, 2, 2), (1, 3, 1))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # amide S(VI)
        #
        #     [OH]                O
        #      |                 //
        #  N = S = N  >>  [NH] - S - [NH]
        #      |                 \\
        #     [OH]                O
        #
        atoms = ({'atom': 'S', 'neighbors': 4}, {'atom': 'O', 'neighbors': 1},
                 {'atom': 'N', 'neighbors': (1, 2), 'hybridization': 2}, {'atom': 'O', 'neighbors': 1},
                 {'atom': 'N', 'neighbors': (1, 2), 'hybridization': 2})
        bonds = ((1, 2, 1), (1, 3, 2), (1, 4, 1), (1, 5, 2))
        atom_fix = {2: {'hybridization': 2}, 3: {'hybridization': 1}, 4: {'hybridization': 2}, 5: {'hybridization': 1}}
        bonds_fix = ((1, 2, 2), (1, 3, 1), (1, 4, 2), (1, 5, 1))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # amide S(VI)
        #
        #     [OH]                O
        #      |                 //
        #  N = S = A  >>  [NH] - S = A
        #      |                 |
        #      A                 A
        #
        atoms = ({'atom': 'S', 'neighbors': 4}, {'atom': 'O', 'neighbors': 1},
                 {'atom': 'N', 'neighbors': (1, 2), 'hybridization': 2}, {'atom': 'A'}, {'atom': 'A'})
        bonds = ((1, 2, 1), (1, 3, 2), (1, 4, 2), (1, 5, 1))
        atom_fix = {2: {'hybridization': 2}, 3: {'hybridization': 1}}
        bonds_fix = ((1, 2, 2), (1, 3, 1))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Silicate Selenite
        #
        #           O             O
        #          /             /
        # [O-] - [Si+]  >>  O = Si
        #          \             \
        #           O             O
        #
        atoms = ({'atom': 'Si', 'charge': 1, 'neighbors': 3}, {'atom': 'O', 'charge': -1, 'neighbors': 1},
                 {'atom': 'O'}, {'atom': 'O'})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        atoms = ({'atom': 'Se', 'charge': 1, 'neighbors': 3}, {'atom': 'O', 'charge': -1, 'neighbors': 1},
                 {'atom': 'O'}, {'atom': 'O'})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Carbon Monoxide
        #
        # [CX1] = O  >> [С-] # [O+]
        #
        atoms = ({'atom': 'C', 'neighbors': 1, 'is_radical': True}, {'atom': 'O', 'neighbors': 1})
        bonds = ((1, 2, 2), )
        atom_fix = {1: {'charge': -1, 'hybridization': 3, 'is_radical': False}, 2: {'charge': 1, 'hybridization': 3}}
        bonds_fix = ((1, 2, 3),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        #
        # C # C - [O,NH,S]H  >> C=C=[O,NH,S]
        #
        atoms = ({'atom': ListElement(['O', 'S', 'N']), 'neighbors': 1}, {'atom': 'C', 'neighbors': 2},
                 {'atom': 'C', 'neighbors': (1, 2)})
        bonds = ((1, 2, 1), (2, 3, 3))
        atom_fix = {1: {'hybridization': 2}, 3: {'hybridization': 2}}
        bonds_fix = ((1, 2, 2), (2, 3, 2))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        #
        # C # C - [NH]R  >> C=C=[NH]R
        #
        atoms = ({'atom': 'N', 'neighbors': 2, 'hybridization': 1}, {'atom': 'C', 'neighbors': 2},
                 {'atom': 'C', 'neighbors': (1, 2)})
        bonds = ((1, 2, 1), (2, 3, 3))
        atom_fix = {1: {'hybridization': 2}, 3: {'hybridization': 2}}
        bonds_fix = ((1, 2, 2), (2, 3, 2))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Ozone
        #
        # [O] -- O -- [O]  >>  O == [O+] -- [O-]
        #
        atoms = ({'atom': 'O', 'neighbors': 1, 'is_radical': True}, {'atom': 'O', 'neighbors': 2},
                 {'atom': 'O', 'neighbors': 1, 'is_radical': True})
        bonds = ((1, 2, 1), (2, 3, 1))
        atom_fix = {1: {'hybridization': 2, 'is_radical': False}, 2: {'charge': 1, 'hybridization': 2},
                    3: {'charge': -1, 'is_radical': False}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Br-ion + I-ion
        #
        #        A           A
        #        |           |
        # [Br-].[I+] >> Br - I
        #        |           |
        #        A           A
        #
        atoms = ({'atom': 'Br', 'charge': -1, 'neighbors': 0}, {'atom': 'A'}, {'atom': 'A'},
                 {'atom': 'I', 'charge': 1, 'neighbors': 2},)
        bonds = ((2, 4, 1), (3, 4, 1))
        atom_fix = {1: {'charge': 0}, 4: {'charge': 0}}
        bonds_fix = ((1, 4, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # CuCl
        #
        #         R - N - C                  R - N - C
        #            /    ||                    /    ||
        # Cl - Cu = C     || >> [Cl-] --- Cu - C     ||
        #            \    ||                   \\    ||
        #         R - N - C                R - [N+]- C
        #
        atoms = ({'atom': 'Cl', 'charge': 0, 'neighbors': 1}, {'atom': 'Cu', 'neighbors': 2},
                 {'atom': 'C', 'charge': 0, 'neighbors': 3},
                 {'atom': 'N', 'neighbors': 3, 'hybridization': 1}, {'atom': 'N', 'neighbors': 3, 'hybridization': 1},
                 {'atom': 'C', 'hybridization': 2}, {'atom': 'C', 'hybridization': 2})
        bonds = ((1, 2, 1), (2, 3, 2), (3, 4, 1), (3, 5, 1), (4, 6, 1), (5, 7, 1), (6, 7, 2))
        atom_fix = {1: {'charge': -1}, 2: {'hybridization': 1}, 4: {'charge': 1, 'hybridization': 2}}
        bonds_fix = ((1, 2, 8), (2, 3, 1), (3, 4, 2))
        rules.append((atoms, bonds, atom_fix, bonds_fix))
        return rules


__all__ = ['Standardize']
