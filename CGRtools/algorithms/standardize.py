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
from itertools import count
from typing import List, TYPE_CHECKING, Union, Tuple
from ..containers import query  # cyclic imports resolve
from ..containers.bonds import Bond
from ..exceptions import ValenceError
from ..periodictable import ListElement, H


if TYPE_CHECKING:
    from CGRtools import MoleculeContainer, ReactionContainer


class Standardize:
    __slots__ = ()

    def canonicalize(self: 'MoleculeContainer', *,
                     logging=False) -> Union[bool, List[Tuple[Tuple[int, ...], int, str]]]:
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

    def standardize(self: Union['MoleculeContainer', 'Standardize'], *, fix_stereo=True,
                    logging=False) -> Union[bool, List[Tuple[Tuple[int, ...], int, str]]]:
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

    def neutralize(self: Union['MoleculeContainer', 'Standardize'], *, fix_stereo=True,
                   logging=False) -> Union[bool, List[int]]:
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
        Remove explicit hydrogen if possible. Works only with Kekule forms of aromatic structures.

        :return: number of removed hydrogens
        """
        atoms = self._atoms
        charges = self._charges
        radicals = self._radicals
        bonds = self._bonds
        explicit = defaultdict(list)
        for n, atom in atoms.items():
            if atom.atomic_number == 1:
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


class StandardizeReaction:
    __slots__ = ()

    def canonicalize(self: 'ReactionContainer', fix_mapping: bool = True, *,
                     logging=False) -> Union[bool, Tuple[int, Tuple[int, ...], int, str]]:
        """
        Convert molecules to canonical forms of functional groups and aromatic rings without explicit hydrogens.
        Works only for Molecules.
        Return True if in any molecule found not canonical group.

        :param fix_mapping: Search AAM errors of functional groups.
        :param logging: return log from molecules with index of molecule at first position.
            Otherwise return True if these groups found in any molecule.
        """
        if logging:
            total = []
        else:
            total = False
        for n, m in enumerate(self.molecules()):
            if not isinstance(m, Standardize):
                raise TypeError('Only Molecules supported')
            out = m.canonicalize(logging=logging)
            if out:
                if logging:
                    total.extend((n, *x) for x in out)
                else:
                    total = True

        if fix_mapping and self.fix_mapping():
            if logging:
                total.append((-1, (), -1, 'mapping fixed'))
                return total
            return True

        if total:
            self.flush_cache()
        return total

    def standardize(self: 'ReactionContainer', fix_mapping: bool = True, *,
                    logging=False) -> Union[bool, Tuple[int, Tuple[int, ...], int, str]]:
        """
        Standardize functional groups. Works only for Molecules.

        :param fix_mapping: Search AAM errors of functional groups.
        :param logging: return log from molecules with index of molecule at first position.
            Otherwise return True if these groups found in any molecule.
        """
        if logging:
            total = []
        else:
            total = False
        for n, m in enumerate(self.molecules()):
            if not isinstance(m, Standardize):
                raise TypeError('Only Molecules supported')
            out = m.standardize(logging=logging)
            if out:
                if logging:
                    total.extend((n, *x) for x in out)
                else:
                    total = True

        if fix_mapping and self.fix_mapping():
            if logging:
                total.append((-1, (), -1, 'mapping fixed'))
                return total
            return True

        if total:
            self.flush_cache()
        return total

    def neutralize(self: 'ReactionContainer', *, logging=False) -> Union[bool, Tuple[int, Tuple[int, ...]]]:
        """
        Transform biradical or dipole resonance structures into neutral form. Works only for Molecules.

        :param logging: return log from molecules with index of molecule at first position.
            Otherwise return True if these groups found in any molecule.
        """
        if logging:
            total = []
        else:
            total = False
        for n, m in enumerate(self.molecules()):
            if not isinstance(m, Standardize):
                raise TypeError('Only Molecules supported')
            out = m.neutralize(logging=logging)
            if out:
                if logging:
                    total.extend((n, tuple(x)) for x in out)
                else:
                    total = True
        if total:
            self.flush_cache()
        return total

    def fix_mapping(self: Union['ReactionContainer', 'StandardizeReaction']) -> bool:
        """
        Fix atom-to-atom mapping of some functional groups. Return True if found AAM errors.
        """
        seen = set()
        if not (self.reactants and self.products):
            return False
        elif not isinstance(self.reactants[0], Standardize):
            raise TypeError('Only Molecules supported')

        for r_pattern, p_pattern, fix in self.__standardize_compiled_rules:
            found = []
            for m in self.reactants:
                for mapping in r_pattern.get_mapping(m, automorphism_filter=False):
                    if mapping[1] not in seen:
                        found.append(({fix.get(k, k): v for k, v in mapping.items()},
                                      {mapping[k]: mapping[v] for k, v in fix.items()}))

            if not found:
                continue
            for m in self.products:
                for mapping in p_pattern.get_mapping(m, automorphism_filter=False):
                    atom = mapping[1]
                    if atom in seen:
                        continue
                    for n, (k, v) in enumerate(found):
                        if k == mapping:
                            break
                    else:
                        continue

                    del found[n]
                    m.remap(v)
                    seen.add(atom)
        if seen:
            self.flush_cache()
            flag = True
            seen = set()
        else:
            flag = False

        for bad_query, good_query, fix, valid in self.__remapping_compiled_rules:
            cgr = ~self
            first_free = max(cgr) + 1
            free_number = count(first_free)
            cgr_c = set(cgr.center_atoms)
            del self.__dict__['__cached_method_compose']

            for mapping in bad_query.get_mapping(cgr, automorphism_filter=False):
                if not seen.isdisjoint(mapping.values()):  # prevent matching same RC
                    continue
                mapping = {mapping[n]: next(free_number) if m is None else mapping[m] for n, m in fix.items()}

                reverse = {m: n for n, m in mapping.items()}
                for m in self.products:
                    m.remap(mapping)

                check = ~self
                check_c = set(check.center_atoms)
                delta = check_c - cgr_c
                
                for m in good_query.get_mapping(check, automorphism_filter=False):
                    if valid.issubset(m) and delta.issubset(m.values()):
                        seen.update(mapping)
                        break      
                else:
                    # restore old mapping
                    for m in self.products:
                        m.remap(reverse)
                    del self.__dict__['__cached_method_compose']
                    free_number = count(first_free)
                    continue
                break

        if seen:
            self.flush_cache()
            return True
        return flag

    @classmethod
    def load_remapping_rules(cls, reactions):
        """
        Load AAM fixing rules. Required pairs of bad mapped and good mapped reactions.
        Reactants in pairs should be fully equal (equal molecules and equal atom orders).
        Products should be equal but with different atom numbers.
        """
        for bad, good in reactions:
            if str(bad) != str(good):
                raise ValueError('bad and good reaction should be equal')

            cgr_good, cgr_bad = ~good, ~bad
            gc = cgr_good.augmented_substructure(cgr_good.center_atoms, deep=1)
            bc = cgr_bad.augmented_substructure(cgr_bad.center_atoms, deep=1)

            atoms = set(bc.atoms_numbers + gc.atoms_numbers)      

            pr_g, pr_b, re_g, re_b = set(), set(), set(), set()
            for pr in good.products:
                pr_g.update(pr)
            for pr in bad.products:
                pr_b.update(pr) 
            for pr in good.reactants:
                re_g.update(pr)
            for pr in bad.reactants:
                re_b.update(pr)
            atoms.update((re_b.difference(pr_b)).intersection(pr_g))

            strange_atoms = pr_b.difference(pr_g)
            atoms.update(strange_atoms)

            bad_query = cgr_bad.substructure(atoms.intersection(cgr_bad), as_query=True)
            good_query = cgr_good.substructure(atoms.intersection(cgr_good), as_query=True)

            rules = []
            fix = {}
            for mb, mg in zip(bad.products, good.products): 
                fix.update({k: v for k, v in zip(mb, mg) if k != v and k in atoms})
            
            valid = set(fix).difference(strange_atoms)
            rules.append((bad_query, good_query, fix, valid))

        cls.__class_cache__[cls] = {'_StandardizeReaction__remapping_compiled_rules': tuple(rules)}

    @class_cached_property
    def __remapping_compiled_rules(self):
        return ()

    def implicify_hydrogens(self: 'ReactionContainer') -> int:
        """
        Remove explicit hydrogens if possible

        :return: number of removed hydrogens
        """
        total = 0
        for m in self.molecules():
            if not isinstance(m, Standardize):
                raise TypeError('Only Molecules supported')
            total += m.implicify_hydrogens()
        if total:
            self.flush_cache()
        return total

    def explicify_hydrogens(self: 'ReactionContainer') -> int:
        """
        Add explicit hydrogens to atoms

        :return: number of added atoms
        """
        total = 0
        start_map = 0
        for m in self.molecules():
            if not isinstance(m, Standardize):
                raise TypeError('Only Molecules supported')
            map_ = max(m, default=0)
            if map_ > start_map:
                start_map = map_

        mapping = defaultdict(list)
        for m in self.reactants:
            maps = m.explicify_hydrogens(return_maps=True, start_map=start_map + 1)
            if maps:
                for n, h in maps:
                    mapping[n].append(h)
                start_map = maps[-1][1]
                total += len(maps)

        for m in self.reagents:
            maps = m.explicify_hydrogens(return_maps=True, start_map=start_map + 1)
            if maps:
                start_map = maps[-1][1]
                total += len(maps)

        for m in self.products:
            maps = m.explicify_hydrogens(return_maps=True, start_map=start_map + 1)
            if maps:
                total += len(maps)
                remap = {}
                free = []
                for n, h in maps:
                    if n in mapping and mapping[n]:
                        remap[h] = mapping[n].pop()
                        free.append(h)
                    elif free:
                        remap[h] = start_map = free.pop(0)
                    else:
                        start_map = h
                m.remap(remap)

        if total:
            self.flush_cache()
        return total

    def thiele(self: 'ReactionContainer') -> bool:
        """
        Convert structures to aromatic form. Works only for Molecules.
        Return True if in any molecule found kekule ring
        """
        total = False
        for m in self.molecules():
            if not isinstance(m, Standardize):
                raise TypeError('Only Molecules supported')
            if m.thiele() and not total:
                total = True
        if total:
            self.flush_cache()
        return total

    def kekule(self: 'ReactionContainer') -> bool:
        """
        Convert structures to kekule form. Works only for Molecules.
        Return True if in any molecule found aromatic ring
        """
        total = False
        for m in self.molecules():
            if not isinstance(m, Standardize):
                raise TypeError('Only Molecules supported')
            if m.kekule() and not total:
                total = True
        if total:
            self.flush_cache()
        return total

    def clean_isotopes(self: 'ReactionContainer') -> bool:
        """
        Clean isotope marks for all molecules in reaction.
        Returns True if in any molecule found isotope.
        """
        flag = False
        for m in self.molecules():
            if not isinstance(m, Standardize):
                raise TypeError('Only Molecules supported')
            if m.clean_isotopes() and not flag:
                flag = True

        if flag:
            self.flush_cache()
        return flag

    def clean_stereo(self: 'ReactionContainer'):
        """
        Remove stereo data
        """
        for m in self.molecules():
            if not hasattr(m, 'clean_stereo'):
                raise TypeError('Only Molecules and Queries supported')
            m.clean_stereo()
        self.flush_cache()

    def clean2d(self: 'ReactionContainer', **kwargs):
        """
        Recalculate 2d coordinates
        """
        for m in self.molecules():
            m.clean2d(**kwargs)
        self.fix_positions()

    def fix_positions(self: Union['ReactionContainer', 'StandardizeReaction']):
        """
        Fix coordinates of molecules in reaction
        """
        shift_x = 0
        reactants = self.reactants
        amount = len(reactants) - 1
        signs = []
        for m in reactants:
            max_x = m._fix_plane_mean(shift_x)
            if amount:
                max_x += .2
                signs.append(max_x)
                amount -= 1
            shift_x = max_x + 1
        arrow_min = shift_x

        if self.reagents:
            for m in self.reagents:
                max_x = m._fix_plane_min(shift_x + .4, .5)
                shift_x = max_x + 1
            if shift_x - arrow_min < 3:
                shift_x = arrow_min + 3
        else:
            shift_x += 3
        arrow_max = shift_x - 1

        products = self.products
        amount = len(products) - 1
        for m in products:
            max_x = m._fix_plane_mean(shift_x)
            if amount:
                max_x += .2
                signs.append(max_x)
                amount -= 1
            shift_x = max_x + 1
        self._arrow = (arrow_min, arrow_max)
        self._signs = tuple(signs)
        self.flush_cache()

    def check_valence(self: 'ReactionContainer') -> List[Tuple[int, Tuple[int, ...]]]:
        """
        Check valences of all atoms of all molecules.

        Works only on molecules with aromatic rings in Kekule form.
        :return: list of invalid molecules with invalid atoms lists
        """
        out = []
        for n, m in enumerate(self.molecules()):
            if not isinstance(m, Standardize):
                raise TypeError('Only Molecules supported')
            c = m.check_valence()
            if c:
                out.append((n, tuple(c)))
        return out

    @class_cached_property
    def __standardize_compiled_rules(self):
        rules = []
        for (r_atoms, r_bonds), (p_atoms, p_bonds), fix in self.__standardize_rules():
            r_q = query.QueryContainer()
            p_q = query.QueryContainer()
            for a in r_atoms:
                r_q.add_atom(**a)
            for n, m, b in r_bonds:
                r_q.add_bond(n, m, b)
            for a in p_atoms:
                p_q.add_atom(**a)
            for n, m, b in p_bonds:
                p_q.add_bond(n, m, b)
            rules.append((r_q, p_q, fix))
        return rules

    @staticmethod
    def __standardize_rules():
        rules = []

        # Nitro
        #
        #      O
        #     //
        # * - N+
        #      \
        #       O-
        #
        atoms = ({'atom': 'N', 'neighbors': 3, 'hybridization': 2, 'charge': 1},
                 {'atom': 'O', 'neighbors': 1, 'charge': -1}, {'atom': 'O', 'neighbors': 1})
        bonds = ((1, 2, 1), (1, 3, 2))
        fix = {2: 3, 3: 2}
        rules.append(((atoms, bonds), (atoms, bonds), fix))

        # Carbonate
        #
        #      O
        #     //
        # * - C
        #      \
        #       O-
        #
        atoms = ({'atom': 'C', 'neighbors': 3, 'hybridization': 2},
                 {'atom': 'O', 'neighbors': 1, 'charge': -1}, {'atom': 'O', 'neighbors': 1})
        bonds = ((1, 2, 1), (1, 3, 2))
        fix = {2: 3, 3: 2}
        rules.append(((atoms, bonds), (atoms, bonds), fix))

        # Carbon Acid
        #
        #      O
        #     //
        # * - C
        #      \
        #       OH
        #
        atoms = ({'atom': 'C', 'neighbors': 3, 'hybridization': 2},
                 {'atom': 'O', 'neighbors': 1}, {'atom': 'O', 'neighbors': 1})
        bonds = ((1, 2, 1), (1, 3, 2))
        fix = {2: 3, 3: 2}
        rules.append(((atoms, bonds), (atoms, bonds), fix))

        # Phosphate
        #
        #      *
        #      |
        #  * - P = O
        #      |
        #      OH
        #
        atoms = ({'atom': 'P', 'neighbors': 4, 'hybridization': 2},
                 {'atom': 'O', 'neighbors': 1}, {'atom': 'O', 'neighbors': 1})
        bonds = ((1, 2, 1), (1, 3, 2))
        fix = {2: 3, 3: 2}
        rules.append(((atoms, bonds), (atoms, bonds), fix))

        # Nitro addition
        #
        #      O             O -- *
        #     //            /
        # * - N+   >>  * = N+
        #      \            \
        #       O-           O-
        #
        atoms = ({'atom': 'N', 'neighbors': 3, 'charge': 1, 'hybridization': 2},
                 {'atom': 'O', 'neighbors': 1, 'charge': -1}, {'atom': 'O', 'neighbors': 1})
        bonds = ((1, 2, 1), (1, 3, 2))
        p_atoms = ({'atom': 'N', 'neighbors': 3, 'charge': 1, 'hybridization': 2},
                   {'atom': 'O', 'neighbors': 1, 'charge': -1}, {'atom': 'O', 'neighbors': 2})
        p_bonds = ((1, 2, 1), (1, 3, 1))
        fix = {2: 3, 3: 2}
        rules.append(((atoms, bonds), (p_atoms, p_bonds), fix))

        # Sulphate addition
        #
        #      O [3]            O -- * [2]
        #     //               /
        # * = S - *   >>  * = S - *
        #     |               \\
        #     O- [2]           O [3]
        #
        atoms = ({'atom': 'S', 'neighbors': 4, 'hybridization': 3},
                 {'atom': 'O', 'neighbors': 1, 'charge': -1}, {'atom': 'O', 'neighbors': 1})
        bonds = ((1, 2, 1), (1, 3, 2))
        p_atoms = ({'atom': 'S', 'neighbors': 4, 'hybridization': 3},
                   {'atom': 'O', 'neighbors': 2}, {'atom': 'O', 'neighbors': 1})
        p_bonds = ((1, 2, 1), (1, 3, 2))
        fix = {2: 3, 3: 2}
        rules.append(((atoms, bonds), (p_atoms, p_bonds), fix))

        return rules


__all__ = ['Standardize', 'StandardizeReaction']
