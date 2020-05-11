# -*- coding: utf-8 -*-
#
#  Copyright 2018-2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ..containers import query  # cyclic imports resolve
from ..exceptions import ValenceError


class Standardize:
    __slots__ = ()
    __class_cache__ = {}

    def standardize(self) -> bool:
        """
        Standardize functional groups. Return True if any not canonical group found.
        """
        neutralized = self.neutralize()
        atom_map = {'charge': self._charges, 'is_radical': self._radicals, 'hybridization': self._hybridizations}
        bonds = self._bonds
        hs = set()
        for pattern, atom_fix, bonds_fix in self.__standardize_compiled_rules:
            for mapping in pattern.get_mapping(self):
                hs.update(mapping.values())

                for n, fix in atom_fix.items():
                    n = mapping[n]
                    for key, value in fix.items():
                        atom_map[key][n] = value
                for n, m, b in bonds_fix:
                    bonds[mapping[n]][mapping[m]]._Bond__order = b
        if hs:
            if not neutralized:
                self.flush_cache()
            for n in hs:
                self._calc_implicit(n)
            return True
        return False

    def neutralize(self) -> bool:
        """
        Transform biradical or dipole resonance structures into neutral form. Return True if structure form changed.
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
                charges[n] += 1
                charges[m] -= 1
                exits.discard(m)
                hs.add(n)
                hs.update(x for _, x, _ in path)
        if hs:
            self.flush_cache()
            for n in hs:
                self._calc_implicit(n)
                self._calc_hybridization(n)
            return True
        return False

    def clean_isotopes(self) -> bool:
        """
        Clean isotope marks from molecule.
        Return True if any isotope found.
        """
        isotopes = [x for x in self._atoms if x.isotope]
        if isotopes:
            for i in isotopes:
                i._Core__isotope = None
            self.flush_cache()
            return True
        return False

    def __patch_path(self, path):
        bonds = self._bonds
        for n, m, b in path:
            bonds[n][m]._Bond__order = b

    def __find_delocalize_path(self, start, finish, constrains):
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

    def __entries(self):
        charges = self._charges
        radicals = self._radicals
        atoms = self._atoms

        transfer = set()
        entries = set()
        exits = set()
        rads = set()
        for n, a in atoms.items():
            if a.atomic_number not in {5, 6, 7, 8, 14, 15, 16, 33, 34, 52}:  # filter non-organic set and halogens
                continue
            transfer.add(n)
            if charges[n] == -1:
                entries.add(n)
            elif charges[n] == 1:
                exits.add(n)
            elif radicals[n]:
                rads.add(n)
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

        # Nitrone
        #
        #      O          O-
        #     //         /
        # C = N  >> C = N+
        #      \         \
        #       C         C
        #
        atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'O'}, {'atom': 'C', 'hybridization': 2}, {'atom': 'C'})
        bonds = ((1, 2, 2), (1, 3, 2), (1, 4, 1))
        atom_fix = {1: {'charge': 1, 'hybridization': 2}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitrone
        #
        # N = N = O >> N = N+ - O-
        #     |            |
        #     C            C
        #
        atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'O'}, {'atom': 'N', 'hybridization': 2}, {'atom': 'C'})
        bonds = ((1, 2, 2), (1, 3, 2), (1, 4, 1))
        atom_fix = {1: {'charge': 1, 'hybridization': 2}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitrone
        #
        # N- - N+ = O >> N = N+ - O-
        #      |             |
        #      C             C
        #
        atoms = ({'atom': 'N', 'charge': 1, 'neighbors': 3}, {'atom': 'O'},
                 {'atom': 'N', 'charge': -1, 'hybridization': 1}, {'atom': 'C'})
        bonds = ((1, 2, 2), (1, 3, 1), (1, 4, 1))
        atom_fix = {2: {'charge': -1, 'hybridization': 1}, 3: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 1), (1, 3, 2))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitro
        #
        #      O         O
        #     //        //
        # C = N  >> C - N+
        #      \         \
        #       OH        O-
        #
        atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'O', 'neighbors': 1}, {'atom': 'O'},
                 {'atom': 'C', 'hybridization': 2})
        bonds = ((1, 2, 1), (1, 3, 2), (1, 4, 2))
        atom_fix = {1: {'charge': 1, 'hybridization': 2}, 2: {'charge': -1}, 4: {'hybridization': 1}}
        bonds_fix = ((1, 4, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitro
        #
        #       O-         O-
        #      /          /
        # C = N+  >>  C - N+
        #      \          \\
        #       OH         O
        #
        atoms = ({'atom': 'N', 'charge': 1, 'neighbors': 3}, {'atom': 'O', 'charge': -1}, {'atom': 'O', 'neighbors': 1},
                 {'atom': 'C', 'hybridization': 2})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 2))
        atom_fix = {3: {'hybridization': 2}, 4: {'hybridization': 1}}
        bonds_fix = ((1, 3, 2), (1, 4, 1))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitro
        #
        #   O          O-
        #  //         /
        #  N  >> C - N+
        #  \\        \\
        #   O         O
        #
        atoms = ({'atom': 'N', 'neighbors': 3, 'hybridization': 3}, {'atom': 'O'}, {'atom': 'O'})
        bonds = ((1, 2, 2), (1, 3, 2))
        atom_fix = {1: {'charge': 1, 'hybridization': 2}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Azide
        #
        #  N- - N+ # N  >>  N = N+ = N-
        #
        atoms = ({'atom': 'N', 'charge': 1, 'neighbors': 2},
                 {'atom': 'N', 'charge': -1, 'neighbors': 2, 'hybridization': 1}, {'atom': 'N', 'neighbors': 1})
        bonds = ((1, 2, 1), (1, 3, 3))
        atom_fix = {2: {'charge': 0, 'hybridization': 2}, 3: {'charge': -1, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2), (1, 3, 2))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Azide
        #
        #  N+ # N = N-  >>  N = N+ = N-
        #
        atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'N', 'charge': 1, 'neighbors': 2},
                 {'atom': 'N', 'charge': -1, 'neighbors': 1})
        bonds = ((1, 2, 3), (1, 3, 2))
        atom_fix = {1: {'charge': 1}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Azide
        #
        # N = N # N >> N = N+ = N-
        #
        atoms = ({'atom': 'N', 'neighbors': 2},
                 {'atom': 'N', 'neighbors': 2, 'hybridization': 2}, {'atom': 'N', 'neighbors': 1})
        bonds = ((1, 2, 2), (1, 3, 3))
        atom_fix = {1: {'charge': 1}, 3: {'charge': -1, 'hybridization': 2}}
        bonds_fix = ((1, 3, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Azide
        #
        # - N = N = N >> - N = N+ = N-
        #
        atoms = ({'atom': 'N', 'neighbors': 2},
                 {'atom': 'N', 'neighbors': 2, 'hybridization': 2}, {'atom': 'N', 'neighbors': 1})
        bonds = ((1, 2, 2), (1, 3, 2))
        atom_fix = {1: {'charge': 1}, 3: {'charge': -1}}
        bonds_fix = ()
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Azide
        #
        # - NH - N # N >> - N = N+ = N-
        #
        atoms = ({'atom': 'N', 'neighbors': 2},
                 {'atom': 'N', 'neighbors': 2, 'hybridization': 1}, {'atom': 'N', 'neighbors': 1})
        bonds = ((1, 2, 1), (1, 3, 3))
        atom_fix = {1: {'charge': 1}, 2: {'hybridization': 2}, 3: {'charge': -1, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2), (1, 3, 2))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitrile oxide
        #
        # - C # N = O >> - C # N+ - O-
        #
        atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'O'}, {'atom': 'C'})
        bonds = ((1, 2, 2), (1, 3, 3))
        atom_fix = {1: {'charge': 1}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Tetriary N-oxide
        #
        #    C          C
        #    |          |
        #  - N -  >>  - N+ -
        #    \\         |
        #     O         O-
        #
        atoms = ({'atom': 'N', 'neighbors': 4, 'hybridization': 2}, {'atom': 'O'}, {'atom': 'C'})
        bonds = ((1, 2, 2), (1, 3, 1))
        atom_fix = {1: {'charge': 1, 'hybridization': 1}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Diazo
        #
        #  C- - N+ # N  >>  C = N+ = N-
        #
        atoms = ({'atom': 'N', 'charge': 1, 'neighbors': 2}, {'atom': 'N', 'neighbors': 1},
                 {'atom': 'C', 'charge': -1, 'hybridization': 1})
        bonds = ((1, 2, 3), (1, 3, 1))
        atom_fix = {2: {'charge': -1, 'hybridization': 2}, 3: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2), (1, 3, 2))
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Diazo
        #
        # C = N # N >> C = N+ = N-
        #
        atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'N', 'neighbors': 1}, {'atom': 'C', 'hybridization': 2})
        bonds = ((1, 2, 3), (1, 3, 2))
        atom_fix = {1: {'charge': 1}, 2: {'charge': -1, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Diazonium
        #
        #  C - N = N+  >>  C - N+ # N
        #
        atoms = ({'atom': 'N', 'neighbors': 2}, {'atom': 'N', 'charge': 1, 'neighbors': 1},
                 {'atom': 'C', 'hybridization': 1})
        bonds = ((1, 2, 2), (1, 3, 1))
        atom_fix = {1: {'charge': 1, 'hybridization': 3}, 2: {'charge': 0, 'hybridization': 3}}
        bonds_fix = ((1, 2, 3),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Isocyanate
        #
        #  N+ - C- = O  >>  N = C = O
        #
        atoms = ({'atom': 'C', 'charge': -1, 'neighbors': 2},
                 {'atom': 'N', 'charge': 1, 'neighbors': (1, 2), 'hybridization': 1}, {'atom': 'O'})
        bonds = ((1, 2, 1), (1, 3, 2))
        atom_fix = {1: {'charge': 0, 'hybridization': 3}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Aromatic N-Oxide
        #
        #  : N :  >>  : N+ :
        #    \\          \
        #     O           O-
        #
        atoms = ({'atom': 'N', 'neighbors': 3, 'hybridization': 4}, {'atom': 'O'})
        bonds = ((1, 2, 2),)
        atom_fix = {1: {'charge': 1}, 2: {'charge': -1, 'hybridization': 1}}
        bonds_fix = ((1, 2, 1),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitroso
        #
        # - N+ - O-  >>  - N = O
        #
        atoms = ({'atom': 'N', 'charge': 1, 'neighbors': 2, 'hybridization': 1}, {'atom': 'O', 'charge': -1})
        bonds = ((1, 2, 1),)
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Iminium
        #
        #  C+ - N  >> C = N+
        #
        atoms = ({'atom': 'N', 'neighbors': 3, 'hybridization': 1}, {'atom': 'C', 'charge': 1, 'hybridization': 1})
        bonds = ((1, 2, 1),)
        atom_fix = {1: {'charge': 1, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Nitrilium
        #
        #  C+ = N  >>  C # N+
        #
        atoms = ({'atom': 'N', 'neighbors': (1, 2), 'hybridization': 2}, {'atom': 'C', 'charge': 1, 'neighbors': 2})
        bonds = ((1, 2, 2),)
        atom_fix = {1: {'charge': 1, 'hybridization': 3}, 2: {'charge': 0, 'hybridization': 3}}
        bonds_fix = ((1, 2, 3),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Phosphonic
        #
        #      O                O
        #      |                |
        #  C - P+ - O-  >>  C - P = O
        #      |                |
        #      O                O
        #
        atoms = ({'atom': 'P', 'charge': 1, 'neighbors': 4}, {'atom': 'O', 'charge': -1}, {'atom': 'O'}, {'atom': 'O'},
                 {'atom': 'C'})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 5, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Phosphonium ylide
        #
        #      C                C
        #      |                |
        #  C - P- - C+  >>  C - P = C
        #      |                |
        #      C                C
        #
        atoms = ({'atom': 'P', 'charge': -1, 'neighbors': 4}, {'atom': 'C', 'charge': 1}, {'atom': 'C'}, {'atom': 'C'},
                 {'atom': 'C'})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1), (1, 5, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Sulfoxide
        #
        #      C                C
        #      |                |
        #  O = S+ - O-  >>  O = S = O
        #      |                |
        #      C                C
        #
        atoms = ({'atom': 'S', 'charge': 1, 'neighbors': 4},
                 {'atom': 'O', 'charge': -1}, {'atom': 'O'}, {'atom': 'C'}, {'atom': 'C'})
        bonds = ((1, 2, 1), (1, 3, 2), (1, 4, 1), (1, 5, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Sulfoxonium ylide
        #
        #      C                C
        #      |                |
        #  C = S+ - O-  >>  C = S = O
        #      |                |
        #      C                C
        #
        atoms = ({'atom': 'S', 'charge': 1, 'neighbors': 4},
                 {'atom': 'O', 'charge': -1}, {'atom': 'C'}, {'atom': 'C'}, {'atom': 'C'})
        bonds = ((1, 2, 1), (1, 3, 2), (1, 4, 1), (1, 5, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Sulfon
        #
        #        C           C
        #       /           /
        # O- - S+  >>  O = S
        #       \           \
        #        C           C
        #
        atoms = ({'atom': 'S', 'charge': 1, 'neighbors': 3},
                 {'atom': 'O', 'charge': -1}, {'atom': 'C'}, {'atom': 'C'})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Sulfonium ylide
        #
        #  C - S- - C+  >>  C - S = C
        #      |                |
        #      C                C
        #
        atoms = ({'atom': 'S', 'charge': -1, 'neighbors': 3},
                 {'atom': 'C', 'charge': 1, 'hybridization': 1}, {'atom': 'C'}, {'atom': 'C'})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Sulfine
        #
        #  O- - S+ = C  >>  O = S = C
        #
        atoms = ({'atom': 'S', 'charge': 1, 'neighbors': 2},
                 {'atom': 'O', 'charge': -1}, {'atom': 'C', 'hybridization': 2})
        bonds = ((1, 2, 1), (1, 3, 2))
        atom_fix = {1: {'charge': 0, 'hybridization': 3}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Silicate Selenite
        #
        #        O            O
        #       /            /
        # O- - Si+  >>  O = Si
        #       \            \
        #        O            O
        #
        atoms = ({'atom': 'Si', 'charge': 1, 'neighbors': 3}, {'atom': 'O', 'charge': -1}, {'atom': 'O'}, {'atom': 'O'})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        atoms = ({'atom': 'Se', 'charge': 1, 'neighbors': 3}, {'atom': 'O', 'charge': -1}, {'atom': 'O'}, {'atom': 'O'})
        bonds = ((1, 2, 1), (1, 3, 1), (1, 4, 1))
        atom_fix = {1: {'charge': 0, 'hybridization': 2}, 2: {'charge': 0, 'hybridization': 2}}
        bonds_fix = ((1, 2, 2),)
        rules.append((atoms, bonds, atom_fix, bonds_fix))

        # Carbon Monoxide
        #
        # [CX1] = O  >> С- # O+
        #

        atoms = ({'atom': 'C', 'neighbors': 1, 'is_radical': True}, {'atom': 'O', 'neighbors': 1})
        bonds = ((1, 2, 2), )
        atom_fix = {1: {'charge': -1, 'hybridization': 3, 'is_radical': False}, 2: {'charge': 1, 'hybridization': 3}}
        bonds_fix = ((1, 2, 3),)
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

        return rules


class StandardizeReaction:
    __slots__ = ()
    __class_cache__ = {}

    def standardize(self, fix_mapping: bool = True) -> bool:
        """
        Standardize functional groups. Works only for Molecules.
        Return True if in any molecule found not canonical group.

        :param fix_mapping: Search AAM errors of functional groups.
        """
        total = False
        for m in self.molecules():
            if hasattr(m, 'standardize'):
                if m.standardize() and not total:
                    total = True

        if fix_mapping and self.fix_mapping():
            return True

        if total:
            self.flush_cache()
        return total

    def fix_mapping(self) -> bool:
        """
        Fix atom-to-atom mapping of some functional groups. Return True if found AAM errors.
        """
        seen = set()
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
            return True
        return False

    def clean_isotopes(self) -> bool:
        """
        Clean isotope marks for all molecules in reaction.
        Returns True if in any molecule found isotope.
        """
        flag = False
        for m in self.molecules():
            if hasattr(m, 'clean_isotopes'):
                if m.slean_isotopes() and not flag:
                    flag = True

        if flag:
            self.flush_cache()
        return flag

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
