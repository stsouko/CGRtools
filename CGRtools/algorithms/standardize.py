# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
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


class Standardize:
    __slots__ = ()

    def standardize(self):
        """
        standardize functional groups
        """
        atom_map = {'charge': self._charges, 'radical': self._radicals, 'hybridization': self._hybridizations}
        bonds = self._bonds
        shg = self._hydrogens
        seen = set()
        hs = set()
        for pattern, atom_fix, bonds_fix in self._standardize_compiled_rules:
            for mapping in pattern.get_mapping(self):
                n = mapping[1]
                if n in seen:  # automorphism filter
                    continue
                seen.add(n)
                hs.update(mapping.values())

                for n, fix in atom_fix.items():
                    n = mapping[n]
                    for key, value in fix.items():
                        atom_map[key][n] = value
                for n, m, b in bonds_fix:
                    bonds[mapping[n]][mapping[m]]._Bond__order = b
        if seen:
            self.flush_cache()
            for n in hs:
                shg[n] = self._calc_implicit(n)

    @staticmethod
    def _standardize_rules():
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
        atoms = ({'atom': 'N', 'neighbors': 3}, {'atom': 'O'}, {'atom': 'O'})
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
        atoms = ({'atom': 'N', 'neighbors': 4}, {'atom': 'O'}, {'atom': 'C'})
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

        return rules


__all__ = ['Standardize']
