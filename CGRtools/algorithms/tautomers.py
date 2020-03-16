# -*- coding: utf-8 -*-
#
#  Copyright 2020 Nail Samikaev <samikaevn@yandex.ru>
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


class Tautomers:
    __slots__ = ()

    def enumerate_tautomers(self):
        atoms = self._atoms
        bonds = self._bonds
        charges = self._charges
        radicals = self._radicals
        hydrogens = self._hydrogens
        yield self


    def input_atoms(self):
        atoms = self._atoms
        bonds = self._bonds
        charges = self._charges
        hydrogens = self._hydrogens

        ent_atoms = []
        for i, atom in atoms.items():
            if bonds[i]:
                if charges[i] < 0:         # [C, N, O, Si, P,  S,  As, Se, Te]
                    if atom.atomic_number in {6, 7, 8, 14, 15, 16, 33, 34, 52}:
                        ent_atoms.append((atom, i, 1))
                elif hydrogens[i]:         # [N, O, Si, P,  S,  As, Se, Te]
                    if atom.atomic_number in {7, 8, 14, 15, 16, 33, 34, 52}:
                        ent_atoms.append((atom, i, 2))
        return ent_atoms


    def get_paths(self, ent_atoms):
        bonds = self._bonds

        for a in ent_atoms:
            path = [a[1]]
            stack = [(i, n.order, 1) for i, n in bonds[a[1]].items() if n.order < 3]

            while stack:
                cur, b, d = stack.pop()

                if len(path) > d:
                    path = path[:d]

                if set(bonds[cur].keys()).difference(path):
                    for x in bonds[cur].values():
                        if (x.order < 3) and (x.order != b):
                            path.append(cur)

                            if len(path) % 2:
                                yield path
                            break

                elif not len(path) % 2:
                    path.append(cur)
                    yield path

                d += 1
                nbg = [(x, y.order, d) for x, y in bonds[cur].items() if (y.order < 3) and (y.order != b)]
                stack.extend(nbg)



__all__ = ['Tautomers']
