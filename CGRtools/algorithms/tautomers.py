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
        radicals = self._radicals
        hydrogens = self._hydrogens

        ent_atoms = []
        for i, atom in atoms.items():
            if bonds[i]:
                if radicals[i]:             # [B, C, N, O, Si, P, S, Cl, As, Se, Br, Te, I, At]
                    if atom.atomic_number in {5, 6, 7, 8, 14, 15, 16, 17, 33, 34, 35, 52, 53, 85}: # ???
                        ent_atoms.append((atom, i, 0)) # 0 - radical
                elif charges[i] < 0:        # [B, C, N, O, Si, P, S, Cl, As, Se, Br, Te, I, At]
                    if atom.atomic_number in {5, 6, 7, 8, 14, 15, 16, 17, 33, 34, 35, 52, 53, 85}: # ???
                        ent_atoms.append((atom, i, 1)) # 1 - has a charge
                elif hydrogens[i]:          # [B, C, N, O, Si, P, S, Cl, As, Se, Br, Te, I, At]
                    if atom.atomic_number in {5, 6, 7, 8, 14, 15, 16, 17, 33, 34, 35, 52, 53, 85}: # ???
                        ent_atoms.append((atom, i, 2)) # 2 - has a hydrogens
        return ent_atoms, bonds

    def get_paths(self, ent_atoms):
        bonds = self._bonds
        paths = []
        for z in ent_atoms:
            a = z[1]
            path = [a]
            stack = []
            for i, n in bonds[a].items():
                if n.order < 3:
                    stack.append((i, 1, n.order))
            while stack:
                cur, d, b = stack.pop()
                if len(path) > d:
                    path = path[:d]
                ngbb = bonds[cur]
                ngb = set(ngbb.keys()).difference(path)
                for c in ngb:
                    if (ngbb[c].order < 3) and (ngbb[c].order != b):
                        stack.append((c, d + 1, ngbb[c].order))
                        path.append(cur)
                    else:
                        path.append(cur)
                        if len(path) > 2:
                            paths.append(path)
                if not ngb:
                    path.append(cur)
                    if len(path) > 2:
                        paths.append(path)
        return paths


__all__ = ['Tautomers']
