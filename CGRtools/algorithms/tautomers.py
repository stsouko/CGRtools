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

    def input_atom(self):
        atoms = self._atoms
        bonds = self._bonds
        charges = self._charges
        radicals = self._radicals
        hydrogens = self._hydrogens

        ent_atoms = []
        for i, a in atoms.items():
            if a.atomic_number in {5, 6, 7, 8, 14, 15, 16, 33, 34, 35, 52, 53, 85}:
                if radicals[i]:
                    ent_atoms.append((a, 0))
                elif charges[i] < 0:
                    ent_atoms.append((a, 1))
                elif hydrogens[i]:
                    ent_atoms.append((a, 2))
            return ent_atoms



__all__ = ['Tautomers']
