# -*- coding: utf-8 -*-
#
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
from collections import defaultdict
from .rw import _MDLWrite
from ...containers import MoleculeContainer


class EMDLWrite(_MDLWrite):
    def _convert_structure(self, g):
        if not isinstance(g, MoleculeContainer):
            raise TypeError('MoleculeContainer expected')

        '''M  V30 1 C 1.0912 0.1892 0.6426 1
            M  V30 2 C 1.5071 -0.9383 -0.3461 0 CFG=2'''
        mol = [f'M  V30 BEGIN CTAB\nM  V30 COUNTS {g.atoms_count} {g.bonds_count} 0 0 0\nM  V30 BEGIN ATOM']

        if self._write3d and g._conformers:
            if self._write3d == 2:
                mol.extend(self.__convert_atoms3d(g, xyz) for xyz in g._conformers)
            else:
                mol.append(self.__convert_atoms3d(g, g._conformers[0]))
        else:
            mol.append(self.__convert_atoms2d(g))
        mol.append('M  V30 END ATOM\nM  V30 BEGIN BOND')
        mapping = {m: n for n, m in enumerate(g, start=1)}

        wedge = defaultdict(set)
        bonds = g._bonds
        i = 0
        for i, (n, m, s) in enumerate(g._wedge_map, start=1):
            mol.append(f'M  V30 {i} {bonds[n][m].order} {mapping[n]} {mapping[m]} CFG={s == 1 and "1" or "3"}')
            wedge[n].add(m)
            wedge[m].add(n)

        for i, (n, m, b) in enumerate(g.bonds(), start=i + 1):
            if m not in wedge[n]:
                mol.append(f'M  V30 {i} {b.order} {mapping[n]} {mapping[m]}')
        mol.append('M  V30 END BOND\nM  V30 END CTAB\n')
        return '\n'.join(mol)

    def __convert_atoms2d(self, g):
        gc = g._charges
        gr = g._radicals
        gp = g._plane

        out = []
        for n, (m, a) in enumerate(g._atoms.items(), start=1):
            x, y = gp[m]
            c = gc[m]
            c = f' CHG={c}' if c else ''
            r = f' RAD=2' if gr[m] else ''

            if self._mapping:
                out.append(f'M  V30 {n} {a.atomic_symbol} {x} {y} 0 {m}{c}{r}')
            else:
                out.append(f'M  V30 {n} {a.atomic_symbol} {x} {y} 0 0{c}{r}')
        return '\n'.join(out)

    def __convert_atoms3d(self, g, xyz):
        gc = g._charges
        gr = g._radicals

        out = []
        for n, (m, a) in enumerate(g._atoms.items(), start=1):
            x, y, z = xyz[m]
            c = gc[m]
            c = f' CHG={c}' if c else ''
            r = f' RAD=2' if gr[m] else ''

            if self._mapping:
                out.append(f'M  V30 {n} {a.atomic_symbol} {x} {y} {z} {m}{c}{r}')
            else:
                out.append(f'M  V30 {n} {a.atomic_symbol} {x} {y} {z} 0{c}{r}')
        return '\n'.join(out)


__all__ = ['EMDLWrite']
