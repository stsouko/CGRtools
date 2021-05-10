# -*- coding: utf-8 -*-
#
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
from collections import defaultdict
from itertools import chain, islice
from .rw import _MDLWrite
from ...containers import MoleculeContainer, CGRContainer, QueryContainer, QueryCGRContainer


class MDLWrite(_MDLWrite):
    def _convert_structure(self, g):
        if isinstance(g, MoleculeContainer):
            bonds = self.__convert_molecule(g)
        elif isinstance(g, CGRContainer):
            bonds = self.__convert_cgr(g)
        elif isinstance(g, QueryContainer):
            bonds = self.__convert_query(g)
        elif isinstance(g, QueryCGRContainer):
            raise TypeError('CGR queries not supported')
        else:
            raise TypeError('Graph expected')
        if max(g) > 999:
            raise ValueError('MOL file support only small molecules')

        head = f'{g.name}\n\n\n{g.atoms_count:3d}{g.bonds_count:3d}  0  0  0  0            999 V2000\n'

        gc = g._charges
        gr = g._radicals
        out = [bonds]
        for n, (m, a) in enumerate(g._atoms.items(), start=1):
            if a.isotope:
                out.append(f'M  ISO  1 {n:3d} {a.isotope:3d}\n')
            if gr[m]:
                out.append(f'M  RAD  1 {n:3d}   2\n')  # invalid for carbenes
            c = gc[m]
            if c in (-4, 4):
                out.append(f'M  CHG  1 {n:3d} {c:3d}\n')
        out.append('M  END\n')

        if self._write3d and isinstance(g, MoleculeContainer) and g._conformers:
            if self._write3d == 2 and len(g._conformers) > 1:
                return [''.join((head, self.__convert_atoms3d(g, xyz), *out)) for xyz in g._conformers]
            else:
                return ''.join((head, self.__convert_atoms3d(g, g._conformers[0]), *out))
        else:
            return ''.join((head, self.__convert_atoms2d(g), *out))

    def __convert_atoms2d(self, g):
        gc = g._charges
        gp = g._plane

        out = []
        for n, (m, a) in enumerate(g._atoms.items(), start=1):
            x, y = gp[m]
            c = gc[m]
            if c in (-4, 4):
                if self._mapping:
                    out.append(f'{x:10.4f}{y:10.4f}    0.0000 {a.atomic_symbol:3s} '
                               f'0  0  0  0  0  0  0  0  0{m:3d}  0  0\n')
                else:
                    out.append(f'{x:10.4f}{y:10.4f}    0.0000 {a.atomic_symbol:3s} '
                               f'0  0  0  0  0  0  0  0  0  0  0  0\n')
            else:
                if self._mapping:
                    out.append(f'{x:10.4f}{y:10.4f}    0.0000 {a.atomic_symbol:3s} 0{self.__charge_map[c]}  0  0  0  0'
                               f'  0  0  0{m:3d}  0  0\n')
                else:
                    out.append(f'{x:10.4f}{y:10.4f}    0.0000 {a.atomic_symbol:3s} 0{self.__charge_map[c]}  0  0  0  0'
                               f'  0  0  0  0  0  0\n')
        return ''.join(out)

    def __convert_atoms3d(self, g, xyz):
        gc = g._charges

        out = []
        for n, (m, a) in enumerate(g._atoms.items(), start=1):
            x, y, z = xyz[m]
            c = gc[m]
            if c in (-4, 4):
                if self._mapping:
                    out.append(f'{x:10.4f}{y:10.4f}{z:10.4f} {a.atomic_symbol:3s} '
                               f'0  0  0  0  0  0  0  0  0{m:3d}  0  0\n')
                else:
                    out.append(f'{x:10.4f}{y:10.4f}{z:10.4f} {a.atomic_symbol:3s} 0  0  0  0  0  0  0  0  0  0  0  0\n')
            else:
                if self._mapping:
                    out.append(f'{x:10.4f}{y:10.4f}{z:10.4f} {a.atomic_symbol:3s} 0{self.__charge_map[c]}  0  0  0  0'
                               f'  0  0  0{m:3d}  0  0\n')
                else:
                    out.append(f'{x:10.4f}{y:10.4f}{z:10.4f} {a.atomic_symbol:3s} 0{self.__charge_map[c]}  0  0  0  0'
                               f'  0  0  0  0  0  0\n')
        return ''.join(out)

    @classmethod
    def __convert_molecule(cls, g):
        bonds = g._bonds
        atoms = {m: n for n, m in enumerate(g._atoms, start=1)}
        wedge = defaultdict(set)
        out = []
        for n, m, s in g._wedge_map:
            out.append(f'{atoms[n]:3d}{atoms[m]:3d}  {bonds[n][m].order}  {s == 1 and "1" or "6"}  0  0  0\n')
            wedge[n].add(m)
            wedge[m].add(n)
        for n, m, b in g.bonds():
            if m not in wedge[n]:
                out.append(f'{atoms[n]:3d}{atoms[m]:3d}  {b.order}  0  0  0  0\n')
        return ''.join(out)

    @staticmethod
    def __convert_cgr(g):
        gpc = g._p_charges
        gpr = g._p_radicals

        atoms = {m: n for n, m in enumerate(g._atoms, start=1)}
        bonds = []
        props = []
        for n, c in g._charges.items():
            pc = gpc[n]
            if c != pc:
                i = len(props) + 1
                props.append(f'M  SAL {i:3d}  1 {atoms[n]:3d}\n'
                             f'M  SDT {i:3d} dynatom\n'
                             f'M  SDD {i:3d}     0.0000{i / 3:10.4f}    DAU   ALL  0       0\n'
                             f'M  SED {i:3d} c{pc - c:+d}\n')
        for n, r in g._radicals.items():
            if r != gpr[n]:
                i = len(props) + 1
                props.append(f'M  SAL {i:3d}  1 {atoms[n]:3d}\n'
                             f'M  SDT {i:3d} dynatom\n'
                             f'M  SDD {i:3d}     0.0000{i / 3:10.4f}    DAU   ALL  0       0\n'
                             f'M  SED {i:3d} r1\n')

        for n, m, b in g.bonds():
            if b.order != b.p_order:
                i = len(props) + 1
                props.append(f'M  SAL {i:3d}  2 {atoms[n]:3d} {atoms[m]:3d}\n'
                             f'M  SDT {i:3d} dynbond\n'
                             f'M  SDD {i:3d}     0.0000{i / 3:10.4f}    DAU   ALL  0       0\n'
                             f'M  SED {i:3d} {b.order or 0}>{b.p_order or 0}\n')
                bonds.append(f'{atoms[n]:3d}{atoms[m]:3d}  8  0  0  0  0\n')
            else:
                bonds.append(f'{atoms[n]:3d}{atoms[m]:3d}  {b.order}  0  0  0  0\n')

        iterator = iter(range(1, len(props) + 1))
        for first in iterator:
            dat = list(chain((first,), islice(iterator, 7)))
            bonds.append(f'M  STY  {len(dat)}')
            bonds.extend(f'{x:4d} DAT' for x in dat)
            bonds.append('\n')

        bonds.extend(props)
        return ''.join(bonds)

    @classmethod
    def __convert_query(cls, g):
        atoms = {m: n for n, m in enumerate(g._atoms, start=1)}
        out = []
        for n, m, b in g.bonds():
            if len(b.order) > 1:
                raise ValueError('supported only simple QueryBond')
            out.append(f'{atoms[n]:3d}{atoms[m]:3d}  {b.order[0]}  0  0  0  0\n')

        props = []
        for n, m in g._neighbors.items():
            if m:
                i = len(props) + 1
                props.append(f'M  SAL {i:3d}  1 {atoms[n]:3d}\n'
                             f'M  SDT {i:3d} neighbors\n'
                             f'M  SDD {i:3d}     0.0000{i / 3:10.4f}    DAU   ALL  0       0\n'
                             f'M  SED {i:3d} {",".join(str(x) for x in m)}\n')
        for n, h in g._hybridizations.items():
            if h:
                i = len(props) + 1
                props.append(f'M  SAL {i:3d}  1 {atoms[n]:3d}\n'
                             f'M  SDT {i:3d} hybridization\n'
                             f'M  SDD {i:3d}     0.0000{i / 3:10.4f}    DAU   ALL  0       0\n'
                             f'M  SED {i:3d} {",".join(str(x) for x in h)}\n')

        iterator = iter(range(1, len(props) + 1))
        for first in iterator:
            dat = list(chain((first,), islice(iterator, 7)))
            out.append(f'M  STY  {len(dat)}')
            out.extend(f'{x:4d} DAT' for x in dat)
            out.append('\n')

        out.extend(props)
        return ''.join(out)

    __charge_map = {-3: '  7', -2: '  6', -1: '  5', 0: '  0', 1: '  3', 2: '  2', 3: '  1'}


__all__ = ['MDLWrite']
