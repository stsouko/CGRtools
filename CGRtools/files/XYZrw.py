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
from importlib.util import find_spec
from itertools import combinations, product
from io import StringIO, TextIOWrapper
from logging import warning
from math import sqrt
from pathlib import Path
from random import shuffle
from traceback import format_exc
from typing import List, Iterable, Tuple, Optional
from warnings import warn
from ._mdl import parse_error
from ..containers import MoleculeContainer


if find_spec('numpy') and find_spec('numba'):  # try to load numba jit
    from numpy import array, uint16, empty
    from numba import njit, f8, u2, u4
    from numba.core.types import Tuple as nTuple

    @njit(nTuple((u2[:, :], f8[:]))(f8[:, :], f8[:], f8),
          {'size': u2, 'max_bonds': u4, 'c': u4, 'n': u2, 'm': u2, 'rn': f8, 'r': f8, 'd': f8,
           'nx': f8, 'ny': f8, 'nz': f8, 'mx': f8, 'my': f8, 'mz': f8}, cache=True)
    def _get_possible_bonds(xyz, radii, multiplier):
        size = len(xyz)
        max_bonds = size * 10  # each atom has less then 10 neighbors approximately
        nm = empty((max_bonds, 2), dtype=uint16)
        ds = empty(max_bonds)
        c = 0
        for n in range(size - 1):
            nx, ny, nz = xyz[n]
            rn = radii[n]
            for m in range(n + 1, size):
                mx, my, mz = xyz[m]
                d = sqrt((nx - mx) ** 2 + (ny - my) ** 2 + (nz - mz) ** 2)
                r = (rn + radii[m]) * multiplier
                if d <= r:
                    nm[c] = n + 1, m + 1
                    ds[c] = d
                    c += 1
        return nm[:c], ds[:c]

    def get_possible_bonds(atoms, conformer, multiplier):
        possible_bonds = {n: {} for n in atoms}  # distance matrix
        radii = array([a.atomic_radius for a in atoms.values()])
        xyz = array(list(conformer.values()))
        nm, ds = _get_possible_bonds(xyz, radii, multiplier)
        for (n, m), d in zip(nm.tolist(), ds.tolist()):
            possible_bonds[n][m] = possible_bonds[m][n] = d
        return possible_bonds
else:
    def get_possible_bonds(atoms, conformer, multiplier):
        possible_bonds = {n: {} for n in atoms}  # distance matrix
        radii = {n: a.atomic_radius for n, a in atoms.items()}
        for (n, (nx, ny, nz)), (m, (mx, my, mz)) in combinations(conformer.items(), 2):
            d = sqrt((nx - mx) ** 2 + (ny - my) ** 2 + (nz - mz) ** 2)
            r = (radii[n] + radii[m]) * multiplier
            if d <= r:
                possible_bonds[n][m] = possible_bonds[m][n] = d
        return possible_bonds


charge_priority = {0: 0, -1: 1, 1: 2, 2: 3, 3: 4, -2: 5, -3: 6, 4: 7, -4: 8}


class XYZ:
    """
    Override class below then inheritance used.
    """
    MoleculeContainer = MoleculeContainer

    def __init__(self, radius_multiplier=1.25, store_log=False):
        """
        :param radius_multiplier: Multiplier of sum of covalent radii of atoms which has bonds
        :param store_log: Store parser log if exists messages to `.meta` by key `CGRtoolsParserLog`.
        """
        self.__radius = radius_multiplier
        self._store_log = store_log
        self._log_buffer = []

    def _info(self, msg):
        self._log_buffer.append(msg)

    def _flush_log(self):
        self._log_buffer.clear()

    def _format_log(self):
        return '\n'.join(self._log_buffer)

    def close(self, force=False):
        """
        Close opened file

        :param force: Force closing of externally opened file or buffer
        """
        if not self._is_buffer or force:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    def read(self) -> List[MoleculeContainer]:
        """
        Parse whole file

        :return: list of parsed molecules
        """
        return list(iter(self))

    def __iter__(self):
        return (x for x in self._data if not isinstance(x, parse_error))

    def __next__(self):
        return next(iter(self))

    def _convert_structure(self, matrix: Iterable[Tuple[str, Optional[int], float, float, float]], charge=0, radical=0):
        mol = self.MoleculeContainer()
        atoms = mol._atoms
        charges = mol._charges
        radicals = mol._radicals

        conformer = {}
        defined_charges = {}
        for a, c, x, y, z in matrix:
            n = mol.add_atom(a, xy=(x, y))
            conformer[n] = (x, y, z)
            defined_charges[n] = c

        if all(x is not None for x in defined_charges.values()):
            charge = sum(defined_charges.values())

        bonds = get_possible_bonds(atoms, conformer, self.__radius)
        saturation, bonds = self.__get_atom_states_and_bonds(atoms, bonds, defined_charges)

        # set single bonds in molecule. collect unsaturated atoms
        seen = set()
        unsaturated = {}
        for n, env in bonds.items():
            s = saturation[n]
            if len(s) == 1:
                c, r, h = s.pop()
                if not h:
                    seen.add(n)
                    if c:
                        charges[n] = c
                    if r:
                        radicals[n] = True
                    for m in env:
                        if m not in seen:
                            mol.add_bond(n, m, 1)
                else:
                    unsaturated[n] = [(c, r, h)]
            else:
                # radicals is the lowest priority
                # hydrogens is the highest priority
                unsaturated[n] = sorted(s, key=lambda x: (x[1], -x[2] + x[0] if x[0] > 0 else -x[2],
                                                          charge_priority[x[0]]))

        # create graph of unsaturated atoms
        bonds_graph = {n: {m for m in env if m in unsaturated} for n, env in bonds.items() if n in unsaturated}
        ua, sb, sa = self.__saturate(bonds_graph, unsaturated)
        for n, m, b in sb:
            mol.add_bond(n, m, b)
        for n, c, r in sa:
            if c:
                charges[n] = c
            if r:
                radicals[n] = True

        combo_ua = []
        for n, s in ua.items():
            if len(s) == 1:
                c, r = s[0]
                if c:
                    charges[n] = c
                if r:
                    radicals[n] = True
            elif s:
                combo_ua.append([(n, c, r) for c, r in s])

        # try randomly set charges and radicals.
        # first pick required radical states.
        # second try to minimize charge delta.
        if combo_ua:
            need_radical = radical - sum(radicals.values())
            for attempt in range(1, len(combo_ua) + 1):
                shuffle(combo_ua)
                rad = []
                chg = []
                for atom in combo_ua:
                    if len(rad) < need_radical:  # pick radicals
                        r = next((x for x in atom if x[2]), None)
                        if r:  # pick random radical states
                            rad.append(r)
                        else:  # not radical
                            chg.append(atom)
                    else:  # pick not radical states
                        c = [x for x in atom if not x[2]]
                        if len(c) > 1:
                            chg.append(c)
                        elif c:
                            n, c, r = c[0]
                            charges[n] = c
                            radicals[n] = r
                        elif attempt == len(combo_ua):  # all states has radical. balancing impossible
                            chg.append(atom)  # fuck it horse. we in last attempt
                            self._info('Radical state not balanced.')
                        else:  # do next attempt
                            break
                else:
                    for n, c, r in rad:
                        charges[n] = c
                        radicals[n] = r

                    current_charge = sum(charges.values())
                    for x in chg:
                        n, c, r = min(x, key=lambda x: abs(current_charge + x[1]))
                        charges[n] = c
                        radicals[n] = r
                        current_charge += c
                    if sum(charges.values()) == charge:
                        break
            else:
                self._info('Charge state not balanced.')
        elif sum(radicals.values()) != radical or sum(charges.values()) != charge:
            self._info('Charge or radical state not balanced.')

        for n in unsaturated:
            mol._calc_implicit(n)
        mol.neutralize()
        mol._conformers.append(conformer)
        return mol

    @staticmethod
    def __get_atom_states_and_bonds(atoms, possible_bonds, charges):
        possible_bonds = {n: md.copy() for n, md in possible_bonds.items()}
        while True:
            saturation = defaultdict(set)
            for n, env in possible_bonds.items():
                el = len(env)
                dc = charges[n]
                for (charge, is_radical, valence), rules in atoms[n]._compiled_valence_rules.items():
                    if valence < el or dc is not None and dc != charge:
                        continue  # skip impossible rules
                    for _, d, h in rules:
                        if d:  # todo: optimize
                            env_atoms = defaultdict(int)
                            for m in env:
                                env_atoms[atoms[m].atomic_number] += 1
                            bonds = 0
                            for (b, a), c in d.items():  # stage 1
                                bonds += b
                                if a in env_atoms:
                                    if env_atoms[a] < c:
                                        break  # rule not matched
                                    env_atoms[a] -= c
                                else:  # rule not matched
                                    break
                            else:  # stage 2. found possible valence
                                unmatched = sum(env_atoms.values())  # atoms outside rule
                                implicit = valence - bonds + h  # implicit H in rule
                                if unmatched:
                                    if implicit >= unmatched:
                                        # number of implicit H should be greater or equal to number of neighbors
                                        # excess of implicit H saved as unsaturated atom
                                        saturation[n].add((charge, is_radical, implicit - unmatched))
                                else:  # pattern fully matched. save implicit H count as unsaturated atom.
                                    saturation[n].add((charge, is_radical, implicit))
                        elif el == valence:   # unspecific rule. found possible valence
                            saturation[n].add((charge, is_radical, h))
                if n not in saturation:  # valence not found
                    break
            else:  # all atoms passed
                break
            out = max(env.items(), key=lambda x: x[1])[0]
            del possible_bonds[out][n]
            del possible_bonds[n][out]
        return saturation, possible_bonds

    @staticmethod
    def __saturate(bonds, atoms):
        dots = {}
        saturation = []
        electrons = []
        while True:
            # get isolated atoms. atoms should be charged or radical
            new_dots = {n: [(c, r) for c, r, h in atoms[n] if not h] for n, env in bonds.items() if not env}
            for n in new_dots:
                del bonds[n]
            dots.update(new_dots)
            if not bonds:
                break

            try:  # get terminal atom
                n = next(n for n, ms in bonds.items() if len(ms) == 1)
            except StopIteration:
                # get ring or linker atom
                n, _ = min(bonds.items(), key=lambda x: len(x[1]))
                m = bonds[n].pop()
                bonds[m].discard(n)

                for (nc, nr, nh), (i, (mc, mr, mh)) in product(atoms[n], enumerate(atoms[m])):
                    if nh == mh:
                        saturation.append((n, m, nh + 1))
                        electrons.append((n, nc, nr))
                        electrons.append((m, mc, mr))

                        for x in bonds.pop(n):
                            saturation.append((n, x, 1))
                            bonds[x].discard(n)
                        for x in bonds.pop(m):
                            saturation.append((m, x, 1))
                            bonds[x].discard(m)
                        break
                    elif nh < mh:
                        electrons.append((n, nc, nr))
                        saturation.append((n, m, nh + 1))
                        atoms[m].pop(i)
                        atoms[m].insert(i, (mc, mr, mh - nh))

                        for x in bonds.pop(n):
                            saturation.append((n, x, 1))
                            bonds[x].discard(n)
                        break
                    elif nh > mh:
                        electrons.append((m, mc, mr))
                        saturation.append((n, m, mh + 1))
                        atoms[n].pop(i)
                        atoms[n].insert(i, (nc, nr, nh - mh))

                        for x in bonds.pop(m):
                            saturation.append((m, x, 1))
                            bonds[x].discard(m)
                        break
            else:
                m = bonds.pop(n).pop()
                bonds[m].discard(n)

                for (nc, nr, nh), (i, (mc, mr, mh)) in product(atoms[n], enumerate(atoms[m])):
                    if nh == mh:
                        saturation.append((n, m, nh + 1))
                        electrons.append((n, nc, nr))
                        electrons.append((m, mc, mr))
                        for x in bonds.pop(m):
                            saturation.append((m, x, 1))
                            bonds[x].discard(m)
                        break
                    elif nh < mh and bonds[m]:
                        electrons.append((n, nc, nr))
                        saturation.append((n, m, nh + 1))
                        atoms[m].pop(i)
                        atoms[m].insert(i, (mc, mr, mh - nh))
                        break
                else:
                    saturation.append((n, m, 1))
                    if not bonds[m]:
                        del bonds[m]
        return dots, saturation, electrons


class XYZRead(XYZ):
    """XYZ files reader. Works similar to opened file object. Support `with` context manager.
    On initialization accept opened in text mode file, string path to file,
    pathlib.Path object or another buffered reader object.

    Supported multiple structures in same file. In second line possible to store total charge of system. Example::

        2
        charge=-1
        O 0.0 0.0 0.0
        H 1.0 0.0 0.0

    """
    def __init__(self, file, **kwargs):
        if isinstance(file, str):
            self._file = open(file)
            self._is_buffer = False
        elif isinstance(file, Path):
            self._file = file.open()
            self._is_buffer = False
        elif isinstance(file, (TextIOWrapper, StringIO)):
            self._file = file
            self._is_buffer = True
        else:
            raise TypeError('invalid file. TextIOWrapper, StringIO subclasses possible')
        super().__init__(**kwargs)
        self.__file = iter(self._file.readline, '')
        self._data = self.__reader()

    def __reader(self):
        failkey = True
        meta = False
        xyz = charge = size = radical = None
        file = self._file
        seekable = file.seekable()
        pos = None
        count = -1
        for n, line in enumerate(self.__file):
            if failkey:
                try:
                    size = int(line)
                except ValueError:
                    pass
                else:
                    if seekable:
                        pos = file.tell() - len(line)
                    count += 1
                    meta = True
                    failkey = False
                    xyz = []
            elif meta:  # second line
                charge = 0
                radical = 0
                for x in line.split():
                    if x.startswith('charge='):
                        try:
                            charge = int(line[7:])
                        except ValueError:
                            failkey = True
                            self._info(f'Line [{n}] {line}: consist errors:\n{format_exc()}')
                            yield parse_error(count, pos, self._format_log(), {})
                            self._flush_log()
                            break
                    elif x.startswith('radical='):
                        try:
                            radical = int(line[8:])
                        except ValueError:
                            failkey = True
                            self._info(f'Line [{n}] {line}: consist errors:\n{format_exc()}')
                            yield parse_error(count, pos, self._format_log(), {})
                            self._flush_log()
                            break
                else:
                    meta = False
            elif len(xyz) < size:  # XYZ block
                try:
                    symbol, x, y, z = line.split()
                    xyz.append((symbol, float(x), float(y), float(z)))
                except ValueError:
                    failkey = True
                    self._info(f'Line [{n}] {line}: consist errors:\n{format_exc()}')
                    yield parse_error(count, pos, self._format_log(), {})
                    self._flush_log()
                else:
                    if len(xyz) == size:
                        try:
                            container = self._convert_structure(xyz, charge, radical)
                        except ValueError:
                            self._info(f'record consist errors:\n{format_exc()}')
                            yield parse_error(count, pos, self._format_log(), {})
                        else:
                            if self._store_log:
                                log = self._format_log()
                                if log:
                                    container.meta['CGRtoolsParserLog'] = log
                            yield container
                        self._flush_log()
                        failkey = True  # trigger end of XYZ
        if not failkey:  # cut XYZ
            self._info('Last structure not finished')
            yield parse_error(count, pos, self._format_log(), {})
            self._flush_log()

    def _convert_structure(self, matrix: Iterable[Tuple[str, float, float, float]], charge=0, radical=0):
        return super()._convert_structure([(e, None, x, y, z) for e, x, y, z in matrix], charge, radical)

    def parse(self, matrix: Iterable[Tuple[str, float, float, float]], charge: int = 0, radical: int = 0):
        try:
            container = self._convert_structure(matrix, charge, radical)
        except ValueError:
            self._flush_log()
        else:
            if self._store_log:
                log = self._format_log()
                if log:
                    container.meta['CGRtoolsParserLog'] = log
            return container

    def from_xyz(self, matrix, charge=0, radical=0):
        warn('.from_xyz() deprecated. Use `CGRtools.xyz` or `XYZRead.create_parser` instead', DeprecationWarning)
        warning('.from_xyz() deprecated. Use `CGRtools.xyz` or `XYZRead.create_parser` instead')
        return self.parse(matrix, charge, radical)

    @classmethod
    def create_parser(cls, *args, **kwargs):
        """
        Create XYZ parser function configured same as XYZRead object.
        """
        obj = object.__new__(cls)
        super(XYZRead, obj).__init__(*args, **kwargs)
        return obj.parse


__all__ = ['XYZRead']
