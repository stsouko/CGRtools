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
from itertools import combinations, product
from io import StringIO, TextIOWrapper
from logging import warning
from math import sqrt
from pathlib import Path
from traceback import format_exc
from typing import List, Iterable, Tuple, Optional
from warnings import warn
from ..containers import MoleculeContainer


class XYZ:
    def __init__(self, radius_multiplier=1.25):
        """
        :param radius_multiplier: Multiplier of sum of covalent radii of atoms which has bonds
        """
        self.__radius = radius_multiplier

    def _convert_structure(self, matrix: Iterable[Tuple[str, Optional[int], float, float, float]], charge=0, radical=0):
        mol = MoleculeContainer()
        atoms = mol._atoms
        charges = mol._charges
        radicals = mol._radicals

        conformer = {}
        defined_charges = {}
        for a, c, x, y, z in matrix:
            n = mol.add_atom(a, xy=(x, y))
            conformer[n] = (x, y, z)
            defined_charges[n] = c

        if any(x is not None for x in defined_charges.values()):
            charge = sum(defined_charges.values())

        bonds = self.__get_neighbors(atoms, conformer)
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
                unsaturated[n] = sorted(s, key=lambda x: (x[2], x[0] if x[0] else 5, x[1]), reverse=True)

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

        for s in product(*([(n, c, r) for c, r in s] for n, s in ua.items() if s)):
            for n, c, r in s:
                charges[n] = c
                radicals[n] = r
            if sum(charges.values()) == charge and sum(radicals.values()) == radical:
                break
        else:
            warning('Charge and radical state not balanced.')

        for n in unsaturated:
            mol._calc_implicit(n)
        mol.neutralize()
        mol._conformers.append(conformer)
        return mol

    def __get_neighbors(self, atoms, conformer):
        possible_bonds = defaultdict(dict)  # distance matrix
        for (n, (nx, ny, nz)), (m, (mx, my, mz)) in combinations(conformer.items(), 2):
            d = sqrt((nx - mx) ** 2 + (ny - my) ** 2 + (nz - mz) ** 2)
            r = (atoms[n].atomic_radius + atoms[m].atomic_radius) * self.__radius
            if d <= r:
                possible_bonds[n][m] = possible_bonds[m][n] = d
        return possible_bonds

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
        # get isolated atoms. atoms should be charged or radical
        dots = {}
        saturation = []
        electrons = []
        while True:
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
            self.__file = open(file)
            self.__is_buffer = False
        elif isinstance(file, Path):
            self.__file = file.open()
            self.__is_buffer = False
        elif isinstance(file, (TextIOWrapper, StringIO)):
            self.__file = file
            self.__is_buffer = True
        else:
            raise TypeError('invalid file. TextIOWrapper, StringIO subclasses possible')
        super().__init__(**kwargs)
        self._data = self.__reader()

    def close(self, force=False):
        """
        Close opened file

        :param force: Force closing of externally opened file or buffer
        """
        if not self.__is_buffer or force:
            self.__file.close()

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
        return (x for x in self._data if x is not None)

    def __next__(self):
        return next(iter(self))

    def __reader(self):
        failkey = True
        meta = False
        xyz = charge = size = radical = None
        for n, line in enumerate(self.__file):
            if failkey:
                try:
                    size = int(line)
                except ValueError:
                    warning(f'Line [{n}] {line}: consist errors:\n{format_exc()}')
                else:
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
                            warning(f'Line [{n}] {line}: consist errors:\n{format_exc()}')
                            yield None
                            break
                    elif x.startswith('radical='):
                        try:
                            radical = int(line[8:])
                        except ValueError:
                            failkey = True
                            warning(f'Line [{n}] {line}: consist errors:\n{format_exc()}')
                            yield None
                            break
                else:
                    meta = False
            elif len(xyz) < size:  # XYZ block
                try:
                    symbol, x, y, z = line.split()
                    xyz.append((symbol, float(x), float(y), float(z)))
                except ValueError:
                    failkey = True
                    warning(f'Line [{n}] {line}: consist errors:\n{format_exc()}')
                    yield None
                else:
                    if len(xyz) == size:
                        yield self._convert_structure(xyz, charge, radical)
                        failkey = True  # trigger end of XYZ
        if not failkey:  # cut XYZ
            warning('Last structure not finished')
            yield None

    def _convert_structure(self, matrix: Iterable[Tuple[str, float, float, float]], charge=0, radical=0):
        return super()._convert_structure([(e, None, x, y, z) for e, x, y, z in matrix], charge, radical)

    def parse(self, matrix: Iterable[Tuple[str, float, float, float]], charge: int = 0, radical: int = 0):
        return self._convert_structure(matrix, charge, radical)

    def from_xyz(self, matrix, charge=0, radical=0):
        warn('.from_xyz() deprecated. Use .parse() instead', DeprecationWarning)
        return self._convert_structure(matrix, charge, radical)


__all__ = ['XYZRead', 'XYZ']
