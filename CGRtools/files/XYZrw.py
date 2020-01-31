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
from itertools import combinations
from io import StringIO, TextIOWrapper
from logging import warning
from math import sqrt
from pathlib import Path
from traceback import format_exc
from typing import List, Iterable, Tuple
from ._CGRrw import CGRRead
from ..containers import MoleculeContainer


class XYZRead(CGRRead):
    """XYZ files reader. Works similar to opened file object. Support `with` context manager.
    On initialization accept opened in text mode file, string path to file,
    pathlib.Path object or another buffered reader object.

    Supported multiple structures in same file. In second line possible to store total charge of system. Example::

        2
        charge=-1
        O 0.0 0.0 0.0
        H 1.0 0.0 0.0

    """
    def __init__(self, file, *args, radius_multiplier=1.25, **kwargs):
        """
        :param radius_multiplier: Multiplier of sum of covalent radii of atoms which has bonds
        """
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
        super().__init__(*args, **kwargs)
        self.__radius = radius_multiplier
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
        xyz = charge = size = None
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
                if line.startswith('charge='):
                    try:
                        charge = int(line[7:])
                    except ValueError:
                        failkey = True
                        warning(f'Line [{n}] {line}: consist errors:\n{format_exc()}')
                        yield None
                        continue
                else:
                    charge = 0
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
                        yield self.from_xyz(xyz, charge)
                        failkey = True  # trigger end of XYZ
        if not failkey:  # cut XYZ
            warning('Last structure not finished')
            yield None

    def from_xyz(self, matrix: Iterable[Tuple[str, float, float, float]], charge=0):
        mol = MoleculeContainer()
        atoms = mol._atoms
        conformer = {}
        mol._conformers.append(conformer)

        for a, x, y, z in matrix:
            conformer[mol.add_atom(a, xy=(x, y))] = (x, y, z)

        dd = defaultdict(dict)  # distance matrix
        for (n, (nx, ny, nz)), (m, (mx, my, mz)) in combinations(conformer.items(), 2):
            d = sqrt((nx - mx) ** 2 + (ny - my) ** 2 + (nz - mz) ** 2)
            r = (atoms[n].atomic_radius + atoms[m].atomic_radius) * self.__radius
            if d < r:
                mol.add_bond(n, m, 1)
                dd[n][m] = dd[m][n] = d

        return mol


__all__ = ['XYZRead']
