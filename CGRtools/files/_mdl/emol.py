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
from csv import reader
from ...exceptions import EmptyMolecule


class EMOLRead:
    def __init__(self, log_buffer=None):
        self.__atoms = []
        self.__bonds = []
        self.__atom_map = {}
        self.__stereo = []
        if log_buffer is None:
            log_buffer = []
        self.__log_buffer = log_buffer

    def getvalue(self):
        if self.__in_mol or self.__in_mol is None:
            raise ValueError('molecule not complete')
        return {'atoms': self.__atoms, 'bonds': self.__bonds, 'stereo': self.__stereo}

    def __call__(self, line, lineu=None):
        if lineu is None:
            lineu = line.upper()
        if self.__in_mol:
            if lineu.startswith('M  V30 END CTAB'):
                self.__in_mol = False
                return True
            elif self.__atoms_count:
                if lineu.startswith('M  V30 END'):
                    x = lineu[11:].strip()
                    cp = self.__parser
                    self.__parser = None

                    if x == 'ATOM':
                        if cp == self.__atom_parser and len(self.__atoms) == self.__atoms_count:
                            return
                    elif x == 'BOND':
                        if cp == self.__bond_parser and len(self.__bonds) == self.__bonds_count:
                            return
                    else:
                        return
                    raise ValueError('invalid number of %s records or invalid CTAB' % x)

                elif self.__parser:
                    collected = self.__record_collector(line)
                    if collected:
                        self.__parser(collected)

                elif lineu.startswith('M  V30 BEGIN ATOM'):
                    self.__parser = self.__atom_parser
                elif lineu.startswith('M  V30 BEGIN BOND'):
                    self.__parser = self.__bond_parser
                elif lineu.startswith(('M  V30 BEGIN SGROUP', 'M  V30 BEGIN COLLECTION')):
                    self.__parser = self.__ignored_block_parser
                else:
                    raise ValueError('invalid CTAB')

            else:  # M  V30 COUNTS line expected
                a, b, *_ = line[13:].split()
                atom_count = int(a)
                if not atom_count:
                    raise EmptyMolecule
                self.__bonds_count = int(b)
                self.__atoms_count = atom_count

        elif self.__in_mol is not None:
            raise SyntaxError('invalid usage')
        elif not lineu.startswith('M  V30 BEGIN CTAB'):
            raise ValueError('invalid CTAB')
        else:
            self.__in_mol = True

    def __record_collector(self, line):
        if not line.endswith('-\n'):
            line = line[7:]
            if self.__record:
                line = self.__record + line
                self.__record = None

            return next(reader([line], delimiter=' ', quotechar='"', skipinitialspace=True))

        line = line[7:-2]
        if not self.__record:
            self.__record = line
        else:
            self.__record += line

    def __bond_parser(self, line):
        _, t, a1, a2, *kvs = line
        try:
            self.__bonds.append((self.__atom_map[a1], self.__atom_map[a2], int(t)))
        except KeyError:
            raise ValueError('invalid atoms numbers')
        for kv in kvs:
            k, v = kv.split('=')
            if k in ('CFG', 'cfg', 'Cfg'):
                if v == '1':
                    self.__stereo.append((self.__atom_map[a1], self.__atom_map[a2], 1))
                elif v == '3':
                    self.__stereo.append((self.__atom_map[a1], self.__atom_map[a2], -1))
                else:
                    self.__log_buffer.append('invalid or unsupported stereo')
                break

    def __atom_parser(self, line):
        n, a, x, y, z, m, *kvs = line
        i = None
        c = 0
        r = False
        for kv in kvs:
            k, v = kv.split('=', 1)
            k = k.upper()
            if k == 'CHG':
                c = int(v)
            elif k == 'MASS':
                i = int(v)
            elif k == 'RAD':
                r = True

        self.__atom_map[n] = n = len(self.__atoms)
        if a.startswith(('[', 'NOT', 'not', 'Not')):
            raise ValueError('list of atoms not supported')
        elif a == 'D':
            if i:
                raise ValueError('isotope on deuterium atom')
            a = 'H'
            i = 2

        self.__atoms.append({'element': a, 'isotope': i, 'charge': c, 'is_radical': r,
                             'x': float(x), 'y': float(y), 'z': float(z), 'mapping': int(m)})

    def __ignored_block_parser(self, line):
        return

    __record = __atoms_count = __in_mol = __parser = None


__all__ = ['EMOLRead']
