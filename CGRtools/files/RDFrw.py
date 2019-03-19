# -*- coding: utf-8 -*-
#
#  Copyright 2014-2019 Ramil Nugmanov <stsouko@live.ru>
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
from bisect import bisect_right, bisect_left
from collections import defaultdict
from os.path import getsize
from itertools import chain
from logging import warning
from subprocess import check_output
from sys import platform
from time import strftime
from traceback import format_exc
from _CGRrw import WithMixin, CGRread, CGRwrite
from _MDLrw import MOLwrite, MOLread, EMOLread, RXNread, ERXNread, prepare_meta
from CGRtools.containers.common import BaseContainer


class RDFread(CGRread, WithMixin):
    def __init__(self, file, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(CGRread, self).__init__(file)
        self.__data = self.__reader()
        self._first = next(self.__data)
        if platform == 'linux':
            self.__gen = True
        if not self._is_buffer and self._first and self.__gen:
            self._size = getsize(self._file.name)
            self._shifts = [int(x.split(':', 1)[0])
                            for x in check_output(["grep", "-bE", "\$[RM]FMT", self._file.name]).decode().split()]

    def read(self):
        return list(self.__data)

    def __iter__(self):
        return self.__data

    def __len__(self):
        if not self._is_buffer:
            if self._first and self.__gen:
                return len(self._shifts)
            else:
                raise Exception
        else:
            raise NotImplementedError

    def __next__(self):
        return next(self.__data)

    def __getitem__(self, item):
        pass

    def seek(self, offset):
        if 0 <= offset < len(self._shifts):
            self._file.seek(self._shifts[offset])
            # self.__data.send(self._shifts[offset])
        else:
            raise IndexError('invalid offset')

    def tell(self):
        t = self._file.tell()
        if t == 0:
            return 0
        else:
            return bisect_left(self._shifts, t)

    def __reader(self):
        record = parser = mkey = None
        failed = False
        if not self._is_buffer:
            if next(self._file).startswith('$RXN'):
                is_reaction = True
                ir = 3
                meta = defaultdict(list)
                yield False
            elif next(self._file).startswith('$DATM'):  # skip header
                ir = 0
                is_reaction = meta = None
                yield True
            else:
                raise Exception('Not valid file')
        else:
            ir = 0
            is_reaction = meta = None

        for line in self._file:
            if failed and not line.startswith(('$RFMT', '$MFMT')):
                continue
            elif parser:
                try:
                    if parser(line):
                        record = parser.getvalue()
                        parser = None
                except ValueError:
                    failed = True
                    parser = None
                    warning(f'line:\n{line}\nconsist errors:\n{format_exc()}')
            elif line.startswith('$RFMT'):
                if record:
                    record['meta'] = prepare_meta(meta)
                    try:
                        seek = yield self._convert_reaction(record) if is_reaction else self._convert_structure(record)
                        if seek:
                            yield
                            continue
                    except Exception:
                        warning(f'record consist errors:\n{format_exc()}')
                    record = None
                is_reaction = True
                ir = 4
                failed = False
                mkey = None
                meta = defaultdict(list)
            elif line.startswith('$MFMT'):
                if record:
                    record['meta'] = prepare_meta(meta)
                    try:
                        seek = yield self._convert_reaction(record) if is_reaction else self._convert_structure(record)
                        if seek:
                            yield
                            continue
                    except ValueError:
                        warning(f'record consist errors:\n{format_exc()}')
                    record = None
                ir = 3
                failed = is_reaction = False
                mkey = None
                meta = defaultdict(list)
            elif record:
                if line.startswith('$DTYPE'):
                    mkey = line[7:].strip()
                    if not mkey:
                        warning(f'invalid metadata entry: {line}')
                elif mkey:
                    data = line.lstrip("$DATUM").strip()
                    if data:
                        meta[mkey].append(data)
            elif ir:
                ir -= 1
            elif not ir:
                try:
                    if is_reaction:
                        if line.startswith('M  V30 COUNTS'):
                            parser = ERXNread(line, self._ignore)
                        else:
                            parser = RXNread(line, self._ignore)
                    else:
                        if 'V2000' in line:
                            parser = MOLread(line)
                        elif 'V3000' in line:
                            parser = EMOLread()
                        else:
                            raise ValueError('invalid MOL entry')
                except ValueError:
                    failed = True
                    warning(f'line:\n{line}\nconsist errors:\n{format_exc()}')
        if record:
            record['meta'] = prepare_meta(meta)
            try:
                yield self._convert_reaction(record) if is_reaction else self._convert_structure(record)
            except ValueError:
                warning(f'record consist errors:\n{format_exc()}')
    __gen = False


cgr = RDFread('example.rdf')
# print('1', next(cgr))
# print('2', next(cgr))
# print('t2', cgr.tell())
cgr.seek(4)
print('3', next(cgr))
print('len', len(cgr))
print('t3', cgr.tell())


class RDFwrite(MOLwrite, WithMixin):
    def __init__(self, file, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(CGRwrite, self).__init__(file, 'w')

    def write(self, data):
        self._file.write(strftime('$RDFILE 1\n$DATM    %m/%d/%y %H:%M\n'))
        self.__write(data)
        self.write = self.__write

    def __write(self, data):
        if isinstance(data, BaseContainer):
            m = self._convert_structure(data)
            self._file.write('$MFMT\n')
            self._file.write(self._format_mol(*m['structure']))
            self._file.write('M  END\n')
        else:
            self._file.write('$RFMT\n$RXN\n\n\n\n'
                             f'{len(data.reactants):3d}{len(data.products):3d}{len(data.reagents):3d}\n')
            for m in chain(data.reactants, data.products, data.reagents):
                m = self._convert_structure(m)
                self._file.write('$MOL\n')
                self._file.write(self._format_mol(*m))
                self._file.write('M  END\n')

        for k, v in data.meta.items():
            self._file.write(f'$DTYPE {k}\n$DATUM {v}\n')


__all__ = ['RDFread', 'RDFwrite']
