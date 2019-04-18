# -*- coding: utf-8 -*-
#
#  Copyright 2014-2019 Ramil Nugmanov <stsouko@live.ru>
#  Copyright 2019 Dinar Batyrshin <batyrshin-dinar@mail.ru>
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
from bisect import bisect_left
from collections import defaultdict
from itertools import chain
from logging import warning
from os.path import getsize
from subprocess import check_output
from sys import platform
from time import strftime
from traceback import format_exc
from ._CGRrw import WithMixin, CGRread, CGRwrite
from ._MDLrw import MOLwrite, MOLread, MDLread, EMOLread, RXNread, ERXNread, prepare_meta
from ..containers.common import BaseContainer
from ..exceptions import InvalidFileType


class RDFread(CGRread, WithMixin, MDLread):
    """
    MDL RDF files reader. works similar to opened file object. support `with` context manager.
    on initialization accept opened in text mode file, string path to file,
    pathlib.Path object or another buffered reader object
    """
    def __init__(self, file, *args, indexable=False, **kwargs):
        """
        :param indexable: if True:
            supported methods seek, tell, object size and subscription, it only works when dealing with a real file
            (the path to the file is specified) because the external grep utility is used, supporting in unix-like OS
            the object behaves like a normal open file
                        if False:
            works like generator converting a record into ReactionContainer and returning each object in order,
            records with errors are skipped
        """
        super().__init__(*args, **kwargs)
        super(CGRread, self).__init__(file)
        self._data = self.__reader()

        if indexable and platform != 'win32' and not self._is_buffer:
            self.__file = iter(self._file.readline, '')
            if next(self._data):
                self._shifts = self._load_cache()
                if self._shifts is None:
                    self._shifts = [int(x.split(b':', 1)[0]) for x in
                                    check_output(['grep', '-boE', r'^\$[RM]FMT', self._file.name]).split()]
                    self._shifts.append(getsize(self._file.name))
                    self._dump_cache(self._shifts)
        else:
            self.__file = self._file
            next(self._data)

    def seek(self, offset):
        """
        shifts on a given number of record in the original file
        :param offset: number of record
        """
        if self._shifts:
            if 0 <= offset < len(self._shifts):
                current_pos = self._file.tell()
                new_pos = self._shifts[offset]
                if current_pos != new_pos:
                    if current_pos == self._shifts[-1]:  # reached the end of the file
                        self._data = self.__reader()
                        self.__file = iter(self._file.readline, '')
                        self._file.seek(0)
                        next(self._data)
                        if offset:  # move not to the beginning of the file
                            self._file.seek(new_pos)
                    else:
                        if not self.__already_seeked:
                            if self._shifts[0] < current_pos:  # in the middle of the file
                                self._data.send(True)
                            self.__already_seeked = True
                        self._file.seek(new_pos)
            else:
                raise IndexError('invalid offset')
        else:
            raise self._implement_error

    def tell(self):
        """
        :return: number of records processed from the original file
        """
        if self._shifts:
            t = self._file.tell()
            if t == self._shifts[0]:
                return 0
            elif t == self._shifts[-1]:
                return len(self._shifts) - 1
            elif t in self._shifts:
                return bisect_left(self._shifts, t)
            else:
                return bisect_left(self._shifts, t) - 1
        raise self._implement_error

    def __reader(self):
        record = parser = mkey = None
        failed = False

        if next(self.__file).startswith('$RXN'):  # parse RXN file
            is_reaction = True
            ir = 3
            meta = defaultdict(list)
            yield False
        elif next(self.__file).startswith('$DATM'):  # skip header
            ir = 0
            is_reaction = meta = None
            yield True
        else:
            raise InvalidFileType

        for line in self.__file:
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
                    yield None
            elif line.startswith('$RFMT'):
                if record:
                    record['meta'] = prepare_meta(meta)
                    try:
                        seek = yield self._convert_reaction(record) if is_reaction else self._convert_structure(record)
                    except ValueError:
                        warning(f'record consist errors:\n{format_exc()}')
                        seek = yield None

                    record = None
                    if seek:
                        yield
                        self.__already_seeked = False
                        continue

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
                    except ValueError:
                        warning(f'record consist errors:\n{format_exc()}')
                        seek = yield None

                    record = None
                    if seek:
                        yield
                        self.__already_seeked = False
                        continue

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
                    yield None
        if record:
            record['meta'] = prepare_meta(meta)
            try:
                yield self._convert_reaction(record) if is_reaction else self._convert_structure(record)
            except ValueError:
                warning(f'record consist errors:\n{format_exc()}')
                yield None

    __already_seeked = False


class RDFwrite(MOLwrite, WithMixin):
    """
    MDL RDF files writer. works similar to opened for writing file object. support `with` context manager.
    on initialization accept opened for writing in text mode file, string path to file,
    pathlib.Path object or another buffered writer object
    """
    def __init__(self, file, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(CGRwrite, self).__init__(file, 'w')

    def write(self, data):
        """
        write single molecule or reaction into file
        """
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
