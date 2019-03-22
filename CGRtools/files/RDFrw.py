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
from itertools import chain, islice
from logging import warning
from os.path import getsize
from subprocess import check_output
from sys import platform
from time import strftime
from traceback import format_exc
from ._CGRrw import WithMixin, CGRread, CGRwrite
from ._MDLrw import MOLwrite, MOLread, EMOLread, RXNread, ERXNread, prepare_meta
from ..containers.common import BaseContainer
from ..exceptions import InvalidFileType


class RDFread(CGRread, WithMixin):
    """
    MDL RDF files reader. works similar to opened file object. support `with` context manager.
    on initialization accept opened in text mode file, string path to file,
    pathlib.Path object or another buffered reader object
    """
    def __init__(self, file, *args, indexable=False, **kwargs):
        super().__init__(*args, **kwargs)
        super(CGRread, self).__init__(file)
        self.__data = self.__reader()

        if indexable and platform != 'win32' and not self._is_buffer:
            self.__file = iter(self._file.readline, '')
            if next(self.__data):
                self.__shifts = [int(x.split(':', 1)[0]) for x in
                                 check_output(["grep", "-boE", "^\$[RM]FMT", self._file.name]).decode().split()]
                self.__shifts.append(getsize(self._file.name))
        else:
            self.__file = self._file
            next(self.__data)

    def read(self):
        """
        parse whole file

        :return: list of parsed molecules or reactions
        """
        return list(self)

    def __iter__(self):
        return (x for x in self.__data if x is not None)

    def __len__(self):
        if self.__shifts:
            return len(self.__shifts) - 1
        raise self.__implement_error

    def __next__(self):
        return next(iter(self))

    def seek(self, offset):
        if self.__shifts:
            if 0 <= offset < len(self.__shifts):
                current_pos = self._file.tell()
                new_pos = self.__shifts[offset]
                if current_pos != new_pos:
                    if current_pos == self.__shifts[-1]:  # reached the end of the file
                        self.__data = self.__reader()
                        self.__file = iter(self._file.readline, '')
                        self._file.seek(0)
                        next(self.__data)
                        if offset:  # move not to the beginning of the file
                            self._file.seek(new_pos)
                    else:
                        if not self.__already_seeked:
                            if self.__shifts[0] < current_pos:  # in the middle of the file
                                self.__data.send(True)
                            self.__already_seeked = True
                        self._file.seek(new_pos)
            else:
                raise IndexError('invalid offset')
        else:
            raise self.__implement_error

    def tell(self):
        if self.__shifts:
            t = self._file.tell()
            if t == self.__shifts[0]:
                return 0
            elif t == self.__shifts[-1]:
                return len(self.__shifts) - 1
            elif t in self.__shifts:
                return bisect_left(self.__shifts, t)
            else:
                return bisect_left(self.__shifts, t) - 1
        raise self.__implement_error

    def __getitem__(self, item):
        """
        getting the item by index from the original file,
        if the required block of the file with an error,
        then only the correct blocks are returned
        :param item: int or slice
        :return: ReactionContainer or list of ReactionContainers
        """
        if self.__shifts:
            _len = len(self.__shifts) - 1
            _current_pos = self.tell()

            if isinstance(item, int):
                if item >= _len or item < -_len:
                    raise IndexError('List index out of range')
                if item < 0:
                    item += _len
                self.seek(item)
                records = next(self.__data)
            elif isinstance(item, slice):
                start, stop, step = item.indices(_len)
                if start == stop:
                    return []

                if step == 1:
                    self.seek(start)
                    records = [x for x in islice(self.__data, 0, stop - start) if x is not None]
                else:
                    records = []
                    for index in range(start, stop, step):
                        self.seek(index)
                        record = next(self.__data)
                        if record:
                            records.append(record)
            else:
                raise TypeError('Indices must be integers or slices')

            self.seek(_current_pos)
            if records is None:
                raise self.__index_error
            return records
        raise self.__implement_error

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

    __shifts = None
    __implement_error = NotImplementedError('Indexable supported in unix-like o.s. and for files stored on disk')
    __already_seeked = False
    __index_error = IndexError('Data block with requested index contain errors')


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
