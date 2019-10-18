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
from bisect import bisect_left
from collections import defaultdict
from io import BytesIO
from logging import warning
from subprocess import check_output
from sys import platform
from traceback import format_exc
from warnings import warn
from ._MDLrw import MDLRead, MDLWrite, MOLRead, EMOLRead


class SDFRead(MDLRead):
    """
    MDL SDF files reader. works similar to opened file object. support `with` context manager.
    on initialization accept opened in text mode file, string path to file,
    pathlib.Path object or another buffered reader object
    """
    def __init__(self, *args, indexable=False, **kwargs):
        """
        :param indexable: if True:
            supported methods seek, tell, object size and subscription, it only works when dealing with a real file
            (the path to the file is specified) because the external grep utility is used, supporting in unix-like OS
            the object behaves like a normal open file
                        if False:
            works like generator converting a record into MoleculeContainer and returning each object in order,
            records with errors are skipped
        """
        super().__init__(*args, **kwargs)
        self._data = self.__reader()

        if indexable and platform != 'win32' and not self._is_buffer:
            self.__file = iter(self._file.readline, '')
            self._shifts = self._load_cache()
            if self._shifts is None:
                self._shifts = [0]
                for x in BytesIO(check_output(['grep', '-bE', r'\$\$\$\$', self._file.name])):
                    _pos, _line = x.split(b':', 1)
                    self._shifts.append(int(_pos) + len(_line))
                self._dump_cache(self._shifts)
        else:
            self.__file = self._file

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
            return bisect_left(self._shifts, t)
        raise self._implement_error

    def __reader(self):
        im = 3
        failkey = False
        mkey = parser = record = None
        meta = defaultdict(list)
        for line in self.__file:
            if failkey and not line.startswith("$$$$"):
                continue
            elif parser:
                try:
                    if parser(line):
                        record = parser.getvalue()
                        parser = None
                except ValueError:
                    failkey = True
                    parser = None
                    warning(f'line:\n{line}\nconsist errors:\n{format_exc()}')
                    yield None

            elif line.startswith("$$$$"):
                if record:
                    record['meta'] = self._prepare_meta(meta)
                    try:
                        yield self._convert_structure(record)
                    except ValueError:
                        warning(f'record consist errors:\n{format_exc()}')
                        yield None
                    record = None

                im = 3
                failkey = False
                mkey = None
                meta = defaultdict(list)
            elif record:
                if line.startswith('>  <'):
                    mkey = line.rstrip()[4:-1].strip()
                    if not mkey:
                        warning(f'invalid metadata entry: {line}')
                elif mkey:
                    data = line.strip()
                    if data:
                        meta[mkey].append(data)
            elif im:
                im -= 1
            elif not im:
                try:
                    if 'V2000' in line:
                        parser = MOLRead(line)
                    elif 'V3000' in line:
                        parser = EMOLRead()
                    else:
                        raise ValueError('invalid MOL entry')
                except ValueError:
                    failkey = True
                    warning(f'line:\n{line}\nconsist errors:\n{format_exc()}')
                    yield None

        if record:  # True for MOL file only.
            record['meta'] = self._prepare_meta(meta)
            try:
                yield self._convert_structure(record)
            except ValueError:
                warning(f'record consist errors:\n{format_exc()}')
                yield None


class SDFWrite(MDLWrite):
    """
    MDL SDF files writer. works similar to opened for writing file object. support `with` context manager.
    on initialization accept opened for writing in text mode file, string path to file,
    pathlib.Path object or another buffered writer object
    """
    def write(self, data):
        """
        write single molecule into file
        """
        self._file.write(self._convert_structure(data))

        for k, v in data.meta.items():
            self._file.write(f'>  <{k}>\n{v}\n')
        self._file.write('$$$$\n')


class SDFread:
    def __init__(self, *args, **kwargs):
        warn('SDFread deprecated. Use SDFRead instead', DeprecationWarning)
        warning('SDFread deprecated. Use SDFRead instead')
        self.__obj = SDFRead(*args, **kwargs)

    def __getattr__(self, item):
        return getattr(self.__obj, item)

    def __iter__(self):
        return iter(self.__obj)

    def __next__(self):
        return next(self.__obj)

    def __getitem__(self, item):
        return self.__obj[item]

    def __enter__(self):
        return self.__obj.__enter__()

    def __exit__(self, _type, value, traceback):
        return self.__obj.__exit__(_type, value, traceback)

    def __len__(self):
        return len(self.__obj)


class SDFwrite:
    def __init__(self, *args, **kwargs):
        warn('SDFwrite deprecated. Use SDFWrite instead', DeprecationWarning)
        warning('SDFwrite deprecated. Use SDFWrite instead')
        self.__obj = SDFWrite(*args, **kwargs)

    def __getattr__(self, item):
        return getattr(self.__obj, item)

    def __enter__(self):
        return self.__obj.__enter__()

    def __exit__(self, _type, value, traceback):
        return self.__obj.__exit__(_type, value, traceback)


__all__ = ['SDFRead', 'SDFWrite', 'SDFread', 'SDFwrite']
