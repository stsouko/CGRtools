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
from logging import warning
from subprocess import check_output
from sys import platform
from traceback import format_exc
from ._CGRrw import WithMixin, CGRread, CGRwrite
from ._MDLrw import MOLwrite, MOLread, EMOLread, prepare_meta


class SDFread(CGRread, WithMixin):
    """
    MDL SDF files reader. works similar to opened file object. support `with` context manager.
    on initialization accept opened in text mode file, string path to file,
    pathlib.Path object or another buffered reader object
    """
    def __init__(self, file, *args, indexable=False, **kwargs):
        super().__init__(*args, **kwargs)
        super(CGRread, self).__init__(file)
        self.__data = self.__reader()

        if indexable and platform != 'win32' and not self._is_buffer:
            self.__file = iter(self._file.readline, '')
            self.__shifts = [0]
            for x in check_output(['grep', '-bE', r'\$\$\$\$', self._file.name]).decode().split():
                _pos, _len = x.split(':', 1)
                self.__shifts.append(int(_pos) + len(_len) + 1)
            print(self.__shifts)
        else:
            self.__file = self._file

    def read(self):
        """
        parse whole file

        :return: list of parsed molecules
        """
        return list(self.__data)

    def __iter__(self):
        return self.__data

    def __next__(self):
        return next(self.__data)

    def __len__(self):
        if self.__shifts:
            return len(self.__shifts) - 1
        raise self.__implement_error

    def seek(self, offset):
        if self.__shifts:
            if 0 <= offset < len(self.__shifts):
                current_pos = self._file.tell()
                print('current_pos', current_pos)
                new_pos = self.__shifts[offset]
                print('new_pos', new_pos)
                if current_pos != new_pos:
                    if current_pos == self.__shifts[-1]:  # reached the end of the file
                        self.__data = self.__reader()
                        self.__file = iter(self._file.readline, '')
                        self._file.seek(0)
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

            elif line.startswith("$$$$"):
                if record:
                    record['meta'] = prepare_meta(meta)
                    try:
                        seek = yield self._convert_structure(record)
                        if seek:
                            yield
                            self.__already_seeked = False
                            continue
                    except ValueError:
                        warning(f'record consist errors:\n{format_exc()}')
                    finally:
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
                        parser = MOLread(line)
                    elif 'V3000' in line:
                        parser = EMOLread()
                    else:
                        raise ValueError('invalid MOL entry')
                except ValueError:
                    failkey = True
                    warning(f'line:\n{line}\nconsist errors:\n{format_exc()}')

        if record:  # True for MOL file only.
            record['meta'] = prepare_meta(meta)
            try:
                yield self._convert_structure(record)
            except ValueError:
                warning(f'record consist errors:\n{format_exc()}')

    __already_seeked = False
    __shifts = None
    __implement_error = NotImplementedError('Indexable supported in unix-like o.s. and for files stored on disk')


class SDFwrite(MOLwrite, WithMixin):
    """
    MDL SDF files writer. works similar to opened for writing file object. support `with` context manager.
    on initialization accept opened for writing in text mode file, string path to file,
    pathlib.Path object or another buffered writer object
    """
    def __init__(self, file, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(CGRwrite, self).__init__(file, 'w')

    def write(self, data):
        """
        write single molecule into file
        """
        m = self._convert_structure(data)
        self._file.write(self._format_mol(*m))
        self._file.write('M  END\n')

        for k, v in data.meta.items():
            self._file.write(f'>  <{k}>\n{v}\n')
        self._file.write('$$$$\n')


__all__ = ['SDFread', 'SDFwrite']
