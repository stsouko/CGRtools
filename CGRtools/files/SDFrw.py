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
from collections import defaultdict
from logging import warning
from traceback import format_exc
from ._CGRrw import WithMixin, CGRread, CGRwrite
from ._MDLrw import MOLwrite, MOLread, EMOLread, prepare_meta


class SDFread(CGRread, WithMixin):
    def __init__(self, file, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(CGRread, self).__init__(file)
        self.__data = self.__reader()

    def read(self):
        return list(self.__data)

    def __iter__(self):
        return self.__data

    def __next__(self):
        return next(self.__data)

    def __reader(self):
        im = 3
        failkey = False
        mkey = parser = record = None
        meta = defaultdict(list)
        for line in self._file:
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
                        yield self._convert_structure(record)
                    except ValueError:
                        warning(f'record consist errors:\n{format_exc()}')
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


class SDFwrite(MOLwrite, WithMixin):
    def __init__(self, file, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(CGRwrite, self).__init__(file, 'w')

    def write(self, data):
        m = self._convert_structure(data)
        self._file.write(self._format_mol(*m))
        self._file.write('M  END\n')

        for k, v in data.meta.items():
            self._file.write(f'>  <{k}>\n{v}\n')
        self._file.write('$$$$\n')


__all__ = ['SDFread', 'SDFwrite']
