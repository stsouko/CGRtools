# -*- coding: utf-8 -*-
#
#  Copyright 2014-2018 Ramil Nugmanov <stsouko@live.ru>
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
from itertools import chain
from logging import warning
from time import strftime
from traceback import format_exc
from ._CGRrw import WithMixin, CGRread, CGRwrite
from ._MDLrw import MOLwrite, MOLread, EMOLread, RXNread, ERXNread, prepare_meta
from ..containers.common import BaseContainer


class RDFread(CGRread, WithMixin):
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
        ir = 0
        isreaction = record = parser = mkey = meta = None
        failkey = True  # skip header
        for line in self._file:
            if failkey and not line.startswith(('$RFMT', '$MFMT', '$RXN')):
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

            elif line.startswith('$RFMT'):
                if record:
                    record['meta'] = prepare_meta(meta)
                    try:
                        yield self._convert_reaction(record) if isreaction else self._convert_structure(record)
                    except Exception:
                        warning(f'record consist errors:\n{format_exc()}')
                    record = None
                isreaction = True
                ir = 4
                failkey = False
                mkey = None
                meta = defaultdict(list)
            elif line.startswith('$MFMT'):
                if record:
                    record['meta'] = prepare_meta(meta)
                    try:
                        yield self._convert_reaction(record) if isreaction else self._convert_structure(record)
                    except ValueError:
                        warning(f'record consist errors:\n{format_exc()}')
                    record = None
                ir = 3
                failkey = isreaction = False
                mkey = None
                meta = defaultdict(list)
            elif line.startswith('$RXN'):  # parse RXN file
                isreaction = True
                ir = 3
                failkey = False
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
                    if isreaction:
                        if line.startswith('M  V30 COUNTS'):
                            parser = ERXNread(line, self._ignore)
                        else:
                            parser = RXNread(line, self._ignore)
                    else:
                        if 'V2000' in line:
                            parser = MOLread(line, self._ignore)
                        elif 'V3000' in line:
                            parser = EMOLread(line, self._ignore)
                        else:
                            raise ValueError('invalid MOL entry')
                except ValueError:
                    failkey = True
                    warning(f'line:\n{line}\nconsist errors:\n{format_exc()}')
        if record:
            record['meta'] = prepare_meta(meta)
            try:
                yield self._convert_reaction(record) if isreaction else self._convert_structure(record)
            except ValueError:
                warning(f'record consist errors:\n{format_exc()}')


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
            colors = m['colors']
        else:
            self._file.write(f'$RFMT\n$RXN\n\n\n\n{len(data.reagents):3d}{len(data.products):3d}\n')
            colors = {}
            s = 0
            rl = len(data.reagents)
            for n, m in enumerate(chain(data.reagents, data.products), start=1):
                m = self._convert_structure(m, s)
                if self._fix_position:
                    s = m['max_x'] + (3 if n == rl else 1)
                self._file.write('$MOL\n')
                self._file.write(self._format_mol(*m['structure']))
                self._file.write('M  END\n')

        for k, v in chain(colors.items(), data.meta.items()):
            self._file.write(f'$DTYPE {k}\n$DATUM {v}\n')


__all__ = ['RDFread', 'RDFwrite']
