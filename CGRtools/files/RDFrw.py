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
from itertools import chain
from time import strftime
from traceback import format_exc
from warnings import warn
from ._CGRrw import WithMixin, CGRread, CGRwrite
from ._MDLrw import MOLwrite, MOLread, EMOLread, RXNread, ERXNread
from ..containers import MoleculeContainer
from ..exceptions import EmptyMolecule


class RDFread(CGRread, WithMixin):
    def __init__(self, file, *args, ignore=False, **kwargs):
        super().__init__(*args, ignore=ignore, **kwargs)
        super(CGRread, self).__init__(file)
        self.__data = self.__reader()
        self.__ignore = ignore

    def read(self):
        return list(self.__data)

    def __iter__(self):
        return self.__data

    def __next__(self):
        return next(self.__data)

    def __reader(self):
        ir = 0
        isreaction = reaction = parser = mkey = None
        failkey = True  # skip header
        for line in self._file:
            if failkey and not line.startswith(('$RFMT', '$MFMT', '$RXN')):
                continue
            elif parser:
                try:
                    if parser(line):
                        reaction = parser.getvalue()
                        parser = None
                except EmptyMolecule:
                    if not (isreaction and self.__ignore):
                        failkey = True
                        parser = None
                    warn('line: \n%s\nconsist errors:\n%s' % (line, format_exc()), ResourceWarning)
                except ValueError:
                    failkey = True
                    parser = None
                    warn('line: \n%s\nconsist errors:\n%s' % (line, format_exc()), ResourceWarning)

            elif line.startswith('$RFMT'):
                if reaction:
                    try:
                        yield self._get_reaction(reaction) if isreaction else self._get_molecule(reaction)
                    except Exception:
                        warn('record consist errors:\n%s' % format_exc(), ResourceWarning)
                    reaction = None
                isreaction = True
                ir = 4
                failkey = False
                mkey = None
            elif line.startswith('$MFMT'):
                if reaction:
                    try:
                        yield self._get_reaction(reaction) if isreaction else self._get_molecule(reaction)
                    except Exception:
                        warn('record consist errors:\n%s' % format_exc(), ResourceWarning)
                    reaction = None
                ir = 3
                failkey = isreaction = False
                mkey = None
            elif line.startswith('$RXN'):  # parse RXN file
                isreaction = True
                ir = 3
                failkey = False
            elif reaction:
                if line.startswith('$DTYPE'):
                    mkey = line[7:].strip()
                    if not mkey:
                        continue
                    if mkey.split('.')[0] in ('PHTYP', 'FFTYP', 'PCTYP', 'EPTYP', 'HBONDCHG', 'CNECHG', 'dynPHTYP',
                                              'dynFFTYP', 'dynPCTYP', 'dynEPTYP', 'dynHBONDCHG', 'dynCNECHG'):
                        target = 'colors'
                    else:
                        target = 'meta'
                    reaction[target][mkey] = []
                elif mkey:
                    data = line.lstrip("$DATUM").strip()
                    if data:
                        reaction[target][mkey].append(data)

            elif ir:
                ir -= 1
            elif not ir:
                try:
                    if isreaction:
                        if line.startswith('M  V30 COUNTS'):
                            parser = ERXNread(line)
                        else:
                            parser = RXNread(line)
                    else:
                        if 'V2000' in line:
                            parser = MOLread(line)
                        elif 'V3000' in line:
                            parser = EMOLread(line)
                        else:
                            raise ValueError('invalid MOL')
                except (EmptyMolecule, ValueError):
                    failkey = True
                    warn('line: \n%s\nconsist errors:\n%s' % (line, format_exc()), ResourceWarning)

        if reaction:
            try:
                yield self._get_reaction(reaction) if isreaction else self._get_molecule(reaction)
            except Exception:
                warn('record consist errors:\n%s' % format_exc(), ResourceWarning)


class RDFwrite(MOLwrite, WithMixin):
    def __init__(self, file, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(CGRwrite, self).__init__(file, 'w')
        self.write = self.__init_write

    def __init_write(self, data):
        self._file.write(strftime("$RDFILE 1\n$DATM    %m/%d/%y %H:%M\n"))
        self.__write(data)
        self.write = self.__write

    def __write(self, data):
        if isinstance(data, MoleculeContainer):
            m = self.get_formatted_cgr(data)
            self._file.write('$MFMT\n')
            self._file.write(m['CGR'])
            self._file.write("M  END\n")
            colors = m['colors']
        else:
            self._file.write('$RFMT\n$RXN\n\n  CGRtools. (c) Dr. Ramil I. Nugmanov\n\n%3d%3d\n' %
                             (len(data.reagents), len(data.products)))
            colors = {}
            s = 0
            rl = len(data.reagents)
            for cnext, m in enumerate(chain(data.reagents + data.products), start=1):
                m = self.get_formatted_cgr(m, s)
                if self._fix_position:
                    s = m['max_x'] + (3 if cnext == rl else 1)
                self._file.write('$MOL\n')
                self._file.write(m['CGR'])
                self._file.write("M  END\n")
                colors.update({'%s.%d' % (k, cnext): v for k, v in m['colors'].items()})

        for p in chain(colors.items(), data.meta.items()):
            self._file.write('$DTYPE %s\n$DATUM %s\n' % p)


__all__ = [RDFread.__name__, RDFwrite.__name__]
