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
from traceback import format_exc
from warnings import warn
from ._CGRrw import WithMixin, CGRread, CGRwrite
from ._MDLrw import MOLwrite, MOLread, EMOLread
from ..exceptions import EmptyMolecule


class SDFread(CGRread, WithMixin):
    def __init__(self, file, *args, is_template=None, **kwargs):
        assert not is_template, 'is_tepmlate works only for reactions'
        super().__init__(*args, is_template=False, **kwargs)
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
        mkey = parser = molecule = None
        for line in self._file:
            if failkey and not line.startswith("$$$$"):
                continue
            elif parser:
                try:
                    if parser(line):
                        molecule = parser.getvalue()
                        parser = None
                except (EmptyMolecule, ValueError):
                    failkey = True
                    parser = None
                    warn('line: \n%s\nconsist errors:\n%s' % (line, format_exc()), ResourceWarning)

            elif line.startswith("$$$$"):
                if molecule:
                    try:
                        yield self._get_molecule(molecule)
                    except Exception:
                        warn('record consist errors:\n%s' % format_exc(), ResourceWarning)
                    molecule = None

                im = 3
                failkey = False
                mkey = None

            elif molecule:
                if line.startswith('>  <'):
                    mkey = line.rstrip()[4:-1].strip()
                    if not mkey:
                        continue
                    if mkey in ('PHTYP', 'FFTYP', 'PCTYP', 'EPTYP', 'HBONDCHG', 'CNECHG',
                                'dynPHTYP', 'dynFFTYP', 'dynPCTYP', 'dynEPTYP', 'dynHBONDCHG', 'dynCNECHG'):
                        target = 'colors'
                    else:
                        target = 'meta'

                    molecule[target][mkey] = []
                elif mkey:
                    data = line.strip()
                    if data:
                        molecule[target][mkey].append(data)

            elif im:
                im -= 1
            elif not im:
                try:
                    if 'V2000' in line:
                        parser = MOLread(line)
                    elif 'V3000' in line:
                        parser = EMOLread(line)
                    else:
                        raise ValueError('invalid MOL')
                except (EmptyMolecule, ValueError):
                    failkey = True
                    warn('line: \n%s\nconsist errors:\n%s' % (line, format_exc()), ResourceWarning)

        if molecule:  # True for MOL file only.
            try:
                yield self._get_molecule(molecule)
            except Exception:
                warn('record consist errors:\n%s' % format_exc(), ResourceWarning)


class SDFwrite(MOLwrite, WithMixin):
    def __init__(self, file, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(CGRwrite, self).__init__(file, 'w')
        self.write = self.__write

    def __write(self, data):
        m = self.get_formatted_cgr(data)
        self._file.write(m['CGR'])
        self._file.write("M  END\n")

        for i in chain(m['colors'].items(), m['meta'].items()):
            self._file.write(">  <%s>\n%s\n" % i)
        self._file.write("$$$$\n")


__all__ = [SDFread.__name__, SDFwrite.__name__]
