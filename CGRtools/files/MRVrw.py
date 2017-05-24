# -*- coding: utf-8 -*-
#
#  Copyright 2017 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
from itertools import chain, repeat, count
from sys import stderr
from traceback import format_exc
from .CGRrw import CGRread, CGRwrite, fromMDL, EmptyMolecule, FinalizedFile
from ..containers import MoleculeContainer


class MRVread(CGRread):
    def __init__(self, file, remap=True):
        self.__file = file
        self.__data = self.__reader()
        CGRread.__init__(self, remap)
        raise Exception('NOT IMPLEMENTED')

    def read(self):
        return list(self.__data)

    def __iter__(self):
        return self.__data

    def __next__(self):
        return next(self.__data)

    def __reader(self):
        for n, structure in enumerate(magic(self.__MRWfile)):
            isreaction = True
            reaction = {'substrats': [], 'products': [], 'meta': {}, 'colors': {}}
            molecule = {'atoms': [], 'bonds': [], 'CGR_DAT': {}}
            try:
                if not atoms:
                    raise EmptyMolecule('Molecule without atoms')
                atomcount = atoms + im
                bondcount = int(line[3:6]) + atomcount
            except (EmptyMolecule, ValueError):
                print('line %d\n\n%s\n consist errors: %s' % (n, line, format_exc()), file=stderr)

            molecule['bonds'].append((int(line[:3]), int(line[3:6]), int(line[6:9])))
            molecule['CGR_DAT'] = self._get_collected()
            molecule['atoms'].append(dict(element=line[31:34].strip(), isotope=int(line[34:36]),
                                          charge=fromMDL[int(line[38:39])],
                                          map=int(line[60:63]), mark=line[54:57].strip(),
                                          x=float(line[:10]), y=float(line[10:20]), z=float(line[20:30])))
            if not mend:
                self._collect(line)
            elif line.startswith('$DTYPE'):
                mkey = line[7:].strip()
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
            try:
                yield self._get_reaction(reaction, stereo=next(self.__stereo)) if isreaction \
                    else self._get_molecule(reaction, stereo=next(self.__stereo))
            except:
                print('line %d\n previous record consist errors: %s' % (n, format_exc()), file=stderr)


class MRVwrite(CGRwrite):
    def __init__(self, file, extralabels=False, mark_to_map=False):
        CGRwrite.__init__(self, extralabels=extralabels, mark_to_map=mark_to_map, _format='mrv')
        self.__file = file
        self.write = self.__init_write

    __finalized = False

    def close(self):
        if not self.__finalized:
            self.finalize()
        self.__file.close()

    def finalize(self):
        self.__file.write('</cml>')
        self.write = self.__write_adhoc
        self.__finalized = True

    @staticmethod
    def __write_adhoc(_):
        raise FinalizedFile('Writer closed')

    def __init_write(self, data):
        self.__file.write('<cml>')
        self.__write(data)
        self.write = self.__write

    def __write(self, data):
        self.__file.write('<MDocument><MChemicalStruct>')

        if isinstance(data, MoleculeContainer):
            m = self.get_formatted_cgr(data)
            self.__file.write('<molecule><propertyList>')
            for k, v in chain(m['colors'].items(), data.meta.items()):
                if '\n' in v:
                    v = '<![CDATA[%s]]>' % v
                self.__file.write('<property title="%s"><scalar>%s</scalar></property>' % (k, v))

            self.__file.write('</propertyList>')
            self.__file.write(m['CGR'])
            self.__file.write('</molecule>')
        else:
            colors = {}
            c = count(1)
            self.__file.write('<reaction>')
            for i, j in (('substrats', 'reactantList'), ('products', 'productList')):
                self.__file.write('<%s>' % j)
                for cnext, m in zip(c, data[i]):
                    m = self.get_formatted_cgr(m)
                    self.__file.write('<molecule>')
                    self.__file.write(m['CGR'])
                    self.__file.write('</molecule>')
                    colors.update({'%s.%d' % (k, cnext): v for k, v in m['colors'].items()})
                self.__file.write('</%s>' % j)

            self.__file.write('<propertyList>')
            for k, v in chain(colors.items(), data.meta.items()):
                if '\n' in v:
                    v = '<![CDATA[%s]]>' % v
                    self.__file.write('<property title="%s"><scalar>%s</scalar></property>' % (k, v))

            self.__file.write('</propertyList></reaction>')

        self.__file.write('</MChemicalStruct></MDocument>')
