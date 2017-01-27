# -*- coding: utf-8 -*-
#
#  Copyright 2014-2016 Ramil Nugmanov <stsouko@live.ru>
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
import periodictable as pt
from itertools import count, chain, repeat
from .CGRrw import CGRread, CGRwrite, fromMDL
from . import MoleculeContainer


class SDFread(CGRread):
    def __init__(self, file, remap=True, stereo=None):
        self.__stereo = stereo and iter(stereo) or repeat(None)
        self.__SDFfile = file
        self.__data = self.__reader()
        CGRread.__init__(self, remap)

    def read(self):
        return list(self.__data)

    def __iter__(self):
        return self.__data

    def __next__(self):
        return next(self.__data)

    def __reader(self):
        im = 3
        atomcount = -1
        bondcount = -1
        failkey = False
        mkey = None
        molecule = None
        mend = False
        for n, line in enumerate(self.__SDFfile):
            if failkey and not line.startswith("$$$$"):
                continue
            elif line.startswith("$$$$"):
                if molecule:
                    try:
                        yield self.get_molecule(molecule, stereo=next(self.__stereo))
                    except:
                        pass

                mkey = None
                im = n + 4
                failkey = False
                mend = False
                molecule = None

            elif n == im:
                try:
                    atomcount = int(line[0:3]) + n
                    bondcount = int(line[3:6]) + atomcount
                    molecule = {'atoms': [], 'bonds': [], 'CGR_DAT': {}, 'meta': {}, 'colors': {}}
                except:
                    atomcount = bondcount = -1
                    failkey = True
                    molecule = None

            elif n <= atomcount:
                molecule['atoms'].append(dict(element=line[31:34].strip(), isotope=int(line[34:36]),
                                              charge=fromMDL.get(int(line[38:39]), 0),
                                              map=int(line[60:63]), mark=line[54:57].strip(),
                                              x=float(line[0:10]), y=float(line[10:20]), z=float(line[20:30])))

            elif n <= bondcount:
                try:
                    molecule['bonds'].append((int(line[0:3]), int(line[3:6]), int(line[6:9])))
                except:
                    failkey = True
                    molecule = None

            elif line.startswith("M  END"):
                mend = True
                molecule['CGR_DAT'] = self.getdata()

            elif molecule and n > bondcount:
                try:
                    if not mend:
                        self.collect(line)
                    elif line.startswith('>  <'):
                        mkey = line.strip()[4:-1]
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
                except:
                    failkey = True
                    molecule = None
        else:
            if molecule:  # True for MOL file only.
                try:
                    yield self.get_molecule(molecule, stereo=next(self.__stereo))
                except:
                    pass


class SDFwrite(CGRwrite):
    def __init__(self, output, extralabels=False, mark_to_map=False):
        CGRwrite.__init__(self, extralabels=extralabels, mark_to_map=mark_to_map)
        self.__file = output

    def close(self):
        self.__file.close()

    def write(self, data):
        m = self.getformattedcgr(data)
        self.__file.write(m['CGR'])
        self.__file.write("M  END\n")

        for i in chain(m['colors'].items(), m['meta'].items()):
            self.__file.write(">  <%s>\n%s\n" % i)
        self.__file.write("$$$$\n")
