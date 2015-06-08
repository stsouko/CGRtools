#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2014, 2015 Ramil Nugmanov <stsouko@live.ru>
# This file is part of FEAR (C) Ramil Nugmanov <stsouko@live.ru>.
#
#  fragger is free software; you can redistribute it and/or modify
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
from .CGRrw import CGRRead


class SDFread(CGRRead):
    def __init__(self, file):
        super().__init__()
        self.__SDFfile = file

    def readprop(self):
        with open(self.__SDFfile) as f:
            prop = set()
            for line in f:
                if '>  <' in line[:5]:
                    prop.add(line.strip()[4:-1].replace(',', '_'))
        return prop

    def readdata(self):
        '''парсер SDF файлов. возвращает пакет данных вида
        {'substrats':substrats, 'products':products, 'molecules':{'номер':{'atomlist':atomlist, 'bondmatrix':matrix}}}
        '''
        with open(self.__SDFfile) as f:
            im = 3
            atomcount = -1
            bondcount = -1
            failkey = False
            meta = None
            molecule = None
            mend = False
            for n, line in enumerate(f):
                if failkey and "$$$$" not in line[0:4]:
                    continue
                elif "$$$$" in line[0:4]:
                    if molecule:
                        yield molecule

                    meta = None
                    im = n + 4
                    failkey = False
                    mend = False

                elif n == im:
                    atoms, bonds, prop = [], [], dict(STEREO_DAT={})

                    try:
                        atomcount = int(line[0:3]) + n
                        bondcount = int(line[3:6]) + atomcount
                    except:
                        failkey = True
                        molecule = None

                elif n <= atomcount:
                    atoms.append(dict(element=line[31:34].strip(), isotop=line[34:36].strip(),
                                      charge=line[38:39].strip(), map=int(line[60:63]),
                                      mark=line[51:54].strip(),
                                      x=float(line[0:10]),
                                      y=float(line[10:20]),
                                      z=float(line[20:30])))

                elif n <= bondcount:
                    try:
                        bonds.append((int(line[0:3]) - 1, int(line[3:6]) - 1, line[6:9].strip(),
                                      line[9:12].strip(), line[15:18].strip()))
                    except:
                        failkey = True
                        molecule = None

                elif "M  END" in line:
                    mend = True
                    molecule = dict(bonds=bonds, atoms=atoms, prop=prop, cgr=self.getdata())

                elif n > bondcount and not mend:
                    try:
                        self.collect(line)
                    except:
                        failkey = True
                        molecule = None

                elif n > bondcount:
                    try:
                        if '>  <' in line:
                            meta = n
                            mkey = line.strip()[4:-1]
                            prop[mkey] = ''
                        elif meta:
                            prop[mkey] += line.strip()
                    except:
                        failkey = True
                        molecule = None
            else:
                raise StopIteration

