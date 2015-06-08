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
from itertools import count

toMDL = {-3: 7, -2: 6, -1: 5, 0: 0, 1: 3, 2: 2, 3: 1}
fromMDL = {0: 0, 1: 3, 2: 2, 3: 1, 4: 0, 5: -1, 6: -2, 7: -3}
bondlabels = {'n': 0, '0': None, '1': 1, '2': 2, '3': 3, 'a': 4, '4': 4, '9': 9, 's': 9}


class CGRRead:
    def __init__(self):
        self.__prop = {}

    def collect(self, line):
        if 'M  STY' in line:
            for i in range(int(line[8])):
                if 'DAT' in line[10 + 8 * i:17 + 8 * i]:
                    self.__prop[int(line[10 + 8 * i:13 + 8 * i])] = {}
        elif 'M  SAL' in line:
            if int(line[7:10]) in self.__prop:
                key = []
                for i in range(int(line[10:13])):
                    key.append(int(line[14 + 4 * i:17 + 4 * i]) - 1)
                self.__prop[int(line[7:10])]['atoms'] = sorted(key)
        elif 'M  SDT' in line and int(line[7:10]) in self.__prop:
            key = line.split()[-1].lower()
            if key not in self.__cgrkeys:
                self.__prop.pop(int(line[7:10]))
            else:
                self.__prop[int(line[7:10])]['type'] = key
        elif 'M  SED' in line and int(line[7:10]) in self.__prop:
            self.__prop[int(line[7:10])]['value'] = line[10:].strip().replace('/', '').lower()

    __cgrkeys = dict(dynatom=1, dynbond=2, dynatomstereo=1, dynbondstereo=2,
                     atomstereo=1, bondstereo=2, extrabond=2)

    def getdata(self):
        prop = []
        #todo: запилить проверки значений от мудаков.
        for i in self.__prop.values():
            if len(i['atoms']) == self.__cgrkeys[i['type']]:
                prop.append(i)

        self.__prop = {}
        return prop


class CGRWrite:
    def getformattedtext(self, data):
        text = []
        for j in count():
            sty = data['CGR_DAT'][j * 8:j * 8 + 8]
            if sty:
                stydat = ' '.join(['%3d DAT' % (x + 1 + j * 8) for x in range(len(sty))])
                text.append('M  STY  %d %s\n' % (len(sty), stydat))
            else:
                break
        for i, j in enumerate(data['CGR_DAT'], start=1):
            cx, cy = self.__getposition(j['atoms'], data['atoms'])
            text.append('M  SAL %3d%3d %s\n' % (i, len(j['atoms']), ' '.join(['%3d' % x for x in j['atoms']])))
            text.append('M  SDT %3d %s\n' % (i, j['type']))
            text.append('M  SDD %3d %10.4f%10.4f    DAU   ALL  0       0\n' % (i, cx, cy))
            text.append('M  SED %3d %s\n' % (i, j['value']))
        return ''.join(text)

    @staticmethod
    def __getposition(inp, atoms):
        cord = []
        for i in inp:
            cord.append(atoms[i - 1])
        if len(cord) > 1:
            x = (cord[-1]['x'] + cord[0]['x']) / 2 + .2
            y = (cord[-1]['y'] + cord[0]['y']) / 2
            dy = cord[-1]['y'] - cord[0]['y']
            dx = cord[-1]['x'] - cord[0]['x']
            if dx > 0:
                if dy > 0:
                    y -= .2
                else:
                    y += .2
            elif dx < 0:
                if dy < 0:
                    y -= .2
                else:
                    y += .2
        else:
            x, y = cord[0]['x'] + .25, cord[0]['y']

        return x, y
