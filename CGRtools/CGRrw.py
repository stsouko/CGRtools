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
        self.__mendeleytable = self.__mendeley
        self.__mendeleyset = set(self.__mendeleytable)

    def collect(self, line):
        if 'M  ALS' in line:
            self.__prop[line[3:10]] = dict(atoms=[int(line[7:10])],
                                           type='atomlist' if line[14] == 'F' else 'atomnotlist',
                                           value=[line[16 + x*4: 20 + x*4].strip() for x in range(int(line[10:13]))])
        elif 'M  STY' in line:
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
                     atomstereo=1, bondstereo=2, extrabond=2, atomnotlist=1, atomlist=1)

    def getdata(self):
        prop = []
        for i in self.__prop.values():
            if len(i['atoms']) == self.__cgrkeys[i['type']]:
                prop.append(i)

        self.__prop = {}
        return prop

    def cgr_dat(self, g, k, atom1, atom2):
        if k['type'] == 'dynatomstereo':
            g.node[atom1]['s_stereo'], *_, g.node[atom1]['p_stereo'] = k['value'].split('>')

        elif k['type'] == 'atomstereo':
            g.node[atom1]['s_stereo'] = g.node[atom1]['p_stereo'] = k['value']

        elif k['type'] == 'dynatom':
            key = k['value'][0]
            diff = int(k['value'][1:])
            if key == 'c':  # update atom charges from CGR
                s_charge = fromMDL.get(g.node[atom1]['s_charge'], 0)
                g.node[atom1]['p_charge'] = toMDL.get(s_charge + diff, 0)

            elif key == '*':
                pass  # not implemented

        elif k['type'] == 'dynbond':
            val = k['value'].split('>')
            g.edge[atom1][atom2]['s_bond'] = bondlabels.get(val[0])
            g.edge[atom1][atom2]['p_bond'] = bondlabels.get(val[-1])

        elif k['type'] == 'bondstereo':
            g.edge[atom1][atom2]['s_stereo'] = g.edge[atom1][atom2]['p_stereo'] = k['value']

        elif k['type'] == 'dynbondstereo':
            val = k['value'].split('>')
            g.edge[atom1][atom2]['s_stereo'], *_, g.edge[atom1][atom2]['p_stereo'] = val

        elif k['type'] == 'extrabond':
            g.edge[atom1][atom2]['s_bond'], g.edge[atom1][atom2]['p_bond'] = k['value']

        elif k['type'] == 'atomlist':
            g.node[atom1]['element'] = k['value']

        elif k['type'] == 'atomnotlist':
            g.node[atom1]['element'] = list(self.__mendeleyset.difference(k['value']))

    __mendeley = {j: i for i, j in enumerate(['H', 'He',
                                              'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
                                              'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
                                              'K', 'Ca',
                                                    'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn',
                                                          'Ga', 'Ge', 'As', 'Se', 'Br', 'Kr',
                                              'Rb', 'Sr',
                                                    'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd',
                                                         'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
                                              'Cs', 'Ba',
                                                    'La', 'Ce', 'Pr', 'Nd', 'Pm',
                                                    'Sm', 'Eu', 'Gd', 'Tb', 'Dy',
                                                    'Ho', 'Er', 'Tm', 'Yb', 'Lu',
                                                          'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg',
                                                         'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
                                              'Fr', 'Ra',
                                                    'Ac', 'Th', 'Pa', 'U', 'Np',
                                                    'Pu', 'Am', 'Cm', 'Bk', 'Cf',
                                                    'Es', 'Fm', 'Md', 'No', 'Lr',
                                                          'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn'])}


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
