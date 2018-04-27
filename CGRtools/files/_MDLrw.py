# -*- coding: utf-8 -*-
#
#  Copyright 2017, 2018 Ramil Nugmanov <stsouko@live.ru>
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
from itertools import count, chain
from ._CGRrw import CGRread, CGRwrite, elements_set


class MOLread(CGRread):
    def __init__(self, *args, **kwargs):
        self.__prop = {}
        super().__init__(*args, **kwargs)

    def _collect(self, line):
        if line.startswith('M  ALS'):
            self.__prop[line[3:10]] = dict(atoms=(int(line[7:10]),),
                                           type='atomlist' if line[14] == 'F' else 'atomnotlist',
                                           value=[line[16 + x*4: 20 + x*4].strip() for x in range(int(line[10:13]))])

        elif line.startswith(('M  ISO', 'M  RAD', 'M  CHG')):
            _type = self.__ctf_data[line[3]]
            for i in range(int(line[6:9])):
                atom = int(line[10 + i * 8:13 + i * 8])
                self.__prop['%s_%d' % (_type, atom)] = dict(atoms=(atom,), type=_type,
                                                            value=line[14 + i * 8:17 + i * 8].strip())

        elif line.startswith('M  STY'):
            for i in range(int(line[8])):
                if 'DAT' in line[10 + 8 * i:17 + 8 * i]:
                    self.__prop[int(line[10 + 8 * i:13 + 8 * i])] = {}
        elif line.startswith('M  SAL') and int(line[7:10]) in self.__prop:
            self.__prop[int(line[7:10])]['atoms'] = tuple(int(line[14 + 4 * i:17 + 4 * i])
                                                          for i in range(int(line[10:13])))
        elif line.startswith('M  SDT') and int(line[7:10]) in self.__prop:
            key = line.split()[-1].lower()
            if key not in self._cgr_keys:
                self.__prop.pop(int(line[7:10]))
            else:
                self.__prop[int(line[7:10])]['type'] = key
        elif line.startswith('M  SED') and int(line[7:10]) in self.__prop:
            self.__prop[int(line[7:10])]['value'] = line[10:].strip().replace('/', '').lower()

    def _flush_collected(self):
        self.__prop.clear()

    def _get_collected(self):
        prop = []
        for i in self.__prop.values():
            if len(i['atoms']) == self._cgr_keys[i['type']]:
                prop.append(i)

        self.__prop.clear()
        return prop

    __ctf_data = {'R': 'radical', 'C': 'charge', 'I': 'isotope'}


class MOLwrite(CGRwrite):
    @classmethod
    def _format_mol(cls, atoms, bonds, extended, cgr_dat):
        mol_prop = []
        for i in extended:
            it, iv, ia = i['type'], i['value'], i['atom']
            if it == 'isotope':
                mol_prop.append('M  ISO  1 %3d %3d\n' % (ia, iv))
            elif it == 'atomlist':
                atomslist, _type = (elements_set.difference(iv), 'T') if len(iv) > cls._half_table else (iv, 'F')
                mol_prop.append('M  ALS %3d%3d %s %s\n' % (ia, len(atomslist), _type,
                                                           ''.join('%-4s' % x for x in atomslist)))
            elif it == 'radical':
                mol_prop.append('M  RAD  1 %3d %3d\n' % (ia, iv))

        for j in count():
            sty = len(cgr_dat[j * 8:j * 8 + 8])
            if sty:
                stydat = ' '.join(['%3d DAT' % (x + j * 8) for x in range(1, 1 + sty)])
                mol_prop.append('M  STY  %d %s\n' % (sty, stydat))
            else:
                break

        for i, j in enumerate(cgr_dat, start=1):
            cx, cy = cls._get_position([atoms[x - 1] for x in j['atoms']])
            mol_prop.append('M  SAL %3d%3d %s\n' % (i, len(j['atoms']), ' '.join(['%3d' % x for x in j['atoms']])))
            mol_prop.append('M  SDT %3d %s\n' % (i, j['type']))
            mol_prop.append('M  SDD %3d %10.4f%10.4f    DAU   ALL  0       0\n' % (i, cx, cy))
            mol_prop.append('M  SED %3d %s\n' % (i, j['value']))

        return ''.join(chain(("\n  CGRtools. (c) Dr. Ramil I. Nugmanov\n"
                             "\n%3s%3s  0  0  0  0            999 V2000\n" % (len(atoms), len(bonds)),),
                             ("%(x)10.4f%(y)10.4f%(z)10.4f %(element)-3s 0%(charge)3s  0  0  0  0  0"
                              "%(mark)3s  0%(map)3s  0  0\n" % i for i in atoms),
                             ("%3d%3d%3s%3d  0  0  0\n" % i for i in bonds), mol_prop))

    _stereo_map = {-1: 6, 0: 0, 1: 1, None: 0}
    _charge_map = {-3: 7, -2: 6, -1: 5, 0: 0, 1: 3, 2: 2, 3: 1}
    _radical_map = {2: 2, 1: 1, 3: 3}
