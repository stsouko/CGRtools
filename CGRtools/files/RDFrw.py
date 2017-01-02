#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2014-2017 Ramil Nugmanov <stsouko@live.ru>
# This file is part of FEAR (Fix Errors in Automapped Reactions).
#
# FEAR is free software; you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
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
import time
import periodictable as pt
from collections import defaultdict
from itertools import count, chain
from .CGRrw import CGRread, CGRwrite, fromMDL
from .SDFrw import SDFread
from . import ReactionContainer, MoleculeContainer


class RDFread(CGRread):
    def __init__(self, file, remap=True):
        self.__remap = remap
        self.__RDFfile = file
        self.__data = self.__reader()

    def read(self):
        return list(self.__data)

    def __iter__(self):
        return self.__data

    def __next__(self):
        return next(self.__data)

    def __reader(self):
        ir = -1
        im = -1
        atomcount = -1
        bondcount = -1
        failkey = True
        reaction = None
        isreaction = True
        mkey = None
        mend = False
        for n, line in enumerate(self.__RDFfile):
            if failkey and not line.startswith(("$RFMT", "$MFMT")):
                continue
            elif line.startswith("$RFMT"):
                if reaction:
                    try:
                        yield self.get_reaction(reaction, self.__remap) if isreaction \
                            else self.get_molecule(reaction, self.__remap)
                    except:
                        pass
                reaction = {'substrats': [], 'products': [], 'meta': {}, 'colors': {}}
                isreaction = True
                mkey = None
                ir = n + 5
                failkey = False
            elif line.startswith("$MFMT"):
                if reaction:
                    try:
                        yield self.get_reaction(reaction, self.__remap) if isreaction \
                            else self.get_molecule(reaction, self.__remap)
                    except:
                        pass
                reaction = {'substrats': [], 'products': [], 'meta': {}, 'colors': {}}
                molecule = {'atoms': [], 'bonds': [], 'CGR_DAT': {}}
                isreaction = False
                substrats, products = 1, 0
                mkey = None
                mend = False
                im = n + 4
                ir = -1
                failkey = False
            elif n == ir:
                try:
                    substrats, products = int(line[:3]), int(line[3:6])
                except ValueError:
                    failkey = True
                    reaction = None
            elif line.startswith("$MOL"):
                molecule = {'atoms': [], 'bonds': [], 'CGR_DAT': {}}
                im = n + 4
                mend = False
            elif n == im:
                try:
                    atomcount = int(line[:3]) + im
                    bondcount = int(line[3:6]) + atomcount
                except ValueError:
                    failkey = True
                    reaction = None
            elif n <= atomcount:
                molecule['atoms'].append(dict(element=line[31:34].strip(), isotop=int(line[34:36]),
                                              charge=fromMDL.get(int(line[38:39]), 0),
                                              map=int(line[60:63]), mark=line[54:57].strip(),
                                              x=float(line[:10]), y=float(line[10:20]), z=float(line[20:30])))
            elif n <= bondcount:
                try:
                    molecule['bonds'].append((int(line[:3]), int(line[3:6]), int(line[6:9])))
                except:
                    failkey = True
                    reaction = None

            elif line.startswith("M  END"):
                mend = True
                molecule['CGR_DAT'] = self.getdata()
                if len(reaction['substrats']) < substrats:
                    reaction['substrats'].append(molecule)
                elif len(reaction['products']) < products:
                    reaction['products'].append(molecule)

            elif n > bondcount:
                try:
                    if not mend:
                        self.collect(line)
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
                except:
                    failkey = True
                    reaction = None
        else:
            if reaction:
                try:
                    yield self.get_reaction(reaction, self.__remap) if isreaction \
                        else self.get_molecule(reaction, self.__remap)
                except:
                    pass

    @staticmethod
    def get_molecule(reaction, remap):
        molecule = reaction['substrats'][0]
        molecule['meta'] = reaction['meta']
        molecule['colors'] = reaction['colors']
        return SDFread.get_molecule(molecule, remap)

    @staticmethod
    def get_reaction(reaction, remap):
        maps = {'substrats': [y['map'] for x in reaction['substrats'] for y in x['atoms']],
                'products': [y['map'] for x in reaction['products'] for y in x['atoms']]}
        length = count(max((max(maps['products'], default=0), max(maps['substrats'], default=0))) + 1)

        ''' map unmapped atoms
        '''
        for i in ('substrats', 'products'):
            if 0 in maps[i]:
                maps[i] = [j or next(length) for j in maps[i]]

            if len(maps[i]) != len(set(maps[i])):
                raise Exception("MapError")
        ''' end
        '''
        ''' find breaks in map. e.g. 1,2,5,6. 3,4 - skipped
        '''
        if remap:
            lose = sorted(set(range(1, next(length))).difference(maps['products']).difference(maps['substrats']),
                          reverse=True)
            if lose:
                for k in ('substrats', 'products'):
                    for i in lose:
                        maps[k] = [j if j < i else j - 1 for j in maps[k]]
        ''' end
        '''
        greaction = ReactionContainer(meta={x: '\n'.join(y) for x, y in reaction['meta'].items()})

        colors = defaultdict(dict)
        for k, v in reaction['colors'].items():
            color_type, mol_num = k.split('.')
            colors[int(mol_num)][color_type] = v

        counter = count(1)
        for i in ('substrats', 'products'):
            shift = -1
            for j in reaction[i]:
                g = MoleculeContainer()
                for k, l in enumerate(j['atoms'], start=1):
                    ks = k + shift
                    g.add_node(maps[i][ks], mark=l['mark'], x=l['x'], y=l['y'], z=l['z'],
                               s_charge=l['charge'], p_charge=l['charge'], sp_charge=l['charge'], map=l['map'])
                    if l['element'] not in ('A', '*'):
                        g.node[maps[i][ks]]['element'] = l['element']
                    if l['isotop']:
                        a = pt.elements.symbol(l['element'])
                        g.node[maps[i][ks]]['isotop'] = max((a[x].abundance, x) for x in a.isotopes)[1] + l['isotop']
                for k, l, m in j['bonds']:
                    g.add_edge(maps[i][k + shift], maps[i][l + shift], s_bond=m, p_bond=m, sp_bond=m)

                for k in j['CGR_DAT']:
                    atom1 = maps[i][k['atoms'][0] + shift]
                    atom2 = maps[i][k['atoms'][-1] + shift]

                    CGRread.cgr_dat(g, k, atom1, atom2)

                for k, v in colors[next(counter)].items():
                    CGRread.parsecolors(g, k, v,
                                        {x: y for x, y in enumerate(maps[i][shift + 1: shift + 1 + len(j['atoms'])],
                                                                    start=1)})

                shift += len(j['atoms'])
                greaction[i].append(g)

        return greaction


class RDFwrite(CGRwrite):
    def __init__(self, file, extralabels=False, mark_to_map=False):
        CGRwrite.__init__(self, extralabels=extralabels, mark_to_map=mark_to_map)
        self.__file = file
        self.write = self.__initwrite

    def close(self):
        self.__file.close()

    def __initwrite(self, data):
        self.__file.write(time.strftime("$RDFILE 1\n$DATM    %m/%d/%y %H:%M\n"))
        self.__writedata(data)
        self.write = self.__writedata

    def __writedata(self, data):
        if isinstance(data, MoleculeContainer):
            m = self.getformattedcgr(data)
            self.__file.write('$MFMT\n')
            self.__file.write(m['CGR'])
            self.__file.write("M  END\n")
            colors = m['colors']
        else:
            self.__file.write('$RFMT\n$RXN\n\n  CGRtools. (c) Dr. Ramil I. Nugmanov\n\n%3d%3d\n' %
                              (len(data.substrats), len(data.products)))
            colors = {}
            for cnext, m in enumerate(chain(data.substrats + data.products), start=1):
                m = self.getformattedcgr(m)
                self.__file.write('$MOL\n')
                self.__file.write(m['CGR'])
                self.__file.write("M  END\n")
                colors.update({'%s.%d' % (k, cnext): v for k, v in m['colors'].items()})

        for p in chain(colors.items(), data.meta.items()):
            self.__file.write('$DTYPE %s\n$DATUM %s\n' % p)
