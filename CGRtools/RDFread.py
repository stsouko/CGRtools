#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2014-2016 Ramil Nugmanov <stsouko@live.ru>
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
from collections import defaultdict
from itertools import count
import networkx as nx
from CGRtools.CGRrw import CGRRead
import periodictable as pt


class RDFread(CGRRead):
    def __init__(self, file):
        CGRRead.__init__(self)
        self.__RDFfile = file

    def readdata(self, remap=True):
        """ парсер RDF файлов
        """
        ir = -1
        im = -1
        atomcount = -1
        bondcount = -1
        failkey = False
        reaction = None
        mkey = None
        mend = False
        for n, line in enumerate(self.__RDFfile):
            if failkey and "$RXN" not in line[0:4]:
                continue
            elif "$RXN" in line[0:4]:
                if reaction:
                    try:
                        yield self.__getGraphs(reaction, remap)
                    except:
                        pass
                reaction = {'substrats': [], 'products': [], 'meta': {}, 'colors': {}}
                mkey = None
                ir = n + 4
                failkey = False
            elif n == ir:
                try:
                    substrats, products = int(line[0:3]), int(line[3:6])
                except:
                    failkey = True
                    reaction = None
            elif "$MOL" in line[0:4]:
                molecule = {'atoms': [], 'bonds': [], 'CGR_DAT': {}}
                im = n + 4
                mend = False
            elif n == im:
                try:
                    atomcount = int(line[0:3]) + im
                    bondcount = int(line[3:6]) + atomcount
                except:
                    failkey = True
                    reaction = None
            elif n <= atomcount:
                molecule['atoms'].append(dict(element=line[31:34].strip(), isotop=int(line[34:36]),
                                              charge=int(line[38:39]),
                                              map=int(line[60:63]), mark=line[54:57].strip(),
                                              x=float(line[0:10]), y=float(line[10:20]), z=float(line[20:30])))
            elif n <= bondcount:
                try:
                    molecule['bonds'].append((int(line[0:3]) - 1, int(line[3:6]) - 1, int(line[6:9])))
                except:
                    failkey = True
                    reaction = None

            elif "M  END" in line:
                mend = True
                molecule['CGR_DAT'] = self.getdata()
                try:
                    if len(reaction['substrats']) < substrats:
                        reaction['substrats'].append(molecule)
                    else:
                        reaction['products'].append(molecule)
                except:
                    failkey = True
                    reaction = None

            elif n > bondcount and not mend:
                try:
                    self.collect(line)
                except:
                    failkey = True
                    reaction = None

            elif n > bondcount:
                try:
                    if '$DTYPE' in line:
                        mkey = line[7:].strip()
                        if mkey.split('.')[0] in ('PHTYP', 'FFTYP', 'PCTYP', 'EPTYP', 'HBONDCHG', 'CNECHG',
                                    'dynPHTYP', 'dynFFTYP', 'dynPCTYP', 'dynEPTYP', 'dynHBONDCHG', 'dynCNECHG'):
                            target = 'colors'
                        else:
                            target = 'meta'
                        reaction[target][mkey] = []
                    elif '$RFMT' not in line and mkey:
                        data = line.lstrip("$DATUM").rstrip()
                        if data:
                            reaction[target][mkey].append(data)
                except:
                    failkey = True
                    reaction = None
        else:
            if reaction:
                try:
                    yield self.__getGraphs(reaction, remap)
                except:
                    pass

    def __getGraphs(self, reaction, remap):
        maps = {'substrats': [y['map'] for x in reaction['substrats'] for y in x['atoms']],
                'products': [y['map'] for x in reaction['products'] for y in x['atoms']]}
        length = count(max((max(maps['products']), max(maps['substrats']))) + 1)

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
        lose = sorted(set(range(1, next(length))).difference(maps['products']).difference(maps['substrats']),
                      reverse=True)
        if lose and remap:
            for k in maps:
                for i in lose:
                    maps[k] = [j if j < i else j - 1 for j in maps[k]]
        ''' end
        '''
        greaction = dict(substrats=[], products=[],
                         meta={x: '\n'.join(z.strip() for z in y) for x, y in reaction['meta'].items()})

        colors = defaultdict(dict)
        for k, v in reaction['colors'].items():
            color_type, mol_num = k.split('.')
            colors[int(mol_num)][color_type] = v

        counter = count(1)
        for i in ('substrats', 'products'):
            shift = 0
            for j in reaction[i]:
                g = nx.Graph()
                for k, l in enumerate(j['atoms']):
                    g.add_node(maps[i][k + shift], element=l['element'], mark=l['mark'], x=l['x'], y=l['y'], z=l['z'],
                               s_charge=l['charge'],  # s_stereo=None, s_hyb=None, s_neighbors=None,
                               p_charge=l['charge'],  # p_stereo=None, p_hyb=None, p_neighbors=None,
                               sp_charge=l['charge'],
                               map=l['map'])
                    if l['isotop']:
                        a = pt.elements.symbol(l['element'])
                        g.node[maps[i][k + shift]]['isotop'] = max((a[x].abundance, x) for x in a.isotopes)[1] + l['isotop']

                for k, l, m in j['bonds']:
                    g.add_edge(maps[i][k + shift], maps[i][l + shift],
                               s_bond=m,  # s_stereo=None,
                               p_bond=m,  # p_stereo=None
                               sp_bond=m)

                for k in j['CGR_DAT']:
                    atom1 = maps[i][k['atoms'][0] + shift - 1]
                    atom2 = maps[i][k['atoms'][-1] + shift - 1]

                    self.cgr_dat(g, k, atom1, atom2)

                for k, v in colors[next(counter)].items():
                    self.parsecolors(g, k, v,
                                     {x: y for x, y in enumerate(maps[i][shift: shift + len(j['atoms'])], start=1)})

                shift += len(j['atoms'])
                greaction[i].append(g)

        return greaction

