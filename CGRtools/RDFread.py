#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2014 Ramil Nugmanov <stsouko@live.ru>
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
from .CGRrw import CGRRead
import networkx as nx
from .CGRrw import fromMDL, toMDL, bondlabels


class RDFread(CGRRead):
    def __init__(self, file):
        CGRRead.__init__(self)
        self.__RDFfile = file

    def readdata(self):
        """парсер RDF файлов
        """
        if "$RDFILE 1" not in self.__RDFfile.readline():
            raise StopIteration
        ir = -1
        im = -1
        atomcount = -1
        bondcount = -1
        failkey = False
        reaction = None
        meta = None
        mend = False
        for n, line in enumerate(self.__RDFfile):
            if failkey and "$RXN" not in line[0:4]:
                continue
            elif "$RXN" in line[0:4]:
                if reaction:
                    try:
                        yield self.__getGraphs(reaction)
                    except:
                        pass
                reaction = {'substrats': [], 'products': [], 'meta': defaultdict(str)}
                meta = None
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
                molecule['atoms'].append(dict(element=line[31:34].strip(), isotop=line[34:36].strip(),
                                              charge=int(line[38:39]),
                                              map=int(line[60:63]), mark=line[51:54].strip(),
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
                        meta = line[7:].strip()
                    elif '$RFMT' not in line and meta:
                        reaction['meta'][meta] += line.strip("$DATUM").strip() + ' '
                except:
                    failkey = True
                    reaction = None
        else:
            if reaction:
                try:
                    yield self.__getGraphs(reaction)
                except:
                    pass

    def __getGraphs(self, reaction):
        maps = {'substrats': [int(y['map']) for x in reaction['substrats'] for y in x['atoms']],
                'products': [int(y['map']) for x in reaction['products'] for y in x['atoms']]}
        length = max((max(maps['products']), max(maps['substrats'])))

        ''' map unmapped atoms
        '''
        for i in ('substrats', 'products'):
            if 0 in maps[i]:
                nmaps = []
                for j in maps[i]:
                    if j:
                        nmaps.append(j)
                    else:
                        length += 1
                        nmaps.append(length)
                maps[i] = nmaps
            if len(maps[i]) != len(set(maps[i])):
                raise Exception("MapError")
        ''' end
        '''
        ''' find breaks in map. e.g. 1,2,5,6. 3,4 - skipped
        '''
        lose = sorted(set(range(1, length + 1)).difference(maps['products']).difference(maps['substrats']),
                      reverse=True)
        if lose:
            for k in maps:
                for i in lose:
                    newmaps = []
                    for j in maps[k]:
                        newmaps.append(j if j < i else j - 1)
                    maps[k] = newmaps
            length -= len(lose)
        ''' end
        '''
        greaction = dict(substrats=[], products=[], meta=reaction['meta'])

        for i in ('substrats', 'products'):
            shift = 0
            for j in reaction[i]:
                g = nx.Graph()
                for k, l in enumerate(j['atoms']):
                    g.add_node(maps[i][k + shift], element=l['element'], mark=l['mark'],  x=l['x'], y=l['y'], z=l['z'],
                               s_charge=l['charge'], s_stereo=None,
                               p_charge=l['charge'], p_stereo=None,
                               map=maps[i][k + shift], isotop=l['isotop'])

                for k, l, m in j['bonds']:
                    g.add_edge(maps[i][k + shift], maps[i][l + shift],
                               s_bond=m, s_stereo=None,
                               p_bond=m, p_stereo=None)

                for k in j['CGR_DAT']:
                    atom1 = maps[i][k['atoms'][0] + shift]
                    atom2 = maps[i][k['atoms'][-1] + shift]

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

                shift += len(j['atoms'])
                greaction[i].append(g)

        return greaction

