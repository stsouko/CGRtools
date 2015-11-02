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
from CGRtools.CGRrw import CGRRead, fromMDL, toMDL, bondlabels
import networkx as nx


class SDFread(CGRRead):
    def __init__(self, file):
        CGRRead.__init__(self)
        self.__SDFfile = file

    def readprop(self):
        with open(self.__SDFfile) as f:
            prop = set()
            for line in f:
                if '>  <' in line[:5]:
                    prop.add(line.strip()[4:-1].replace(',', '_'))
        return prop

    def readdata(self):
        im = 3
        atomcount = -1
        bondcount = -1
        failkey = False
        meta = None
        molecule = None
        mend = False
        for n, line in enumerate(self.__SDFfile):
            if failkey and "$$$$" not in line[0:4]:
                continue
            elif "$$$$" in line[0:4]:
                if molecule:
                    try:
                        yield self.__getGraphs(molecule)
                    except:
                        pass

                meta = None
                im = n + 4
                failkey = False
                mend = False

            elif n == im:
                molecule = {'atoms': [], 'bonds': [], 'CGR_DAT': {}, 'meta': {}}

                try:
                    atomcount = int(line[0:3]) + n
                    bondcount = int(line[3:6]) + atomcount
                except:
                    failkey = True
                    molecule = None

            elif n <= atomcount:
                molecule['atoms'].append(dict(element=line[31:34].strip(), isotop=line[34:36].strip(),
                                              charge=int(line[38:39]),
                                              map=int(line[60:63]), mark=line[51:54].strip(),
                                              x=float(line[0:10]), y=float(line[10:20]), z=float(line[20:30])))

            elif n <= bondcount:
                try:
                    molecule['bonds'].append((int(line[0:3]), int(line[3:6]), int(line[6:9])))
                except:
                    failkey = True
                    molecule = None

            elif "M  END" in line:
                mend = True
                molecule['CGR_DAT'] = self.getdata()

            elif n > bondcount and not mend:
                try:
                    self.collect(line)
                except:
                    failkey = True
                    molecule = None

            elif n > bondcount:
                try:
                    if '>  <' in line:
                        meta = True
                        mkey = line.strip()[4:-1]
                        molecule['meta'][mkey] = ''
                    elif meta:
                        molecule['meta'][mkey] += line.strip() + ' '
                except:
                    failkey = True
                    molecule = None
        else:
            raise StopIteration

    def __getGraphs(self, molecule):
        g = nx.Graph()
        for k, l in enumerate(molecule['atoms'], start=1):
            g.add_node(k, element=l['element'], mark=l['mark'], x=l['x'], y=l['y'], z=l['z'],
                       s_charge=l['charge'], s_stereo=None,
                       p_charge=l['charge'], p_stereo=None,
                       map=l['map'], isotop=l['isotop'])

        for k, l, m in molecule['bonds']:
            g.add_edge(k, l,
                       s_bond=m, s_stereo=None,
                       p_bond=m, p_stereo=None)

        for k in molecule['CGR_DAT']:
            atom1 = k['atoms'][0] + 1
            atom2 = k['atoms'][-1] + 1

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

        return dict(meta=molecule['meta'], structure=g)
