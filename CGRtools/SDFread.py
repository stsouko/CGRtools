#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2014-2016 Ramil Nugmanov <stsouko@live.ru>
# This file is part of CGRTools (C) Ramil Nugmanov <stsouko@live.ru>.
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
from CGRtools.CGRrw import CGRRead
import periodictable as pt
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
                molecule['atoms'].append(dict(element=line[31:34].strip(), isotop=int(line[34:36]),
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
                        molecule['meta'][mkey] += line.strip()
                except:
                    failkey = True
                    molecule = None
        else:
            raise StopIteration

    def __getGraphs(self, molecule):
        g = nx.Graph()
        for k, l in enumerate(molecule['atoms'], start=1):
            g.add_node(k, element=l['element'], mark=l['mark'], x=l['x'], y=l['y'], z=l['z'],
                       s_charge=l['charge'], s_stereo=None, s_hyb=None, s_neighbors=None,
                       p_charge=l['charge'], p_stereo=None, p_hyb=None, p_neighbors=None,
                       map=l['map'])
            if l['isotop']:
                a = pt.elements.symbol(l['element'])
                g.node[k]['isotop'] = max((a[x].abundance, x) for x in a.isotopes)[1] + l['isotop']

        for k, l, m in molecule['bonds']:
            g.add_edge(k, l,
                       s_bond=m, s_stereo=None,
                       p_bond=m, p_stereo=None)

        for k in molecule['CGR_DAT']:
            atom1 = k['atoms'][0] + 1
            atom2 = k['atoms'][-1] + 1

            self.cgr_dat(g, k, atom1, atom2)

        return dict(meta=molecule['meta'], structure=g)
