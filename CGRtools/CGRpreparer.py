# -*- coding: utf-8 -*-
#
# Copyright 2014 Ramil Nugmanov <stsouko@live.ru>
# This file is part of cgrtools.
#
# cgrtools is free software; you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
from collections import defaultdict
import copy
from itertools import product
import networkx as nx
from .CGRrw import fromMDL


class CGRPreparer(object):
    def __init__(self, **kwargs):
        self.__cgrtype = self.__getcgrtype(kwargs['type'])
        self.__cgr = self.__tocgr()
        self.__stereo = dict(atom=self.__getatomstereo, bond=self.__getbondstereo) if kwargs['stereo'] else \
            dict(atom=lambda *_: None, bond=lambda *_: None)

    def acceptrepair(self):
        return self.__cgrtype not in (1, 2, 3, 4, 5, 6)

    def prepare(self, data):
        matrix = self.__creatematrix(data)
        return matrix

    def preparetemplate(self, data):
        return dict(substrats=nx.union_all(data['substrats']), products=nx.union_all(data['products']))

    def getformattedcgr(self, graph):
        data = dict(atoms=[], bonds=[], CGR_DAT=[])
        renum = {}
        for n, (i, j) in enumerate(graph.nodes(data=True), start=1):
            renum[i] = n
            meta = self.__charge(j['s_charge'], j['p_charge'])
            if meta:
                meta['atoms'] = (renum[i],)
                data['CGR_DAT'].append(meta)

            meta = self.__stereo['atom'](j['s_stereo'], j['p_stereo'])
            if meta:
                meta['atoms'] = (renum[i],)
                data['CGR_DAT'].append(meta)

            data['atoms'].append(dict(map=i, charge=j['s_charge'], isotop=j['isotop'], element=j['element'],
                                      x=j['x'], y=j['y'], z=j['z'], mark=j['mark']))
        for i, l, j in graph.edges(data=True):
            bond, cbond, btype = self.__cgr[j['s_bond']][j['p_bond']]
            if btype:
                data['CGR_DAT'].append({'value': cbond, 'type': btype, 'atoms': (renum[i], renum[l])})

            data['bonds'].append((renum[i], renum[l], bond))

            meta = self.__stereo['bond'](j['s_stereo'], j['p_stereo'])
            if meta:
                meta['atoms'] = (renum[i], renum[l])
                data['CGR_DAT'].append(meta)

        return data

    def __getcgrtype(self, cgrtype):
        needed = [int(x) for x in cgrtype.split(',')]
        if needed[0] == 0:
            t = 0  # CGR
        elif needed[0] == 1:
            t = 1  # all reagents
        elif needed[0] == 2:
            t = 2  # all products
        elif not any(True for x in needed if -200 < x < -100) and not any(True for x in needed if -300 < x < -200) and \
                any(True for x in needed if 100 < x < 200) and any(True for x in needed if 200 < x < 300):
            t = 7  # CGR on included parts of reagents and products
        elif any(True for x in needed if -200 < x < -100) and any(True for x in needed if -300 < x < -200):
            t = 8  # CGR on excluded parts of reagents and products
        elif any(True for x in needed if -200 < x < -100) and not any(True for x in needed if -300 < x < -200) and \
                any(True for x in needed if 200 < x < 300):
            t = 9  # CGR on excluded part of reagents and included part of products
        elif not any(True for x in needed if -200 < x < -100) and any(True for x in needed if -300 < x < -200) and \
                any(True for x in needed if 100 < x < 200):
            t = 10  # CGR on excluded part of products and included part of reagents
        elif 100 < needed[0] < 200:
            t = 3  # only included part of reagents
        elif 200 < needed[0] < 300:
            t = 4  # only included part of products
        elif -200 < needed[0] < -100:
            t = 5  # only excluded part of reagents
        elif -300 < needed[0] < -200:
            t = 6  # only excluded part of products
        else:
            t = 0

        if t > 2:
            self.__needed = dict(substrats=sorted([abs(x) - 101 for x in needed if 100 < abs(x) < 200], reverse=True),
                                 products=sorted([abs(x) - 201 for x in needed if 200 < abs(x) < 300], reverse=True))
        return t

    def __creatematrix(self, data):
        def getmols(t):
            mols = []
            for x in self.__needed[t]:
                try:
                    mols.append(data[t][x])
                except IndexError:
                    pass
            return mols

        def excmols(t):
            mols = copy.deepcopy(data[t])
            for x in self.__needed[t]:
                try:
                    mols.pop(x)
                except IndexError:
                    pass
            return mols

        matrix = {}

        if self.__cgrtype == 0:
            matrix['substrats'] = nx.union_all(data['substrats'])
            matrix['products'] = nx.union_all(data['products'])
        elif self.__cgrtype == 1:
            matrix['products'] = matrix['substrats'] = nx.union_all(data['substrats'])
        elif self.__cgrtype == 2:
            matrix['products'] = matrix['substrats'] = nx.union_all(data['products'])
        elif self.__cgrtype == 3:
            matrix['products'] = matrix['substrats'] = nx.union_all(getmols('substrats'))
        elif self.__cgrtype == 4:
            matrix['products'] = matrix['substrats'] = nx.union_all(getmols('products'))
        elif self.__cgrtype == 5:
            matrix['products'] = matrix['substrats'] = nx.union_all(excmols('substrats'))
        elif self.__cgrtype == 6:
            matrix['products'] = matrix['substrats'] = nx.union_all(excmols('products'))
        elif self.__cgrtype == 7:
            matrix['substrats'] = nx.union_all(getmols('substrats'))
            matrix['products'] = nx.union_all(getmols('products'))
        elif self.__cgrtype == 8:
            matrix['substrats'] = nx.union_all(excmols('substrats'))
            matrix['products'] = nx.union_all(excmols('products'))
        elif self.__cgrtype == 9:
            matrix['substrats'] = nx.union_all(excmols('substrats'))
            matrix['products'] = nx.union_all(getmols('products'))
        elif self.__cgrtype == 10:
            matrix['substrats'] = nx.union_all(getmols('substrats'))
            matrix['products'] = nx.union_all(excmols('products'))

        for i, k in (('substrats', 's_part'), ('products', 'p_part')):
            for j in matrix[i]:
                matrix[i].node[j][k] = True

        return matrix

    def __propertyswitcher(self, g, s, t):

        pass

    @staticmethod
    def __charge(s, p):
        ss = fromMDL.get(s)
        pp = fromMDL.get(p)
        diff = pp - ss
        if diff:
            meta = {'value': 'c%+d' % diff, 'type': 'dynatom'}
        else:
            meta = None
        return meta

    @staticmethod
    def __tocgr():
        ways = [None, 0, 1, 2, 3, 4, 9]
        rep = {None: 0}
        cgrdict = defaultdict(dict)
        for x, y in product(ways, repeat=2):
            cgrdict[x][y] = ('8', '%s>%s' % (rep.get(x, x), rep.get(y, y)), 'dynbond') if x != y else \
                ('8', 's', 'extrabond') if x == 9 else (str(rep.get(x, x)), None, False)
        return cgrdict

    def __getbondstereo(self, s, p):
        return self.__getstereo(s, p, 'bondstereo', 'dynbondstereo')

    def __getatomstereo(self, s, p):
        return self.__getstereo(s, p, 'atomstereo', 'dynatomstereo')

    def __getstereo(self, s, p, t1, t2):
        rep = {None: 'n'}
        s = rep.get(s, s)
        p = rep.get(p, p)
        if s == p:
            stereo = s
            stype = t1
        else:
            stereo = '%s>%s' % (s, p)
            stype = t2
        if stereo != 'n':
            return {'value': stereo, 'type': stype}
        else:
            return None
