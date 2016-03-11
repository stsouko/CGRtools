# -*- coding: utf-8 -*-
#
# Copyright 2014-2016 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.CGRrw import fromMDL
from CGRtools.RDFread import RDFread
from CGRtools.SDFread import SDFread


class CGRPreparer(object):
    def __init__(self, cgrtype):
        self.__cgrtype = self.__getcgrtype(cgrtype)
        self.__cgr = self.__tocgr()

    def acceptrepair(self):
        return self.__cgrtype not in (1, 2, 3, 4, 5, 6)

    def gettemplates(self, templates, isreaction=True):
        if templates:
            source = RDFread(templates) if isreaction else SDFread(templates)

            templates = []
            for template in source.readdata():
                if isreaction:
                    matrix = self.preparetemplate(template)
                    nx.relabel_nodes(matrix['substrats'], {x: x + 1000 for x in matrix['substrats']}, copy=False)
                    nx.relabel_nodes(matrix['products'], {x: x + 1000 for x in matrix['products']}, copy=False)
                else:
                    matrix = dict(substrats=template['structure'], meta=template['meta'])
                templates.append(matrix)
            return templates

        return None

    @staticmethod
    def preparetemplate(data):
        res = dict(meta=data['meta'])
        for i in ('products', 'substrats'):
            x = nx.union_all(data[i])
            for n, d in x.nodes(data=True):
                if d['element'] in ('A', '*'):
                    x.node[n].pop('element')
            res[i] = x

        for i in set(res['products']).intersection(res['substrats']):
            for j in {'sp_charge', 'sp_hyb', 'sp_neighbors', 'sp_stereo'}.intersection(res['products'].node[i]):
                if isinstance(res['products'].node[i][j], list):
                    res['products'].node[i][j] = {x: y for x, y in zip(res['substrats'].node[i][j],
                                                                       res['products'].node[i][j])}
        print(res['products'].node)
        return res

    def getformattedcgr(self, graph):
        data = dict(atoms=[], bonds=[], CGR_DAT=[], extended=[])
        renum = {}
        for n, (i, j) in enumerate(graph.nodes(data=True), start=1):
            renum[i] = n
            meta = self.__charge(j['s_charge'], j['p_charge'])
            if meta:
                meta['atoms'] = (renum[i],)
                data['CGR_DAT'].append(meta)

            meta = self.__getstate(j.get('s_stereo'), j.get('p_stereo'), 'atomstereo', 'dynatomstereo')
            if meta:
                meta['atoms'] = (renum[i],)
                data['CGR_DAT'].append(meta)

            meta = self.__getstate(j.get('s_hyb'), j.get('p_hyb'), 'atomhyb', 'dynatomhyb')
            if meta:
                meta['atoms'] = (renum[i],)
                data['CGR_DAT'].append(meta)

            meta = self.__getstate(j.get('s_neighbors'), j.get('p_neighbors'), 'atomneighbors', 'dynatomneighbors')
            if meta:
                meta['atoms'] = (renum[i],)
                data['CGR_DAT'].append(meta)

            if j.get('isotop'):
                data['extended'].append(dict(atoms=(renum[i],), value=j.get('isotop'), type='isotop'))

            data['atoms'].append(dict(map=i, charge=j['s_charge'], element=j.get('element', 'A'),
                                      x=j['x'], y=j['y'], z=j['z'], mark=j['mark']))
        for i, l, j in graph.edges(data=True):
            bond, cbond, btype = self.__cgr[j['s_bond']][j['p_bond']]
            if btype:
                data['CGR_DAT'].append({'value': cbond, 'type': btype, 'atoms': (renum[i], renum[l])})

            data['bonds'].append((renum[i], renum[l], bond))

            meta = self.__getstate(j.get('s_stereo'), j.get('p_stereo'), 'bondstereo', 'dynbondstereo')
            if meta:
                meta['atoms'] = (renum[i], renum[l])
                data['CGR_DAT'].append(meta)

        return data

    def prepare(self, data, setlabels=False):
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

        if self.__cgrtype == 0:
            substrats = nx.union_all(data['substrats'])
            products = nx.union_all(data['products'])
        elif self.__cgrtype == 1:
            products = substrats = nx.union_all(data['substrats'])
        elif self.__cgrtype == 2:
            products = substrats = nx.union_all(data['products'])
        elif self.__cgrtype == 3:
            products = substrats = nx.union_all(getmols('substrats'))
        elif self.__cgrtype == 4:
            products = substrats = nx.union_all(getmols('products'))
        elif self.__cgrtype == 5:
            products = substrats = nx.union_all(excmols('substrats'))
        elif self.__cgrtype == 6:
            products = substrats = nx.union_all(excmols('products'))
        elif self.__cgrtype == 7:
            substrats = nx.union_all(getmols('substrats'))
            products = nx.union_all(getmols('products'))
        elif self.__cgrtype == 8:
            substrats = nx.union_all(excmols('substrats'))
            products = nx.union_all(excmols('products'))
        elif self.__cgrtype == 9:
            substrats = nx.union_all(excmols('substrats'))
            products = nx.union_all(getmols('products'))
        elif self.__cgrtype == 10:
            substrats = nx.union_all(getmols('substrats'))
            products = nx.union_all(excmols('products'))
        res = dict(substrats=self.__setlabels(substrats), products=self.__setlabels(products)) if setlabels else \
            dict(substrats=substrats, products=products)
        return res

    @staticmethod
    def __setlabels(g):
        tmp = g.copy()
        for i in tmp.nodes():
            label = {'s_hyb': 1, 'p_hyb': 1, 'sp_hyb': 1, 's_neighbors': 0, 'p_neighbors': 0, 'sp_neighbors': 0}
            #  hyb 1- sp3; 2- sp2; 3- sp1; 4- aromatic
            for b, h, n in (('s_bond', 's_hyb', 's_neighbors'), ('p_bond', 'p_hyb', 'p_neighbors')):
                for node, bond in tmp[i].items():
                    if tmp.node[node]['element'] != 'H' and bond[b]:
                        label[n] += 1

                    if bond[b] in (1, None):
                        pass
                    elif bond[b] == 4:
                        label[h] = 4
                    elif bond[b] == 3 or (bond[b] == 2 and label[h] == 2):  # Если есть 3-я или две 2-х связи, то sp1
                        label[h] = 3
                    elif bond[b] == 2:  # Если есть 2-я связь, но до этого не было найдено другой 2-й, 3-й, или аром.
                        label[h] = 2
            for n, m, h in (('s_hyb', 'p_hyb', 'sp_hyb'), ('s_neighbors', 'p_neighbors', 'sp_neighbors')):
                label[h] = (label[n], label[m]) if label[n] != label[m] else label[n]

            for k in list(label):
                if tmp.node[i].get(k) is not None:
                    label.pop(k)

            tmp.node[i].update(label)
        return tmp

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

    @staticmethod
    def __getstate(s, p, t1, t2):
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
