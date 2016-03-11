#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2014, 2015 Ramil Nugmanov <stsouko@live.ru>
# This file is part of FEAR (Find Errors in Automapped Reactions).
#
# FEAR is free software; you can redistribute it and/or modify
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
from collections import Counter, defaultdict
from itertools import chain, count
import networkx as nx
from CGRtools.weightable import mendeley, atommass
from CGRtools.CGRrw import fromMDL


class FEAR(object):
    def __init__(self, deep=0, isotop=False):
        self.__deep = list(range(deep))
        self.__isotop = isotop

    def sethashlib(self, data):
        self.__whash = {x['FEAR_WHASH'] for x in data}
        self.__lhash = {x['FEAR_LHASH'] for x in data}
        self.__shash = {x['FEAR_SHASH'] for x in data}

    def chkreaction(self, g, full=True):
        chkd = set()
        report = set()

        def dochk(x):
            weights = self.__getMorgan(x)
            whash = '.'.join(str(x) for x in sorted(weights.values()))
            if whash in self.__whash:
                levels = self.__getLevels(x)
                lhash = '.'.join(str(x) for x in sorted(levels.values()))
                if lhash in self.__lhash:
                    shash = self.getsmarts(x, weights, levels)
                    if shash in self.__shash:
                        chkd.update(x)
                        report.add(shash)
                        return True
            return False

        centers = self.getcenters(g)
        for i in sorted(centers, reverse=True)[:-1]:
            for j in centers[i]:
                if not chkd.issuperset(j):
                    dochk(j)
        tmp = [chkd.issuperset(i) or dochk(i) for i in centers[0]]
        return all(tmp) if full else any(tmp), report

    def getcenters(self, g):
        nodes = [set()]
        """ get nodes of reaction center (dynamic bonds, stereo or charges).
        """
        for n, node_attr in g.nodes(data=True):
            for k, l in (('s_charge', 'p_charge'), ('s_stereo', 'p_stereo')):
                if node_attr.get(k) != node_attr.get(l):
                    nodes[0].add(n)

        for *n, node_attr in g.edges(data=True):
            for k, l in (('s_bond', 'p_bond'), ('s_stereo', 'p_stereo')):
                if node_attr.get(k) != node_attr.get(l):
                    nodes[0].update(n)

        """ get nodes of reaction center and neighbors
        """
        for i in self.__deep:
            nodes.append(set(chain.from_iterable(g.edges(nodes[i]))))

        centers = defaultdict(list)
        for n, i in enumerate(nodes):
            for j in nx.connected_component_subgraphs(g.subgraph(i)):
                centers[n].append(j)

        return centers

    def getsmarts(self, g, weights, levels):
        newmaps = dict()
        countmap = count(1)
        countcyc = count(1)

        def getnextatom(atoms):
            if len(atoms) == 1:
                return [i for i in atoms][0]
            maxw = max(weights[i] for i in atoms)
            morgans = [i for i in atoms if weights[i] == maxw]
            if len(morgans) == 1:
                return morgans[0]
            else:
                maxl = max(levels[i] for i in morgans)
                return next(i for i in morgans if levels[i] == maxl)

        def dosmarts(trace, inter, prev):
            smi = [(lambda x: lambda ch, el=True, ms=False, st=False, hb=False: '[%s%s%s%s:%s]' %
                                                                                (x['isotop'] if ms else '',
                                                                                 x['element'] if el else '*',
                                                                                 ';%s;' % ','.join(
                                                                                     [x['s_stereo'] if st == 's' else
                                                                                      x['p_stereo'] if st == 'p' else
                                                                                      '',
                                                                                      x['s_hyb'] if hb == 's' else
                                                                                      x['p_hyb'] if hb == 'p' else '']),
                                                                                 x['s_chagre'] if ch == 's' else
                                                                                 x['p_charge'] if ch == 'p' else '',
                                                                                 x['map']))(
                dict(isotop=g.node[inter].get('isotop', ''), element=g.node[inter]['element'],
                     s_stereo=self.__stereotypes[g.node[inter].get('s_stereo')],
                     p_stereo=self.__stereotypes[g.node[inter].get('p_stereo')],
                     s_hyb=self.__hybtypes[g.node[inter].get('s_hyb')],
                     p_hyb=self.__hybtypes[g.node[inter].get('p_hyb')],
                     map=newmaps.get(inter) or newmaps.setdefault(inter, next(countmap)),
                     s_charge='%+d' % fromMDL.get(g.node[inter].get('s_charge', 0)),
                     p_charge='%+d' % fromMDL.get(g.node[inter].get('p_charge', 0))
                     ))]
            smis = ['[%s:%d]' % (g.node[inter]['element'],
                                 newmaps.get(inter) or newmaps.setdefault(inter, next(countmap)))]
            smip = smis.copy()
            concat = []
            stoplist = []
            iterlist = set(g.neighbors(inter)).difference([prev])
            while iterlist:
                i = getnextatom(iterlist)
                iterlist.discard(i)
                if i in trace:
                    if i not in stoplist:  # костыль для циклов. чтоб не было 2х проходов.
                        cyc = next(countcyc)
                        concat.append((i, cyc, inter))
                        smis.append('%s%d' % (self.__tosmiles[g[inter][i].get('s_bond')], cyc))
                        smip.append('%s%d' % (self.__tosmiles[g[inter][i].get('p_bond')], cyc))
                    continue

                deep = dosmarts(set(chain(trace, [i])), i, inter)
                trace.update(deep[0])
                if deep[3]:
                    concat.extend(deep[3])
                    for j in deep[3]:
                        if j[0] == inter:
                            stoplist.append(j[2])
                            smis.append('%s%d' % (self.__tosmiles[g[inter][i].get('s_bond')], j[1]))
                            smip.append('%s%d' % (self.__tosmiles[g[inter][i].get('p_bond')], j[1]))
                smis.extend(['(' if iterlist else ''] + ['%s' % self.__tosmiles[g[inter][i].get('s_bond')]] + deep[1] +
                            [')' if iterlist else ''])
                smip.extend(['(' if iterlist else ''] + ['%s' % self.__tosmiles[g[inter][i].get('p_bond')]] + deep[2] +
                            [')' if iterlist else ''])
            return trace, smis, smip, concat

        smirks = dosmarts(set(), getnextatom(g), 0)
        return '%s>>%s' % (''.join(smirks[1]), ''.join(smirks[2]))

    __tosmiles = {1: '-', 2: '=', 3: '#', 4: ':', None: '.', 9: '~'}

    __hybtypes = {4: 'Ha', 3: 'Ht', 2: 'Hd', 1: 'Hs', None: ''}

    __stereotypes = {None: '', 1: '@s', 2: '@d'}

    def __getLevels(self, g):
        return {n: mendeley[attr['element']] +
                   (atommass[attr['element']] - attr.get('isotop', atommass[attr['element']]) if self.__isotop else 0) -
                   (attr['s_charge'] or 4) * (attr['p_charge'] or 4) -
                   max(10 * (eattr.get('s_bond') or 0) + (eattr.get('p_bond') or 0) for eattr in g[n].values())
                for n, attr in g.nodes(data=True)}

    @staticmethod
    def __getMorgan(g):
        weights = dict.fromkeys(g, 1)
        oldnumb = numb = len(g)
        maxcount = 0
        stab = 0
        while oldnumb >= numb and maxcount != 1 and stab < 3:
            oldnumb = numb
            tmp = dict.fromkeys(g, 0)
            for n, m in g.edges():
                tmp[n] += weights[m]
                tmp[m] += weights[n]

            numb = len(set(tmp.values()))
            if numb == oldnumb:
                x = Counter(tmp.values())
                stab += 1
                maxcount = x[max(x)]
            else:
                stab = 0
                maxcount = 0
            weights = tmp

        return weights
