#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2014-2016 Ramil Nugmanov <stsouko@live.ru>
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
import operator
from collections import Counter, defaultdict
from itertools import chain, count
from functools import reduce
import networkx as nx
import periodictable as pt


def eratosthenes():
    """Yields the sequence of prime numbers via the Sieve of Eratosthenes."""
    D = {}  # map each composite integer to its first-found prime factor
    for q in count(2):  # q gets 2, 3, 4, 5, ... ad infinitum
        p = D.pop(q, None)
        if p is None:
            # q not a key in D, so q is prime, therefore, yield it
            yield q
            # mark q squared as not-prime (with q as first-found prime factor)
            D[q * q] = q
        else:
            # let x <- smallest (N*p)+q which wasn't yet known to be composite
            # we just learned x is composite, with p first-found prime factor,
            # since p is the first-found prime factor of q -- find and mark it
            x = p + q
            while x in D:
                x += p
            D[x] = p


class FEAR(object):
    def __init__(self, isotop=False, stereo=False, hyb=False, element=True, deep=0):
        self.__primes = tuple(x for _, x in zip(range(1000), eratosthenes()))
        self.__isotop = isotop
        self.__stereo = stereo
        self.__deep = deep
        self.__hyb = hyb
        self.__element = element
        self.__whash = set()
        self.__lhash = set()
        self.__shash = set()

    def sethashlib(self, data):
        self.__whash = {x['CGR_FEAR_WHASH'] for x in data}
        self.__lhash = {x['CGR_FEAR_LHASH'] for x in data}
        self.__shash = {x['CGR_FEAR_SHASH'] for x in data}

    def chkreaction(self, g, full=True, gennew=False):
        chkd = set()
        report = set()
        newhash = set()

        def dochk(x):
            weights = self.__getMorgan(x)
            whash = '.'.join(str(x) for x in sorted(weights.values()))
            if whash in self.__whash or gennew:
                shash = self.__getsmarts(x, weights)
                if shash in self.__shash:
                    chkd.update(x)
                    report.add((whash, shash))
                    return True
                elif gennew:
                    newhash.add((whash, shash))
            return False

        centers = self.__getcenters(g)
        ''' check for very specific centers
        '''
        for i in sorted(centers, reverse=True)[:-1]:
            for j in centers[i]:
                if not chkd.issuperset(j):
                    dochk(j)
        ''' final check for common centers
        '''
        tmp = (chkd.issuperset(i) or dochk(i) for i in centers[0])
        if gennew:
            tmp = list(tmp)
        return all(tmp) if full else any(tmp), report, newhash

    def getreactionhash(self, g):
        return self.__getsmarts(g, self.__getMorgan(g))

    def __getcenters(self, g):
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
        for i in range(self.__deep):
            nodes.append(set(chain.from_iterable(g.edges(nodes[i]))))

        centers = defaultdict(list)
        for n, i in enumerate(nodes):
            for j in nx.connected_component_subgraphs(g.subgraph(i)):
                centers[n].append(j)

        return centers

    def __getsmarts(self, g, weights):
        newmaps = dict()
        countmap = count(1)
        countcyc = count(1)
        visited = set()

        def getnextatom(atoms):
            if len(atoms) == 1:
                nextatom = list(atoms)[0]
            else:
                nextatom = sorted((weights[i], i) for i in atoms)[-1][1]

            visited.add(nextatom)
            return nextatom

        def dosmarts(trace, inter, prev):
            smis = ['[%s%s%s%s:%d]' % (g.node[inter].get('isotop', '') if self.__isotop else '',
                                       g.node[inter]['element'] if self.__element else '*',
                                       ';%s;' % ','.join(
                                           [self.__stereotypes[g.node[inter].get('s_stereo')] if self.__stereo else '',
                                            self.__hybtypes[g.node[inter].get('s_hyb')] if self.__hyb else '']),
                                       '%+d' % g.node[inter].get('s_charge', 0) if self.__element else '',
                                       newmaps.get(inter) or newmaps.setdefault(inter, next(countmap)))]
            smip = ['[%s%s%s%s:%d]' % (g.node[inter].get('isotop', '') if self.__isotop else '',
                                       g.node[inter]['element'] if self.__element else '*',
                                       ';%s;' % ','.join(
                                           [self.__stereotypes[g.node[inter].get('p_stereo')] if self.__stereo else '',
                                            self.__hybtypes[g.node[inter].get('p_hyb')] if self.__hyb else '']),
                                       '%+d' % g.node[inter].get('p_charge', 0) if self.__element else '',
                                       newmaps.get(inter))]
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
                            smis.append('%s%d' % (self.__tosmiles[g[inter][j[2]].get('s_bond')], j[1]))
                            smip.append('%s%d' % (self.__tosmiles[g[inter][j[2]].get('p_bond')], j[1]))
                smis.extend(['(' if iterlist else ''] + ['%s' % self.__tosmiles[g[inter][i].get('s_bond')]] + deep[1] +
                            [')' if iterlist else ''])
                smip.extend(['(' if iterlist else ''] + ['%s' % self.__tosmiles[g[inter][i].get('p_bond')]] + deep[2] +
                            [')' if iterlist else ''])
            return trace, smis, smip, concat

        has_next = g
        ssmiles, psmiles = [], []
        while has_next:
            firstatom = getnextatom(has_next)
            smirks = dosmarts({firstatom}, firstatom, firstatom)
            ssmiles.append(''.join(smirks[1]))
            psmiles.append(''.join(smirks[2]))
            has_next = set(g).difference(visited)
        return '%s>>%s' % ('.'.join(ssmiles), '.'.join(psmiles))

    __tosmiles = {1: '-', 2: '=', 3: '#', 4: ':', None: '.', 9: '~'}

    __hybtypes = {4: 'Ha', 3: 'Ht', 2: 'Hd', 1: 'Hs', None: ''}

    __stereotypes = {None: '', 1: '@s', 2: '@r'}

    def __getMorgan(self, g):
        newlevels = {}
        countprime = iter(self.__primes)

        params = {n: (self.__primes[pt.elements.symbol(attr['element']).number] if self.__element else 1,
                      self.__primes[10 * attr['s_charge'] + attr['p_charge']] if self.__element else 1,
                      reduce(operator.mul,
                             (self.__primes[10 * (eattr.get('s_bond') or 0) + (eattr.get('p_bond') or 0)]
                              for eattr in g[n].values()), 1),
                      self.__primes[attr['isotop']] if self.__isotop and 'isotop' in attr else 1,
                      self.__primes[10 * (attr.get('s_stereo') or 0) + (attr.get('p_stereo') or 0)]
                      if self.__stereo else 1,
                      reduce(operator.mul,
                             (self.__primes[10 * (eattr.get('s_stereo') or 0) + (eattr.get('p_stereo') or 0)]
                              for eattr in g[n].values()), 1) if self.__stereo else 1)
                  for n, attr in g.nodes(data=True)}

        weights = {x: (newlevels.get(y) or newlevels.setdefault(y, next(countprime)))
                   for x, y in sorted(params.items(), key=operator.itemgetter(1))}

        oldnumb = numb = len(g)
        maxcount = 0
        stab = 0

        scaf = {}
        for n, m in g.edge.items():
            scaf[n] = tuple(m)

        while oldnumb >= numb and maxcount != 1 and stab < 3:
            oldnumb = numb
            neweights = {}
            countprime = iter(self.__primes)

            tmp = {}
            for n, m in scaf.items():
                """ if don't have neighbors use self weight
                """
                tmp[n] = reduce(operator.mul, (weights[x] for x in m), weights[n]**2)

            numb = len(set(tmp.values()))
            if numb == oldnumb:
                x = Counter(tmp.values())
                stab += 1
                maxcount = x[max(x)]
            else:
                stab = 0
                maxcount = 0

            weights = {x: (neweights.get(y) or neweights.setdefault(y, next(countprime)))
                       for x, y in sorted(tmp.items(), key=operator.itemgetter(1))}

        return weights
