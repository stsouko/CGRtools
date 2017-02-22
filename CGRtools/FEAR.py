# -*- coding: utf-8 -*-
#
# Copyright 2014-2017 Ramil Nugmanov <stsouko@live.ru>
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
from operator import mul, itemgetter
from collections import Counter
from itertools import chain, count
from functools import reduce
from networkx import connected_component_subgraphs
from periodictable import elements


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
    def __init__(self, isotope=False, stereo=False, hyb=False, element=True, deep=0):
        self.__primes = tuple(x for _, x in zip(range(1000), eratosthenes()))
        self.__isotope = isotope
        self.__stereo = stereo
        self.__deep = deep
        self.__hyb = hyb
        self.__element = element
        if stereo:
            self.__rc_node.append(('s_stereo', 'p_stereo'))
            self.__rc_edge.append(('s_stereo', 'p_stereo'))

    __whash = set()
    __shash = set()
    __rc_node = [('s_charge', 'p_charge')]
    __rc_edge = [('s_bond', 'p_bond')]

    def set_hashlib(self, data):
        self.__whash = {x['CGR_FEAR_WHASH'] for x in data}
        self.__shash = {x['CGR_FEAR_SHASH'] for x in data}

    def check_cgr(self, g, full=True, gennew=False):
        chkd = set()
        report = set()

        def dochk(x):
            weights = self.get_morgan(x)
            whash = '.'.join(str(x) for x in sorted(weights.values()))
            if whash in self.__whash or gennew:
                shash = self.__get_smarts(x, weights)
                if shash in self.__shash:
                    chkd.update(x)
                    report.add((whash, shash))
                    return True
                elif gennew:
                    report.add((whash, shash))
            return False

        centers = self.get_environment(g, self.get_center_atoms(g), dante=True)
        ''' check for very specific centers
        '''
        for s in centers[:0:-1]:
            for c in connected_component_subgraphs(s):
                if not chkd.issuperset(c):
                    dochk(c)
        ''' final check for common centers
        '''
        tmp = (chkd.issuperset(c) or dochk(c) for c in connected_component_subgraphs(centers[0]))
        if gennew:
            tmp = list(tmp)
        return all(tmp) if full else any(tmp), report

    def get_cgr_string(self, g):
        return self.__get_smarts(g, self.get_morgan(g))

    def get_center_atoms(self, g):
        """ get atoms of reaction center (dynamic bonds, stereo or charges).
        """
        nodes = set()
        for n, node_attr in g.nodes(data=True):
            if any(node_attr.get(k) != node_attr.get(l) for k, l in self.__rc_node):
                nodes.add(n)

        for *n, node_attr in g.edges(data=True):
            if any(node_attr.get(k) != node_attr.get(l) for k, l in self.__rc_edge):
                nodes.update(n)

        return list(nodes)

    def get_environment(self, g, atoms, dante=False):
        """ get subgraph with atoms and their neighbors
        """
        nodes = [atoms]
        for i in range(self.__deep):
            nodes.append(set(chain.from_iterable(g.edges(nodes[i]))))

        if dante:
            centers = [g.subgraph(a) for a in nodes]
        else:
            centers = g.subgraph(nodes[-1])

        return centers

    def __get_smarts(self, g, weights):
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
            s_sh = '%s%s' % (self.__stereo and g.node[inter].get('s_stereo') or '',
                             self.__hyb_types[g.node[inter].get('s_hyb')] if self.__hyb else '')
            p_sh = '%s%s' % (self.__stereo and g.node[inter].get('p_stereo') or '',
                             self.__hyb_types[g.node[inter].get('p_hyb')] if self.__hyb else '')

            smis = ['[%s%s%s%s:%d]' %
                    (self.__isotope and g.node[inter].get('isotope') or '',
                     self.__element and g.node[inter]['element'] or '*',
                     s_sh and ';%s;' % s_sh or '',
                     self.__element and g.node[inter]['s_charge'] and '%+d' % g.node[inter]['s_charge'] or '',
                     newmaps.get(inter) or newmaps.setdefault(inter, next(countmap)))]
            smip = ['[%s%s%s%s:%d]' %
                    (self.__isotope and g.node[inter].get('isotope') or '',
                     self.__element and g.node[inter]['element'] or '*',
                     p_sh and ';%s;' % p_sh or '',
                     self.__element and g.node[inter]['p_charge'] and '%+d' % g.node[inter]['p_charge'] or '',
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
                        smis.append('%s%s%d' % (self.__to_smiles[g[inter][i].get('s_bond')],
                                                self.__stereo and g[inter][i].get('s_stereo', '') or '', cyc))
                        smip.append('%s%s%d' % (self.__to_smiles[g[inter][i].get('p_bond')],
                                                self.__stereo and g[inter][i].get('p_stereo', '') or '', cyc))
                    continue

                deep = dosmarts(set(chain(trace, [i])), i, inter)
                trace.update(deep[0])
                if deep[3]:
                    concat.extend(deep[3])
                    for j in deep[3]:
                        if j[0] == inter:
                            stoplist.append(j[2])
                            smis.append('%s%s%d' % (self.__to_smiles[g[inter][j[2]].get('s_bond')],
                                                    self.__stereo and g[inter][j[2]].get('s_stereo', '') or '', j[1]))
                            smip.append('%s%s%d' % (self.__to_smiles[g[inter][j[2]].get('p_bond')],
                                                    self.__stereo and g[inter][j[2]].get('p_stereo', '') or '', j[1]))
                smis.extend(['(' if iterlist else ''] +
                            ['%s%s' % (self.__to_smiles[g[inter][i].get('s_bond')],
                                       self.__stereo and g[inter][i].get('s_stereo', '') or '')] + deep[1] +
                            [')' if iterlist else ''])
                smip.extend(['(' if iterlist else ''] +
                            ['%s%s' % (self.__to_smiles[g[inter][i].get('p_bond')],
                                       self.__stereo and g[inter][i].get('p_stereo', '') or '')] + deep[2] +
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
        jssmiles = '.'.join(ssmiles)
        jpsmiles = '.'.join(psmiles)
        return '%s>>%s' % (jssmiles, jpsmiles) if jssmiles != jpsmiles else jssmiles

    __to_smiles = {1: '-', 2: '=', 3: '#', 4: ':', None: '.', 9: '~'}

    __hyb_types = {4: ',a', 3: ',t', 2: ',d', 1: ',s', None: ''}

    __stereo_types = {None: 0, 'u': 1, 'e': 2, 'z': 3, 'r': 4, 's': 5, 're': 6, 'si': 7}

    def get_morgan(self, g):
        newlevels = {}
        countprime = iter(self.__primes)

        params = {n: (self.__primes[elements.symbol(attr['element']).number] if self.__element else 1,
                      self.__primes[10 * attr['s_charge'] + attr['p_charge']] if self.__element else 1,
                      reduce(mul,
                             (self.__primes[10 * (eattr.get('s_bond') or 0) + (eattr.get('p_bond') or 0)]
                              for eattr in g[n].values()), 1),
                      self.__primes[attr['isotope']] if self.__isotope and 'isotope' in attr else 1,
                      self.__primes[10 * self.__stereo_types[attr.get('s_stereo')] +
                                    self.__stereo_types[attr.get('p_stereo')]]
                      if self.__stereo else 1,
                      reduce(mul,
                             (self.__primes[10 * self.__stereo_types[eattr.get('s_stereo')] +
                                            self.__stereo_types[eattr.get('p_stereo')]]
                              for eattr in g[n].values()), 1) if self.__stereo else 1)
                  for n, attr in g.nodes(data=True)}

        weights = {x: (newlevels.get(y) or newlevels.setdefault(y, next(countprime)))
                   for x, y in sorted(params.items(), key=itemgetter(1))}

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
                tmp[n] = reduce(mul, (weights[x] for x in m), weights[n]**2)

            numb = len(set(tmp.values()))
            if numb == oldnumb:
                x = Counter(tmp.values())
                stab += 1
                maxcount = x[max(x)]
            else:
                stab = 0
                maxcount = 0

            weights = {x: (neweights.get(y) or neweights.setdefault(y, next(countprime)))
                       for x, y in sorted(tmp.items(), key=itemgetter(1))}

        return weights
