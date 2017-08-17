# -*- coding: utf-8 -*-
#
#  Copyright 2017 Ramil Nugmanov <stsouko@live.ru>
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
from collections import Counter
from hashlib import md5, sha256
from itertools import chain, count
from functools import reduce
from operator import mul, itemgetter
from periodictable import elements


def eratosthenes():
    """Yields the sequence of prime numbers via the Sieve of Eratosthenes."""
    d = {}  # map each composite integer to its first-found prime factor
    for q in count(2):  # q gets 2, 3, 4, 5, ... ad infinitum
        p = d.pop(q, None)
        if p is None:
            # q not a key in D, so q is prime, therefore, yield it
            yield q
            # mark q squared as not-prime (with q as first-found prime factor)
            d[q * q] = q
        else:
            # let x <- smallest (N*p)+q which wasn't yet known to be composite
            # we just learned x is composite, with p first-found prime factor,
            # since p is the first-found prime factor of q -- find and mark it
            x = p + q
            while x in d:
                x += p
            d[x] = p


def hash_cgr_string(string):
    """
    concatenated md5 and sha256 hashes of cgr string 
    :param string: 
    :return: 48 bytes length string  
    """
    bs = string.encode()
    return md5(bs).digest() + sha256(bs).digest()

primes = tuple(x for _, x in zip(range(1000), eratosthenes()))
to_smiles = {1: '-', 2: '=', 3: '#', 4: ':', None: '.', 9: '~'}
hyb_types = {4: ',a', 3: ',t', 2: ',d', 1: ',s', None: ''}


def get_morgan(g, isotope=False, element=True):
    newlevels = {}
    countprime = iter(primes)

    params = {n: (primes[elements.symbol(attr['element']).number] if element else 1,
                  primes[attr['isotope']] if isotope and 'isotope' in attr else 1,
                  primes[10 * attr['s_charge'] + attr['p_charge']] if element else 1,
                  reduce(mul, (primes[10 * (eattr.get('s_bond') or 0) + (eattr.get('p_bond') or 0)]
                               for eattr in g[n].values()), 1))
              for n, attr in g.nodes(data=True)}

    weights = {x: newlevels.get(y) or newlevels.setdefault(y, next(countprime))
               for x, y in sorted(params.items(), key=itemgetter(1))}

    numb = len(set(weights.values()))
    stab = 0

    scaf = {}
    for n, m in g.adjacency():
        scaf[n] = tuple(m)

    while True:
        oldnumb = numb
        neweights = {}
        countprime = iter(primes)

        tmp = {}
        for n, m in scaf.items():
            """ if don't have neighbors use self weight
            """
            tmp[n] = reduce(mul, (weights[x] for x in m), weights[n]**2)

        weights = {x: (neweights.get(y) or neweights.setdefault(y, next(countprime)))
                   for x, y in sorted(tmp.items(), key=itemgetter(1))}

        numb = len(set(weights.values()))
        if numb == oldnumb:
            x = Counter(weights.values())
            if x[max(x)] > 1:
                if stab == 3:
                    break
            elif stab >= 2:
                break

            stab += 1
        elif stab:
            stab = 0

    return weights


def get_cgr_string(g, weights, isotope=False, stereo=False, hyb=False, element=True):
    newmaps = dict()
    countmap = count(1)
    countcyc = count(1)
    visited = set()
    stereo = False  # disable stereo. not implemented!

    def getnextatom(atoms):
        if len(atoms) == 1:
            nextatom = list(atoms)[0]
        else:
            nextatom = sorted((weights[i], i) for i in atoms)[-1][1]

        visited.add(nextatom)
        return nextatom

    def dosmarts(trace, inter, prev):
        s_sh = '%s%s' % (stereo and g.nodes[inter].get('s_stereo') or '',
                         hyb_types[g.nodes[inter].get('s_hyb')] if hyb else '')
        p_sh = '%s%s' % (stereo and g.nodes[inter].get('p_stereo') or '',
                         hyb_types[g.nodes[inter].get('p_hyb')] if hyb else '')

        smis = ['[%s%s%s%s:%d]' %
                (isotope and g.nodes[inter].get('isotope') or '', element and g.nodes[inter]['element'] or '*',
                 s_sh and ';%s;' % s_sh or '',
                 element and g.nodes[inter]['s_charge'] and '%+d' % g.nodes[inter]['s_charge'] or '',
                 newmaps.get(inter) or newmaps.setdefault(inter, next(countmap)))]
        smip = ['[%s%s%s%s:%d]' %
                (isotope and g.nodes[inter].get('isotope') or '', element and g.nodes[inter]['element'] or '*',
                 p_sh and ';%s;' % p_sh or '',
                 element and g.nodes[inter]['p_charge'] and '%+d' % g.nodes[inter]['p_charge'] or '',
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
                    smis.append('%s%s%d' % (to_smiles[g[inter][i].get('s_bond')],
                                            stereo and g[inter][i].get('s_stereo', '') or '', cyc))
                    smip.append('%s%s%d' % (to_smiles[g[inter][i].get('p_bond')],
                                            stereo and g[inter][i].get('p_stereo', '') or '', cyc))
                continue

            deep = dosmarts(set(chain(trace, [i])), i, inter)
            trace.update(deep[0])
            if deep[3]:
                concat.extend(deep[3])
                for j in deep[3]:
                    if j[0] == inter:
                        stoplist.append(j[2])
                        smis.append('%s%s%d' % (to_smiles[g[inter][j[2]].get('s_bond')],
                                                stereo and g[inter][j[2]].get('s_stereo', '') or '', j[1]))
                        smip.append('%s%s%d' % (to_smiles[g[inter][j[2]].get('p_bond')],
                                                stereo and g[inter][j[2]].get('p_stereo', '') or '', j[1]))
            smis.extend(['(' if iterlist else ''] +
                        ['%s%s' % (to_smiles[g[inter][i].get('s_bond')],
                                   stereo and g[inter][i].get('s_stereo', '') or '')] + deep[1] +
                        [')' if iterlist else ''])
            smip.extend(['(' if iterlist else ''] +
                        ['%s%s' % (to_smiles[g[inter][i].get('p_bond')],
                                   stereo and g[inter][i].get('p_stereo', '') or '')] + deep[2] +
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
