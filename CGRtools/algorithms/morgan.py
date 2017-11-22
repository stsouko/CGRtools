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
from functools import reduce
from itertools import count
from operator import mul, itemgetter
from ..periodictable import elements


def get_morgan(g, isotope=False, element=True, stereo=False):
    newlevels = {}
    countprime = iter(primes)

    params = {n: (primes[elements.index(attr['element'])] if element else 1,
                  primes[attr['isotope']] if isotope and 'isotope' in attr else 1,
                  primes[10 * attr['s_charge'] + attr.get('p_charge', 0)] if element else 1,
                  primes[10 * (attr.get('s_radical') or 0) + (attr.get('p_radical') or 0)] if element else 1,
                  primes[10 * (attr.get('s_stereo') or 0) + (attr.get('p_stereo') or 0)] if stereo else 1,
                  reduce(mul, (primes[10 * (eattr.get('s_bond') or 0) + (eattr.get('p_bond') or 0)]
                               for eattr in g[n].values()), 1),
                  reduce(mul, (primes[10 * (eattr.get('s_stereo') or 0) + (eattr.get('p_stereo') or 0)]
                               for eattr in g[n].values()), 1) if stereo else 1)
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
            tmp[n] = reduce(mul, (weights[x] for x in m), weights[n])

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


def _eratosthenes():
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


primes = tuple(x for _, x in zip(range(1000), _eratosthenes()))
