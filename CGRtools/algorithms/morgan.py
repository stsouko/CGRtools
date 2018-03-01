# -*- coding: utf-8 -*-
#
#  Copyright 2017, 2018 Ramil Nugmanov <stsouko@live.ru>
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
from sys import stderr
from ..exceptions import InvalidConfig
from ..periodictable import elements


def get_morgan(g, isotope=False, element=True, stereo=False, hybridization=False, neighbors=False, labels=('s', 'p')):
    """
    Morgan like algorithm for graph nodes ordering

    :param g: (CGR|Molecule)Container
    :param isotope: differentiate isotopes
    :param element: differentiate elements and charges
    :param stereo: differentiate stereo atoms and bonds
    :param hybridization: differentiate hybridization of atoms
    :param neighbors: differentiate neighbors of atoms. useful for queries structures
    :param labels: for MoleculeContainer usable only default value, None, ('s',) or 's'.
                   for CGRContainer labels allow to control ordering.

                   if labels = ('s',) or 's':
                       ordering is equal to reagents ordering of reaction from which CGR is created.
                   if labels = ('p',) or 'p':
                       ordering is equal to products ordering.
                   if labels = ('s', 'p') or 'sp' or None:
                       ordering is default for given CGR
                   labels = ('p', 's') or 'ps':
                       ordering equal to reversed CGR. CGR of back reaction.
    :return: dict of atom: weights
    """
    if labels is None:
        s_charge = 's_charge'
        p_charge = 'p_charge'
        s_radical = 's_radical'
        p_radical = 'p_radical'
        s_stereo = 's_stereo'
        p_stereo = 'p_stereo'
        s_bond = 's_bond'
        p_bond = 'p_bond'
        s_hyb = 's_hyb'
        p_hyb = 'p_hyb'
        s_neighbors = 's_neighbors'
        p_neighbors = 'p_neighbors'
    elif len(labels) == 2:
        s, p = labels
        if not (s == 's' and p == 'p' or s == 'p' and p == 's'):
            raise InvalidConfig('invalid labels')
        s_charge = '%s_charge' % s
        p_charge = '%s_charge' % p
        s_radical = '%s_radical' % s
        p_radical = '%s_radical' % p
        s_stereo = '%s_stereo' % s
        p_stereo = '%s_stereo' % p
        s_bond = '%s_bond' % s
        p_bond = '%s_bond' % p
        s_hyb = '%s_hyb' % s
        p_hyb = '%s_hyb' % p
        s_neighbors = '%s_neighbors' % s
        p_neighbors = '%s_neighbors' % p
    else:
        s = labels[0]
        if s not in 'sp':
            raise InvalidConfig('invalid labels')
        s_charge = '%s_charge' % s
        s_radical = '%s_radical' % s
        s_stereo = '%s_stereo' % s
        s_bond = '%s_bond' % s
        s_hyb = '%s_hyb' % s
        s_neighbors = '%s_neighbors' % s
        p_charge = p_radical = p_stereo = p_bond = p_hyb = p_neighbors = None

    params = {n: (elements.index(attr['element']) if element else 1,
                  attr.get('isotope', 1) if isotope else 1,
                  10 * attr[s_charge] + attr.get(p_charge, 0) if element else 1,
                  10 * (attr.get(s_radical) or 0) + (attr.get(p_radical) or 0) if element else 1,
                  10 * (attr.get(s_stereo) or 0) + (attr.get(p_stereo) or 0) if stereo else 1,
                  10 * (attr.get(s_hyb) or 0) + (attr.get(p_hyb) or 0) if hybridization else 1,
                  10 * (attr.get(s_neighbors) or 0) + (attr.get(p_neighbors) or 0) if neighbors else 1,
                  reduce(mul, (primes[10 * (eattr.get(s_bond) or 0) + (eattr.get(p_bond) or 0)]
                               for eattr in g[n].values() if p_bond or eattr.get(s_bond)), 1),
                  reduce(mul, (primes[10 * (eattr.get(s_stereo) or 0) + (eattr.get(p_stereo) or 0)]
                               for eattr in g[n].values() if p_bond or eattr.get(s_bond)), 1) if stereo else 1)
              for n, attr in g.nodes(data=True)}

    newlevels = {}
    countprime = iter(primes)
    weights = {x: newlevels.get(y) or newlevels.setdefault(y, next(countprime))
               for x, y in sorted(params.items(), key=itemgetter(1))}

    numb = len(set(weights.values()))
    stab = 0

    scaf = {}
    for n, m in g.adjacency():
        scaf[n] = tuple(i for i, j in m.items() if p_bond or j.get(s_bond))

    tries = len(g) * 4  # limit for searching
    while tries:
        oldnumb = numb
        neweights = {}
        countprime = iter(primes)

        tmp = {n: reduce(mul, (weights[x] for x in m), weights[n]) for n, m in scaf.items()}

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

        tries -= 1
        if not tries and numb < oldnumb:
            print('morgan. number of attempts exceeded. uniqueness has decreased. last attempt will be made',
                  file=stderr)
            tries = 1
    else:
        print('morgan. number of attempts exceeded', file=stderr)

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
