# -*- coding: utf-8 -*-
#
#  Copyright 2017, 2018 Ramil Nugmanov <stsouko@live.ru>
#  Copyright 2017 Timur Madzhidov <tmadzhidov@gmail.com>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from collections import Counter
from functools import reduce
from itertools import count, repeat
from operator import mul, itemgetter
from warnings import warn


def initial_weights(g, isotope=False, element=True, stereo=False, hybridization=False, neighbors=False):
    return {n: (atom.number if element else 0,
                atom.isotope if isotope else 0,  # isotope precedence over electron state
                (atom.charge, atom.multiplicity or 0) if element else 0,
                (atom.stereo or 0) if stereo else 0,
                (atom.hybridization or 0) if hybridization else 0,
                (atom.neighbors or 0) if neighbors else 0,
                tuple(sorted(bond.order or 0 for bond in g._adj[n].values())),
                tuple(sorted(bond.stereo or 0 for bond in g._adj[n].values())) if stereo else 0)
            for n, atom in g._node.items()}, {n: tuple(m) for n, m in g._adj.items()}


def initial_weights_cgr(g, isotope=False, element=True, stereo=False, hybridization=False, neighbors=False):
    return {n: (atom.number if element else 0,
                atom.isotope if isotope else 0,
                (atom.charge, atom.multiplicity or 0, atom.p_charge, atom.p_multiplicity or 0) if element else 0,
                (atom.stereo or 0, atom.p_stereo or 0) if stereo else 0,
                (atom.hybridization or 0, atom.p_hybridization or 0) if hybridization else 0,
                (atom.neighbors or 0, atom.p_neighbors or 0) if neighbors else 0,
                tuple(sorted((bond.order or 0, bond.p_order or 0) for bond in g._adj[n].values())),
                tuple(sorted((bond.stereo or 0, bond.p_stereo or 0) for bond in g._adj[n].values())) if stereo else 0)
            for n, atom in g._node.items()}, {n: tuple(m) for n, m in g._adj.items()}


def initial_weights_cgr_reversed(g, isotope=False, element=True, stereo=False, hybridization=False, neighbors=False):
    return {n: (atom.number if element else 0,
                atom.isotope if isotope else 0,
                (atom.p_charge, atom.p_multiplicity or 0, atom.charge, atom.multiplicity or 0) if element else 0,
                (atom.p_stereo or 0, atom.stereo or 0) if stereo else 0,
                (atom.p_hybridization or 0, atom.hybridization or 0) if hybridization else 0,
                (atom.p_neighbors or 0, atom.neighbors or 0) if neighbors else 0,
                tuple(sorted((bond.p_order or 0, bond.order or 0) for bond in g._adj[n].values())),
                tuple(sorted((bond.p_stereo or 0, bond.stereo or 0) for bond in g._adj[n].values())) if stereo else 0)
            for n, atom in g._node.items()}, {n: tuple(m) for n, m in g._adj.items()}


def initial_weights_cgr_reagents(g, isotope=False, element=True, stereo=False, hybridization=False, neighbors=False):
    return {n: (atom.number if element else 0,
                atom.isotope if isotope else 0,  # isotope precedence over electron state
                (atom.charge, atom.multiplicity or 0) if element else 0,
                (atom.stereo or 0) if stereo else 0,
                (atom.hybridization or 0) if hybridization else 0,
                (atom.neighbors or 0) if neighbors else 0,
                tuple(sorted(bond.order or 0 for bond in g._adj[n].values())),
                tuple(sorted(bond.stereo or 0 for bond in g._adj[n].values() if bond.order)) if stereo else 0)
            for n, atom in g._node.items()}, {n: tuple(i for i, j in m.items() if j.order) for n, m in g._adj.items()}


def initial_weights_cgr_products(g, isotope=False, element=True, stereo=False, hybridization=False, neighbors=False):
    return {n: (atom.number if element else 0,
                atom.isotope if isotope else 0,  # isotope precedence over electron state
                (atom.p_charge, atom.p_multiplicity or 0) if element else 0,
                (atom.p_stereo or 0) if stereo else 0,
                (atom.p_hybridization or 0) if hybridization else 0,
                (atom.p_neighbors or 0) if neighbors else 0,
                tuple(sorted(bond.p_order or 0 for bond in g._adj[n].values())),
                tuple(sorted(bond.p_stereo or 0 for bond in g._adj[n].values() if bond.p_order)) if stereo else 0)
            for n, atom in g._node.items()}, {n: tuple(i for i, j in m.items() if j.p_order) for n, m in g._adj.items()}


def get_morgan(g, initial, *args, **kwargs):
    """
    Morgan like algorithm for graph nodes ordering

    :param g: (CGR|Molecule)Container
    :param initial: callable which returns dict with keys from graph nodes and orderable values and
                    dict with keys from graph nodes and values as lists of neighbors atoms
                   for MoleculeContainer initial_weights function available.

                   for CGRContainer allowed to control ordering.
                   if initial_weights_cgr_reagents passed:
                       ordering is equal to reagents ordering of reaction from which CGR is created.
                   if initial_weights_cgr_products passed:
                       ordering is equal to products ordering.
                   if initial_weights_cgr passed:
                       ordering is default for given CGR
                   if initial_weights_cgr_reversed passed:
                       ordering equal to reversed CGR. CGR of back reaction.
    :param isotope: differentiate isotopes
    :param element: differentiate elements and charges
    :param stereo: differentiate stereo atoms and bonds
    :param hybridization: differentiate hybridization of atoms
    :param neighbors: differentiate neighbors of atoms. useful for queries structures
    :return: dict of atom-weight pairs
    """
    if not len(g):  # for empty containers
        return {}

    params, scaf = initial(g, *args, **kwargs)

    newlevels = {}
    countprime = iter(primes)
    weights = {x: newlevels.get(y) or newlevels.setdefault(y, next(countprime))
               for x, y in sorted(params.items(), key=itemgetter(1))}

    numb = len(set(weights.values()))
    stab = 0

    tries = len(g) * 4  # limit for searching
    while tries:
        oldnumb = numb
        neweights = {}
        countprime = iter(primes)

        # weights[n] ** 2 NEED for differentiation of molecules like A-B or any other complete graphs.
        tmp = {n: reduce(mul, (weights[x] for x in m), weights[n] ** 2) for n, m in scaf.items()}

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
            warn('morgan. number of attempts exceeded. uniqueness has decreased. last attempt will be made')
            tries = 1
    else:
        warn('morgan. number of attempts exceeded')

    return weights


def _tupled(a):
    if isinstance(a, list):
        return tuple(sorted(a))
    return a,


def _tupled_sp(s, p):
    if isinstance(s, list):
        if None in s:
            s = (0 if x is None else x for x in s)
        if p is None:
            p = repeat(0)
        elif None in p:
            p = (x or 0 for x in p)
        return tuple(10 * x + y for x, y in sorted(zip(s, p)))
    if s is None:
        s = 0
    if p is None:
        p = 0
    return 10 * s + p,


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
