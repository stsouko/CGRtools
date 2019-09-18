# -*- coding: utf-8 -*-
#
#  Copyright 2017-2019 Ramil Nugmanov <stsouko@live.ru>
#  Copyright 2017-2019 Timur Madzhidov <tmadzhidov@gmail.com> atom ordering algorithm
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
from CachedMethods import cached_property
from collections import Counter
from functools import reduce
from itertools import count
from logging import warning
from operator import mul, itemgetter
from typing import Dict


class Morgan:
    __slots__ = ()

    @cached_property
    def atoms_order(self) -> Dict[int, int]:
        """
        Morgan like algorithm for graph nodes ordering

        :return: dict of atom-weight pairs
        """
        atoms = self._atoms
        bonds = self._bonds
        if not len(atoms):  # for empty containers
            return {}
        elif len(atoms) == 1:  # optimize single atom containers
            return dict.fromkeys(atoms, 2)

        params = {n: (int(a), tuple(sorted(int(b) for b in bonds[n].values()))) for n, a in atoms.items()}
        return self._morgan(self._sorted_primed(params))

    @staticmethod
    def _sorted_primed(params: Dict[int, int]) -> Dict[int, int]:
        levels = {}
        iter_primes = iter(primes)
        primed = {}
        for x, y in sorted(params.items(), key=itemgetter(1)):
            try:
                primed[x] = levels[y]
            except KeyError:
                primed[x] = levels[y] = next(iter_primes)
        return primed

    def _morgan(self, weights: Dict[int, int]) -> Dict[int, int]:
        atoms = self._atoms
        bonds = self._bonds

        tries = len(atoms) * 4
        numb = len(set(weights.values()))
        stab = 0

        while tries:
            oldnumb = numb

            # weights[n] ** 2 NEED for differentiation of molecules like A-B or any other complete graphs.
            tmp = {n: reduce(mul, (weights[x] for x in m), weights[n] ** 2) for n, m in bonds.items()}
            weights = self._sorted_primed(tmp)

            numb = len(set(weights.values()))
            if numb == len(atoms):  # each atom now unique
                break
            elif numb == oldnumb:
                x = Counter(weights.values())
                if x[min(x)] > 1:
                    if stab == 3:
                        break
                elif stab >= 2:
                    break

                stab += 1
            elif stab:
                stab = 0

            tries -= 1
            if not tries and numb < oldnumb:
                warning('morgan. number of attempts exceeded. uniqueness has decreased. next attempt will be made')
                tries = 1
        else:
            warning('morgan. number of attempts exceeded')
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


primes = {x: n for n, x in zip(range(1000), _eratosthenes())}


__all__ = ['Morgan']
