# -*- coding: utf-8 -*-
#
#  Copyright 2017-2019 Ramil Nugmanov <stsouko@live.ru>
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
from itertools import count
from logging import warning
from operator import mul, itemgetter
from ..cache import cached_property


class Morgan:
    @cached_property
    def atoms_order(self):
        """
        Morgan like algorithm for graph nodes ordering

        :return: dict of atom-weight pairs
        """
        if not len(self):  # for empty containers
            return {}
        elif len(self) == 1:  # optimize single atom containers
            return dict.fromkeys(self, 2)

        params = {n: (int(node), tuple(sorted(int(edge) for edge in self._adj[n].values())))
                  for n, node in self.atoms()}
        newlevels = {}
        countprime = iter(primes)
        weights = {x: newlevels.get(y) or newlevels.setdefault(y, next(countprime))
                   for x, y in sorted(params.items(), key=itemgetter(1))}

        tries = len(self) * 4

        numb = len(set(weights.values()))
        stab = 0

        while tries:
            oldnumb = numb
            neweights = {}
            countprime = iter(primes)

            # weights[n] ** 2 NEED for differentiation of molecules like A-B or any other complete graphs.
            tmp = {n: reduce(mul, (weights[x] for x in m), weights[n] ** 2) for n, m in self._adj.items()}

            weights = {x: (neweights.get(y) or neweights.setdefault(y, next(countprime)))
                       for x, y in sorted(tmp.items(), key=itemgetter(1))}

            numb = len(set(weights.values()))
            if numb == len(self):  # each atom now unique
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
