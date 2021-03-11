# -*- coding: utf-8 -*-
#
#  Copyright 2017-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import groupby
from logging import warning
from operator import itemgetter
from typing import Dict
from .._functions import tuple_hash


class Morgan:
    __slots__ = ()

    @cached_property
    def atoms_order(self) -> Dict[int, int]:
        """
        Morgan like algorithm for graph nodes ordering

        :return: dict of atom-order pairs
        """
        atoms = self._atoms
        if not atoms:  # for empty containers
            return {}
        elif len(atoms) == 1:  # optimize single atom containers
            return dict.fromkeys(atoms, 1)
        ring = self.ring_atoms
        return self._morgan({n: tuple_hash((hash(a), n in ring)) for n, a in atoms.items()})

    def _morgan(self, weights: Dict[int, int]) -> Dict[int, int]:
        atoms = self._atoms
        bonds = self._bonds

        tries = len(atoms) - 1
        numb = len(set(weights.values()))
        stab = old_numb = 0

        for _ in range(tries):
            weights = {n: tuple_hash((weights[n],
                                      *(x for x in sorted((weights[m], int(b)) for m, b in ms.items()) for x in x)))
                       for n, ms in bonds.items()}
            old_numb, numb = numb, len(set(weights.values()))
            if numb == len(atoms):  # each atom now unique
                break
            elif numb == old_numb:  # not changed. molecules like benzene
                if stab == 3:
                    break
                stab += 1
            elif stab:  # changed unique atoms number. reset stability check.
                stab = 0
        else:
            if numb < old_numb:
                warning('morgan. number of attempts exceeded. uniqueness has decreased.')

        return {n: i for i, (_, g) in enumerate(groupby(sorted(weights.items(), key=itemgetter(1)), key=itemgetter(1)),
                                                start=1) for n, _ in g}


__all__ = ['Morgan']
