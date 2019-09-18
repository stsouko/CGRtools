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
from CachedMethods import cached_property
from itertools import chain
from typing import Set, List, Dict, Union, Any, Tuple


class SSSR:
    """ SSSR calculation. based on idea of PID matrices from:
        Lee, C. J., Kang, Y.-M., Cho, K.-H., & No, K. T. (2009).
        A robust method for searching the smallest set of smallest rings with a path-included distance matrix.
        Proceedings of the National Academy of Sciences of the United States of America, 106(41), 17355–17358.
        http://doi.org/10.1073/pnas.0813040106
    """
    __slots__ = ()

    @cached_property
    def skin_atoms(self) -> Tuple[int, ...]:
        """
        atoms of rings and rings linkers [without terminal atoms]
        """
        return tuple(self.__shave(self._bonds))

    @cached_property
    def sssr(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Smallest Set of Smallest Rings

        :return rings atoms numbers
        """
        # ignore isolated atoms. optimization.
        return self._sssr(self._bonds)

    @classmethod
    def _sssr(cls, bonds: Dict[int, Union[List[int], Set[int], Tuple[int, ...], Dict[int, Any]]]) -> \
            Tuple[Tuple[int, ...], ...]:
        """
        Smallest Set of Smallest Rings of any adjacency matrix
        """
        bonds = cls.__shave(bonds)
        if bonds:
            terminated, n_sssr = cls.__bfs(bonds)
            if n_sssr:
                return cls.__pid(terminated, n_sssr)
        return ()

    @staticmethod
    def __shave(bonds):
        bonds = {n: set(ms) for n, ms in bonds.items() if ms}
        while True:  # skip not-cycle chains
            try:
                n = next(n for n, ms in bonds.items() if len(ms) <= 1)
            except StopIteration:
                break
            for m in bonds.pop(n):
                bonds[m].discard(n)
        return bonds

    @staticmethod
    def __bfs(bonds):
        n_sssr = sum(len(x) for x in bonds.values()) // 2 - len(bonds) + 1
        atoms = set(bonds)
        terminated = {}
        tail = atoms.pop()
        next_stack = {x: [((tail, x), ())] for x in bonds[tail] & atoms}

        while True:
            next_front = set()
            found_odd = set()
            stack, next_stack = next_stack, {}
            for tail, broom in stack.items():
                next_front.add(tail)
                neighbors = bonds[tail] & atoms
                if len(neighbors) == 1:
                    n = neighbors.pop()
                    if n in found_odd:
                        continue
                    next_broom = [((*path, n), ticks) for path, ticks in broom]
                    if n in stack:  # odd rings
                        found_odd.add(tail)
                        if n in next_stack:
                            next_stack[n].extend(next_broom)
                        else:
                            stack[n].extend(next_broom)  # not visited
                            terminated[n] = stack[n]
                    elif n in next_stack:  # even rings
                        next_stack[n].extend(next_broom)
                        if n not in terminated:
                            terminated[n] = next_stack[n]
                    else:
                        next_stack[n] = next_broom
                elif neighbors:
                    for n in neighbors:
                        if n in found_odd:
                            continue
                        next_broom = [((*path, n), (*ticks, len(path) - 1)) for path, ticks in broom]
                        if n in stack:  # odd rings
                            found_odd.add(tail)
                            if n in next_stack:
                                next_stack[n].extend(next_broom)
                            else:
                                stack[n].extend(next_broom)  # not visited
                                terminated[n] = stack[n]
                        elif n in next_stack:  # even rings
                            next_stack[n].extend(next_broom)
                            if n not in terminated:
                                terminated[n] = next_stack[n]
                        else:
                            next_stack[n] = next_broom

            atoms.difference_update(next_front)
            if not atoms:
                break
            elif not next_stack:
                n_sssr += 1
                tail = atoms.pop()
                next_stack = {x: [((tail, x), ())] for x in bonds[tail] & atoms}
        return terminated, n_sssr

    @staticmethod
    def __pid(terminated, n_sssr):
        pid1 = {}
        pid2 = {}
        for j, paths_ticks in terminated.items():
            for paths, ticks in paths_ticks:
                for path in chain((paths,), (paths[x:] for x in ticks)):
                    i = path[0]
                    k = (i, j)
                    if k in pid1:
                        ls = len(pid1[k][0])
                        lp = len(path)
                        if lp == ls:
                            pid1[k].append(path)
                        elif ls - lp == 1:
                            pid2[k], pid1[k] = pid1[k], [path]
                        elif lp - ls == 1:
                            pid2[k].append(path)
                        elif lp < ls:
                            pid1[k] = [path]
                            pid2[k] = []
                    else:
                        pid1[k] = [path]
                        pid2[k] = []

        c_set = []
        for k, p1ij in pid1.items():
            dij = len(p1ij[0]) * 2 - 2
            p2ij = pid2[k]
            if len(p1ij) == 1:  # one shortest
                if not p2ij:  # need shortest + 1 path
                    continue
                c_set.append((dij + 1, p1ij, p2ij))
            elif not p2ij:  # one or more odd rings
                c_set.append((dij, p1ij, None))
            else:  # odd and even rings found (e.g. bicycle)
                c_set.append((dij, p1ij, None))
                c_set.append((dij + 1, p1ij, p2ij))

        c_sssr = {}
        for c_num, p1ij, p2ij in sorted(c_set):
            if c_num % 2:  # odd rings
                c1 = p1ij[0]  # any shortest acceptable. sssr is not a unique set of rings
                c11 = c1[1]
                c12 = c1[-2]
                for c2 in p2ij:
                    if c11 == c2[1] or c12 == c2[-2]:
                        continue
                    c = c1 + c2[-2:0:-1]
                    ck = tuple(sorted(c))
                    if ck not in c_sssr:
                        c_sssr[ck] = c
                        if len(c_sssr) == n_sssr:
                            return tuple(c_sssr.values())
            else:
                for c1, c2 in zip(p1ij, p1ij[1:]):
                    if c1[1] == c2[1] or c1[-2] == c2[-2]:
                        continue
                    c = c1 + c2[-2:0:-1]
                    ck = tuple(sorted(c))
                    if ck not in c_sssr:
                        c_sssr[ck] = c
                        if len(c_sssr) == n_sssr:
                            return tuple(c_sssr.values())


__all__ = ['SSSR']
