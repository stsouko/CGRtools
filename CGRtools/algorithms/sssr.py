# -*- coding: utf-8 -*-
#
#  Copyright 2017-2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from collections import defaultdict
from itertools import chain, combinations
from typing import Set, Dict, Union, Any, Tuple


class SSSR:
    """ SSSR calculation. based on idea of PID matrices from:
        Lee, C. J., Kang, Y.-M., Cho, K.-H., & No, K. T. (2009).
        A robust method for searching the smallest set of smallest rings with a path-included distance matrix.
        Proceedings of the National Academy of Sciences of the United States of America, 106(41), 17355–17358.
        http://doi.org/10.1073/pnas.0813040106
    """
    __slots__ = ()

    @cached_property
    def sssr(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Smallest Set of Smallest Rings

        :return rings atoms numbers
        """
        # ignore isolated atoms. optimization.
        return self._sssr(self._bonds)

    @classmethod
    def _sssr(cls, bonds: Dict[int, Union[Set[int], Dict[int, Any]]]) -> Tuple[Tuple[int, ...], ...]:
        """
        Smallest Set of Smallest Rings of any adjacency matrix
        """
        bonds = cls._skin_graph(bonds)
        if bonds:
            terminated, n_sssr = cls.__bfs(bonds)
            if n_sssr:
                return cls.__rings_filter(cls.__pid(terminated), n_sssr, bonds)
        return ()

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
    def __pid(terminated):
        # collect shortest and +1 paths
        pid1 = {}
        pid2 = {}
        pid1l = {}
        for j, paths_ticks in terminated.items():
            for paths, ticks in paths_ticks:
                for path in chain((paths,), (paths[x:] for x in ticks)):
                    k = (path[0], j)
                    k2 = (path[1], path[-2])
                    lp = len(path)
                    if k in pid1:
                        ls = pid1l[k]
                        if lp == ls:
                            pid1[k][k2] = path
                        elif ls - lp == 1:
                            pid2[k], pid1[k] = pid1[k], {k2: path}
                            pid1l[k] = lp
                        elif lp - ls == 1:
                            pid2[k][k2] = path
                        elif lp < ls:
                            pid1[k] = {k2: path}
                            pid2[k] = {}
                            pid1l[k] = lp
                    else:
                        pid1[k] = {k2: path}
                        pid2[k] = {}
                        pid1l[k] = lp

        c_set = []
        for k, p1ij in pid1.items():
            dij = pid1l[k] * 2 - 2
            p1ij = list(p1ij.values())
            p2ij = list(pid2[k].values())
            if len(p1ij) == 1:  # one shortest
                if not p2ij:  # need shortest + 1 path
                    continue
                c_set.append((dij + 1, p1ij, p2ij))
            elif not p2ij:  # one or more odd rings
                c_set.append((dij, p1ij, None))
            else:  # odd and even rings found (e.g. bicycle)
                p1ij = list(p1ij)
                c_set.append((dij, p1ij, None))
                c_set.append((dij + 1, p1ij, p2ij))

        for c_num, p1ij, p2ij in sorted(c_set):
            if c_num % 2:  # odd rings
                for c1 in p1ij:
                    for c2 in p2ij:
                        c = c1 + c2[-2:0:-1]
                        if len(set(c)) == len(c):
                            yield c
            else:
                for c1, c2 in zip(p1ij, p1ij[1:]):
                    c = c1 + c2[-2:0:-1]
                    if len(set(c)) == len(c):
                        yield c

    @staticmethod
    def __rings_filter(rings, n_sssr, bonds):
        hold_rings = {}  # rings with neighbours

        # step 1: collect isolated rings
        c = next(rings)
        if n_sssr == 1:
            return c,

        c_rings = {frozenset(c): c}
        for c in rings:
            ck = frozenset(c)

            if any(len(eck.intersection(ck)) > 1 for eck in c_rings):  # hold potential bicycle rings
                hold_rings[ck] = c
            else:
                c_rings[ck] = c
                if len(c_rings) == n_sssr:
                    return tuple(c_rings.values())

        # filter hold rings
        for ck, c in hold_rings.items():
            flag = False
            neighbors = {x: set() for x in c_rings if len(x.intersection(ck)) > 1}

            # check for bicycle with dead ends.
            # for example bfs found (124653) and (12478953) with lengths 6 and 8, but (478956) is smaller.
            #
            #   2 -- 4 -- 7
            #  /      \    \
            # 1       6     8
            #  \      /    /
            #   3 -- 5 -- 9
            for eck in neighbors:
                mc = ck ^ eck
                if len(mc) < len(c) - 2:
                    mb = {n for n in eck & ck if not mc.isdisjoint(bonds[n])}
                    if len(mb) == 2:
                        ck = mc | mb
                        c = tuple(n for n in c if n in ck)
                        ec = tuple(n for n in c_rings[eck] if n in ck)
                        if len(ec) < 3:
                            continue
                        elif ec[0] == c[0]:
                            if ec[-1] != c[-1]:
                                continue
                            c += ec[-2: 0: -1]
                        else:
                            if ec[0] != c[-1] or ec[-1] != c[0]:
                                continue
                            c += ec[1:-1]
                        if not flag:
                            flag = True
            if flag:
                c_rings[ck] = c
                if len(c_rings) == n_sssr:
                    return tuple(sorted(c_rings.values(), key=len))
                continue

            # create graph of connected neighbour rings
            for i, j in combinations(neighbors, 2):
                if len(i.intersection(j)) > 1:
                    neighbors[i].add(j)
                    neighbors[j].add(i)
            # check if hold rings is combination of existing. (123654) is combo of (1254) and (2365)
            #
            # 1--2--3
            # |  |  |
            # 4--5--6
            #
            # modified NX.dfs_labeled_edges
            # https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.traver\
            # sal.depth_first_search.dfs_labeled_edges.html
            depth_limit = len(neighbors) - 1
            for start, nbrs in neighbors.items():
                if not nbrs:
                    continue
                seen = {start}
                stack = [(start, depth_limit, iter(neighbors[start]))]
                while stack:
                    parent, depth_now, children = stack[-1]
                    try:
                        child = next(children)
                    except StopIteration:
                        stack.pop()
                    else:
                        if child not in seen:
                            seen.add(child)
                            mc = parent ^ child
                            mb = {n for n in parent & child if not mc.isdisjoint(bonds[n])}
                            if len(mb) == 2:
                                mc |= mb
                                if ck == mc:  # macrocycle found
                                    break
                                if depth_now and len(mc) < len(c):
                                    stack.append((mc, depth_now - 1, iter(neighbors[child])))
                else:
                    continue
                break
            else:
                c_rings[ck] = c
                if len(c_rings) == n_sssr:
                    return tuple(sorted(c_rings.values(), key=len))


__all__ = ['SSSR']
