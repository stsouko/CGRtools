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
from itertools import combinations, product
from operator import itemgetter
from typing import Any, Dict, Set, Tuple, Union


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
        Smallest Set of Smallest Rings.

        :return rings atoms numbers
        """
        return self._sssr(self._bonds, self.rings_count)

    @classmethod
    def _sssr(cls, bonds: Dict[int, Union[Set[int], Dict[int, Any]]], n_sssr: int) -> Tuple[Tuple[int, ...], ...]:
        """
        Smallest Set of Smallest Rings of any adjacency matrix.
        Number of rings required.
        """
        bonds = cls._skin_graph(bonds)
        return cls.__rings_filter(cls.__c_set(*cls.__make_pid(bonds)), n_sssr, bonds)

    @staticmethod
    def __make_pid(bonds):
        lb = len(bonds)
        pid1 = defaultdict(lambda: defaultdict(dict))
        pid2 = defaultdict(lambda: defaultdict(dict))
        distances = defaultdict(lambda: defaultdict(lambda: lb))
        for n, ms in bonds.items():
            dn = distances[n]
            pn = pid1[n]
            for m in ms:
                pn[m][(m, n)] = (n, m)
                dn[m] = 1

        for k in bonds:
            new_distances = defaultdict(dict)
            dk = distances[k]
            ndk = new_distances[k]
            for i in bonds:
                if i == k:
                    continue
                di = distances[i]
                ndi = new_distances[i]
                ndk[i] = ndi[k] = di[k]
                for j in bonds:
                    if j == k or j == i:
                        continue
                    ij = di[j]
                    ikj = di[k] + dk[j]
                    if ij - ikj == 1:  # A new shortest path == previous shortest path - 1
                        pid2[i][j] = pid1[i][j]
                        pid1[i][j] = {(ni, mj): ip[:-1] + jp for ((ni, _), ip), ((_, mj), jp) in
                                      product(pid1[i][k].items(), pid1[k][j].items())}
                        ndi[j] = ikj
                    elif ij > ikj:  # A new shortest path
                        pid2[i][j] = {}
                        pid1[i][j] = {(ni, mj): ip[:-1] + jp for ((ni, _), ip), ((_, mj), jp) in
                                      product(pid1[i][k].items(), pid1[k][j].items())}
                        ndi[j] = ikj
                    elif ij == ikj:  # Another shortest path
                        pid1[i][j].update({(ni, mj): ip[:-1] + jp for ((ni, _), ip), ((_, mj), jp) in
                                           product(pid1[i][k].items(), pid1[k][j].items())})
                        ndi[j] = ij
                    elif ikj - ij == 1:  # Shortest+1 path
                        pid2[i][j].update({(ni, mj): ip[:-1] + jp for ((ni, _), ip), ((_, mj), jp) in
                                           product(pid1[i][k].items(), pid1[k][j].items())})
                        ndi[j] = ij
                    else:
                        ndi[j] = ij
            distances = new_distances
        return pid1, pid2, distances

    @staticmethod
    def __c_set(pid1, pid2, pid1l):
        c_set = []
        seen = set()
        for i, p1i in pid1.items():
            seen.add(i)
            di = pid1l[i]
            p2i = pid2[i]

            for j, p1ij in p1i.items():
                if j in seen:
                    continue
                p1ij = list(p1ij.values())
                p2ij = list(p2i[j].values())
                dij = di[j] * 2

                if len(p1ij) == 1:  # one shortest
                    if not p2ij:  # need shortest + 1 path
                        continue
                    c_set.append((dij + 1, p1ij, p2ij))
                elif not p2ij:  # one or more odd rings
                    c_set.append((dij, p1ij, None))
                else:  # odd and even rings found (e.g. bicycle)
                    c_set.append((dij, p1ij, None))
                    c_set.append((dij + 1, p1ij, p2ij))

        for c_num, p1ij, p2ij in sorted(c_set, key=itemgetter(0)):
            if c_num % 2:  # odd rings
                for c1 in p1ij:
                    for c2 in p2ij:
                        yield c1 + c2[-2:0:-1]
            else:
                for c1, c2 in zip(p1ij, p1ij[1:]):
                    yield c1 + c2[-2:0:-1]

    @staticmethod
    def __rings_filter(rings, n_sssr, bonds):
        # step 1: collect isolated rings
        c = next(rings)
        if n_sssr == 1:
            return c,

        ck = frozenset(c)
        seen_rings = {ck}
        sssr = {ck: c}
        for c in rings:
            ck = frozenset(c)
            if len(ck) != len(c) or ck in seen_rings:
                continue

            # create graph of connected neighbour rings
            neighbors = {x: set() for x in sssr if len(x.intersection(ck)) > 1}
            seen_rings.add(ck)
            if neighbors:
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
                    stack = [(start, depth_limit, iter(neighbors[start]), {start})]
                    while stack:
                        parent, depth_now, children, seen = stack[-1]
                        try:
                            child = next(children)
                        except StopIteration:
                            stack.pop()
                        else:
                            if child not in seen:
                                unique = parent ^ child
                                common = parent & child
                                up = parent - common
                                uc = child - common

                                border = set()
                                for n in common:
                                    ms = bonds[n]
                                    if not up.isdisjoint(ms) and not uc.isdisjoint(ms):
                                        border.add(n)
                                if len(border) == 2:
                                    unique |= border
                                    if ck == unique:  # macrocycle found
                                        break
                                    if depth_now and len(unique) < len(c):
                                        stack.append((unique, depth_now - 1, iter(neighbors[child]), {child} | seen))
                    else:
                        continue
                    break
                else:
                    sssr[ck] = c
                    if len(sssr) == n_sssr:
                        return tuple(sssr.values())
            else:
                sssr[ck] = c
                if len(sssr) == n_sssr:
                    return tuple(sssr.values())


__all__ = ['SSSR']
