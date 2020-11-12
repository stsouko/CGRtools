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
from itertools import chain, combinations, product
from logging import warning
from operator import itemgetter
from typing import Any, Dict, Set, Tuple, Union
from ..exceptions import ImplementationError


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
        if self.rings_count:
            return self._sssr(self._bonds, self.rings_count)
        return ()

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
        condensed_rings = {ck}  # collection of contours of condensed rings
        hold = []
        for c in rings:
            ck = frozenset(c)
            if len(ck) != len(c) or ck in seen_rings:
                continue
            seen_rings.add(ck)

            # create graph of connected neighbour rings
            neighbors = {x: set() for x in chain(sssr, condensed_rings) if len(x & ck) > 1}
            if neighbors:
                for i, j in combinations(neighbors, 2):
                    if len(i & j) > 1:
                        neighbors[i].add(j)
                        neighbors[j].add(i)
                # check if hold rings is combination of existing. (123654) is combo of (1254) and (2365)
                #
                # 1--2--3
                # |  |  |
                # 4--5--6
                #
                # modified NX.dfs_labeled_edges
                # https://networkx.github.io/documentation/stable/reference/algorithms/generated/networkx.algorithms.\
                # traversal.depth_first_search.dfs_labeled_edges.html
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
                                common = parent & child
                                if len(common) > 2:  # only terminal common atoms required
                                    mc = parent ^ child | {n for n in common if len(common.intersection(bonds[n])) == 1}
                                else:  # neighbors at least have 2 common atoms
                                    mc = parent | child
                                if ck == mc:  # macrocycle found
                                    break
                                elif depth_now and 2 < len(mc) <= len(c):
                                    stack.append((mc, depth_now - 1, iter(neighbors[child]), {child} | seen))
                    else:
                        continue
                    break
                else:
                    # update condensed rings. need for cuban and fullerene like structures.
                    tmp = set()
                    ckc = ck
                    for r in condensed_rings:
                        common = r & ckc
                        lc = len(common)
                        if lc == 2:
                            ckc |= r
                        elif lc == len(ckc):
                            unique = r - ckc
                            ckc = unique | {n for n in common if not unique.isdisjoint(bonds[n])}
                        # elif lc == len(r):  # impossible?
                        #     ckc -= r
                        elif lc > 2:
                            ckc = r ^ ckc | {n for n in common if len(common.intersection(bonds[n])) == 1}
                            while True:
                                term = {n for n in ckc if len(ckc.intersection(bonds[n])) == 1}
                                if not term:
                                    break
                                ckc -= term
                        else:
                            tmp.add(r)
                    condensed_rings = tmp

                    if ckc != ck and ckc in seen_rings:
                        # check ring for full surrounding by other rings
                        # reduced to existing ring. finish reached?
                        neighbors = set()  # bonds of neighbors
                        for r in sssr:
                            if len(r & ck) > 1:
                                for n in r:
                                    for m in bonds[n]:
                                        if m in r:
                                            neighbors.add((n, m))
                        if (c[0], c[-1]) in neighbors and all(x in neighbors for x in zip(c, c[1:])):
                            condensed_rings.add(ck)  # add ring to condensed. required for combined rings detection.
                            hold.append(c)
                            continue

                    condensed_rings.add(ckc)
                    seen_rings.add(ckc)
                    sssr[ck] = c
                    if len(sssr) == n_sssr:
                        return tuple(sssr.values())
            else:
                condensed_rings.add(ck)
                sssr[ck] = c
                if len(sssr) == n_sssr:
                    return tuple(sssr.values())

        if len(hold) < n_sssr - len(sssr):
            raise ImplementationError
        elif len(hold) > n_sssr - len(sssr):
            warning('Number of closing rings more than SSSR count. Possible errors.')
        return (*sssr.values(), *hold[:n_sssr - len(sssr)])


__all__ = ['SSSR']
