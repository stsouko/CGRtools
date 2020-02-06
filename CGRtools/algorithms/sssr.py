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
        pid1 = {}
        pid2 = {}
        pid1l = {}
        for j, paths_ticks in terminated.items():
            for paths, ticks in paths_ticks:
                for path in chain((paths,), (paths[x:] for x in ticks)):
                    k = (path[0], j)
                    if k in pid1:
                        ls = pid1l[k]
                        lp = len(path)
                        if lp == ls:
                            pid1[k].add(path)
                        elif ls - lp == 1:
                            pid2[k], pid1[k] = pid1[k], {path}
                            pid1l[k] = lp
                        elif lp - ls == 1:
                            pid2[k].add(path)
                        elif lp < ls:
                            pid1[k] = {path}
                            pid2[k] = set()
                            pid1l[k] = lp
                    else:
                        pid1[k] = {path}
                        pid2[k] = set()
                        pid1l[k] = len(path)

        pidk = defaultdict(list)
        for k in pid1:
            pidk[k[0]].append(k)

        for k in list(pid1):
            i, j = k
            for path in chain(pid1[k], pid2[k]):
                path = path[::-1]
                for fk in pidk[i]:
                    if fk == k:
                        continue
                    for fpath in chain(pid1[fk], pid2[fk]):
                        m_path = path + fpath[1:]
                        if len(set(m_path)) != len(m_path):
                            continue  # skip noose
                        mk = (j, fk[1])
                        if mk in pid1:
                            ls = pid1l[mk]
                            lp = len(m_path)
                            if lp == ls:
                                pid1[mk].add(m_path)
                            elif ls - lp == 1:
                                pid2[mk], pid1[mk] = pid1[mk], {m_path}
                                pid1l[mk] = lp
                            elif lp - ls == 1:
                                pid2[mk].add(m_path)
                            elif lp < ls:
                                pid1[mk] = {m_path}
                                pid2[mk] = set()
                                pid1l[mk] = lp
                        else:
                            pid1[mk] = {m_path}
                            pid2[mk] = set()
                            pid1l[mk] = len(m_path)

        c_set = []
        for k, p1ij in pid1.items():
            dij = pid1l[k] * 2 - 2
            p2ij = pid2[k]
            if len(p1ij) == 1:  # one shortest
                if not p2ij:  # need shortest + 1 path
                    continue
                c_set.append((dij + 1, list(p1ij), p2ij))
            elif not p2ij:  # one or more odd rings
                c_set.append((dij, list(p1ij), None))
            else:  # odd and even rings found (e.g. bicycle)
                p1ij = list(p1ij)
                c_set.append((dij, p1ij, None))
                c_set.append((dij + 1, p1ij, p2ij))

        for c_num, p1ij, p2ij in sorted(c_set):
            if c_num % 2:  # odd rings
                c1 = p1ij[0]  # any shortest acceptable. sssr is not a unique set of rings
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
        c_rings = {}
        ck_filter = set()
        hold_rings = {}
        for c in rings:
            ck = frozenset(c)
            if ck in ck_filter:
                continue
            ck_filter.add(ck)

            neighbors = [eck for eck in c_rings if not eck.isdisjoint(ck)]
            if len(neighbors) > 1:
                n_atoms = set(chain.from_iterable(neighbors))
                if ck < n_atoms:
                    hold_rings[ck] = c
                else:
                    c_rings[ck] = c
                    if len(c_rings) == n_sssr:
                        return tuple(c_rings.values())
            else:
                c_rings[ck] = c
                if len(c_rings) == n_sssr:
                    return tuple(c_rings.values())

        for ck, c in hold_rings.items():
            lc = len(c)
            neighbors = {x: set() for x in c_rings if not x.isdisjoint(ck)}
            for i, j in combinations(neighbors, 2):
                if not i.isdisjoint(j):
                    neighbors[i].add(j)
                    neighbors[j].add(i)

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
                                if depth_now and len(mc) < lc:
                                    stack.append((mc, depth_now - 1, iter(neighbors[child])))
                else:
                    continue
                break
            else:
                c_rings[ck] = c
                if len(c_rings) == n_sssr:
                    return tuple(sorted(c_rings.values(), key=len))


__all__ = ['SSSR']
