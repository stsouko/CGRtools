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
from itertools import combinations, product
from networkx import shortest_simple_paths, NetworkXNoPath, number_connected_components
from ..cache import cached_property


class SSSR:
    """ SSSR calculation. based on idea of PID matrices from:
        Lee, C. J., Kang, Y.-M., Cho, K.-H., & No, K. T. (2009).
        A robust method for searching the smallest set of smallest rings with a path-included distance matrix.
        Proceedings of the National Academy of Sciences of the United States of America, 106(41), 17355–17358.
        http://doi.org/10.1073/pnas.0813040106
    """
    @cached_property
    def sssr(self):
        adj = self._adj
        if len(adj) < 3:
            return []

        atoms = {x for x, y in adj.items() if len(y)}  # ignore isolated atoms
        terminals = {x for x, y in adj.items() if len(y) == 1}
        if terminals:
            bubble = terminals
            while True:
                bubble = {y for x in bubble for y in adj[x].keys() - terminals if len(adj[y].keys() - terminals) < 2}
                if not bubble:
                    break
                terminals.update(bubble)
            atoms.difference_update(terminals)  # skip not-cycle chains

        if not atoms:
            return []

        n_sssr = self.bonds_count - len(self) + 1
        terminated = {}
        tail = atoms.pop()
        next_stack = {x: [[tail, x]] for x in adj[tail].keys() & atoms}

        while True:
            next_front = set()
            found_odd = set()
            stack, next_stack = next_stack, {}
            for broom in stack.values():
                tail = broom[0][-1]
                next_front.add(tail)
                neighbors = adj[tail].keys() & atoms
                if len(neighbors) == 1:
                    n = neighbors.pop()
                    if n in found_odd:
                        continue
                    next_broom = [branch + [n] for branch in broom]
                    if n in stack:  # odd rings
                        found_odd.add(tail)
                        if n in next_stack:
                            next_stack[n].extend(next_broom)
                        else:
                            next_broom.extend(stack[n])
                            terminated[n] = next_stack[n] = next_broom
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
                        next_broom = [[tail, n]]
                        for branch in broom:
                            next_broom.append(branch + [n])
                        if n in stack:  # odd rings
                            found_odd.add(tail)
                            if n in next_stack:
                                next_stack[n].extend(next_broom)
                            else:
                                next_broom.extend(stack[n])
                                terminated[n] = next_stack[n] = next_broom
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
                next_stack = {x: [[tail, x]] for x in adj[tail].keys() & atoms}

        if not n_sssr:
            return []

        pid1 = {}
        pid2 = {}
        for j, paths in terminated.items():
            for path in paths:
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
            dij = len(p1ij[0]) - 1
            p2ij = pid2[k]

            if not p2ij and len(p1ij) == 1:
                continue

            c_num = 2 * dij
            if p2ij:
                c_num += 1

            c_set.append((c_num, p1ij, p2ij))

        n_ringidx, c_sssr = 0, {}
        for c_num, p1ij, p2ij in sorted(c_set):
            if c_num % 2:
                c1 = p1ij[0]
                c11 = c1[1]
                for c2 in p2ij:
                    if c11 != c2[1]:
                        c = c1 + c2[-2:0:-1]
                        ck = tuple(sorted(c))
                        if ck not in c_sssr:
                            c_sssr[ck] = c
                            n_ringidx += 1
                        if n_ringidx == n_sssr:
                            return list(c_sssr.values())
            else:
                for c1, c2 in zip(p1ij, p1ij[1:]):
                    if c1[1] != c2[1]:
                        c = c1 + c2[-2:0:-1]
                        ck = tuple(sorted(c))
                        if ck not in c_sssr:
                            c_sssr[ck] = c
                            n_ringidx += 1
                        if n_ringidx == n_sssr:
                            return list(c_sssr.values())


__all__ = ['SSSR']
