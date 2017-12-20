# -*- coding: utf-8 -*-
#
#  Copyright 2017 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
""" SSSR calculation. based on:
    Lee, C. J., Kang, Y.-M., Cho, K.-H., & No, K. T. (2009).
    A robust method for searching the smallest set of smallest rings with a path-included distance matrix.
    Proceedings of the National Academy of Sciences of the United States of America, 106(41), 17355–17358.
    http://doi.org/10.1073/pnas.0813040106
"""
from collections import defaultdict
from itertools import product, permutations


def find_sssr(g):
    """
    SSSR search.
    
    
    :param g: Molecule Container 
    :return: list of sets of edges or None. edges coded as nodes labels.
    """
    n_sssr = g.number_of_edges() - len(g) + 1
    if not n_sssr:
        return None

    n_ringidx, c_sssr = 0, []
    for c_num, p1ij, p2ij in _make_c_set(*_make_dist_pids_matrices(g)):
        if c_num % 2:
            c1 = p1ij[0]
            lc1 = len(c1)
            for c2 in p2ij:
                c = c1.symmetric_difference(c2)
                if len(c) == lc1 + len(c2):
                    c_sssr.append(c)
                    n_ringidx += 1
                if n_ringidx == n_sssr:
                    return c_sssr
        else:
            for c1, c2 in zip(p1ij, p1ij[1:]):
                c = c1.union(c2)
                if c not in c_sssr:
                    c_sssr.append(c)
                    n_ringidx += 1
                if n_ringidx == n_sssr:
                    return c_sssr


def _make_dist_pids_matrices(g):
    dist = defaultdict(lambda: defaultdict(lambda: float('inf')))
    pid1 = defaultdict(lambda: defaultdict(list))
    pid2 = defaultdict(lambda: defaultdict(list))
    for i in g:
        dist[i][i] = 0
    for i, j in g.edges():
        dist[i][j] = dist[j][i] = 1
        pid1[i][j] = pid1[j][i] = [{(i, j) if i < j else (j, i)}]

    for k, i, j in permutations(g, 3):
        pid_sum = [x.union(y) for x, y in product(pid1[i][k], pid1[k][j])] or pid1[i][k] or pid1[k][j]
        if pid_sum:
            dij = dist[i][j]
            dist_sum = dist[i][k] + dist[k][j]
            if dij > dist_sum:
                if dij == dist_sum + 1:
                    pid2[i][j] = pid1[i][j]
                else:
                    pid2[i][j] = []

                dist[i][j] = dist_sum
                pid1[i][j] = pid_sum
            elif dij == dist_sum:
                    pid1[i][j].extend(pid_sum)
            elif dij == dist_sum - 1:
                pid2[i][j].extend(pid_sum)

    return dist, pid1, pid2


def _make_c_set(dist, pid1, pid2):
    c_set = []
    for i, j in permutations(dist, 2):
        dij = dist[i][j]
        p1ij = pid1[i][j]
        p2ij = pid2[i][j]
        if not p2ij and len(p1ij) == 1:
            continue

        c_num = 2 * dij
        if p2ij:
            c_num += 1

        c_set.append((c_num, p1ij, p2ij))

    return sorted(c_set)
