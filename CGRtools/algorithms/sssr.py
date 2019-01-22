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
from itertools import combinations
from networkx import shortest_simple_paths, NetworkXNoPath, number_connected_components
from ..cache import cached_property


class SSSR:
    """ SSSR calculation. based on idea from:
        Lee, C. J., Kang, Y.-M., Cho, K.-H., & No, K. T. (2009).
        A robust method for searching the smallest set of smallest rings with a path-included distance matrix.
        Proceedings of the National Academy of Sciences of the United States of America, 106(41), 17355–17358.
        http://doi.org/10.1073/pnas.0813040106
    """
    @cached_property
    def sssr(self):
        """
        SSSR search.

        :return: list of lists of rings nodes
        """
        n_sssr = self.bonds_count - len(self) + number_connected_components(self)
        if not n_sssr:
            return []

        pid1 = {}
        pid2 = {}

        for ij in combinations(self, 2):
            paths = shortest_simple_paths(self, *ij)  # slowest part of algorithm
            try:  # stop if path not reachable
                path = next(paths)
            except NetworkXNoPath:
                break

            pid1[ij] = [path]
            ls = len(path)
            lsp = ls + 1

            for path in paths:
                lp = len(path)
                if lp == ls:
                    pid1[ij].append(path)
                elif lp == lsp:
                    if ij not in pid2:
                        pid2[ij] = [path]
                    else:
                        pid2[ij].append(path)
                else:
                    break

        c_set = []
        for ij, p1ij in pid1.items():
            dij = len(p1ij[0]) - 1
            p2ij = pid2.get(ij)

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
                cs1 = set(c1)
                for c2 in p2ij:
                    if len(cs1.intersection(c2)) == 2:
                        c = c1 + c2[-2:0:-1]
                        ck = tuple(sorted(c))
                        if ck not in c_sssr:
                            c_sssr[ck] = c
                            n_ringidx += 1
                        if n_ringidx == n_sssr:
                            return list(c_sssr.values())
            else:
                for c1, c2 in zip(p1ij, p1ij[1:]):
                    if len(set(c1).intersection(c2)) == 2:
                        c = c1 + c2[-2:0:-1]
                        ck = tuple(sorted(c))
                        if ck not in c_sssr:
                            c_sssr[ck] = c
                            n_ringidx += 1
                        if n_ringidx == n_sssr:
                            return list(c_sssr.values())


__all__ = ['SSSR']
