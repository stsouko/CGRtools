# -*- coding: utf-8 -*-
#
#  Copyright 2019 Ramil Nugmanov <stsouko@live.ru>
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
from collections import defaultdict


def _pyramid_volume(n, u, v, w):
    nx, ny, nz = n
    ux, uy, uz = u
    vx, vy, vz = v
    wx, wy, wz = w
    ux -= nx
    uy -= ny
    uz -= nz
    vx -= nx
    vy -= ny
    vz -= nz
    wx -= nx
    wy -= ny
    wz -= nz
    return ux * (vy * wz - vz * wy) + uy * (vz * wx - vx * wz) + uz * (vx * wy - vy * wx)


class Stereo:
    def chiral(self):
        """
        list of tuples of chiral atom and their neighbors
        """
        weights = self.atoms_order
        chiral = []

        # detect allenes and cis/trans double bonds.
        adj = defaultdict(set)  # carbon double bonds adjacency matrix
        for n, m, bond in self._bonds():
            if bond.order == 2 and self._node[n].element == self._node[m].element == 'C':
                adj[n].add(m)
                adj[m].add(n)

        terminals = set(x for x, y in adj.items() if len(y) == 1)
        cis_trans = []
        allene = []

        while terminals:
            n = terminals.pop()
            path = list(plain_bfs(adj, n))
            terminals.difference_update(path)
            lp = len(path)
            if lp == 2:
                # ignore double bonds in small rings
                sp = set(path)
                for ring in self.sssr:
                    if not sp.issubset(ring) or len(ring) > 7:
                        cis_trans.append(path)
            elif lp == 3:  # simple allene
                allene.append(path)
            elif lp % 2:
                allene.append(path[::lp // 2])
            else:
                cis_trans.append(path[::lp - 1])

        # detect chiral allenes
        for x, n, y in allene:
            if len(self._adj[x]) > 1 and len(self._adj[y]) > 1:  # fast test
                []

        # tetrahedral should be single bonded and contain zero or one H
        for n, m_bond in self._adj.items():
            if n in adj:  # fast filter
                continue
            if not all(x.order == 1 for x in m_bond.values()):
                continue
            if len({weights[x] for x in m_bond}) != len(m_bond):  # chiral atom all times contain unique neighbors
                continue
            if self.atom_total_h(n) > 1:  # more then one H impossible
                continue

            if self.atom_implicit_h(n) + len(m_bond) == 4:
                chiral.append((n, *m_bond))
        return chiral


__all__ = ['Stereo']
