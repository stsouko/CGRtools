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
from .morgan import primes
from ..cache import cached_property


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
    @cached_property
    def chiral_atoms(self):
        """
        dict of chiral atoms valued with ordered neighbors
        """
        chiral = self.__tetrahedron()
        return chiral

    def __allenes(self):
        # detect allenes and cis/trans double bonds.
        adj = defaultdict(list)  # carbon double bonds adjacency matrix
        for n, m, bond in self.bonds():
            if bond.order == 2 and self.atom(n).element == self.atom(m).element == 'C':
                adj[n].append(m)
                adj[m].append(n)

        terminals = {x for x, y in adj.items() if len(y) == 1}
        cis_trans = []
        allene = []

        while terminals:
            n = terminals.pop()
            m = adj[n][0]  # terminal atom has 1 neighbor always
            path = [n, m]
            while m not in terminals:
                n, m = m, next(x for x in adj[m] if x != n)
                path.append(m)
            terminals.discard(m)
            print(path)
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

    @cached_property
    def __potentially_tetrahedron(self):
        # tetrahedral should be single bonded and contain zero or one H
        chiral = []
        for n, m_bond in self._adj.items():
            if len(m_bond) < 3:
                continue
            if not all(x.order == 1 for x in m_bond.values()):
                continue
            if self.atom_total_h(n) > 1:  # more then one H impossible
                continue
            if self.atom_implicit_h(n) + len(m_bond) != 4:
                continue
            chiral.append(n)
        return chiral

    def __tetrahedron(self):
        chiral = {}
        weights = self.atoms_order
        for n in self.__potentially_tetrahedron:
            m_bond = self._adj[n]
            # chiral atom all times contain unique neighbors
            if len(m_bond) == len({weights[x] for x in m_bond}):
                chiral[n] = tuple(m_bond)  # neighbors order
        """
        if chiral:  # generate new weights
            weights = weights.copy()
            countprime = iter(list(primes)[len(self):])
            for n in chiral:
                weights[n] = next(countprime)
            weights = self._morgan(weights)
            for n in self.__potentially_tetrahedron:
                if n in chiral:
                    continue
                m_bond = self._adj[n]
                if len(m_bond) == len({weights[x] for x in m_bond}):
                    chiral[n] = tuple(m_bond)  # neighbors order
        """
        return chiral


__all__ = ['Stereo']
