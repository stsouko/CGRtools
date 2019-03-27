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
from networkx.drawing.layout import kamada_kawai_layout
from random import uniform


class Calculate2D:
    def calculate2d(self, force=False, scale=1):
        """
        recalculate 2d coordinates. currently rings can be calculated badly.

        :param scale: rescale calculated positions.
        :param force: ignore existing coordinates of atoms
        """
        dist = {}
        # length forces
        for n, m_bond in self._adj.items():
            dist[n] = {}
            for m in m_bond:
                dist[n][m] = .825

        # angle forces
        for n, m_bond in self._adj.items():
            if len(m_bond) == 2:  # single-single or single-double bonds has angle = 120, other 180
                (m1, b1), (m2, b2) = m_bond.items()
                dist[m1][m2] = dist[m2][m1] = 1.43 if b1.order + b2.order in (2, 3) else 1.7  # +.05
            elif len(m_bond) == 3:
                m1, m2, m3 = m_bond
                dist[m1][m2] = dist[m1][m3] = dist[m2][m3] = dist[m3][m2] = dist[m2][m1] = dist[m3][m1] = 1.43
            elif len(m_bond) == 4:
                #    1
                #
                # 2  X  4
                #
                #    3
                m1, m2, m3, m4 = m_bond
                dist[m1][m2] = dist[m1][m4] = dist[m2][m1] = dist[m2][m3] = 1.17
                dist[m3][m2] = dist[m3][m4] = dist[m4][m1] = dist[m4][m3] = 1.17
                dist[m1][m3] = dist[m3][m1] = dist[m2][m4] = dist[m4][m2] = 1.7  # +.05

        # cycle forces
        for r in self.sssr:
            if len(r) == 6:
                #    6
                #
                # 1     5
                #
                # 2     4
                #
                #    3
                m1, m2, m3, m4, m5, m6 = r
                dist[m1][m4] = dist[m4][m1] = dist[m2][m5] = dist[m5][m2] = dist[m3][m6] = dist[m6][m3] = 1.7  # +.05

        if force:
            pos = None
        else:
            pos = {n: (atom.x or uniform(0, .01), atom.y or uniform(0, .01)) for n, atom in self.atoms()}

        for n, xy in kamada_kawai_layout(self, dist=dict(dist), pos=pos, scale=scale).items():
            atom = self._node[n]
            atom.x, atom.y = xy

        self.flush_cache()


__all__ = ['Calculate2D']
