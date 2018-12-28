# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ravil Mukhametgaleev <sonic-mc@mail.ru>
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
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


class Centers:
    def get_centers_list(self, stereo=False):
        """ get a list of lists of atoms of reaction centers

        :param stereo: take into account stereo changes
        """
        center = set()
        adj = defaultdict(set)
        for n, a in self._node.items():
            if a._reagent != a._product or stereo and a.sterep != a.p_stereo:
                center.add(n)

        for n, m, b in self._bonds():
            if b._reagent != b._product or stereo and b.sterep != b.p_stereo:
                adj[n].add(m)
                adj[m].add(n)
                center.add(n)
                center.add(m)

        out = []
        while center:
            n = center.pop()
            if n in adj:
                c = set(self.__plain_bfs(adj, n))
                out.append(list(c))
                center.difference_update(c)
            else:
                out.append([n])

        return out

    def get_center_atoms(self, stereo=False):
        """ get list of atoms of reaction center (atoms with dynamic: bonds, stereo, charges, radicals).
        """
        nodes = set()
        for n, atom in self._node.items():
            if atom._reagent != atom._product or stereo and atom.stereo != atom.p_stereo:
                nodes.add(n)

        for n, m, bond in self._bonds():
            if bond._reagent != bond._product or stereo and bond.stereo != bond.p_stereo:
                nodes.add(n)
                nodes.add(m)

        return list(nodes)

    @staticmethod
    def __plain_bfs(adj, source):
        """A fast BFS node generator"""
        seen = set()
        nextlevel = {source}
        while nextlevel:
            thislevel = nextlevel
            nextlevel = set()
            for v in thislevel:
                if v not in seen:
                    yield v
                    seen.add(v)
                    nextlevel.update(adj[v])


__all__ = ['Centers']
