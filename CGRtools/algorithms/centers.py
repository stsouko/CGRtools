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
        """ get a list of lists of atoms, each of them is a reaction center
        
        :param stereo:  
        """
        center = set()
        adj = defaultdict(set)
        for n, a in self._node.items():
            if stereo:
                if not a._reagent != a._product:
                    center.add(n)
            elif not a._reagent == a._product:
                center.add(n)

        for n, m, b in self._bonds():
            if stereo:
                if not b._reagent != b._product:
                    adj[n].add(m)
                    adj[m].add(n)
                    center.add(n)
                    center.add(m)
            elif not b._reagent == b._product:
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
