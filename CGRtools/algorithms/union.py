# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#


class Union:
    def __or__(self, other):
        """
        G | H is union of graphs
        """
        return self.union(other)

    def union(self, other):
        if not isinstance(other, Union):
            raise TypeError('BaseContainer subclass expected')
        if self._node.keys() & set(other):
            raise KeyError('mapping of graphs is not disjoint')

        # dynamic container resolving
        qc = self._get_subclass('QueryCGRContainer')
        qq = self._get_subclass('QueryContainer')
        cc = self._get_subclass('CGRContainer')

        if isinstance(self, qc):
            u = type(self)()
        elif isinstance(other, qc):
            u = type(other)()
        elif isinstance(self, cc):
            if isinstance(other, qq):  # force QueryCGRContainer
                u = qc()
            else:
                u = type(self)()
        elif isinstance(other, cc):
            if isinstance(self, qq):
                u = qc()
            else:
                u = type(other)()
        elif isinstance(self, qq):  # self has precedence
            u = type(self)()
        elif isinstance(other, qq):
            u = type(other)()
        else:
            u = type(self)()

        for n, a in self.atoms():
            u.add_atom(a, n)
        for n, a in other.atoms():
            u.add_atom(a, n)

        for n, m, b in self.bonds():
            u.add_bond(n, m, b)
        for n, m, b in other.bonds():
            u.add_bond(n, m, b)
        return u


__all__ = ['Union']
