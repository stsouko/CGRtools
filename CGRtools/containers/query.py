# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
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
from . import cgr, molecule  # cyclic imports resolve
from .common import Graph


class QueryContainer(Graph):
    def _matcher(self, other):
        """
        QueryContainer < MoleculeContainer
        QueryContainer < QueryContainer[more general]
        QueryContainer < QueryCGRContainer[more general]
        """
        if isinstance(other, molecule.MoleculeContainer):
            return GraphMatcher(other, self, lambda x, y: y == x, lambda x, y: y == x)
        elif isinstance(other, (QueryContainer, QueryCGRContainer)):
            return GraphMatcher(other, self, lambda x, y: x == y, lambda x, y: x == y)
        raise TypeError('only query-molecule, query-query or query-cgr_query possible')


class QueryCGRContainer(Graph):
    def _matcher(self, other):
        """
        QueryCGRContainer < CGRContainer
        QueryContainer < QueryCGRContainer[more general]
        """
        if isinstance(other, CGRContainer):
            return GraphMatcher(other, self, lambda x, y: y == x, lambda x, y: y == x)
        elif isinstance(other, QueryCGRContainer):
            return GraphMatcher(other, self, lambda x, y: x == y, lambda x, y: x == y)
        raise TypeError('only cgr_query-cgr or cgr_query-cgr_query possible')


__all__ = []
