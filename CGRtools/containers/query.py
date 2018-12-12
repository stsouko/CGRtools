# -*- coding: utf-8 -*-
#
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
from networkx.algorithms.isomorphism import GraphMatcher
from .cgr import CGRContainer
from .common import BaseContainer
from .molecule import MoleculeContainer
from ..algorithms import StringMolecule, StringCGR
from ..attributes import QueryAtom, QueryBond, DynQueryAtom, DynQueryBond


class QueryContainer(StringMolecule, BaseContainer):
    node_attr_dict_factory = QueryAtom
    edge_attr_dict_factory = QueryBond

    def _matcher(self, other):
        """
        QueryContainer < MoleculeContainer
        QueryContainer < QueryContainer[more general]
        QueryContainer < QueryCGRContainer[more general]
        """
        if isinstance(other, MoleculeContainer):
            return GraphMatcher(other, self, lambda x, y: y == x, lambda x, y: y == x)
        elif isinstance(other, (QueryContainer, QueryCGRContainer)):
            return GraphMatcher(other, self, lambda x, y: x == y, lambda x, y: x == y)
        raise TypeError('only query-molecule, query-query or query-cgr_query possible')

    def __str__(self):
        return self.get_signature(isotope=True, stereo=True, hybridization=True, neighbors=True)


class QueryCGRContainer(StringCGR, BaseContainer):
    node_attr_dict_factory = DynQueryAtom
    edge_attr_dict_factory = DynQueryBond

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

    def __str__(self):
        return self.get_signature(isotope=True, stereo=True, hybridization=True, neighbors=True)


__all__ = ['QueryContainer', 'QueryCGRContainer']
