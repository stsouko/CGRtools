# -*- coding: utf-8 -*-
#
#  Copyright 2017, 2018 Ramil Nugmanov <stsouko@live.ru>
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
"""
implements all internal structures, which represents: molecules, reactions, CGR and over
"""
from collections import namedtuple
from .cgr import CGRContainer
from .molecule import MoleculeContainer
from .query import QueryContainer
from .reaction import ReactionContainer, MergedReaction


CGRTemplate = namedtuple('CGRTemplate', ['pattern', 'patch', 'meta'])
MatchContainer = namedtuple('MatchContainer', ['mapping', 'patch', 'meta'])


CGRTemplate.__doc__ = '''container for [sub]structure queries. 
                         contains query structure and [sub]structure for replacement of found atoms and bonds'''
CGRTemplate.pattern.__doc__ = 'query structure. CGRContainer'
CGRTemplate.patch.__doc__ = '''replacement structure. CGRContainer.
                               Atom-to-atom mapping can be intersect with query at least in one atom.
                               replacement example for ketones:

                               * pattern = C[C:1](=[O:2])C, patch = [C:1]=[N:2], result = C[C:1](=[N:2])C
                               * pattern = C[C:1](=[O:2])C, patch = [C:1]=N, result = C[C:1](=[O:2])(=N)C
                            '''

MatchContainer.__doc__ = '''container with [sub]structure query result'''
MatchContainer.patch.__doc__ = '''replacement structure. CGRContainer.
                                  remapped to queried structure patch from CGRTemplate'''
MatchContainer.mapping.__doc__ = '''dictionary of queried structure atoms (keys) mapped to query atoms (values)'''

CGRTemplate.meta.__doc__ = MatchContainer.meta.__doc__ = 'dictionary of metadata. like DTYPE-DATUM in RDF'
