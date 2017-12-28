# -*- coding: utf-8 -*-
#
#  Copyright 2017 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
"""
implements all internal structures, which represents: molecules, reactions, CGR and over
"""
from collections import namedtuple
from .cgr import CGRContainer
from .molecule import MoleculeContainer
from .reaction import ReactionContainer
from ..algorithms import hash_cgr_string


CGRTemplate = namedtuple('CGRTemplate', ['pattern', 'patch', 'meta'])
MatchContainer = namedtuple('MatchContainer', ['mapping', 'patch', 'meta'])


class MergedReaction(namedtuple('MergedReaction', ['reagents', 'products', 'meta'])):
    """represent reactions as single disjointed reagents and single disjointed products graphs"""
    __slots__ = ()

    def copy(self):
        """
        get copy of object

        :return: MergedReaction
        """
        return self.__class__(self.reagents.copy(), self.products.copy(), self.meta.copy())

    def get_fear_hash(self, isotope=False, stereo=False, hyb=False, element=True, flush_cache=False):
        """
        get 40bytes hash of fear string. see get_fear

        :return: bytes
        """
        return hash_cgr_string(self.get_fear(isotope, stereo, hyb, element, flush_cache))

    def get_fear(self, isotope=False, stereo=False, hyb=False, element=True, flush_cache=False):
        """
        return string representation of reaction with unique atoms and molecules order
        CAUTION: if reaction contains CGRs. fear will be unobvious

        :param isotope: set isotope marks to string
        :param stereo: set stereo marks
        :param hyb: set hybridization mark of atom
        :param element: set elements marks and charges of atoms
        :param flush_cache: recalculate fear if True
        """
        if flush_cache or self.__fear is None:
            r = self.reagents.get_fear(isotope=isotope, stereo=stereo, hyb=hyb, element=element)
            p = self.products.get_fear(isotope=isotope, stereo=stereo, hyb=hyb, element=element)
            self.__fear = '%s>>%s' % ('{%s}' % r if isinstance(self.reagents, CGRContainer) else r,
                                      '{%s}' % p if isinstance(self.products, CGRContainer) else p)
        return self.__fear

    def flush_cache(self):
        """clear cached fears and representation strings. use if structure changed"""
        self.__pickle = self.__fear = None

    def __str__(self):
        return self.get_fear(True, True)

    def __repr__(self):
        if self.__pickle is None:
            self.__pickle = '%s(%s, %s)' % (self.__class__.__name__, repr(self.reagents), repr(self.products))
        return self.__pickle

    __fear = __pickle = None


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

MergedReaction.reagents.__doc__ = 'disjointed reagents graph'
MergedReaction.products.__doc__ = 'disjointed products graph'

CGRTemplate.meta.__doc__ = MatchContainer.meta.__doc__ = MergedReaction.meta.__doc__ = \
    'dictionary of metadata. like DTYPE-DATUM in RDF'
