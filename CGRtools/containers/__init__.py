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
from collections import namedtuple
from .cgr import CGRContainer
from .molecule import MoleculeContainer
from .reaction import ReactionContainer


CGRTemplate = namedtuple('CGRTemplate', ['reagents', 'products', 'meta'])
MatchContainer = namedtuple('MatchContainer', ['mapping', 'meta', 'patch'])


class MergedReaction(namedtuple('MergedReaction', ['reagents', 'products'])):
    def copy(self):
        return MergedReaction(self.reagents.copy(), self.products.copy())

    def get_fear(self, isotope=False, stereo=False, hyb=False, element=True, flush_cache=False):
        """
        return mapless fear of reaction.
        CAUTION: if reaction contains CGRs. fear will be unobvious.

        :param isotope: set isotope marks
        :param stereo: set stereo marks
        :param hyb: set hybridization mark of atom
        :param element: set elements marks
        :param flush_cache: recalculate fear if True
        :return: string representation of Reaction
        """
        if flush_cache or self.__fear is None:
            r = self.reagents.get_fear(isotope=isotope, stereo=stereo, hyb=hyb, element=element)
            p = self.products.get_fear(isotope=isotope, stereo=stereo, hyb=hyb, element=element)
            self.__fear = '%s>>%s' % ('{%s}' % r if isinstance(self.reagents, CGRContainer) else r,
                                      '{%s}' % p if isinstance(self.products, CGRContainer) else p)
        return self.__fear

    def __str__(self):
        return self.get_fear(True, True)

    def __repr__(self):
        return '%s(%s, %s)' % (self.__class__.__name__, repr(self.reagents), repr(self.products))

    __fear = None
