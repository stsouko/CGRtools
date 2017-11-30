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
from warnings import warn
from .cgr import CGRContainer


class ReactionContainer(object):
    __slots__ = ('__reagents', '__products', '__reactants', '__meta')

    def __init__(self, substrats=None, products=None, reactants=None, meta=None, reagents=None):
        if substrats:
            warn('deprecated key. use reagents instead', DeprecationWarning)

        self.__reagents = substrats or reagents or []
        self.__products = products or []
        self.__reactants = reactants or []
        self.__meta = meta or {}

    def __getitem__(self, item):
        if item == 'substrats':
            warn('deprecated key. use reagents instead', DeprecationWarning)
            return self.__reagents
        elif item == 'reagents':
            return self.__reagents
        elif item == 'products':
            return self.__products
        elif item == 'reactants':
            return self.__reactants
        elif item == 'meta':
            return self.__meta
        else:
            raise Exception('Invalid key')

    def pickle(self):
        """ return json serializable reaction
        """
        return dict(reagents=[x.pickle() for x in self.__reagents], meta=self.meta,
                    products=[x.pickle() for x in self.__products], reactants=[x.pickle() for x in self.__reactants])

    @staticmethod
    def unpickle(data):
        """ convert json serializable reaction into ReactionContainer object instance
        """
        return ReactionContainer(reagents=[CGRContainer.unpickle(x) for x in
                                           (data['reagents'] if 'reagents' in data else data['substrats'])],
                                 products=[CGRContainer.unpickle(x) for x in data['products']], meta=data['meta'],
                                 reactants=[CGRContainer.unpickle(x) for x in data.get('reactants', [])])

    @property
    def substrats(self):  # reverse compatibility
        warn('deprecated key. use reagents instead', DeprecationWarning)
        return self.__reagents

    @property
    def reagents(self):
        return self.__reagents

    @property
    def products(self):
        return self.__products

    @property
    def reactants(self):
        return self.__reactants

    @property
    def meta(self):
        return self.__meta

    def copy(self):
        return ReactionContainer(reagents=[x.copy() for x in self.__reagents], meta=self.__meta.copy(),
                                 products=[x.copy() for x in self.__products],
                                 reactants=[x.copy() for x in self.__reactants])
