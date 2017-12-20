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
from ..algorithms import hash_cgr_string


class MindfulList:
    def __init__(self, data=None):
        self.__data = [] if data is None else data
        self.__check = True

    def get_state(self):
        tmp = self.__check
        self.__check = False
        return tmp

    def append(self, obj):
        self.__check = True
        self.__data.append(obj)

    def clear(self):
        self.__check = True
        self.__data.clear()

    def extend(self, iterable):
        self.__check = True
        self.__data.extend(iterable)

    def insert(self, index, obj):
        self.__check = True
        self.__data.insert(index, obj)

    def pop(self, index=None):
        self.__check = True
        return self.__data.pop(index)

    def __delitem__(self, index):
        self.__check = True
        del self.__data[index]

    def __iadd__(self, obj):
        self.__check = True
        self.__data.append(obj)
        return self

    def __getitem__(self, index):
        return self.__data[index]

    def __iter__(self):
        return iter(self.__data)

    def __len__(self):
        return len(self.__data)

    def __setitem__(self, key, value):
        self.__check = True
        self.__data[key] = value

    def __repr__(self):
        return '%s([%s])' % (self.__class__.__name__, ', '.join(repr(x) for x in self.__data))

    def __str__(self):
        return '[%s]' % ', '.join(str(x) for x in self.__data)


class ReactionContainer:
    __slots__ = ('__reagents', '__products', '__reactants', '__meta', '__fear', '__pickle')

    def __init__(self, substrats=None, products=None, reactants=None, meta=None, reagents=None):
        if substrats:
            warn('deprecated key. use reagents instead', DeprecationWarning)

        self.__reagents = MindfulList(substrats or reagents)
        self.__products = MindfulList(products)
        self.__reactants = MindfulList(reactants)
        self.__meta = meta or {}
        self.__fear = {}
        self.__pickle = None

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
            raise Exception('invalid key: %s' % item)

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

    def get_fear_hash(self, isotope=False, stereo=False, hyb=False, element=True, flush_cache=False):
        return hash_cgr_string(self.get_fear(isotope, stereo, hyb, element, flush_cache))

    def get_fear(self, isotope=False, stereo=False, hyb=False, element=True, flush_cache=False):
        """
        return fear of reaction with molecules in same order.
        CAUTION: if reaction contains CGRs. fear will be unobvious.

        :param isotope: set isotope marks
        :param stereo: set stereo marks
        :param hyb: set hybridization mark of atom
        :param element: set elements marks
        :param flush_cache: recalculate fear if True
        :return: string representation of Reaction
        """
        if flush_cache or self.__fear is None or any(x.get_state() for x in
                                                     (self.__reagents, self.__reactants, self.__products)):
            self.__fear = {}

        k = (isotope, element, stereo, hyb)
        if k not in self.__fear:
            self.__fear[k] = '%s>%s>%s' % tuple('.'.join('{%s}' % str(x) if isinstance(x, CGRContainer) else str(x)
                                                         for x in l)
                                                for l in (self.__reagents, self.__reactants, self.__products))

        return self.__fear[k]

    def flush_cache(self):
        self.__pickle = self.__fear = None

    def __str__(self):
        return self.get_fear(True, True)

    def __repr__(self):
        if self.__pickle is None or any(x.get_state() for x in (self.__reagents, self.__reactants, self.__products)):
            self.__pickle = '%s.unpickle(%s)' % (self.__class__.__name__, self.pickle())
        return self.__pickle
