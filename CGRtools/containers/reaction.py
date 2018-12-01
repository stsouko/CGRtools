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
from collections.abc import MutableSequence
from hashlib import md5, sha256
from .cgr import CGRContainer
from .common import BaseContainer
from .query import QueryCGRContainer


class MindfulList(MutableSequence):
    """list with self-checks of modification. need for control of ReactionContainer caches actuality"""
    def __init__(self, data=None):
        self.__data = [] if data is None else list(data)
        self.__check = True

    def get_state(self):
        """return True if structure data list changed from previous checking time"""
        tmp = self.__check
        self.__check = False
        return tmp

    def insert(self, index, obj):
        self.__check = True
        self.__data.insert(index, obj)

    def __delitem__(self, index):
        self.__check = True
        del self.__data[index]

    def __add__(self, other):
        return self.__class__(self.__data + list(other))

    def __getitem__(self, index):
        return self.__data[index]

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
    """reaction storage. contains reagents, products and reactants lists"""
    __slots__ = ('__reagents', '__products', '__reactants', '__meta', '__signatures', '__pickle')

    def __init__(self, reagents=None, products=None, reactants=None, meta=None):
        """
        new empty or filled reaction object creation

        :param reagents: list of MoleculeContainers [or other Structure Containers] in left side of reaction
        :param products: right side of reaction. see reagents
        :param reactants: middle side of reaction: solvents, catalysts, etc. see reagents
        :param meta: dictionary of metadata. like DTYPE-DATUM in RDF

        """
        self.__reagents = MindfulList(reagents)
        self.__products = MindfulList(products)
        self.__reactants = MindfulList(reactants)
        self.__meta = dict(meta) or {}
        self.__signatures = {}
        self.__pickle = None

    def __getitem__(self, item):
        if item == 'reagents':
            return self.__reagents
        elif item == 'products':
            return self.__products
        elif item == 'reactants':
            return self.__reactants
        elif item == 'meta':
            return self.__meta
        raise AttributeError('invalid attribute')

    def __getstate__(self):
        return dict(reagents=list(self.__reagents), products=list(self.__products), reactants=list(self.__reactants),
                    meta=self.meta)

    def __setstate__(self, state):
        self.__init__(**state)

    def pickle(self):
        """return json serializable reaction"""
        return dict(reagents=[x.pickle() for x in self.__reagents], products=[x.pickle() for x in self.__products],
                    reactants=[x.pickle() for x in self.__reactants], meta=self.meta)

    @classmethod
    def unpickle(cls, data):
        """convert json serializable reaction into ReactionContainer object instance"""
        return cls(reagents=[BaseContainer.unpickle(x) for x in data['reagents']],
                   products=[BaseContainer.unpickle(x) for x in data['products']],
                   reactants=[BaseContainer.unpickle(x) for x in data['reactants']],
                   meta=data['meta'])

    @property
    def reagents(self):
        """reagents list. see products"""
        return self.__reagents

    @property
    def products(self):
        """list of CGRs or/and Molecules in products side"""
        return self.__products

    @property
    def reactants(self):
        """reactants list. see products"""
        return self.__reactants

    @property
    def meta(self):
        """dictionary of metadata. like DTYPE-DATUM in RDF"""
        return self.__meta

    def copy(self):
        """
        get copy of object

        :return: ReactionContainer
        """
        return self.__class__(reagents=[x.copy() for x in self.__reagents], meta=self.__meta.copy(),
                              products=[x.copy() for x in self.__products],
                              reactants=[x.copy() for x in self.__reactants])

    def get_signature_hash(self, *args, **kwargs):
        """
        get 40bytes hash of signature string. see get_signature

        :return: bytes
        """
        bs = self.get_signature(*args, **kwargs).encode()
        return md5(bs).digest() + sha256(bs).digest()

    def get_signature(self, atom=True, isotope=False, stereo=False, hybridization=False, neighbors=False,
                      flush_cache=False):
        """
        return string representation of reaction with molecules
        in order same as in lists of reagents, reactants, products.
        CAUTION: if reaction contains CGRs. signature will not be obvious

        :param isotope: set isotope marks
        :param stereo: set stereo marks
        :param hybridization: set hybridization mark of atom
        :param neighbors: set neighbors count mark of atom
        :param atom: set elements marks
        :param flush_cache: recalculate signature if True
        """
        if flush_cache or self.__signatures is None or any(x.get_state() for x in
                                                           (self.__reagents, self.__reactants, self.__products)):
            self.__signatures = {}

        k = (atom, isotope, stereo, hybridization, neighbors)
        if k in self.__signatures:
            return self.__signatures[k]

        sig = []
        for ml in (self.__reagents, self.__reactants, self.__products):
            ms = []
            for m in ml:
                mol = m.get_signature(atom=atom, isotope=isotope, stereo=stereo, hybridization=hybridization,
                                      neighbors=neighbors)
                ms.append(f'{mol}' if isinstance(m, (CGRContainer, QueryCGRContainer)) else mol)
            sig.append('.'.join(ms))
        self.__signatures[k] = out = '>'.join(sig)
        return out

    def implicify_hydrogens(self):
        """
        remove explicit hydrogens if possible

        :return: number of removed hydrogens
        """
        total = 0
        for ml in (self.__reagents, self.__reactants, self.__products):
            for m in ml:
                if hasattr(m, 'implicify_hydrogens'):
                    total += m.implicify_hydrogens()
        if total:
            self.flush_cache()
        return total

    def explicify_hydrogens(self):
        """
        add explicit hydrogens to atoms

        :return: number of added atoms
        """
        total = 0
        for ml in (self.__reagents, self.__reactants, self.__products):
            for m in ml:
                if hasattr(m, 'explicify_hydrogens'):
                    total += m.explicify_hydrogens()
        if total:
            self.flush_cache()
        return total

    def aromatize(self):
        """
        convert structures to aromatic form. works only for Molecules

        :return: number of processed molecules
        """
        total = 0
        for ml in (self.__reagents, self.__reactants, self.__products):
            for m in ml:
                if hasattr(m, 'aromatize'):
                    if m.aromatize():
                        total += 1
        if total:
            self.flush_cache()
        return total

    def flush_cache(self):
        """clear cached signatures and representation strings. use if structures objects in reaction object changed"""
        self.__pickle = self.__signatures = None

    def __str__(self):
        return self.get_signature(isotope=True, stereo=True, hybridization=True, neighbors=True)

    def __repr__(self):
        if self.__pickle is None or any(x.get_state() for x in (self.__reagents, self.__reactants, self.__products)):
            self.__pickle = '%s.unpickle(%s)' % (type(self).__name__, self.pickle())
        return self.__pickle


__all__ = ['ReactionContainer']
