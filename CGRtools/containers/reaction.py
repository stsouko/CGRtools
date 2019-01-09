# -*- coding: utf-8 -*-
#
#  Copyright 2017-2019 Ramil Nugmanov <stsouko@live.ru>
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
from functools import reduce
from hashlib import sha512
from operator import mul, or_
from .cgr import CGRContainer
from .molecule import MoleculeContainer
from .query import QueryCGRContainer
from ..cache import cached_method


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

    def __str__(self):
        return '[%s]' % ', '.join(str(x) for x in self.__data)


class ReactionContainer:
    """reaction storage. contains reagents, products and reactants lists"""
    __slots__ = ('__reagents', '__products', '__reactants', '__meta', '__dict__')

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
        if meta is None:
            self.__meta = {}
        else:
            self.__meta = dict(meta)

    def __getitem__(self, item):
        if item in ('reagents', 0):
            return self.__reagents
        elif item in ('products', 1):
            return self.__products
        elif item in ('reactants', 2):
            return self.__reactants
        elif item == 'meta':
            return self.__meta
        raise AttributeError('invalid attribute')

    def __getstate__(self):
        return dict(reagents=list(self.__reagents), products=list(self.__products), reactants=list(self.__reactants),
                    meta=self.meta)

    def __setstate__(self, state):
        self.__init__(**state)

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

    def standardize(self):
        """
        standardize functional groups and convert structures to aromatic form. works only for Molecules

        :return: number of processed molecules
        """
        total = 0
        for ml in (self.__reagents, self.__reactants, self.__products):
            for m in ml:
                if hasattr(m, 'standardize'):
                    if m.standardize():
                        total += 1
        if total:
            self.flush_cache()
        return total

    def reset_query_marks(self):
        """
        set or reset hyb and neighbors marks to atoms.
        """
        for ml in (self.__reagents, self.__reactants, self.__products):
            for m in ml:
                if hasattr(m, 'reset_query_marks'):
                    m.reset_query_marks()
        self.flush_cache()

    @cached_method
    def compose(self):
        """
        get CGR of reaction

        reactants will be presented as unchanged molecules
        :return: CGRContainer
        """
        rr = self.__reagents + self.__reactants
        if rr:
            if not all(isinstance(x, (MoleculeContainer, CGRContainer)) for x in rr):
                raise TypeError('Queries not composable')
            r = reduce(or_, rr)
        else:
            r = MoleculeContainer()
        if self.__products:
            if not all(isinstance(x, (MoleculeContainer, CGRContainer)) for x in self.__products):
                raise TypeError('Queries not composable')
            p = reduce(or_, self.__products)
        else:
            p = MoleculeContainer()
        return r ^ p

    def __invert__(self):
        """
        get CGR of reaction
        """
        return self.compose()

    @cached_method
    def depict(self):
        pass  # todo: depict components

    def _repr_svg_(self):
        return self.depict()

    @cached_method
    def __str__(self):
        sig = []
        for ml in (self.__reagents, self.__reactants, self.__products):
            ms = []
            for m in sorted(ml, key=lambda x: reduce(mul, x.atoms_order) if hasattr(x, 'atoms_order') else 0):
                ms.append('{%s}' % m if isinstance(m, (CGRContainer, QueryCGRContainer)) else str(m))
            sig.append('.'.join(ms))
        return '>'.join(sig)

    @cached_method
    def __bytes__(self):
        return sha512(str(self).encode()).digest()

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(self) == str(other)

    def flush_cache(self):
        self.__dict__.clear()


__all__ = ['ReactionContainer']
