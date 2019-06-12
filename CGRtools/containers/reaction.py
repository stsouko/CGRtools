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
from CachedMethods import cached_method
from collections.abc import Iterable
from itertools import chain
from functools import reduce
from operator import or_
from .cgr import CGRContainer
from .common import Graph
from .molecule import MoleculeContainer
from .query import QueryContainer


class ReactionContainer:
    """
    reaction storage. contains reactants, products and reagents lists.

    reaction storages hashable and comparable. based on reaction unique signature (SMIRKS).
    for reactions with query containers hash and comparison may give errors due to non-uniqueness.
    query containers itself not support hashing and comparison.
    """
    __slots__ = ('__reactants', '__products', '__reagents', '__meta', '_arrow', '__dict__')

    def __init__(self, reactants=(), products=(), reagents=(), meta=None):
        """
        new empty or filled reaction object creation

        :param reactants: list of MoleculeContainers [or other Structure Containers] in left side of reaction
        :param products: right side of reaction. see reactants
        :param reagents: middle side of reaction: solvents, catalysts, etc. see reactants
        :param meta: dictionary of metadata. like DTYPE-DATUM in RDF

        """
        if not isinstance(reactants, Iterable) or isinstance(reactants, (str, bytes)) or \
                not isinstance(products, Iterable) or isinstance(products, (str, bytes)) or \
                not isinstance(reagents, Iterable) or isinstance(reagents, (str, bytes)):
            raise TypeError('iterator of molecules or CGRs or queries expected')
        reactants = tuple(reactants)
        products = tuple(products)
        reagents = tuple(reagents)
        if any(not isinstance(x, Graph) for x in chain(reactants, products, reagents)):
            raise TypeError('molecule or CGR or query expected')
        self.__reactants = reactants
        self.__products = products
        self.__reagents = reagents
        if meta is None:
            self.__meta = {}
        else:
            self.__meta = dict(meta)
        self._arrow = None

    def __getitem__(self, item):
        if item in ('reactants', 0):
            return self.__reactants
        elif item in ('products', 1):
            return self.__products
        elif item in ('reagents', 2):
            return self.__reagents
        elif item == 'meta':
            return self.__meta
        raise KeyError('invalid attribute')

    def __getstate__(self):
        return dict(reactants=self.__reactants, products=self.__products, reagents=self.__reagents, meta=self.__meta)

    def __setstate__(self, state):
        if next(iter(state)) == 'reagents':  # 3.0 compatibility
            state['reagents'], state['reactants'] = state['reactants'], state['reagents']
        self.__reactants = state['reactants']
        self.__products = state['products']
        self.__reagents = state['reagents']
        self.__meta = state['meta']

    @property
    def reactants(self):
        """reactants list. see products"""
        return self.__reactants

    @property
    def reagents(self):
        """reagents list. see products"""
        return self.__reagents

    @property
    def products(self):
        """list of CGRs or/and Molecules in products side"""
        return self.__products

    @property
    def meta(self):
        """dictionary of metadata. like DTYPE-DATUM in RDF"""
        return self.__meta

    def copy(self):
        """
        get copy of object

        :return: ReactionContainer
        """
        return self.__class__(reagents=(x.copy() for x in self.__reagents), meta=self.__meta.copy(),
                              products=(x.copy() for x in self.__products),
                              reactants=(x.copy() for x in self.__reactants))

    def implicify_hydrogens(self):
        """
        remove explicit hydrogens if possible

        :return: number of removed hydrogens
        """
        total = 0
        for m in chain(self.__reagents, self.__reactants, self.__products):
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
        for m in chain(self.__reagents, self.__reactants, self.__products):
            if hasattr(m, 'explicify_hydrogens'):
                total += m.explicify_hydrogens()
        if total:
            self.flush_cache()
        return total

    def aromatize(self):
        """
        convert structures to aromatic form. works only for Molecules

        :return: number of processed rings
        """
        total = 0
        for m in chain(self.__reagents, self.__reactants, self.__products):
            if hasattr(m, 'aromatize'):
                total += m.aromatize()
        if total:
            self.flush_cache()
        return total

    def standardize(self):
        """
        standardize functional groups and convert structures to aromatic form. works only for Molecules

        :return: number of processed molecules
        """
        total = 0
        for m in chain(self.__reagents, self.__reactants, self.__products):
            if hasattr(m, 'standardize'):
                total += m.standardize()
        if total:
            self.flush_cache()
        return total

    @cached_method
    def compose(self):
        """
        get CGR of reaction

        reagents will be presented as unchanged molecules
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

    def calculate2d(self):
        """
        recalculate 2d coordinates
        """
        for m in chain(self.__reagents, self.__reactants, self.__products):
            m.calculate2d()
        self.fix_positions()

    def fix_positions(self):
        """
        fix coordinates of molecules in reaction
        """
        shift_x = 0
        for m in self.__reactants:
            max_x = self.__fix_positions(m, shift_x, 0)
            shift_x = max_x + 1
        arrow_min = shift_x

        if self.__reagents:
            for m in self.__reagents:
                max_x = self.__fix_positions(m, shift_x, 1.5)
                shift_x = max_x + 1
        else:
            shift_x += 3
        arrow_max = shift_x - 1

        for m in self.__products:
            max_x = self.__fix_positions(m, shift_x, 0)
            shift_x = max_x + 1
        self._arrow = (arrow_min, arrow_max)
        self.flush_cache()

    @staticmethod
    def __fix_positions(molecule, shift_x, shift_y):
        plane = molecule._plane
        min_x = min(x for x, _ in plane.values()) - shift_x
        max_x = max(x for x, _ in plane.values()) - min_x
        min_y = min(y for _, y in plane.values()) - shift_y
        for n, (x, y) in plane.items():
            plane[n] = (x - min_x, y - min_y)
        return max_x

    @cached_method
    def __str__(self):
        """
        SMIRKS of reaction. query and CGR containers in reaction {surrounded by curly braces}
        """
        sig = []
        for ml in (self.__reactants, self.__reagents, self.__products):
            sig.append(self.__get_smiles(ml) if ml else '')
        return '>'.join(sig)

    @staticmethod
    def __get_smiles(molecules):
        mc = []
        cc = []
        qc = []
        qcc = []
        smiles = []
        for m in molecules:
            if isinstance(m, MoleculeContainer):
                mc.append(m)
            elif isinstance(m, CGRContainer):
                cc.append(m)
            elif isinstance(m, QueryContainer):
                qc.append(m)
            else:
                qcc.append(m)

        if mc:
            smiles.append(str(reduce(or_, mc)))
        if cc:
            smiles.append(str(reduce(or_, cc)))
        if qc:
            smiles.append(str(reduce(or_, qc)))
        if qcc:
            smiles.append(str(reduce(or_, qcc)))
        return '.'.join(smiles)

    def flush_cache(self):
        self.__dict__.clear()


__all__ = ['ReactionContainer']
