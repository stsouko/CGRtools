# -*- coding: utf-8 -*-
#
#  Copyright 2017-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from functools import reduce
from hashlib import sha512
from itertools import chain
from operator import or_
from typing import Dict, Iterable as TIterable, Iterator, Optional, Tuple, Union
from .cgr import CGRContainer
from .cgr_query import QueryCGRContainer
from .molecule import MoleculeContainer
from .query import QueryContainer
from ..algorithms.components import ReactionComponents
from ..algorithms.depict import DepictReaction
from ..algorithms.standardize import StandardizeReaction


graphs = Union[MoleculeContainer, QueryContainer, CGRContainer, QueryCGRContainer]


class ReactionContainer(StandardizeReaction, ReactionComponents, DepictReaction):
    """
    Reaction storage. Contains reactants, products and reagents lists.

    Reaction storage hashable and comparable. based on reaction unique signature (SMIRKS).
    """
    __slots__ = ('__reactants', '__products', '__reagents', '__meta', '__name', '_arrow', '_signs', '__dict__')
    __class_cache__ = {}

    def __init__(self, reactants: TIterable[graphs] = (), products: TIterable[graphs] = (),
                 reagents: TIterable[graphs] = (), meta: Optional[Dict] = None, name: Optional[str] = None):
        """
        New reaction object creation

        :param reactants: list of MoleculeContainers [or other Structure Containers] in left side of reaction
        :param products: right side of reaction. see reactants
        :param reagents: middle side of reaction: solvents, catalysts, etc. see reactants
        :param meta: dictionary of metadata. like DTYPE-DATUM in RDF

        """
        if not isinstance(reactants, Iterable) or isinstance(reactants, (str, bytes)) or \
                not isinstance(products, Iterable) or isinstance(products, (str, bytes)) or \
                not isinstance(reagents, Iterable) or isinstance(reagents, (str, bytes)):
            raise TypeError('Iterator of molecules or CGRs or queries expected')
        reactants = tuple(reactants)
        products = tuple(products)
        reagents = tuple(reagents)
        try:
            base_type = next(chain(reactants, products, reagents)).__class__
        except StopIteration:
            raise ValueError('At least one graph object required')
        if not any(issubclass(base_type, x) for x in
                   (MoleculeContainer, QueryContainer, CGRContainer, QueryCGRContainer)):
            raise TypeError('MoleculeContainer or QueryContainer or CGRContainer or QueryCGRContainer expected')
        elif not all(isinstance(x, base_type) for x in chain(reactants, products, reagents)):
            raise TypeError(f'{base_type.__name__} expected for all graphs')

        self.__reactants = reactants
        self.__products = products
        self.__reagents = reagents
        if meta is None:
            self.__meta = {}
        else:
            self.__meta = dict(meta)
        if name is None:
            self.__name = ''
        else:
            self.name = name
        self._arrow = None
        self._signs = None

    @classmethod
    def from_cgr(cls, cgr: CGRContainer) -> 'ReactionContainer':
        """
        Decompose CGR into reaction
        """
        if not isinstance(cgr, CGRContainer):
            raise TypeError('CGR expected')
        r, p = ~cgr
        reaction = object.__new__(cls)
        reaction._ReactionContainer__reactants = tuple(r.split())
        reaction._ReactionContainer__products = tuple(p.split())
        reaction._ReactionContainer__reagents = ()
        reaction._ReactionContainer__meta = cgr.meta.copy()
        reaction._ReactionContainer__name = cgr.name
        reaction._arrow = None
        reaction._signs = None
        return reaction

    def __getitem__(self, item):
        if item in ('reactants', 0):
            return self.__reactants
        elif item in ('products', 1):
            return self.__products
        elif item in ('reagents', 2):
            return self.__reagents
        elif item == 'meta':
            return self.__meta
        elif item == 'name':
            return self.__name
        raise KeyError('invalid attribute')

    def __getstate__(self):
        state = {'reactants': self.__reactants, 'products': self.__products, 'reagents': self.__reagents,
                 'meta': self.__meta, 'name': self.__name, 'arrow': self._arrow, 'signs': self._signs}

        if MoleculeContainer.__class_cache__.get('save_cache', False):
            state['cache'] = self.__dict__
        return state

    def __setstate__(self, state):
        if next(iter(state)) == 'reagents':  # 3.0 compatibility
            state['reagents'], state['reactants'] = state['reactants'], state['reagents']
        self.__reactants = state['reactants']
        self.__products = state['products']
        self.__reagents = state['reagents']
        self.__meta = state['meta']
        self.__name = state.get('name', '')  # 4.0.9 compatibility
        if 'signs' in state:  # >= 4.1.15
            self._arrow = state['arrow']
            self._signs = state['signs']
        else:
            self._arrow = None
            self._signs = None
        if 'cache' in state:  # >= 4.1.15
            self.__dict__.update(state['cache'])

    @classmethod
    def pickle_save_cache(cls, arg: bool):
        """
        Store cache of reaction into pickle for speedup loading
        """
        MoleculeContainer.__class_cache__['save_cache'] = arg

    @property
    def reactants(self) -> Tuple[graphs, ...]:
        return self.__reactants

    @property
    def reagents(self) -> Tuple[graphs, ...]:
        return self.__reagents

    @property
    def products(self) -> Tuple[graphs, ...]:
        return self.__products

    def molecules(self) -> Iterator[graphs]:
        """
        Iterator of all reaction molecules
        """
        return chain(self.__reactants, self.__reagents, self.__products)

    @property
    def meta(self) -> Dict:
        """
        Dictionary of metadata.
        Like DTYPE-DATUM in RDF
        """
        return self.__meta

    @property
    def name(self) -> str:
        return self.__name

    @name.setter
    def name(self, name: str):
        if not isinstance(name, str):
            raise TypeError('name should be string up to 80 symbols')
        self.__name = name

    def copy(self) -> 'ReactionContainer':
        """
        Get copy of object
        """
        copy = object.__new__(self.__class__)
        copy._ReactionContainer__reactants = tuple(x.copy() for x in self.__reactants)
        copy._ReactionContainer__products = tuple(x.copy() for x in self.__products)
        copy._ReactionContainer__reagents = tuple(x.copy() for x in self.__reagents)
        copy._ReactionContainer__meta = self.__meta.copy()
        copy._ReactionContainer__name = self.__name
        copy._arrow = self._arrow
        copy._signs = self._signs
        return copy

    @cached_method
    def compose(self) -> CGRContainer:
        """
        Get CGR of reaction

        Reagents will be presented as unchanged molecules
        :return: CGRContainer
        """
        rr = self.__reagents + self.__reactants
        if rr:
            if not isinstance(rr[0], (MoleculeContainer, CGRContainer)):
                raise TypeError('Queries not composable')
            r = reduce(or_, rr)
        else:
            r = MoleculeContainer()
        if self.__products:
            if not isinstance(self.__products[0], (MoleculeContainer, CGRContainer)):
                raise TypeError('Queries not composable')
            p = reduce(or_, self.__products)
        else:
            p = MoleculeContainer()
        c = r ^ p
        c.meta.update(self.__meta)
        return c

    def __invert__(self):
        """
        Get CGR of reaction
        """
        return self.compose()

    def __eq__(self, other):
        return isinstance(other, ReactionContainer) and str(self) == str(other)

    @cached_method
    def __hash__(self):
        return hash(str(self))

    @cached_method
    def __bytes__(self):
        return sha512(str(self).encode()).digest()

    def __bool__(self):
        """
        Exists both reactants and products
        """
        return bool(self.__reactants and self.__products)

    @cached_method
    def __str__(self):
        sig = []
        for ml in (self.__reactants, self.__reagents, self.__products):
            sig.append('.'.join(sorted(str(x) for x in ml)))
        return '>'.join(sig)

    def __format__(self, format_spec):
        """
        :param format_spec: see specification of nested containers.
            !c - Keep nested containers order
        """
        sig = []
        if '!c' in format_spec:
            for ml in (self.__reactants, self.__reagents, self.__products):
                sig.append('.'.join(format(x, format_spec) for x in ml))
        else:
            for ml in (self.__reactants, self.__reagents, self.__products):
                sig.append('.'.join(sorted(format(x, format_spec) for x in ml)))
        return '>'.join(sig)

    @cached_method
    def __len__(self):
        return len(self.__reactants) + len(self.__products) + len(self.__reagents)

    def flush_cache(self):
        self.__dict__.clear()
        for m in self.molecules():
            m.flush_cache()


__all__ = ['ReactionContainer']
