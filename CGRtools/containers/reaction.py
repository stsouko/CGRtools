# -*- coding: utf-8 -*-
#
#  Copyright 2017-2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from .molecule import MoleculeContainer
from .query import QueryContainer
from ..algorithms.depict import DepictReaction
from ..algorithms.standardize import StandardizeReaction


class ReactionContainer(StandardizeReaction, DepictReaction):
    """
    Reaction storage. Contains reactants, products and reagents lists.

    Reaction storage hashable and comparable. based on reaction unique signature (SMIRKS).
    """
    __slots__ = ('__reactants', '__products', '__reagents', '__meta', '__name', '_arrow', '_signs', '__dict__')

    def __init__(self, reactants: TIterable[Union[MoleculeContainer, QueryContainer, CGRContainer]] = (),
                 products: TIterable[Union[MoleculeContainer, QueryContainer, CGRContainer]] = (),
                 reagents: TIterable[Union[MoleculeContainer, QueryContainer, CGRContainer]] = (),
                 meta: Optional[Dict] = None, name: Optional[str] = None):
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
        if not any(issubclass(base_type, x) for x in (MoleculeContainer, QueryContainer, CGRContainer)):
            raise TypeError('MoleculeContainer or QueryContainer or CGRContainer expected')
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
        reaction._ReactionContainer__meta = cgr._Graph__meta.copy()
        reaction._ReactionContainer__name = cgr._Graph__name
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
        return dict(reactants=self.__reactants, products=self.__products, reagents=self.__reagents, meta=self.__meta,
                    name=self.__name)

    def __setstate__(self, state):
        if next(iter(state)) == 'reagents':  # 3.0 compatibility
            state['reagents'], state['reactants'] = state['reactants'], state['reagents']
        self.__reactants = state['reactants']
        self.__products = state['products']
        self.__reagents = state['reagents']
        self.__meta = state['meta']
        self.__name = state.get('name', '')  # 4.0.9 compatibility
        self._arrow = None
        self._signs = None

    @property
    def reactants(self) -> Tuple[Union[MoleculeContainer, QueryContainer, CGRContainer], ...]:
        return self.__reactants

    @property
    def reagents(self) -> Tuple[Union[MoleculeContainer, QueryContainer, CGRContainer], ...]:
        return self.__reagents

    @property
    def products(self) -> Tuple[Union[MoleculeContainer, QueryContainer, CGRContainer], ...]:
        return self.__products

    def molecules(self) -> Iterator[Union[MoleculeContainer, QueryContainer, CGRContainer]]:
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
        if len(name) > 80:
            raise ValueError('name should be string up to 80 symbols')
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

    @property
    def centers_list(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Union reaction centers by leaving or substitute group
        """
        if not self.__reactants or self.__products:
            return ()  # no rc
        elif not isinstance(self.__reactants[0], MoleculeContainer):
            raise TypeError('Only Molecules supported')

        reactants = reduce(or_, self.__reactants)
        products = reduce(or_, self.__products)
        cgr = reactants ^ products
        all_atoms = set(reactants) ^ set(products)
        all_groups = cgr.substructure(all_atoms).connected_components
        new_centers_list = list(cgr.centers_list)

        for x in all_groups:
            x = set(x)
            intersection = []
            for i, y in enumerate(new_centers_list):
                if not x.isdisjoint(y):
                    intersection.append(i)

            if len(intersection) > 1:
                union = []
                for i in reversed(intersection):
                    union.extend(new_centers_list.pop(i))
                new_centers_list.append(union)

        return tuple(tuple(x) for x in new_centers_list)

    def implicify_hydrogens(self) -> int:
        """
        Remove explicit hydrogens if possible

        :return: number of removed hydrogens
        """
        total = 0
        for m in self.molecules():
            if not isinstance(m, MoleculeContainer):
                raise TypeError('Only Molecules supported')
            total += m.implicify_hydrogens()
        if total:
            self.flush_cache()
        return total

    def explicify_hydrogens(self) -> int:
        """
        Add explicit hydrogens to atoms

        :return: number of added atoms
        """
        total = 0
        for m in self.molecules():
            if not isinstance(m, MoleculeContainer):
                raise TypeError('Only Molecules supported')
            total += m.explicify_hydrogens()
        if total:
            self.flush_cache()
        return total

    def thiele(self) -> bool:
        """
        Convert structures to aromatic form. Works only for Molecules.
        Return True if in any molecule found kekule ring
        """
        total = False
        for m in self.molecules():
            if not isinstance(m, MoleculeContainer):
                raise TypeError('Only Molecules supported')
            if m.thiele() and not total:
                total = True
        if total:
            self.flush_cache()
        return total

    def kekule(self) -> bool:
        """
        Convert structures to kekule form. Works only for Molecules.
        Return True if in any molecule found aromatic ring
        """
        total = False
        for m in self.molecules():
            if not isinstance(m, MoleculeContainer):
                raise TypeError('Only Molecules supported')
            if m.kekule() and not total:
                total = True
        if total:
            self.flush_cache()
        return total

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

    def clean2d(self):
        """
        Recalculate 2d coordinates
        """
        for m in self.molecules():
            m.clean2d()
        self.fix_positions()

    def fix_positions(self):
        """
        Fix coordinates of molecules in reaction
        """
        shift_x = 0
        reactants = self.__reactants
        amount = len(reactants) - 1
        signs = []
        for m in reactants:
            max_x = self.__fix_positions(m, shift_x)
            if amount:
                max_x += .2
                signs.append(max_x)
                amount -= 1
            shift_x = max_x + 1
        arrow_min = shift_x

        if self.__reagents:
            for m in self.__reagents:
                max_x = self.__fix_reagent_positions(m, shift_x)
                shift_x = max_x + 1
            if shift_x - arrow_min < 3:
                shift_x = arrow_min + 3
        else:
            shift_x += 3
        arrow_max = shift_x - 1

        products = self.__products
        amount = len(products) - 1
        for m in products:
            max_x = self.__fix_positions(m, shift_x)
            if amount:
                max_x += .2
                signs.append(max_x)
                amount -= 1
            shift_x = max_x + 1
        self._arrow = (arrow_min, arrow_max)
        self._signs = tuple(signs)
        self.flush_cache()

    @staticmethod
    def __fix_reagent_positions(molecule, shift_x):
        plane = molecule._plane
        shift_y = .5

        values = plane.values()
        min_x = min(x for x, _ in values) - shift_x
        max_x = max(x for x, _ in values) - min_x
        min_y = min(y for _, y in values) - shift_y
        for n, (x, y) in plane.items():
            plane[n] = (x - min_x, y - min_y)
        return max_x

    @staticmethod
    def __fix_positions(molecule, shift_x):
        plane = molecule._plane
        atoms = molecule._atoms

        values = plane.values()
        min_x = min(x for x, _ in values) - shift_x

        right_atom, right_atom_plane = max((x for x in plane.items()), key=lambda x: x[1][0])
        max_x = right_atom_plane[0]
        max_x -= min_x

        min_y = min(y for _, y in values)
        max_y = max(y for _, y in values)
        mean_y = (max_y + min_y) / 2
        for n, (x, y) in plane.items():
            plane[n] = (x - min_x, y - mean_y)

        r_y = plane[right_atom][1]
        if isinstance(molecule, MoleculeContainer) and len(atoms[right_atom].atomic_symbol) == 2 and -.18 <= r_y <= .18:
            factor = molecule._hydrogens[right_atom]
            if factor == 1:
                max_x += .15
            elif factor:
                max_x += .25
        return max_x

    def __eq__(self, other):
        return isinstance(other, ReactionContainer) and str(self) == str(other)

    @cached_method
    def __hash__(self):
        return hash(str(self))

    @cached_method
    def __bytes__(self):
        return sha512(str(self).encode()).digest()

    @cached_method
    def __str__(self):
        sig = []
        for ml in (self.__reactants, self.__reagents, self.__products):
            sig.append(''.join(sorted(str(x) for x in ml)))
        return '>'.join(sig)

    def __format__(self, format_spec):
        sig = []
        for ml in (self.__reactants, self.__reagents, self.__products):
            sig.append(''.join(sorted(format(x, format_spec) for x in ml)))
        return '>'.join(sig)

    def flush_cache(self):
        self.__dict__.clear()
        for m in self.molecules():
            m.flush_cache()


__all__ = ['ReactionContainer']
