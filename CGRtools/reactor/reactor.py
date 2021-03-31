# -*- coding: utf-8 -*-
#
#  Copyright 2019-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Adelia Fatykhova <adelik21979@gmail.com>
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
from functools import reduce
from itertools import count, permutations
from logging import info
from operator import or_
from typing import Iterable
from .base import BaseReactor
from .._functions import lazy_product
from ..containers import QueryContainer, MoleculeContainer, ReactionContainer


class Reactor(BaseReactor):
    """
    Reactor for molecules transformations.
    Generates reaction from input molecules using
    transformation template (CGRtools ReactionContainer).

    Reactor calling transforms reactants to products and
    returns generator of reaction transformations with all
    possible reactions.
    """
    def __init__(self, template, *, delete_atoms: bool = True, one_shot: bool = True,
                 polymerise_limit: int = 0, automorphism_filter: bool = True):
        """
        :param template: CGRtools ReactionContainer
        :param delete_atoms: if True atoms exists in reactants but
                            not exists in products will be removed
        :param one_shot: do only single reaction center then True, else do all possible combinations of reactions.
        :param polymerise_limit: limit of self reactions. Zero by default - prevent polymerization.
        """
        reactants, products = template.reactants, template.products
        if not reactants or not products:
            raise ValueError('empty template')

        if isinstance(reactants[0], QueryContainer):
            self.__patterns = reactants
            products_ = reduce(or_, products)
        else:
            raise TypeError('only Queries supported')

        self.__one_shot = one_shot
        self.__polymerise_limit = polymerise_limit
        self.__products_atoms = tuple(set(m) for m in products)
        self.__automorphism_filter = automorphism_filter
        self.__meta = template.meta.copy()
        super().__init__(reduce(or_, reactants), products_, delete_atoms)

    def __call__(self, structures: Iterable[MoleculeContainer]):
        if any(not isinstance(structure, MoleculeContainer) for structure in structures):
            raise TypeError('only list of Molecules possible')

        structures = self.__remap(structures)
        if self.__one_shot:
            s_nums = set(range(len(structures)))
            for chosen in permutations(s_nums, len(self.__patterns)):
                ignored = [structures[x] for x in s_nums.difference(chosen)]
                chosen = [structures[x] for x in chosen]
                yield from (ReactionContainer(structures, new + ignored, meta=self.__meta)
                            for new in self.__single_stage(chosen, ignored))
        else:
            ...

    def __single_stage(self, chosen, ignored):
        ignored_numbers = {x for x in ignored for x in x}
        max_ignored_number = max(ignored_numbers, default=0)
        united_chosen = reduce(or_, chosen)
        for match in lazy_product(*(x.get_mapping(y, automorphism_filter=self.__automorphism_filter) for x, y in
                                    zip(self.__patterns, chosen))):
            mapping = match[0]
            for m in match[1:]:
                mapping.update(m)

            new = self._patcher(united_chosen, mapping)
            collision = set(new).intersection(ignored_numbers)
            if collision:
                new.remap(dict(zip(collision, count(max(max_ignored_number, max(new.atoms_numbers)) + 1))))
            if len(self.__products_atoms) > 1:
                components = new._connected_components(new._bonds)
                # stack unmatched components of reactants
                # works only for partially matched molecules
                matched = set(mapping.values())
                components = self.__stack_components(components, [set(x) - matched for x in chosen])
                # multi-component molecules in product side of template will be kept.
                components = self.__stack_components(components,
                                                     [{mapping[x] for x in x if x in mapping}
                                                      for x in self.__products_atoms])
                yield [new.substructure(c) for c in components]
            else:
                yield [new]

    @staticmethod
    def __stack_components(components, groups):
        out = []
        for g in groups:
            if not g:
                continue
            common = []
            tmp = []
            for c in components:
                if not g.isdisjoint(c):
                    common.append(c)
                else:
                    tmp.append(c)
            if common:
                components = tmp
                common = reduce(or_, common)
                out.append(common)
        out.extend(components)
        return out

    @staticmethod
    def __remap(structures):
        checked = []
        checked_atoms = set()
        for structure in structures:
            intersection = set(structure).intersection(checked_atoms)
            if intersection:
                mapping = dict(zip(intersection, count(max(max(checked_atoms), max(structure.atoms_numbers)) + 1)))
                structure = structure.remap(mapping, copy=True)
                info('some atoms in input structures had the same numbers.\n'
                     f'atoms {list(mapping)} were remapped to {list(mapping.values())}')
            checked_atoms.update(structure)
            checked.append(structure)
        return checked

    def __getstate__(self):
        return {'patterns': self.__patterns, 'meta': self.__meta, 'products_atoms': self.__products_atoms,
                'polymerise_limit': self.__polymerise_limit, 'one_shot': self.__one_shot,
                'automorphism_filter': self.__automorphism_filter, **super().__getstate__()}

    def __setstate__(self, state):
        if 'split' in state:
            raise ValueError('Reactor pickled with incompatible version.')
        self.__patterns = state['patterns']
        self.__meta = state['meta']
        self.__one_shot = state['one_shot']
        self.__polymerise_limit = state['polymerise_limit']
        self.__products_atoms = state['products_atoms']
        self.__automorphism_filter = state['automorphism_filter']
        super().__setstate__(state)


__all__ = ['Reactor']
