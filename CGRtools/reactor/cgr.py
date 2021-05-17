# -*- coding: utf-8 -*-
#
#  Copyright 2014-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from operator import or_
from typing import Union
from .base import BaseReactor
from ..containers import QueryContainer, QueryCGRContainer, MoleculeContainer, CGRContainer, ReactionContainer


class CGRReactor(BaseReactor):
    """
    Editor for CGRs and molecules.
    generates modified CGR's/molecules from input CGR/molecule using template.
    Template should contain one reactant and one product.
    CGRReactor calling returns generator of all possible replacements.
    """
    def __init__(self, template: ReactionContainer, delete_atoms: bool = True):
        """
        :param template: CGRtools ReactionContainer
        :param delete_atoms: if True atoms exists in reactant but
                            not exists in product will be removed
        """
        reactants, products = template.reactants, template.products
        if not reactants or not products:
            raise ValueError('empty template')
        if isinstance(reactants[0], QueryCGRContainer):
            reactants = reduce(or_, reactants)
            products = reduce(or_, products)
        elif isinstance(reactants[0], QueryContainer):
            reactants = reduce(or_, reactants)
            products = reduce(or_, products)
        else:
            raise TypeError('only QueryCGRContainer or QueryContainer supported')

        self.__pattern = reactants
        self.__meta = template.meta.copy()
        super().__init__(reactants, products, delete_atoms)

    def __call__(self, structure: Union[MoleculeContainer, CGRContainer], automorphism_filter: bool = True):
        if not isinstance(structure, (MoleculeContainer, CGRContainer)):
            raise TypeError('only Molecules and CGRs possible')

        for mapping in self.__pattern.get_mapping(structure, automorphism_filter=automorphism_filter):
            new = self._patcher(structure, mapping)
            new.meta.update(self.__meta)
            yield new

    def __getstate__(self):
        return {'pattern': self.__pattern, 'meta': self.__meta, **super().__getstate__()}

    def __setstate__(self, state):
        self.__pattern = state['pattern']
        self.__meta = state['meta']
        super().__setstate__(state)


__all__ = ['CGRReactor']
