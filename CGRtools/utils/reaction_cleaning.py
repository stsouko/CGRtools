# -*- coding: utf-8 -*-
#
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2021 Timur Gimadiev <timur.gimadiev@gmail.com>
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
from ..containers import ReactionContainer
from ..exceptions import MappingError


def remove_reagents(reaction, keep_reagents=False):
    """
    Preprocess reaction according to mapping, using the following idea: molecules(each separated graph) will be
    placed to reagents if it is not changed in the reaction (no bonds, charges reorders)
    """
    if not isinstance(reaction, ReactionContainer):
        raise TypeError('Reaction container only supported')
    try:  # check if CGR can be build, else raise error
        cgr = ~reaction
    except ValueError:
        raise ValueError('Problem with CGR construction')
    if cgr.center_atoms:
        active = set(cgr.center_atoms)
        reactants = []
        products = []
        reagents = set()
        for i in reaction.reactants:
            if not active.isdisjoint(i):
                reactants.append(i)
            else:
                reagents.add(i)
        for i in reaction.products:
            if not active.isdisjoint(i):
                products.append(i)
            else:
                reagents.add(i)
        if keep_reagents:
            reaction = ReactionContainer(reactants=reactants, reagents=reaction.reagents+tuple(reagents),
                                         products=products, meta=reaction.meta)
            return reaction
        else:
            reaction = ReactionContainer(reactants=reactants,
                                         products=products, meta=reaction.meta)
            return reaction
    raise MappingError("Reaction center is absent according to mapping")


__all__ = ['remove_reagents']
