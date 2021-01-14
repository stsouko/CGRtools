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


def remove_reagents(reaction, keep_reagents=False, cleaning=True, filter_rc=None):
    """
    Preprocess reaction according to mapping, using the following idea: molecules(each separated graph) will be
    placed to reagents if it is not changed in the reaction (no bonds, charges reorders)
    """
    if not isinstance(reaction, ReactionContainer):
        raise TypeError('Reaction container only supported')
    if cleaning:
        try:  # just pass if any error in standardization for now
            reaction.clean_isotopes()
            reaction.clean_stereo()
            reaction.canonicalize()
        except:  # temporary fix
            pass
    try:  # check if CGR can be build, else return None instead of reaction
        cgr = ~reaction
    except ValueError:
        return None
    if cgr.center_atoms:
        active = set()
        for n, i in enumerate(cgr.split()):
            active.update(i.center_atoms)
        reactants = []
        products = []
        reagents = set()
        for i in reaction.reactants:
            if set(i).intersection(active):
                reactants.append(i)
            else:
                reagents.add(i)
        for i in reaction.products:
            if set(i).intersection(active):
                products.append(i)
            else:
                reagents.add(i)
        if keep_reagents:
            reaction = ReactionContainer(reactants=reactants, reagents=reagents,
                                         products=products, meta=reaction.meta)
        else:
            reaction = ReactionContainer(reactants=reactants,
                                         products=products, meta=reaction.meta)
        if not filter_rc:
            return reaction
        elif filter_rc and len(reaction.centers_list) > filter_rc:
            return reaction
    return None

