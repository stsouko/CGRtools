# -*- coding: utf-8 -*-
#
#  Copyright 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2020 Dinar Batyrshin <batyrshin-dinar@mail.ru>
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


def functional_groups(molecule, deep):
    """
    :param molecule: MoleculeContainer
    :param deep: chain length

    return: list of molecule functional groups
    """
    bonds = molecule._bonds

    if deep < 1:
        raise ValueError('Deep should be >= 1')

    response = []
    groups = set()
    stack = [([a], list(n)) for a, n in bonds.items()]
    while stack:
        aug, nei = stack.pop(0)
        for x in nei:
            augx = (*aug, x)
            if augx not in groups:
                groups.add(augx)
                response.append(molecule.substructure(augx, as_query=True))
                nt = nei.copy()
                nt.remove(x)
                nt.extend(list(bonds[x]))
                if len(augx) < deep:
                    stack.append((augx, nt))

    return response


__all__ = ['functional_groups']
