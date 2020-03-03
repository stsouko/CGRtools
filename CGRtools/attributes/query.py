# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ..periodictable import QueryElement


class QueryAtom:
    __slots__ = ('atom', 'charge', 'is_radical', 'neighbors', 'hybridization', 'xy')

    def __setstate__(self, state):
        atom = state['atom']
        if atom['element'] is None or len(atom['element']) > 1:
            raise TypeError('Any element in query not supported')
        self.xy = (atom['x'], atom['y'])
        self.charge = atom['charge']
        self.is_radical = bool(atom['multiplicity'])
        self.neighbors = tuple(sorted(atom['neighbors']))
        self.hybridization = tuple(sorted(atom['hybridization']))
        self.atom = QueryElement.from_symbol(atom['element'][0])(atom['isotope'])


__all__ = ['QueryAtom']
