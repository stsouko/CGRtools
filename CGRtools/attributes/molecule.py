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


class Atom:
    __slots__ = ('atom', 'charge', 'is_radical', 'parsed_mapping', 'xy')

    def __setstate__(self, state):
        cm = state['atom']._backward
        self.xy = (state['x'], state['y'])
        self.parsed_mapping = state['mapping']
        self.charge = cm[0]
        self.is_radical = bool(cm[1])
        self.atom = state['atom']


class Bond:
    __slots__ = ('order',)

    def __setstate__(self, state):
        self.order = state['order']


__all__ = ['Atom', 'Bond']
