# -*- coding: utf-8 -*-
#
#  Copyright 2020 Nail Samikaev <samikaevn@yandex.ru>
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


class Tautomers:
    __slots__ = ()

    def enumerate_tautomers(self):
        atoms = self._atoms
        bonds = self._bonds
        charges = self._charges
        radicals = self._radicals
        hydrogens = self._hydrogens
        yield self


__all__ = ['Tautomers']
