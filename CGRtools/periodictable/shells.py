# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
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
from .tables import classes, elements


_noble_shells = ((0,) + classes['noble'])[::-1]
_wiswesser = sorted([(n, l) for n in range(1, 8) for l in range(n)], key=lambda x: x[0] + x[1] - (x[1] / (x[1] + 1)))
_electrons = [{(1, 0): 1}]
for x in [(n, l) for n, l in _wiswesser for _ in range(4 * l + 2)][1:]:
    new = _electrons[-1].copy()
    new[x] = new.get(x, 0) + 1
    _electrons.append(new)

electrons_configuration = _ = dict(zip(elements, _electrons))
_['Cr'].update({(3, 2): 5, (4, 0): 1})
_['Cu'].update({(3, 2): 10, (4, 0): 1})
_['Nb'].update({(4, 2): 4, (5, 0): 1})
_['Mo'].update({(4, 2): 5, (5, 0): 1})
_['Ru'].update({(4, 2): 7, (5, 0): 1})
_['Rh'].update({(4, 2): 8, (5, 0): 1})
_['Pd'][(4, 2)] = 10; del _['Pd'][(5, 0)]
_['Ag'].update({(4, 2): 10, (5, 0): 1})
_['Pt'].update({(5, 2): 9, (6, 0): 1})
_['Au'].update({(5, 2): 10, (6, 0): 1})
_['La'][(5, 2)] = 1; del _['La'][(4, 3)]
_['Ce'].update({(4, 3): 1, (5, 2): 1})
_['Gd'].update({(4, 3): 7, (5, 2): 1})
_['Ac'][(6, 2)] = 1; del _['Ac'][(5, 3)]
_['Th'][(6, 2)] = 2; del _['Th'][(5, 3)]
_['Pa'].update({(5, 3): 2, (6, 2): 1})
_['U'].update({(5, 3): 3, (6, 2): 1})
_['Np'].update({(5, 3): 4, (6, 2): 1})
_['Cm'].update({(5, 3): 7, (6, 2): 1})

valence_electrons = {s: next(n - x for x in _noble_shells if n > x) for n, s in enumerate(elements, start=1)}
orbitals_names = ['s', 'p', 'd', 'f', 'g']


__all__ = ['valence_electrons', 'electrons_configuration', 'orbitals_names']
