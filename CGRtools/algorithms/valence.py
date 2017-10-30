# -*- coding: utf-8 -*-
#
#  Copyright 2017 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#


def is_valid(g):
    g.reset_query_marks()


# http://onlinelibrarystatic.wiley.com/marvin/help/sci/ValenceCalculator.html
# elements, charge, radical, bonds, implicitH, aromatic
_valence = [
    # 1 bond is accepted
    (('H', 'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr'), 0, 0, 1, False, False),
    # H+, Li+, Na+, K+, Rb+, Cs+, Fr+ are accepted.
    (('H', 'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr'), 1, 0, 0, False, False),
    # H·, Me· is accepted
    (('H', 'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr'), 0, 1, 0, False, False),
    # 2 bonds are accepted.
    (('Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra'), 0, 0, 2, False, False),
    #  Each added charge decreases by one the number of possible bonds.
    (('Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra'), 1, 0, 1, False, False),
    (('Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra'), 2, 0, 0, False, False),
    # Each added radical decreases by one the number of possible bonds
    (('Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra'), 0, 1, 1, False, False),
    (('Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra'), 0, 2, 0, False, False),
    # ·Me+ is accepted
    (('Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra'), 1, 1, 0, False, False),
    # Implicit Hydrogens are added to have 3 bonds.
    ('B', 0, 0, 3, True, True),
    # Every added radical on these atoms decreases by one the number of necessary bonds.
    ('B', 0, 3, 0, True, True),
    ('B', 0, 2, 1, True, True),
    ('B', 0, 1, 2, True, True),
    # Every added charge decreases by one the number of possible bonds. Maximum charge accepted is +3.
    ('B', 3, 0, 0, True, True),
    ('B', 2, 0, 1, True, True),
    ('B', 1, 0, 2, True, True),
    # BH4-, BH32- BH23- BH4- and B5- are accepted.
    ('B', -1, 0, 4, True, True),
    ('B', -2, 0, 3, True, True),
    ('B', -3, 0, 2, True, True),
    ('B', -4, 0, 1, True, True),
    ('B', -5, 0, 0, True, True),
    # 3 bonds around these atoms are accepted.
    (('Al', 'Ga', 'In', 'Tl'), 0, 0, 3, False, False),
    # Negative charge: 1- charge with 4 ligands is accepted.
    (('Al', 'Ga', 'In', 'Tl'), -1, 0, 4, False, False),
    #  Every added charge decreases by one the number of possible bonds. Maximum charge accepted is +3.
    (('Al', 'Ga', 'In', 'Tl'), 1, 0, 2, False, False),
    (('Al', 'Ga', 'In', 'Tl'), 2, 0, 1, False, False),
    (('Al', 'Ga', 'In', 'Tl'), 3, 0, 0, False, False),
    # Al+, Ga+, In+, Tl+  are accepted.
    (('Al', 'Ga', 'In', 'Tl'), 1, 0, 0, False, False),
    # Every added radical on these atoms decreases by one the number of necessary bonds.
    (('Al', 'Ga', 'In', 'Th'), 0, 1, 2, False, False),
    (('Al', 'Ga', 'In', 'Th'), 0, 2, 1, False, False),
    (('Al', 'Ga', 'In', 'Th'), 0, 3, 0, False, False),
    #  ·Me2+ and :̣Me+ are allowed.
    (('Al', 'Ga', 'In', 'Th'), 2, 1, 0, False, False),
    (('Al', 'Ga', 'In', 'Th'), 1, 2, 0, False, False),
    # CH4, SiH4, GeH4 respectively.
    (('C', 'Si', 'Ge'), 0, 0, 4, True, True),
    # Every added charge on these atoms (no matter it is positive or negative)
    # decreases by one the number of possible bonds. Maximum charge accepted is ±4.
    (('C', 'Si', 'Ge'), -4, 0, 0, True, True),
    (('C', 'Si', 'Ge'), -3, 0, 1, True, True),
    (('C', 'Si', 'Ge'), -2, 0, 2, True, True),
    (('C', 'Si', 'Ge'), -1, 0, 3, True, True),
    (('C', 'Si', 'Ge'), 3, 0, 1, True, True),
    (('C', 'Si', 'Ge'), 2, 0, 2, True, True),
    (('C', 'Si', 'Ge'), 1, 0, 3, True, True),
    # Every added radical on these atoms decreases by one the number accepted of bonds.
    (('C', 'Si', 'Ge'), 0, 4, 0, True, True),
    (('C', 'Si', 'Ge'), 0, 3, 1, True, True),
    (('C', 'Si', 'Ge'), 0, 2, 2, True, True),
    (('C', 'Si', 'Ge'), 0, 1, 3, True, True),

]
