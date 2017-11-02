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
# elements, charge, radical, bonds, implicitH
_valence = (
    # elemental Me
    (('Li', 'Na', 'K', 'Rb', 'Cs', 'Fr'), 0, 0, 0, 0),
    # 1 bond is accepted
    (('H', 'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr'), 0, 0, 1, 0),
    # Li+, Na+, K+, Rb+, Cs+, Fr+ are accepted.
    (('H', 'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr'), 1, 0, 0, 0),
    # Me· is accepted
    (('H', 'Li', 'Na', 'K', 'Rb', 'Cs', 'Fr'), 0, 1, 0, 0),

    # elemental
    (('Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra'), 0, 0, 0, 0),
    # 2 bonds are accepted.
    (('Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra'), 0, 0, 2, 0),
    #  Each added charge decreases by one the number of possible bonds.
    (('Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra'), 1, 0, 1, 0),
    (('Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra'), 2, 0, 0, 0),
    # Each added radical decreases by one the number of possible bonds
    (('Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra'), 0, 1, 1, 0),
    (('Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra'), 0, 2, 0, 0),
    # ·Me+ is accepted
    (('Be', 'Mg', 'Ca', 'Sr', 'Ba', 'Ra'), 1, 1, 0, 0),

    # elemental
    (('Al', 'Ga', 'In', 'Tl'), 0, 0, 0, 0),
    # Implicit Hydrogens are added to have 3 bonds.
    ('B', 0, 0, 3, 3),
    (('Al', 'Ga', 'In', 'Tl'), 0, 0, 3, 2),
    # Every added charge decreases by one the number of possible bonds. Maximum charge accepted is +3.
    ('B', 1, 0, 2, 2), (('Al', 'Ga', 'In', 'Tl'), 1, 0, 2, 1),
    ('B', 2, 0, 1, 1), (('Al', 'Ga', 'In', 'Tl'), 2, 0, 1, 0),
    (('B', 'Al', 'Ga', 'In', 'Tl'), 3, 0, 0, 0),
    # [B|Al|...]H4- are accepted.
    (('B', 'Al', 'Ga', 'In', 'Tl'), -1, 0, 4, 4),
    #  BH32- are accepted.
    ('B', -2, 0, 3, 3),
    # Al+, Ga+, In+, Tl+  are accepted.
    (('Al', 'Ga', 'In', 'Tl'), 1, 0, 0, 0),
    # Every added radical on these atoms decreases by one the number of necessary bonds.
    (('B', 'Al', 'Ga', 'In', 'Tl'), 0, 1, 2, 2),
    (('B', 'Al', 'Ga', 'In', 'Tl'), 0, 2, 1, 1),
    #  ·Me2+ and :̣Me+ are allowed.
    (('Al', 'Ga', 'In', 'Th'), 2, 1, 0, 0),
    (('Al', 'Ga', 'In', 'Th'), 1, 2, 0, 0),
    # RMe+·
    (('B', 'Al', 'Ga', 'In', 'Th'), 1, 1, 1, 1),

    # elemental
    (('C', 'Si', 'Ge', 'Sn', 'Pb'), 0, 0, 0, 0),
    # CH4, SiH4, GeH4, SnHxR(4-x) respectively.
    (('C', 'Si', 'Ge'), 0, 0, 4, 4), ('Sn', 0, 0, 4, 3),
    # Sn Pb 5 or 6 bonds accepted.
    (('Sn', 'Pb'), 0, 0, 5, 0),
    (('Sn', 'Pb'), 0, 0, 6, 0),
    # 1 or 2 bonds: 1 or 0 implicit hydrogen is added, respectively.
    ('Pb', 0, 0, 2, 1),
    # 3 or 4 bonds: 1 or 0 implicit hydrogen is added, respectively.
    ('Pb', 0, 0, 4, 1),
    # Every added charge on these atoms (no matter it is positive or negative. except Sn)
    # decreases by one the number of possible bonds. Maximum charge accepted is -2 +4.
    (('C', 'Si', 'Ge'), -1, 0, 3, 3),
    (('C', 'Si', 'Ge'), -2, 0, 2, 2),
    (('C', 'Si', 'Ge', 'Sn'), 1, 0, 3, 3),
    # Pb+[RH]
    ('Pb', 1, 0, 1, 1),
    # Pb+R2[RH]
    ('Pb', 1, 0, 3, 1),
    ('Pb', 1, 0, 4, 0),
    # Sn+R5 Pb+R5
    (('Sn', 'Pb'), 1, 0, 5, 0),
    (('C', 'Si', 'Ge'), 2, 0, 2, 2),
    # Sn2+ , Pb2+ are accepted.
    (('Sn', 'Pb'), 2, 0, 0, 0),
    # Sn2+RX , Pb2+RX are accepted X = R or H.
    (('Sn', 'Pb'), 2, 0, 2, 1),
    (('Ge', 'Sn', 'Pb'), 3, 0, 1, 1),
    # Pb3+R2
    ('Pb', 3, 0, 2, 0),
    # Pb3+R3
    ('Pb', 3, 0, 3, 0),
    (('Ge', 'Sn', 'Pb'), 4, 0, 0, 0),
    # Pb4+R
    ('Pb', 4, 0, 1, 0),
    # Pb4+R2
    ('Pb', 4, 0, 2, 0),
    # Every added radical on these atoms decreases by one the number accepted of bonds.
    (('C', 'Si', 'Ge', 'Sn'), 0, 1, 3, 3),
    (('C', 'Si', 'Ge', 'Sn'), 0, 2, 2, 2),
    ('Sn', 0, 1, 5, 0),
    ('Pb', 0, 1, 1, 1),
    ('Pb', 0, 1, 3, 1),
    ('Pb', 0, 1, 4, 0),
    ('Pb', 0, 1, 5, 0),
    ('Pb', 0, 2, 0, 0),
    ('Pb', 0, 2, 2, 1),
    ('Pb', 0, 2, 3, 0),
    ('Pb', 0, 2, 4, 0),
    # radical-charge combo
    (('C', 'Si', 'Ge'), -1, 1, 2, 2),
    (('C', 'Si', 'Ge', 'Sn'), 1, 1, 2, 2),
    ('Pb', 1, 1, 0, 0),
    ('Pb', 1, 1, 2, 1),
    ('Pb', 1, 1, 3, 0),
    ('Pb', 1, 1, 4, 0),
    (('Ge', 'Sn', 'Pb'), 1, 2, 1, 1),
    ('Pb', 1, 2, 3, 0),
    ('Pb', 1, 2, 2, 0),
    (('Ge', 'Sn', 'Pb'), 2, 1, 1, 1),
    ('Pb', 2, 1, 3, 0),
    ('Pb', 2, 1, 2, 0),
    (('Ge', 'Sn', 'Pb'), 2, 2, 0, 0),
    ('Pb', 2, 2, 1, 0),
    ('Pb', 2, 2, 2, 0),
    (('Ge', 'Sn', 'Pb'), 3, 1, 0, 0),
    ('Pb', 3, 2, 1, 0),
    ('Pb', 3, 1, 1, 0),
    ('Pb', 3, 1, 2, 0),

    # elemental
    ('Bi', 0, 0, 0, 0),
    # hydrogen
    ('Bi', 0, 0, 3, 2),
    (('N', 'P', 'As', 'Sb'), 0, 0, 3, 3),
    (('P', 'As', 'Sb', 'Bi'), 0, 0, 5, 1),
    (('N', 'P', 'As', 'Sb', 'Bi'), 1, 0, 4, 4),
    (('N', 'P', 'As', 'Sb', 'Bi'), 2, 0, 3, 3),
    (('N', 'P', 'As', 'Sb', 'Bi'), -1, 0, 2, 2),
    (('N', 'P', 'As', 'Sb', 'Bi'), -2, 0, 1, 1),
    (('N', 'P', 'As', 'Sb', 'Bi'), -3, 0, 0, 0),
    # radical
    (('N', 'P', 'As', 'Sb', 'Bi'), 0, 1, 2, 2),
    (('N', 'P', 'As', 'Sb', 'Bi'), 0, 2, 1, 1),
    # radical-charge combo
    (('N', 'P', 'As', 'Sb', 'Bi'), 1, 1, 3, 3),
    (('N', 'P', 'As', 'Sb', 'Bi'), 1, 2, 2, 2),
    (('N', 'P', 'As', 'Sb', 'Bi'), 2, 1, 2, 2),
    (('N', 'P', 'As', 'Sb', 'Bi'), 2, 2, 1, 1),
    (('N', 'P', 'As', 'Sb', 'Bi'), -1, 1, 1, 1),
    (('N', 'P', 'As', 'Sb', 'Bi'), -1, 2, 0, 0),
    (('N', 'P', 'As', 'Sb', 'Bi'), -2, 1, 0, 0),

)


_aromatic = ('B', 'C', 'N', 'P')
