# -*- coding: utf-8 -*-
#
#  Copyright 2017-2019 Ramil Nugmanov <stsouko@live.ru>
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
from collections import defaultdict
from itertools import product

# http://onlinelibrarystatic.wiley.com/marvin/help/sci/ValenceCalculator.html
# elements, charge, radical, bonds, [implicitH]
_valence_rules = (
    ('H', -1, 0, 0, 0),
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
    (('Al', 'Ga', 'In', 'Tl'), 2, 1, 0, 0),
    (('Al', 'Ga', 'In', 'Tl'), 1, 2, 0, 0),
    # RMe+·
    (('B', 'Al', 'Ga', 'In', 'Tl'), 1, 1, 1, 1),

    # elemental
    (('Sn', 'Pb'), 0, 0, 0, 0),
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
    (('P', 'As', 'Sb', 'Bi'), -1, 0, 4, 0),
    # radical
    (('N', 'P', 'As', 'Sb', 'Bi'), 0, 1, 2, 2),
    (('N', 'P', 'As', 'Sb', 'Bi'), 0, 2, 1, 1),
    (('P', 'As', 'Sb', 'Bi'), 0, 1, 4, 1),
    (('P', 'As', 'Sb', 'Bi'), 0, 2, 3, 1),
    # radical-charge combo
    (('N', 'P', 'As', 'Sb', 'Bi'), 1, 1, 3, 3),
    (('N', 'P', 'As', 'Sb', 'Bi'), 1, 2, 2, 2),
    (('N', 'P', 'As', 'Sb', 'Bi'), 2, 1, 2, 2),
    (('N', 'P', 'As', 'Sb', 'Bi'), 2, 2, 1, 1),
    (('N', 'P', 'As', 'Sb', 'Bi'), -1, 1, 1, 1),
    (('N', 'P', 'As', 'Sb', 'Bi'), -1, 2, 0, 0),
    (('N', 'P', 'As', 'Sb', 'Bi'), -2, 1, 0, 0),

    # implicit H
    (('O', 'S', 'Se', 'Te', 'Po'), 0, 0, 2, 2),
    (('S', 'Se', 'Te', 'Po'), 0, 0, 4, 0),
    (('S', 'Se', 'Te', 'Po'), 0, 0, 6, 1),
    # charges
    (('O', 'S', 'Se', 'Te', 'Po'), 1, 0, 3, 3),
    (('S', 'Se', 'Te', 'Po'), 1, 0, 5, 1),
    (('Se', 'Te', 'Po'), 1, 0, 1, 1),
    (('O', 'S', 'Se', 'Te', 'Po'), -1, 0, 1, 1),
    (('O', 'S', 'Se', 'Te', 'Po'), -2, 0, 0, 0),
    # radicals
    (('O', 'S', 'Se', 'Te', 'Po'), 0, 1, 1, 1),
    (('O', 'S', 'Se', 'Te', 'Po'), 0, 2, 0, 0),
    (('O', 'S', 'Se', 'Te', 'Po'), 1, 1, 0, 0),
    (('O', 'S', 'Se', 'Te', 'Po'), -1, 1, 0, 0),

    # implicit H
    (('F', 'Cl', 'Br', 'I', 'At'), 0, 0, 1, 1),
    (('F', 'Cl', 'Br', 'I', 'At'), 1, 0, 2, 2),
    (('F', 'Cl', 'Br', 'I', 'At'), -1, 0, 0, 0),
    (('F', 'Cl', 'Br', 'I', 'At'), 0, 1, 0, 0),

    (('He', 'Ne', 'Ar', 'Kr', 'Xe', 'Rn'), 0, 0, 0, 0)
)

_extra_valence_rules = (
    (('Sc', 'Y',
      'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
      'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'), 0, 0, 3),
    (('Sc', 'Y',
      'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
      'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr'), 3, 0, 0),
    (('Ti', 'Zr', 'Hf', 'Rf'), 0, 0, 4)
)

_valence_rules = _valence_rules + _extra_valence_rules

atom_valences = {}
atom_charge_radical = set()
for a, c, r, *vh in _valence_rules:
    h = 0 if len(vh) == 1 else vh[1]
    v = vh[0]
    x = range(v - h, v + 1) if h else [v]
    for i in (a if isinstance(a, tuple) else [a]):
        for j in x:
            atom_valences[(i, c, r, j)] = v
            atom_charge_radical.add((i, c, r))


atom_implicit_h = {('N', 0, 0, 4): 1}
for a, c, r, *vh in _valence_rules:
    h = 0 if len(vh) == 1 else vh[1]
    v = vh[0]
    if h:
        for k in range(1, h + 1):
            for i in (a if isinstance(a, tuple) else [a]):
                atom_implicit_h[(i, c, r, v - k)] = k


# elements charge radical Bond Atom BondAtom
_exceptions = (
    (('S', 'Se', 'Te', 'Po'), -1, 0, ('=', 'O', '=', 'O', '-', ('C', 'O', 'S'))),
    # Pentavalent N is accepted in traditional form by default, but keep in mind that these forms are not correct.
    ('N', 0, 0, ('=', 'O', '=', ('O', 'C'), '-', ('O', 'C', 'N'))),
    ('N', 0, 0, ('=', 'O', '-', 'C', '-', ('O', 'C', 'N')), 5),
    ('N', 0, 0, ('=', 'O', '-', 'C', '-', 'C', '-', ('O', 'C', 'N'))),
    # XeF2, XeF4, XeF6, XeO3, XeO4, XeO2F2, XeOF4,  XeO3F2
    ('Xe', 0, 0, ('-', 'F', '-', 'F')),
    ('Xe', 0, 0, ('-', 'F', '-', 'F', '-', 'F', '-', 'F')),
    ('Xe', 0, 0, ('-', 'F', '-', 'F', '-', 'F', '-', 'F', '-', 'F', '-', 'F')),
    ('Xe', 0, 0, ('=', 'O', '=', 'O', '=', 'O')),
    ('Xe', 0, 0, ('=', 'O', '=', 'O', '=', 'O', '=', 'O')),
    ('Xe', 0, 0, ('=', 'O', '=', 'O', '-', 'F', '-', 'F')),
    ('Xe', 0, 0, ('=', 'O', '-', 'F', '-', 'F', '-', 'F', '-', 'F')),
    ('Xe', 0, 0, ('=', 'O', '=', 'O', '=', 'O', '-', 'F', '-', 'F')),
    # 7 valence: IF7 is accepted
    ('I', 0, 0, ('-', 'F', '-', 'F', '-', 'F', '-', 'F', '-', 'F', '-', 'F', '-', 'F')),
    # HXOy, where X = Cl, Br, I or At; and y = 2, 3 or 4.
    (('Cl', 'Br', 'I', 'At'), 0, 0, ('-', 'O', '=', 'O')),
    (('Cl', 'Br', 'I', 'At'), 0, 0, ('-', 'O', '=', 'O', '=', 'O')),
    (('Cl', 'Br', 'I', 'At'), 0, 0, ('-', 'O', '=', 'O', '=', 'O', '=', 'O')),
    # Two of the ligands should be F, Cl or oxygen; and one more atom from the 14th, 15th, 16th columns.
    (('Cl', 'Br', 'I', 'At'), 0, 0, ('-', ('C', 'Si', 'N', 'P', 'S', 'Se', 'O'), '-', ('O', 'F', 'Cl'),
                                     '-', ('O', 'F', 'Cl'))),
    # Double bonded oxygen and one any other atom from 14th, 15th, 16th columns are also accepted.
    (('Cl', 'Br', 'I', 'At'), 0, 0, ('-', ('C', 'Si', 'N', 'P', 'S', 'Se', 'O'), '=', 'O')),
    # Four of the ligands should be F, Cl or oxygen; and one more atom from the 14th, 15th, 16th columns.
    (('Cl', 'Br', 'I', 'At'), 0, 0, ('-', ('C', 'Si', 'N', 'P', 'S', 'Se', 'O'), '-', ('O', 'F', 'Cl'),
                                     '-', ('O', 'F', 'Cl'), '-', ('O', 'F', 'Cl'), '-', ('O', 'F', 'Cl'))),
    # Double bonded oxygens are accepted.
    (('Cl', 'Br', 'I', 'At'), 0, 0, ('-', ('C', 'Si', 'N', 'P', 'S', 'Se', 'O'), '=', 'O',
                                     '-', ('O', 'F', 'Cl'), '-', ('O', 'F', 'Cl'))),
    (('Cl', 'Br', 'I', 'At'), 0, 0, ('-', ('C', 'Si', 'N', 'P', 'S', 'Se', 'O'), '=', 'O', '=', 'O')),
)

bonds_map = {1: 1, 2: 2, 3: 3, 4: 1.5, 9: 1, '-': 1, '=': 2, '#': 3, ':': 1.5, '.': 0, None: 0}

atom_valences_exceptions = defaultdict(dict)
for a, c, r, b, *bs in _exceptions:
    ap = list(product(*(x if isinstance(x, tuple) else [x] for x in b[1::2])))
    bl = [bonds_map[x] for x in b[::2]]
    bs = bs[0] if bs else int(sum(bl))
    for i in (a if isinstance(a, tuple) else [a]):
        for ape in ap:
            atom_valences_exceptions[(i, c, r, len(bl))][tuple(sorted(zip(ape, bl)))] = bs

atom_valences_exceptions = dict(atom_valences_exceptions)


__all__ = ['atom_valences', 'atom_valences_exceptions', 'atom_implicit_h', 'atom_charge_radical', 'bonds_map']
