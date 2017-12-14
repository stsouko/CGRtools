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


class ValenceError(Exception):
    __message = 'need neighbors atoms for validation in order of bonds connections'

    def __init__(self, message=None):
        super().__init__(message or self.__message)


class Valence:
    def __init__(self):
        self.__valence_rules = self.__prepare_valence_rules()
        self.__implicit_rules = self.__prepare_implicit_h_rules()

    def _check_charge_radical(self, element, charge, radical=None):
        if radical is None:
            radical = 0

        return (element, charge, radical) in self.__valence_rules

    def _check_valence(self, element, charge, bonds, radical=None, neighbors=None):
        if radical is None:
            radical = 0

        sum_bonds = self.__bonds_sum(bonds)
        len_bonds = len(bonds)

        if element in ('S', 'Se', 'Te', 'Po') and charge == -1 and len_bonds == 3:
            # Sulfur (S), Selenium (Se), Tellurium (Te), Polonium (Po)
            # Exception: 2 double bonded oxygens and 1 more ligand connected by a single bond)
            if neighbors is None:
                raise ValenceError()
            return sum_bonds == 5 and all(y == 'O' for x, y in zip(bonds, neighbors) if x == 2)

        if element == 'N' and charge == 0 and sum_bonds in (5, 4):
            # Pentavalent N is accepted in traditional form by default,
            # but keep in mind that these forms are not correct.
            if neighbors is None:
                raise ValenceError()
            return any(y == 'O' for x, y in zip(bonds, neighbors) if x == 2)

        if element == 'Xe' and charge == 0 and len_bonds > 0:
            # XeF2, XeF4, XeF6, XeO3, XeO4, XeOF4, XeO2F2, XeO3F2
            if neighbors is None:
                raise ValenceError()
            if all(y == 'O' and x == 2 or y == 'F' and x == 1 for x, y in zip(bonds, neighbors)):
                ofc = (neighbors.count('O'), neighbors.count('F'))
                return ofc in ((0, 2), (0, 4), (0, 6), (3, 0), (4, 0), (1, 4), (2, 2), (3, 2))
            return False

        if element == 'I' and charge == 0 and len_bonds == 7:
            # 7 valence: IF7 is accepted.
            if neighbors is None:
                raise ValenceError()
            return all(x == 'F' for x in neighbors)

        if element in ('Cl', 'Br', 'I', 'At') and charge == 0 and sum_bonds in (3, 5, 7):
            if neighbors is None:
                raise ValenceError()
            if all(x == 'O' for x in neighbors):
                # HXOy, where X = Cl, Br, I or At; and y = 2, 3 or 4.
                return 2 * len_bonds - 1 == sum_bonds
            if {'C', 'Si', 'N', 'P', 'S', 'Se', 'O', 'F', 'Cl'}.issuperset(neighbors):
                # 3 valence: Two of the ligands should be F, Cl or oxygen;
                # and one more atom from the 14th, 15th, 16th columns.
                # Double bonded oxygen and one any other atom from 14th, 15th, 16th columns are also accepted.
                # 5 valence: Four of the ligands should be F, Cl or oxygen;
                # and one more atom from the 14th, 15th, 16th columns.
                # Double bonded oxygens are accepted.
                return sum(x for x, y in zip(bonds, neighbors) if y in ('C', 'Si', 'N', 'P', 'S', 'Se')) == 1 and \
                       all(y in ('F', 'Cl') and x == 1 or y == 'O' and x in (1, 2)
                           for x, y in zip(bonds, neighbors) if y in ('F', 'Cl', 'O'))
            return False

        return sum_bonds in self.__valence_rules.get((element, charge, radical), [])

    def _get_implicit_h(self, element, charge, bonds, radical=None):
        if radical is None:
            radical = 0

        sum_bonds = self.__bonds_sum(bonds)

        if element == 'N' and charge == 0 and sum_bonds == 4:
            return 1

        return self.__implicit_rules.get((element, charge, radical, sum_bonds), 0)

    @classmethod
    def __bonds_sum(cls, bonds):
        return int(sum(cls.__bonds[x] for x in bonds))

    @classmethod
    def __prepare_valence_rules(cls):
        out = {}
        for a, c, r, v, h in cls.__valence:
            x = range(v - h, v + 1) if h else [v]
            for i in (a if isinstance(a, tuple) else [a]):
                out.setdefault((i, c, r), []).extend(x)
        return out

    @classmethod
    def __prepare_implicit_h_rules(cls):
        out = {}
        for a, c, r, v, h in cls.__valence:
            if h:
                for k in range(1, h + 1):
                    for i in (a if isinstance(a, tuple) else [a]):
                        out[(i, c, r, v - k)] = k
        return out

    # http://onlinelibrarystatic.wiley.com/marvin/help/sci/ValenceCalculator.html
    # elements, charge, radical, bonds, implicitH
    __valence = (
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

    __aromatic = ('B', 'C', 'N', 'P', 'O', 'S')

    __inorganic_molecules = {'transition': ('Mn,=O,=O,=O,-O', 'Mn,=O,=O')}
    __bonds = {1: 1, 2: 2, 3: 4, 4: 1.5, 9: 1}
