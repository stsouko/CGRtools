#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  weighttable.py
#
#  Copyright 2014 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of cgrtools.
#
#  cgrtools is free software; you can redistribute it and/or modify
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


table = [('H', 1), ('He', 4),
         ('Li', 7), ('Be', 9), ('B', 11), ('C', 12), ('N', 14), ('O', 16), ('F', 19), ('Ne', 20),
         ('Na', 23), ('Mg', 24), ('Al', 26), ('Si', 28), ('P', 31), ('S', 32), ('Cl', 35), ('Ar', 40),
         ('K', 39), ('Ca', 40),
         ('Sc', 45), ('Ti', 48), ('V', 51), ('Cr', 52), ('Mn', 55),
         ('Fe', 56), ('Co', 59), ('Ni', 59), ('Cu', 63), ('Zn', 65),
         ('Ga', 70), ('Ge', 73), ('As', 75), ('Se', 79), ('Br', 80), ('Kr', 84),
         ('Rb', 85), ('Sr', 88),
         ('Y', 89), ('Zr', 91), ('Nb', 93), ('Mo', 96), ('Tc', 99),
         ('Ru', 101), ('Rh', 103), ('Pd', 106), ('Ag', 108), ('Cd', 112),
         ('In', 115), ('Sn', 119), ('Sb', 122), ('Te', 128), ('I', 127), ('Xe', 131),
         ('Cs', 133), ('Ba', 137),
         ('La', 139), ('Ce', 140), ('Pr', 141), ('Nd', 144), ('Pm', 145), ('Sm', 150), ('Eu', 152),
         ('Gd', 157), ('Tb', 159), ('Dy', 162), ('Ho', 165), ('Er', 167), ('Tm', 169), ('Yb', 173), ('Lu', 175),
         ('Hf', 178), ('Ta', 181), ('W', 184), ('Re', 186), ('Os', 190),
         ('Ir', 192), ('Pt', 195), ('Au', 197), ('Hg', 201),
         ('Tl', 204), ('Pb', 207), ('Bi', 209), ('Po', 210), ('At', 210), ('Rn', 222),
         ('Fr', 223), ('Ra', 226),
         ('Ac', 227), ('Th', 232), ('Pa', 231), ('U', 238), ('Np', 237), ('Pu', 244), ('Am', 243),
         ('Cm', 247), ('Bk', 247), ('Cf', 251), ('Es', 252), ('Fm', 257), ('Md', 258), ('No', 259), ('Lr', 262),
         ('Rf', 261), ('Db', 262), ('Sg', 263), ('Bh', 264), ('Hs', 277), ('Mt', 268), ('Ds', 271), ('Rg', 272)]

mendeley = {j[0]: i for i, j in enumerate(table)}
atommass = dict(table)
