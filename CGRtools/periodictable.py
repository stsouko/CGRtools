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
"""
contains perodictable, atom masses, electronegativity of elements
"""
from functools import reduce


_table_e = '''H                                                  He
              Li Be                               B  C  N  O  F  Ne
              Na Mg                               Al Si P  S  Cl Ar
              K  Ca Sc Ti V  Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr
              Rb Sr Y  Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I  Xe
              Cs Ba La

                    Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu

                       Hf Ta W  Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn
              Fr Ra Ac

                    Th Pa U  Np Pu Am Cm Bk Cf Es Fm Md No Lr
                     
                       Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og
           '''

# most stable isotope
_table_i = '''  1                                                                   4
                7   9                                          11  12  14  16  19  20
               23  24                                          27  28  31  32  35  40
               39  40  45  48  51  52  55  56  59  58  63  64  69  74  75  80  79  84
               85  88  89  90  93  98 115 102 103 106 107 114 115 120 121 130 127 132
              133 138 139

                      140 141 142 163 152 153 158 159 164 165 166 169 174 175

                          180 181 184 187 192 193 195 197 202 205 208 209 218 223 228
              232 234 236

                      232 231 238 244 247 249 252 254 256 257 259 260 262 263

                          264 265 266 267 277 271 281 272 285 286 289 289 293 294 294
           '''

# Luo Benson electronegativity
_table_p = '''2.7                                                                                  0
              0.75 2.08                                                   3.66 5.19 6.67 8.11 9.92 0
              0.65 1.54                                                   2.4  3.41 4.55 5.77 7.04 0
              0.51 1.15 1.49 1.57 1.65 1.72 1.71 1.72 1.83 1.92 2.3  1.87 2.38 3.24 4.2  5.13 6.13 0
              0.48 1.05 1.31 1.4  1.43 1.46 1.56 1.65 1.69 1.8  1.79 1.56 2.0  2.83 3.62 4.38 5.25 0
              0.43 1.01 1.17

                        None 1.2  None 1.23 None 1.23 None 1.28 None 1.31 None 1.33 1.34 1.36
                             1.41 1.44 1.45 1.46 1.46 1.46 1.49 1.5  1.51 1.91 2.6  3.29 4.03 4.67 0
              None None None
                        None None None None None None None None None None None None None None
                             None None None None None None None None None None None None None None 0
           '''

_groups = ((1, 3, 11, 19, 37, 55, 87),
           (4, 12, 20, 38, 56, 88),
           (21, 39) + tuple(range(57, 72)) + tuple(range(89, 104))) + \
          tuple((22 + x, 40 + x, 72 + x, 104 + x) for x in range(9)) + \
          tuple((5 + x, 13 + x, 31 + x, 49 + x, 81 + x, 113 + x) for x in range(5)) + ((2, 10, 18, 36, 54, 86, 118), )

_periods = ((1, 2), range(3, 11), range(11, 19), range(19, 37), range(37, 55), range(55, 87), range(87, 119))


elements = tuple(_table_e.split())
isotopes = {k: int(v) for k, v in zip(elements, _table_i.split())}
groups = {n: tuple(elements[y - 1] for y in x) for n, x in enumerate(_groups, start=1)}
periods = {n: tuple(elements[y - 1] for y in x) for n, x in enumerate(_periods, start=1)}
negativity = {k: float(v) if v != 'None' else None for k, v in zip(elements, _table_p.split())}
types = {'alkali': groups[1][1:], 'alkaline': groups[2], 'actinide': periods[7][2:17], 'lanthanide': periods[6][2:17],
         'transition': groups[3][:2] + reduce(lambda x, y: x + y, (groups[x] for x in range(4, 9))) +
         reduce(lambda x, y: x + y, (groups[x][:-1] for x in range(9, 12))),
         'post_transition': groups[12] + groups[13][1:-1] + groups[14][3:-1] + groups[15][4:-1] + groups[16][4:-1],
         'metalloid': groups[13][:1] + groups[14][1:3] + groups[15][2:4] + groups[16][3:4] + groups[17][4:5],
         'polyatomic': (groups[14][0], groups[15][1]) + groups[16][1:3],
         'diatomic': (groups[1][0], groups[15][0], groups[16][0]) + groups[17][:4], 'noble': groups[18][:-1],
         'unknown': periods[7][22:25] + periods[7][26:]}


# negativity fill
# lanthanoids fix
for x, l, r in zip(periods[6][3:14:2], periods[6][2:13:2], periods[6][4:15:2]):
    if negativity[x] is None:
        negativity[x] = (negativity[l] + negativity[r]) / 2

# actinoids fix
for l, a in zip(periods[6][3:17], periods[7][3:17]):
    if negativity[a] is None:
        negativity[a] = 2 * negativity[l] - negativity['Y']

# 7th period fix. linear correlation in groups
for g5, g6, g7 in zip(periods[5], periods[6][:3] + periods[6][17:], periods[7][:3] + periods[7][17:]):
    if negativity[g7] is None:
        negativity[g7] = 2 * negativity[g6] - negativity[g5]
