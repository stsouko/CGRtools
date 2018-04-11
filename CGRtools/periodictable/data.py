# -*- coding: utf-8 -*-
#
#  Copyright 2017, 2018 Ramil Nugmanov <stsouko@live.ru>
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
from functools import reduce


table_e = tuple('''H                                                  He
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
              '''.split())

# most stable isotope
table_i = tuple(map(int, '''  1                                                                   4
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
                        '''.split()))

groups = ((1, 3, 11, 19, 37, 55, 87),
          (4, 12, 20, 38, 56, 88),
          (21, 39) + tuple(range(57, 72)) + tuple(range(89, 104))) + \
         tuple((22 + x, 40 + x, 72 + x, 104 + x) for x in range(9)) + \
         tuple((5 + x, 13 + x, 31 + x, 49 + x, 81 + x, 113 + x) for x in range(5)) + ((2, 10, 18, 36, 54, 86, 118), )

periods = ((1, 2), tuple(range(3, 11)), tuple(range(11, 19)), tuple(range(19, 37)), tuple(range(37, 55)),
           tuple(range(55, 87)), tuple(range(87, 119)))


types = {'alkali': groups[0][1:], 'alkaline': groups[1], 'actinide': periods[6][2:17], 'lanthanide': periods[5][2:17],
         'transition': groups[2][:2] + reduce(lambda x, y: x + y, (groups[x] for x in range(3, 8))) +
         reduce(lambda x, y: x + y, (groups[x][:-1] for x in range(8, 11))),
         'post_transition': groups[11] + groups[12][1:-1] + groups[13][3:-1] + groups[14][4:-1] + groups[15][4:-1],
         'metalloid': groups[12][:1] + groups[13][1:3] + groups[14][2:4] + groups[15][3:4] + groups[16][4:5],
         'polyatomic': (groups[13][0], groups[14][1]) + groups[15][1:3],
         'diatomic': (groups[0][0], groups[14][0], groups[15][0]) + groups[16][:4], 'noble': groups[17][:-1],
         'unknown': periods[6][22:25] + periods[6][26:]}

electrons = (0,) + types['noble']
