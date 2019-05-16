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
from functools import reduce


elements = tuple('''H                                                  He
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

# averaged isotopes
isotopes = tuple(map(int, '''  1                                                                   4
                               7   9                                          11  12  14  16  19  20
                              23  24                                          27  28  31  32  35  40
                              39  40  45  48  51  52  55  56  59  59  64  65  70  73  75  79  80  84
                              85  88  89  91  93  96  98 101 103 106 108 112 115 119 122 128 127 131
                             133 137 139

                                     140 141 144 145 150 152 157 159 163 165 167 169 173 175

                                         178 181 184 186 190 192 195 197 201 204 207 209 209 210 222 
                             223 226 227

                                     232 231 238 237 244 243 247 247 251 252 257 258 259 260

                                         261 270 269 270 270 278 281 281 285 278 289 289 293 297 294
                          '''.split()))

weights = tuple(map(float, '''
  1.007825035                                                                                       4.002600000 
  7.016000000   9.012180000  11.009300000  12.000000000  14.003074000  15.994914630  18.998403220  19.992440000 
 22.989770000  23.985000000  26.981540000  27.976927100  30.973762000  31.972070700  34.968852730  39.962400000 
 38.963700000  39.962600000

               44.955910000  47.947950000  50.943960000  51.940500000  54.938050000 
               55.934900000  58.933200000  57.935300000  62.929600000  63.929147000

                             68.925600000  73.921177400  74.921594200  79.916519600  78.918336100  83.911500000 
 84.911800000  87.905600000

               88.905860000  89.904700000  92.906400000  97.905400000  97.907200000 
              101.904300000 102.905500000 105.903500000 106.905100000 113.903400000

                            114.903900000 119.902200000 120.903800000 129.906200000 126.904500000 131.904100000 
132.905430000 137.905200000

              138.906360000 139.905400000 140.907660000 141.907719000 144.912800000 151.919700000 152.921200000 
              157.924099000 158.925350000 163.929200000 164.930300000 165.930300000 168.934230000 173.938900000

              174.940800000 179.946600000 180.948010000 183.951000000 186.955800000
              191.961500000 192.962900000 194.964800000 196.966560000 201.970617000
 
                            204.974400000 207.976627000 208.980390000 208.982400000 209.987100000 222.017500000
223.019700000 226.025410000

              227.027750000 232.038050000 231.035880000 238.050790000 237.048170000 244.064200000 243.061370000 
              247.070300000 247.070300000 251.079600000 252.082800000 257.095100000 258.098600000 259.100900000

              260.105400000 261.108700000 270.131000000 269.129000000 270.133000000
              270.134000000 278.156000000 281.165000000 281.166000000 285.177000000

                            278.000000000 289.190000000 289.000000000 293.204000000 297.000000000 294.000000000
'''.split()))

# H changed from FFFFFF to 909090, C changed from 909090 to
cpk = tuple('''
 #909090                                                                                         #D9FFFF
 #CC80FF #C2FF00                                         #FFB5B5 #000000 #3050F8 #FF0D0D #90E050 #B3E3F5
 #AB5CF2 #8AFF00                                         #BFA6A6 #F0C8A0 #FF8000 #C6C600 #1FF01F #80D1E3
 #8F40D4 #3DFF00 #E6E6E6 #BFC2C7 #A6A6AB #8A99C7 #9C7AC7 
                 #E06633 #F090A0 #50D050 #C88033 #7D80B0 #C28F8F #668F8F #BD80E3 #FFA100 #A62929 #5CB8D1
 #702EB0 #00FF00 #94FFFF #94E0E0 #73C2C9 #54B5B5 #3B9E9E 
                 #248F8F #0A7D8C #006985 #C0C0C0 #FFD98F #A67573 #668080 #9E63B5 #D47A00 #940094 #429EB0
 #57178F #00C900 #70D4FF 
                 #FFFFC7 #D9FFC7 #C7FFC7 #A3FFC7 #8FFFC7 #61FFC7 #45FFC7
                 #30FFC7 #1FFFC7 #00FF9C #00E675 #00D452 #00BF38 #00AB24 
                         #4DC2FF #4DA6FF #2194D6 #267DAB
                 #266696 #175487 #D0D0E0 #FFD123 #B8B8D0 #A6544D #575961 #9E4FB5 #AB5C00 #754F45 #428296
 #420066 #007D00 #70ABFA
                 #00BAFF #00A1FF #008FFF #0080FF #006BFF #545CF2 #785CE3
                 #8A4FE3 #A136D4 #B31FD4 #B31FBA #B30DA6 #BD0D87 #C70066
                         #CC0059 #D1004F #D90045 #E00038 
                 #E6002E #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026
'''.split())

# elements by groups
groups = ((1, 3, 11, 19, 37, 55, 87),
          (4, 12, 20, 38, 56, 88),
          (21, 39) + tuple(range(57, 72)) + tuple(range(89, 104))) + \
         tuple((22 + x, 40 + x, 72 + x, 104 + x) for x in range(9)) + \
         tuple((5 + x, 13 + x, 31 + x, 49 + x, 81 + x, 113 + x) for x in range(5)) + ((2, 10, 18, 36, 54, 86, 118), )

# elements by periods
periods = ((1, 2), tuple(range(3, 11)), tuple(range(11, 19)), tuple(range(19, 37)), tuple(range(37, 55)),
           tuple(range(55, 87)), tuple(range(87, 119)))

# elements by types
classes = {'alkali': groups[0][1:], 'alkaline': groups[1], 'actinide': periods[6][2:17], 'lanthanide': periods[5][2:17],
           'transition': groups[2][:2] + reduce(lambda x, y: x + y, (groups[x] for x in range(3, 8))) +
           reduce(lambda x, y: x + y, (groups[x][:-1] for x in range(8, 11))),
           'post_transition': groups[11] + groups[12][1:-1] + groups[13][3:-1] + groups[14][4:-1] + groups[15][4:-1],
           'metalloid': groups[12][:1] + groups[13][1:3] + groups[14][2:4] + groups[15][3:4] + groups[16][4:5],
           'polyatomic': (groups[13][0], groups[14][1]) + groups[15][1:3],
           'diatomic': (groups[0][0], groups[14][0], groups[15][0]) + groups[16][:4], 'noble': groups[17][:-1],
           'unknown': periods[6][22:25] + periods[6][26:]}


__all__ = ['elements', 'isotopes', 'groups', 'periods', 'classes', 'weights', 'cpk']
