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

_table = '''  1 H                                                                   4 He
              7 Li   9 Be  11 B   12 C   14 N   16 O   19 F                        20 Ne
             23 Na  24 Mg  27 Al  28 Si  31 P   32 S   35 Cl                       40 Ar
             39 K   40 Ca  45 Sc  48 Ti  51 V   52 Cr  55 Mn  56 Fe  59 Co  58 Ni
             63 Cu  64 Zn  69 Ga  74 Ge  75 As  80 Se  79 Br                       84 Kr
             85 Rb  88 Sr  89 Y   90 Zr  93 Nb  98 Mo 115 Tc 102 Ru 103 Rh 106 Pd 
            107 Ag 114 Cd 115 In 120 Sn 121 Sb 130 Te 127 I                       132 Xe
            133 Cs 138 Ba 139 La

                     140 Ce 141 Pr 142 Nd 163 Pm 152 Sm 153 Eu 158 Gd 159 Tb 164 Dy 165 Ho 166 Er 169 Tm 174 Yb 175 Lu

                                 180 Hf 181 Ta 184 W  187 Re 192 Os 193 Ir 195 Pt
            197 Au 202 Hg 205 Tl 208 Pb 209 Bi 218 Po 223 At 228 Rn
            232 Fr 234 Ra 236 Ac

                     232 Th 231 Pa 238 U  244 Np 247 Pu 249 Am 252 Cm 254 Bk 256 Cf 257 Es 259 Fm 260 Md 262 No 263 Lr
                     
                                 264 Rf 265 Db 266 Sg 267 Bh 277 Hs 271 Mt 281 Ds 
            272 Rg 285 Cn 286 Nh 289 Fl 289 Mc 293 Lv 294 Ts 294 Og
         '''

_atoms = _table.split()

elements = tuple(_atoms[1::2])
isotopes = {k: int(v) for k, v in zip(elements, _atoms[::2])}
