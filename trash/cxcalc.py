# -*- coding: utf-8 -*-
#
#  Copyright 2014, 2017 Dr. Ramil Nugmanov <stsouko@live.ru>
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
from subprocess import check_output

rep_dict = dict(UNKNOWN='u', WIGGLY='u', E='e', Z='z', R='r', S='s', r='re', s='si', M='m', P='p')


def stereo(file):
    result = []
    for data in check_output(['cxcalc', 'stereoanalysis', file]).decode().rstrip().split('\n\n'):
        structure = dict(atomstereo={}, bondstereo={}, allene={})
        if data != 'NONE':
            for x in data.split('\n'):
                line = x.split()
                m = rep_dict[line[-1]]
                mark = dict(s_stereo=m, p_stereo=m, sp_stereo=m)
                if line[0] == 'TETRAHEDRAL':
                    structure['atomstereo'][int(line[1][1:-1])] = mark
                elif line[0] == 'CISTRANS':
                    structure['bondstereo'][(int(line[1][1:-1]), int(line[2][:-1]))] = mark
                elif line[0] == 'AXIAL':
                    structure['allene'][(int(line[1][1:-1]), int(line[2][:-1]))] = mark

        result.append(structure)
    return result
