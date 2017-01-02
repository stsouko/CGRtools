#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright 2014-2016 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGR tools.
#
#  CGR tools is free software; you can redistribute it and/or modify
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
import sys
import traceback
from ..files.RDFrw import RDFread
from ..files.SDFrw import SDFwrite
from ..CGRpreparer import CGRcombo


def condenser_core(**kwargs):
    inputdata = RDFread(kwargs['input'])
    outputdata = SDFwrite(kwargs['output'], extralabels=kwargs['save_extralabels'])

    worker = CGRcombo(cgr_type=kwargs['cgr_type'], extralabels=kwargs['extralabels'], speed=kwargs['speed'],
                      b_templates=kwargs['b_templates'], m_templates=kwargs['m_templates'],
                      stereo=kwargs['stereo'], isotop=kwargs['isotop'], element=kwargs['element'], deep=kwargs['deep'])

    err = 0
    num = 0
    for num, data in enumerate(inputdata, start=1):
        if num % 100 == 1:
            print("reaction: %d" % num, file=sys.stderr)
        try:
            a = worker.getCGR(data)
            outputdata.write(a)
        except Exception:
            err += 1
            print('reaction %d consist errors: %s' % (num, traceback.format_exc()), file=sys.stderr)
    print('%d from %d reactions condensed' % (num - err, num), file=sys.stderr)

    return 0 if num and not err else 1 if num - err else 2
