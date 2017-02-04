# -*- coding: utf-8 -*-
#
#  Copyright 2014-2017 Ramil Nugmanov <stsouko@live.ru>
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
import sys
import traceback
from ..files.RDFrw import RDFread, RDFwrite
from ..CGRcore import CGRcore
from ..FEAR import FEAR


def fear_core(**kwargs):
    inputdata = RDFread(kwargs['input'])
    outputdata = RDFwrite(kwargs['output'])

    fear = FEAR(isotope=kwargs['isotope'], stereo=kwargs['stereo'], hyb=kwargs['extralabels'],
                element=kwargs['element'], deep=kwargs['deep'])
    cgr = CGRcore(extralabels=kwargs['extralabels'])
    err = 0
    num = 0
    report = set()

    for num, data in enumerate(inputdata, start=1):
        if num % 100 == 1:
            print("reaction: %d" % num, file=sys.stderr)
        try:
            g = cgr.getCGR(data)
            h = fear.check_cgr(g, gennew=True)[1]
            report.update(x[1] for x in h)
            data.meta['REACTION_HASHES'] = ' : '.join(x[1] for x in h)
            outputdata.write(data)
        except Exception:
            err += 1
            print('reaction %d consist errors: %s' % (num, traceback.format_exc()), file=sys.stderr)
            break
    print(report)
    print('%d from %d reactions checked' % (num - err, num), file=sys.stderr)
