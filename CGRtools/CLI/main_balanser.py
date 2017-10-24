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
from sys import stderr
from traceback import format_exc
from ..files.RDFrw import RDFread, RDFwrite
from ..preparer import CGRpreparer


def balanser_core(**kwargs):
    inputdata = RDFread(kwargs['input'])
    outputdata = RDFwrite(kwargs['output'], extralabels=kwargs['save_extralabels'])

    worker = CGRpreparer(cgr_type=kwargs['cgr_type'], extralabels=kwargs['extralabels'], stereo=kwargs['stereo'],
                         templates=RDFread(kwargs['templates'],
                                           is_template=True).read() if kwargs['b_templates'] else None,
                         balance=kwargs['balance'], element=kwargs['element'], isotope=kwargs['isotope'])

    err = 0
    num = 0
    for num, data in enumerate(inputdata, start=1):
        if num % 100 == 1:
            print("reaction: %d" % num, file=stderr)
        try:
            a = worker.condense(data)
            a = worker.dissociate(a)
            outputdata.write(a)
        except Exception:
            err += 1
            print('reaction %d consist errors: %s' % (num, format_exc()), file=stderr)
    print('%d from %d reactions balanced' % (num - err, num), file=stderr)

    return 0 if num and not err else 1 if num - err else 2
