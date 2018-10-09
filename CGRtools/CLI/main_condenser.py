# -*- coding: utf-8 -*-
#
#  Copyright 2014-2018 Ramil Nugmanov <stsouko@live.ru>
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
from sys import stderr
from traceback import format_exc
from ..files import RDFread, SDFwrite
from ..preparer import CGRpreparer


def condenser_core(**kwargs):
    inputdata = RDFread(kwargs['input'])
    outputdata = SDFwrite(kwargs['output'], extralabels=kwargs['save_extralabels'] or kwargs['extralabels'])

    worker = CGRpreparer(cgr_type=kwargs['cgr_type'], extralabels=kwargs['extralabels'], balance=kwargs['balance'])

    err = 0
    num = 0
    for num, data in enumerate(inputdata, start=1):
        if num % 100 == 1:
            print("reaction: %d" % num, file=stderr)
        try:
            a = worker.condense(data)
            outputdata.write(a)
        except Exception:
            err += 1
            print('reaction %d consist errors: %s' % (num, format_exc()), file=stderr)
    print('%d from %d reactions condensed' % (num - err, num), file=stderr)

    return 0 if num and not err else 1 if num - err else 2
