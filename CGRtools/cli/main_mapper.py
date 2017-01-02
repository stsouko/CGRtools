#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright 2015, 2016 Ramil Nugmanov <stsouko@live.ru>
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
from ..files.RDFrw import RDFread, RDFwrite
from ..Reactmap import ReactMap


def mapper_core(**kwargs):
    inputdata = RDFread(kwargs['input'])
    outputdata = RDFwrite(kwargs['output'])
    mapper = ReactMap(**kwargs)
    err = 0
    num = 0

    for num, data in enumerate(inputdata, start=1):
        if num % 100 == 1:
            print("reaction: %d" % num, file=sys.stderr)
        try:
            a = mapper.map(data)
            outputdata.write(a)
        except Exception:
            err += 1
            print('reaction: %d' % num, file=sys.stderr)
            print('reaction %d consist errors: %s' % (num, traceback.format_exc()), file=sys.stderr)
            break

    dump = RDFwrite(kwargs['dump_templates'])
    for i in mapper.templates:
        dump.write(i)

    print('%d from %d reactions mapped' % (num - err, num), file=sys.stderr)
