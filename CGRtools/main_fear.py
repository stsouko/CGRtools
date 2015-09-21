#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright 2014 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.RDFread import RDFread
from CGRtools.RDFwrite import RDFwrite
from CGRtools.FEAR import FEAR


def fear_core(args):
    options = vars(args)
    inputdata = RDFread(args.input)
    outputdata = RDFwrite(args.output)
    fear = FEAR()
    #con = CGRcore(balance=0, type='0', **options)
    err = 0
    num = -1
    #for num, data in enumerate(inputdata.readdata()):
    #    if num % 100 == 0:
    #        print("reaction: %d" % (num + 1), file=sys.stderr)
        #try:
        #a = con.getFreaction(data)
        #outputdata.writedata(a)
        #except:
        #    err += 1
        #    print('reaction %d consist errors' % (num + 1), file=sys.stderr)
    #print('%d from %d reactions checked' % (num + 1 - err, num + 1), file=sys.stderr)


