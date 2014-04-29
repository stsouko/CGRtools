#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  condenser.py
#
#  Copyright 2013 nougmanoff <stsouko@live.ru>
#  This file is part of condenser.
#
#  condenser is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#

import argparse

from condenserpkg.rdfrw import Rdfrw
from condenserpkg.func import Condenser
from condenserpkg.version import version


def main():
    rawopts = argparse.ArgumentParser(description="CGR generator", epilog="created by stsouko")
    rawopts.add_argument("--version", "-v", action="version", version=version(), default=False)
    rawopts.add_argument("--input", "-i", type=str, default="input.rdf", help="RDF inputfile")
    rawopts.add_argument("--output", "-o", type=str, default="output.sdf", help="SDF outputfile")
    rawopts.add_argument("--type", "-t", type=int, default=0, help="graph type: 0 - CGR, 1 - reagents only, "
                                                                   "2 - products only, 11+ - reagent 1,2 or later, "
                                                                   "21+ - product 1,2… (e.g. 22 - second product) "
                                                                   "-11+ or -21+ - exclude reagent or product")
    rawopts.add_argument("--coords", "-c", type=int, default=1, help="write to SDF coordinates of: 1 - reagents, "
                                                                     " 2 - products, 3 - both (changed spec of MOL)")
    rawopts.add_argument("--charge", "-g", type=int, default=0,
                         help="charge data modifier: 1 - add changed pseudocharges,"
                              " 2 - add charge diff (only changed spec of MOL), "
                              "3 - both 1 and 2, 4 - use products charges (only changed spec of MOL),"
                              " 5 - use products pseudocharges, 6 - both 4 and 5, 7 - 1 and 4, 8 - 2 and 5")
    options = vars(rawopts.parse_args())

    rw = Rdfrw(options['input'], options['output'], options['coords'])
    result = []
    if options['type'] > 29 or options['type'] < -29 or 2 < options['type'] < 11 or -11 < options['type'] < 0 or \
                    options['type'] in (20, -20):
        print ('USE CORRECT TYPE, MAZAFAKA')
        return 0
    if options['charge'] > 8 or options['charge'] < 0:
        print ('USE CORRECT charge parameter, MAZAFAKA')
        return 0

    con = Condenser(options['charge'], options['type'])
    parseddata = rw.readdata()
    if parseddata:
        for num, data in enumerate(parseddata):
            try:
                result += [con.calc(data)]
            except:
                print 'reaction %d consist errors' % (num + 1)
                continue
        print '%d from %d reactions condensed' % (len(result), len(parseddata))
        if rw.writedata(result):
            print ('successfully done')
    return 0


if __name__ == '__main__':
    main()
