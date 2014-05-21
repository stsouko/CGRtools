#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  csvconcat.py
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

import os
import argparse
from condenserpkg.rdfrw import Rdfrw


def filechk(filelist):
    if (os.path.exists(filelist['rdf']) or os.path.exists(filelist['sdf'])) and os.path.exists(filelist['csv']):
        if os.path.exists(filelist['header']):
            return filelist
        elif os.path.exists(filelist['csv'][:-3] + 'hdr'):
            filelist['header'] = filelist['csv'][:-3] + 'hdr'
            return filelist
        else:
            return None
    else:
        return None


def parser(rdfmeta, headlist, csvlist):
    metalist = []
    for i in rdfmeta:
        metalist += i['meta'].keys()
    metalist = list(set(metalist))

    for i in enumerate(metalist):
        print '%2d -- %s' % (i[0], i[1])
    meta = list(set(metalist) - set([metalist[int(x)] for x in raw_input('select params. e.g. 1,2,3… : ').split(',')]))
    headline = [x.replace(' ', '_') for x in meta] + [x[8:].strip() for x in headlist]
    head = [','.join(headline)]
    headlen = len(headline)

    for line in zip(csvlist, rdfmeta):
        dline = []
        for i in meta:
            try:
                try:
                    dval = '%1.20f' % float(line[1]['meta'][i])
                except:
                    dval = line[1]['meta'][i].replace(',', '.')
                dline += [dval]
            except:
                dline += ['']
        dline += line[0][2:].rstrip().split(';')
        head += [','.join(dline + ['0'] * (headlen - len(dline)))]
    return head


def parseSDF(file):
    with open(file) as f:
        meta = []
        a = 0
        for i in f:
            if "M  END" in i:
                buf = {}
            if ">  <" in i:
                a = i.strip()[4:-1]
            elif a:
                buf[a] = i.strip()
                a = 0
            elif "$$$$" in i:
                meta.append({'meta': buf})
    return meta


def main():
    rawopts = argparse.ArgumentParser(description="fragmentor csv refactor", epilog="created by stsouko")
    rawopts.add_argument("--version", "-v", action="version", version="0.1", default=False)
    rawopts.add_argument("--rdf", "-r", type=str, default="input.rdf", help="RDF inputfile")
    rawopts.add_argument("--sdf", "-s", type=str, default="input.sdf", help="SDF inputfile")
    rawopts.add_argument("--csv", "-c", type=str, default="input.csv", help="CSV inputfile")
    rawopts.add_argument("--header", "-d", type=str, default="input.hdr", help="HDR inputfile")
    rawopts.add_argument("--output", "-o", type=str, default="output.csv", help="CSV outputfile")
    options = vars(rawopts.parse_args())

    files = filechk(options)

    if files:
        try:
            out = open(files['output'], 'w')
        except IOError:
            print("error: can't write to outfile")
            return 1
        if os.path.exists(files['sdf']):
            parsedrdf = parseSDF(files['sdf'])
        else:
            rdf = Rdfrw(files['rdf'], '')
            parsedrdf = rdf.readdata()

        lines = parser(parsedrdf, open(files['header']), open(files['csv']))
        out.write('\n'.join(lines))
        return 0


if __name__ == '__main__':
    main()
