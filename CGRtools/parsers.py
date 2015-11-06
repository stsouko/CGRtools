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
import argparse
from CGRtools.main_condenser import condenser_core
from CGRtools.main_fear import fear_core
from CGRtools.main_balanser import balanser_core
from CGRtools.main_mapper import mapper_core


def fear_common(parser):
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--c_rules", "-cr", type=argparse.FileType('r'), help="correct reactions type file (SDF)")
    group.add_argument("--e_rules", "-er", type=argparse.FileType('r'), help="incorrect reactions type file (RDF)")

    parser.add_argument("--map_repair", "-r", action='store_true', help="repair incorrect AAM [experimental]")


def condenser_common(parser):
    parser.add_argument("--type", "-t", type=str, default='0',
                        help="graph type: 0 - CGR, 1 - reagents only, 2 - products only, 101-199 - reagent 1,"
                             "2 or later, 201+ - product 1,2… (e.g. 202 - second product) -101/-199 or -201+ "
                             "- exclude reagent or product. comma-separated list of selected or excluded "
                             "reagents/products also supported (e.g 101,102 - only 1,2 molecules of reagents;"
                             " -101,[-]103 <second [-] no sense> - exclude 1 and 3 molecules of reagents). "
                             "also supported CGR on parts of reagents or/and products molecules (e.g. 101,102,-201 - "
                             "CGR on only first and second reagents molecules with all products molecules excluded "
                             "first).")
    parser.add_argument("--stereo", "-s", action='store_true', help="add stereo data")


def balanser(subparsers):
    parser = subparsers.add_parser('balanser', help='reaction balanser')
    parser.add_argument("--input", "-i", default="input.rdf", type=argparse.FileType('r'),
                        help="RDF inputfile")
    parser.add_argument("--output", "-o", default="output.rdf", type=argparse.FileType('w'),
                        help="RDF outputfile")

    condenser_common(parser)
    fear_common(parser)

    parser.add_argument("--balance", "-b", type=int, default=0, choices=[1, 2],
                        help="repair unbalanced reactions [experimental]. 1 - template based, 2 - template in existing "
                             "inputfile")
    parser.add_argument("--b_templates", "-B", type=argparse.FileType('r'),
                        help="RDF with reactions balancing rules")

    parser.set_defaults(func=balanser_core)


def condenser(subparsers):
    parser = subparsers.add_parser('condenser', help='CGR generator')
    parser.add_argument("--input", "-i", default="input.rdf", type=argparse.FileType('r'),
                        help="RDF inputfile")
    parser.add_argument("--output", "-o", default="output.sdf", type=argparse.FileType('w'),
                        help="SDF outputfile")

    condenser_common(parser)
    fear_common(parser)

    parser.add_argument("--balance", "-b", type=int, default=0, choices=[0, 1, 2],
                        help="repair unbalanced reactions [experimental]. 1 - template based, 2 - template in existing "
                             "inputfile")
    parser.add_argument("--b_templates", "-B", type=argparse.FileType('r'),
                        help="RDF with reactions balancing rules")

    parser.add_argument("--format", "-f", action='store_true', help="generate old format of CGR SDF")
    parser.set_defaults(func=condenser_core)


def fear(subparsers):
    parser = subparsers.add_parser('fear', help='reaction atom-to-atom mapping (AAM) checker')
    parser.add_argument("--input", "-i", default="input.rdf", type=argparse.FileType('r'), help="RDF inputfile")
    parser.add_argument("--output", "-o", default="output.rdf", type=argparse.FileType('w'), help="RDF outputfile")

    fear_common(parser)

    parser.set_defaults(func=fear_core)


def reactmap(subparsers):
    parser = subparsers.add_parser('reactmap', help='reaction atom-to-atom mapper (AAM)')
    parser.add_argument("--input", "-i", default="input.rdf", type=argparse.FileType('r'), help="RDF inputfile")
    parser.add_argument("--output", "-o", default="output.rdf", type=argparse.FileType('w'), help="RDF outputfile")

    parser.add_argument("--templates", "-t", type=argparse.FileType('r'),
                        help="RDF with reactions mapping rules")
    parser.set_defaults(func=mapper_core)
