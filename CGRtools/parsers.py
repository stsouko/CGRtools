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
import argparse
import importlib
from importlib.util import find_spec
from .version import version
from .cli.main_mapper import mapper_core
from .cli.main_balanser import balanser_core
from .cli.main_condenser import condenser_core
from .cli.main_fear import fear_core


def fear_common(parser):
    parser.add_argument("--element", "-e", action='store_false', help="FEAR use element data")
    parser.add_argument("--isotop", "-I", action='store_true', help="FEAR use isotop data")
    parser.add_argument("--deep", "-d", default=0, type=int, help="FEAR deep of scan")


def stereo_common(parser):
    parser.add_argument("--stereo", "-s", action='store_true', help="add stereo data")


def extra_common(parser):
    parser.add_argument("--save_extralabels", "-se", action='store_true', help="save extralabels data")


def condenser_common(parser):
    parser.add_argument("--cgr_type", "-t", type=str, default='0',
                        help="graph type: 0 - CGR, 1 - reagents only, 2 - products only, 101-199 - reagent 1,"
                             "2 or later, 201+ - product 1,2… (e.g. 202 - second product) -101/-199 or -201+ "
                             "- exclude reagent or product. comma-separated list of selected or excluded "
                             "reagents/products also supported (e.g 101,102 - only 1,2 molecules of reagents;"
                             " -101,[-]103 <second [-] no sense> - exclude 1 and 3 molecules of reagents). "
                             "also supported CGR on parts of reagents or/and products molecules (e.g. 101,102,-201 - "
                             "CGR on only first and second reagents molecules with all products molecules excluded "
                             "first).")
    parser.add_argument("--extralabels", "-E", action='store_true',
                        help="generate atom hybridization and neighbors labels")


def balanser_common(parser):
    parser.add_argument("--speed", "-S", action='store_true', help="speed up by FEAR")
    parser.add_argument("--m_templates", "-M", type=argparse.FileType('r'), default=None,
                        help="RDF with reactions remapping rules")
    parser.add_argument("--b_templates", "-B", type=argparse.FileType('r'), default=None,
                        help="RDF with reactions balancing rules")


def balanser(subparsers):
    parser = subparsers.add_parser('balanser', help='reaction balanser',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", "-i", default="input.rdf", type=argparse.FileType('r'),
                        help="RDF inputfile")
    parser.add_argument("--output", "-o", default="output.rdf", type=argparse.FileType('w'),
                        help="RDF outputfile")

    condenser_common(parser)
    extra_common(parser)
    stereo_common(parser)
    balanser_common(parser)
    fear_common(parser)

    parser.set_defaults(func=balanser_core)


def condenser(subparsers):
    parser = subparsers.add_parser('condenser', help='CGR generator',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", "-i", default="input.rdf", type=argparse.FileType('r'),
                        help="RDF inputfile")
    parser.add_argument("--output", "-o", default="output.sdf", type=argparse.FileType('w'),
                        help="SDF outputfile")

    condenser_common(parser)
    extra_common(parser)
    stereo_common(parser)
    balanser_common(parser)
    fear_common(parser)

    parser.set_defaults(func=condenser_core)


def fear(subparsers):
    parser = subparsers.add_parser('fear', help='reaction atom-to-atom mapping (AAM) checker',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", "-i", default="input.rdf", type=argparse.FileType('r'), help="RDF inputfile")
    parser.add_argument("--output", "-o", default="output.rdf", type=argparse.FileType('w'), help="RDF outputfile")

    fear_common(parser)
    stereo_common(parser)
    condenser_common(parser)

    parser.set_defaults(func=fear_core)


def reactmap(subparsers):
    parser = subparsers.add_parser('reactmap', help='reaction atom-to-atom mapper (AAM)',
                                   formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", "-i", default="input.rdf", type=argparse.FileType('r'), help="RDF inputfile")
    parser.add_argument("--output", "-o", default="output.rdf", type=argparse.FileType('w'), help="RDF outputfile")
    parser.add_argument("--templates", "-t", type=argparse.FileType('r'), help="RDF with reactions mapping rules")
    parser.add_argument("--dump_templates", "-d", default="dump.rdf", type=argparse.FileType('w'),
                        help="save RDF with reactions mapping rules")

    stereo_common(parser)

    parser.set_defaults(func=mapper_core)


def parse_args():
    parser = argparse.ArgumentParser(description="CGR tools", epilog="(c) Dr. Ramil Nugmanov", prog='cgrtools')
    parser.add_argument("--version", "-v", action="version", version=version(), default=False)
    subparsers = parser.add_subparsers(title='subcommands', description='available utilities')

    condenser(subparsers)
    balanser(subparsers)
    fear(subparsers)
    reactmap(subparsers)

    if find_spec('argcomplete'):
        argcomplete = importlib.import_module('argcomplete')
        argcomplete.autocomplete(parser)

    return parser


def launcher():
    parser = parse_args()
    args = parser.parse_args()
    if 'func' in args:
        args.func(**vars(args))
    else:
        parser.print_help()
