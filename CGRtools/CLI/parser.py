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
from argparse import ArgumentParser, ArgumentDefaultsHelpFormatter, FileType
from importlib.util import find_spec
from .main_condenser import condenser_core
from ..version import version


def _condenser(subparsers):
    parser = subparsers.add_parser('condenser', help='CGR generator',
                                   formatter_class=ArgumentDefaultsHelpFormatter)
    parser.add_argument("--input", "-i", default="input.rdf", type=FileType(),
                        help="RDF inputfile")
    parser.add_argument("--output", "-o", default="output.sdf", type=FileType('w'),
                        help="SDF outputfile")

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
    parser.add_argument("--save_extralabels", "-se", action='store_true', help="save extralabels data")
    parser.add_argument("--balance", "-B", action='store_true', help="Balance reactions")
    parser.set_defaults(func=condenser_core)


def argparser():
    parser = ArgumentParser(description="CGRtools", epilog="(c) Dr. Ramil Nugmanov", prog='cgrtools')
    parser.add_argument("--version", "-v", action="version", version=version(), default=False)
    subparsers = parser.add_subparsers(title='subcommands', description='available utilities')

    _condenser(subparsers)

    if find_spec('argcomplete'):
        from argcomplete import autocomplete
        autocomplete(parser)

    return parser
