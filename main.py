#!/usr/bin/env python3
# PYTHON_ARGCOMPLETE_OK
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
import importlib
from importlib.util import find_spec
from CGRtools.version import version
from CGRtools.parsers import condenser, balanser, fear, reactmap


def parse_args():
    parser = argparse.ArgumentParser(description="CGR tools", epilog="(c) Dr. Ramil Nugmanov", prog='cgrtools')
    parser.add_argument("--version", "-v", action="version", version=version(), default=False)
    subparsers = parser.add_subparsers(title='subcommands', description='available utilities')

    condenser(subparsers)
    balanser(subparsers)
    reactmap(subparsers)
    fear(subparsers)

    if find_spec('argcomplete'):
        argcomplete = importlib.import_module('argcomplete')
        argcomplete.autocomplete(parser)

    return parser


def main():
    parser = parse_args()
    args = parser.parse_args()
    if 'func' in args:
        args.func(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
