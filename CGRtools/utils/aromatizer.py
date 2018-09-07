# -*- coding: utf-8 -*-
#
#  Copyright 2017, 2018 Ramil Nugmanov <stsouko@live.ru>
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
from io import TextIOWrapper
from pkg_resources import resource_stream
from ..files import RDFread
from ..reactor import CGRreactor


class Aromatizer:
    def __init__(self):
        with resource_stream(__package__, 'aromatize.rdf') as fa, resource_stream(__package__, 'dearomatize.rdf') as fd:
            raw_templates_a = RDFread(TextIOWrapper(fa), is_template=True).read()
            raw_templates_d = RDFread(TextIOWrapper(fd), is_template=True).read()

        self.__reactor = r = CGRreactor()
        self.__searcher_a = r.get_template_searcher(r.prepare_templates(raw_templates_a))
        self.__searcher_d = r.get_template_searcher(r.prepare_templates(raw_templates_d))

    def __call__(self, g):
        flag = False
        patcher = self.__reactor.patcher
        while True:  # dearomatize pyroles (furans, thiophenes) and quinones
            searcher = self.__searcher_d(g)
            match = next(searcher, None)
            if not match:
                break

            flag = True
            g = patcher(g, match.patch)

            for match in searcher:
                g = patcher(g, match.patch)

        while True:  # aromatize benzenes
            searcher = self.__searcher_a(g)
            match = next(searcher, None)
            if not match:
                return g, flag

            flag = True
            g = patcher(g, match.patch)

            for match in searcher:
                g = patcher(g, match.patch)

    __visible = None
