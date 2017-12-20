# -*- coding: utf-8 -*-
#
#  Copyright 2017 Ramil Nugmanov <stsouko@live.ru>
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
from pathlib import Path
from ..files import RDFread
from ..reactor import CGRreactor


class Aromatizer(CGRreactor):
    def __init__(self):
        CGRreactor.__init__(self)
        p = Path(__file__).parent
        with (p / 'aromatize.rdf').open() as f_a, (p / 'dearomatize.rdf').open() as f_d:
            raw_templates_a = RDFread(f_a, is_template=True).read()
            raw_templates_d = RDFread(f_d, is_template=True).read()
        self.__searcher_a = self.get_template_searcher(self.prepare_templates(raw_templates_a))
        self.__searcher_d = self.get_template_searcher(self.prepare_templates(raw_templates_d))

    def __call__(self, g):
        flag = False
        while True:  # dearomatize pyroles (furans, thiophenes) and quinones
            searcher = self.__searcher_d(g)
            match = next(searcher, None)
            if not match:
                break

            flag = True
            g = self.patcher(g, match.patch)

            for match in searcher:
                g = self.patcher(g, match.patch)

        while True:  # aromatize benzenes
            searcher = self.__searcher_a(g)
            match = next(searcher, None)
            if not match:
                return g, flag

            flag = True
            g = self.patcher(g, match.patch)

            for match in searcher:
                g = self.patcher(g, match.patch)
