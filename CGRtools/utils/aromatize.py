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
from ..files.RDFrw import RDFread
from ..reactor import CGRreactor, patcher


class Aromatize(CGRreactor):
    def __init__(self):
        CGRreactor.__init__(self)
        with (Path(__file__).parent / 'aromatize.rdf').open() as f:
            raw_templates = RDFread(f).read()
        self.__searcher = self.get_template_searcher(self.get_templates(raw_templates))
        #todo: Implement optimizations!

    def get(self, g):
        flag = False
        while True:
            patch = next(self.__searcher(g), None)
            if patch:
                g = patcher(patch)
                flag = True
            else:
                break
        return g, flag
