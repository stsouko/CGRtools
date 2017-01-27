# -*- coding: utf-8 -*-
#
#  Copyright 2014-2017 Ramil Nugmanov <stsouko@live.ru>
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
from networkx import Graph


class MoleculeContainer(Graph):
    def __init__(self, meta=None):
        super(MoleculeContainer, self).__init__()
        self.graph['meta'] = meta or {}

    @property
    def meta(self):
        if 'meta' not in self.graph:
            self.graph['meta'] = {}
        return self.graph['meta']


class ReactionContainer(dict):
    def __init__(self, substrats=None, products=None, meta=None):
        super(ReactionContainer, self).__init__(substrats=substrats or [], products=products or [], meta=meta or {})

    @property
    def substrats(self):
        return self['substrats']

    @property
    def products(self):
        return self['products']

    @property
    def meta(self):
        return self['meta']

    def copy(self):
        return ReactionContainer(substrats=[x.copy() for x in self['substrats']], meta=self['meta'].copy(),
                                 products=[x.copy() for x in self['products']])
