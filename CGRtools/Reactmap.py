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
from CGRtools.CGRpreparer import CGRPreparer
from CGRtools.CGRreactor import CGRReactor
from CGRtools.RDFread import RDFread
from CGRtools.CGRcore import patcher
import networkx as nx


class ReactMap(CGRPreparer, CGRReactor):
    def __init__(self, **kwargs):
        CGRPreparer.__init__(self, type='0', stereo=True, **kwargs)
        CGRReactor.__init__(self, self.__gettemplates(kwargs['templates']), deep=True)

    def __gettemplates(self, templates):
        if templates:
            source = RDFread(templates)
            templates = []
            for template in source.readdata():
                matrix = self.preparetemplate(template)
                nx.relabel_nodes(matrix['substrats'], {x: x + 1000 for x in matrix['substrats']}, copy=False)
                nx.relabel_nodes(matrix['products'], {x: x + 1000 for x in matrix['products']}, copy=False)
                templates.append(matrix)
            return templates

        return None

    def getMap(self, data):
        self.__maps = []

        matrix = self.prepare(data)
        self.__reactpath(matrix, first=True, deep=10)
        print(self.__maps)
        goodmap = self.__maps[0]
        data['substrats'] = [self.getformattedcgr(nx.relabel_nodes(mol, goodmap, copy=True)) for mol in data['substrats']]
        data['products'] = [self.getformattedcgr(mol) for mol in data['products']]
        return data

    def __reactpath(self, matrix, deep=10, first=False):
        if deep:
            paths = [dict(substrats=matrix['substrats'], products=nx.Graph())] if first \
                else self.searchpatch(matrix['substrats'])
            for i in paths:
                intermediate = patcher(i)
                gm = self.spgraphmatcher(intermediate, matrix['products'])
                if gm.subgraph_is_isomorphic():
                    self.__maps.append(gm.mapping)
                else:
                    self.__reactpath(dict(substrats=intermediate, products=matrix['products']), deep=deep - 1)
