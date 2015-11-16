#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright 2015, 2016 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.CGRcore import patcher
import networkx as nx


class ReactMap(CGRPreparer, CGRReactor):
    def __init__(self, debug=False, **kwargs):
        CGRPreparer.__init__(self, '0')
        CGRReactor.__init__(self, kwargs['stereo'])

        templates = self.gettemplates(kwargs['templates'])
        self.__searchpatch = self.searchtemplate(templates)
        self.__searcharomatic = self.searchtemplate(self.__aromatize_templates())
        self.__debug = debug

    def getMap(self, data):
        matrix = self.prepare(data)
        """ get aromatized products
        """
        matrix['products'], arflag = self.__aromatize(matrix['products'])

        """ prototype start. тут надо пилить и пилить чтобы на выхое получитьнормальный goodmap.
        """
        res = self.__reactpath(matrix, first=True, deep=3, aromatized=arflag)

        """ end. start normal final.
        """
        if res:
            data['meta']['CGRtools mapping'] = res[1]
        elif self.__debug:
            raise Exception

        data['substrats'] = [self.getformattedcgr((nx.relabel_nodes(mol, res[0], copy=True) if res else mol))
                             for mol in data['substrats']]
        data['products'] = [self.getformattedcgr(mol) for mol in data['products']]
        return data

    def __reactpath(self, matrix, aromatized=False, deep=3, first=False):
        morth = []
        if deep:
            paths = [dict(substrats=matrix['substrats'], products=nx.Graph())] if first \
                else self.__searchpatch(matrix['substrats'])
            for i in paths:
                intermediate = patcher(i)
                forcheck, ar = self.__aromatize(intermediate) if aromatized else (intermediate, False)
                if ar == aromatized:  # microoptimization
                    gm = self.spgraphmatcher(forcheck, matrix['products'])
                    if gm.subgraph_is_isomorphic():
                        return gm.mapping, float(i['meta']['AAMSCORE'])
                if deep > 1:
                    morth.append(intermediate)
            for i in morth:
                tmp = self.__reactpath(dict(substrats=i, products=matrix['products']), aromatized=aromatized, deep=deep - 1)
                if tmp:
                    return tmp
        return None

    def __aromatize(self, g):
        flag = False
        while True:
            patch = next(self.__searcharomatic(g), None)
            if patch:
                g = patcher(patch)
                flag = True
            else:
                break
        return g, flag

    @staticmethod
    def __aromatize_templates():
        templates = []
        benzene = nx.Graph()
        benzene.add_edges_from([(1001, 1002, dict(s_bond=1, p_bond=1)), (1002, 1003, dict(s_bond=2, p_bond=2)),
                                (1003, 1004, dict(s_bond=1, p_bond=1)), (1004, 1005, dict(s_bond=2, p_bond=2)),
                                (1005, 1006, dict(s_bond=1, p_bond=1)), (1006, 1001, dict(s_bond=2, p_bond=2))])
        aromatic = nx.Graph()
        aromatic.add_edges_from([(1001, 1002, dict(s_bond=4, p_bond=4)), (1002, 1003, dict(s_bond=4, p_bond=4)),
                                 (1003, 1004, dict(s_bond=4, p_bond=4)), (1004, 1005, dict(s_bond=4, p_bond=4)),
                                 (1005, 1006, dict(s_bond=4, p_bond=4)), (1006, 1001, dict(s_bond=4, p_bond=4))])

        benzene1 = nx.Graph()
        benzene1.add_edges_from([(1001, 1002, dict(s_bond=1, p_bond=1)), (1002, 1003, dict(s_bond=2, p_bond=2)),
                                 (1003, 1004, dict(s_bond=1, p_bond=1)), (1004, 1005, dict(s_bond=2, p_bond=2)),
                                 (1005, 1006, dict(s_bond=1, p_bond=1)), (1006, 1001, dict(s_bond=4, p_bond=4))])
        benzene2 = nx.Graph()
        benzene2.add_edges_from([(1001, 1002, dict(s_bond=1, p_bond=1)), (1002, 1003, dict(s_bond=2, p_bond=2)),
                                 (1003, 1004, dict(s_bond=1, p_bond=1)), (1004, 1005, dict(s_bond=4, p_bond=4)),
                                 (1005, 1006, dict(s_bond=1, p_bond=1)), (1006, 1001, dict(s_bond=4, p_bond=4))])

        templates.append(dict(substrats=benzene, products=aromatic, meta=None))
        templates.append(dict(substrats=benzene1, products=aromatic, meta=None))
        templates.append(dict(substrats=benzene2, products=aromatic, meta=None))
        return templates
