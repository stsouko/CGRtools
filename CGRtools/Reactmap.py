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
from .CGRpreparer import CGRbalanser
from .CGRreactor import patcher
import networkx as nx


class ReactMap(CGRbalanser):
    def __init__(self, debug=False, **kwargs):
        CGRbalanser.__init__(self, kwargs['templates'], balanse_groups=False, stereo=kwargs['stereo'],
                             extralabels=True, isotop=True)

        self.__coretemplates = self.get_templates(kwargs['templates'])
        self.__preparesearchpatcher()

        self.__searcharomatic = self.get_template_searcher(self.__aromatize_templates())
        self.__debug = debug

    def __preparesearchpatcher(self, templates=None):
        if templates:
            self.templates = templates + self.templates
        else:
            self.templates = self.__coretemplates
        self.__searchpatch = self.get_template_searcher(self.templates)

    def map(self, data):
        matrix = self.prepare(data)
        """ get aromatized products
        """
        matrix['products'], arflag = self.__aromatize(matrix['products'])

        """ prototype start. тут надо пилить и пилить чтобы на выхое получить нормальный map.
        """
        res = self.__reactpath(matrix, first=True, deep=10, aromatized=arflag)

        """ end. start normal final.
        """
        if res:
            data['meta']['CGRtools mapping'] = res['score']
            self.__preparesearchpatcher(templates=[dict(substrats=res.get('substrats', data['substrats']),
                                                        products=res['products'],
                                                        meta=dict(AAMSCORE=res['score'], AAMID=10))])
        elif self.__debug:
            raise Exception

        data['substrats'] = [self.getformattedcgr((nx.relabel_nodes(mol, res['mapping'], copy=True) if res else mol))
                             for mol in data['substrats']]
        data['products'] = [self.getformattedcgr(mol) for mol in data['products']]
        return data

    def __reactpath(self, matrix, aromatized=False, deep=3, first=False, patchlist=None, score=0):
        morth = []
        if first:
            patchlist = []

        if deep:
            paths = [dict(substrats=matrix['substrats'], products=nx.Graph(), meta=dict(AAMSCORE=0, AAMID=0))] \
                if first else self.__searchpatch(matrix['substrats'])
            for i in paths:
                intermediate = patcher(i)
                forcheck, ar = self.__aromatize(intermediate) if aromatized else (intermediate, False)
                if ar == aromatized:  # microoptimization
                    gm = self.get_CGR_matcher(forcheck, matrix['products'])
                    if gm.subgraph_is_isomorphic():
                        # todo: generate new mapping rule
                        base = nx.Graph()
                        print(len(patchlist))
                        for x in (patchlist + [i['products']]):
                            base = patcher(dict(substrats=base, products=x))
                        return dict(mapping=gm.mapping, score=score + float(i['meta']['AAMSCORE']), products=base)
                if deep > 1:
                    morth.append((intermediate, patchlist + [i['products']], score + float(i['meta']['AAMSCORE'])))
            for i, j, s in morth:
                tmp = self.__reactpath(dict(substrats=i, products=matrix['products']), aromatized=aromatized,
                                       deep=deep-1, patchlist=j, score=s)
                if tmp:
                    if first:
                        # todo: slice substrat atoms
                        g = matrix['substrats'].copy()
                        g.remove_nodes_from(set(g).difference(tmp['products']))
                        tmp['substrats'] = g
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
        benzene0 = nx.Graph()
        benzene0.add_edges_from([(1001, 1002, dict(s_bond=1, p_bond=1)), (1002, 1003, dict(s_bond=2, p_bond=2)),
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

        templates.append(dict(substrats=benzene0, products=aromatic, meta=None))
        templates.append(dict(substrats=benzene1, products=aromatic, meta=None))
        templates.append(dict(substrats=benzene2, products=aromatic, meta=None))
        return templates
