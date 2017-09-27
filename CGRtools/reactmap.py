# -*- coding: utf-8 -*-
#
#  Copyright 2015-2017 Ramil Nugmanov <stsouko@live.ru>
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
from networkx import Graph, relabel_nodes
from .preparer import CGRpreparer
from .reactor import CGRreactor, patcher
from .utils.aromatize import Aromatize


aromatize = Aromatize()


class ReactMap(CGRreactor, CGRpreparer):
    def __init__(self, templates, stereo=False):
        CGRreactor.__init__(self, stereo=stereo, extralabels=True, isotope=True)
        CGRpreparer.__init__(self, extralabels=True)

        self.__core_templates = self.get_templates(templates)
        self.__prepare_search_patcher()

    def __prepare_search_patcher(self, templates=None):
        if templates:
            self.templates = templates + self.templates
        else:
            self.templates = self.__core_templates
        self.__search_patch = self.get_template_searcher(self.templates)

    def map(self, data):
        res = self.merge_mols(data)
        self.set_labels(res['substrats'], copy=False)
        self.set_labels(res['products'], copy=False)

        """ get aromatized products
        """
        res['products'], arflag = aromatize.get(res['products'])

        """ prototype start. тут надо пилить и пилить чтобы на выхое получить нормальный map.
        """
        res = self.__reactpath(res, first=True, deep=10, aromatized=arflag)

        """ end. start normal final.
        """
        if res:
            data['meta']['CGRtools mapping'] = res['score']
            self.__prepare_search_patcher(templates=[dict(substrats=res.get('substrats', data['substrats']),
                                                          products=res['products'],
                                                          meta=dict(AAMSCORE=res['score'], AAMID=10))])

        if res:
            data['substrats'] = [relabel_nodes(mol, res['mapping']) for mol in data['substrats']]

        return data

    def __reactpath(self, matrix, aromatized=False, deep=3, first=False, patchlist=None, score=0):
        morth = []
        if first:
            patchlist = []

        if deep:
            paths = [dict(substrats=matrix['substrats'], products=Graph(), meta=dict(AAMSCORE=0, AAMID=0))] \
                if first else self.__search_patch(matrix['substrats'])
            for i in paths:
                intermediate = patcher(i)
                forcheck, ar = aromatize.get(intermediate) if aromatized else (intermediate, False)
                if ar == aromatized:  # microoptimization
                    gm = self.get_cgr_matcher(forcheck, matrix['products'])
                    if gm.subgraph_is_isomorphic():
                        # todo: generate new mapping rule
                        base = Graph()
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
