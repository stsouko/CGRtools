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
from .CGRcore import CGRcore
from .files import MoleculeContainer
from .CGRreactor import CGRreactor, patcher


class ReactMap(CGRreactor, CGRcore):
    def __init__(self, templates, stereo=False):
        CGRreactor.__init__(self, stereo=stereo, hyb=True, neighbors=True, isotope=True)
        CGRcore.__init__(self, extralabels=True)

        self.__core_templates = self.get_templates(templates)
        self.__prepare_search_patcher()

        self.__search_aromatic = self.get_template_searcher(self.__aromatize_templates())

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
        res['products'], arflag = self.__aromatize(res['products'])

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
                forcheck, ar = self.__aromatize(intermediate) if aromatized else (intermediate, False)
                if ar == aromatized:  # microoptimization
                    gm = self.get_CGR_matcher(forcheck, matrix['products'])
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

    def __aromatize(self, g):
        flag = False
        while True:
            patch = next(self.__search_aromatic(g), None)
            if patch:
                g = patcher(patch)
                flag = True
            else:
                break
        return g, flag

    def __aromatize_templates(self):
        relabel = {x: 1000 + x for x in range(1, 7)}

        aromatic = Graph()
        aromatic.add_edges_from([(1, 2, dict(s_bond=4, p_bond=4)), (2, 3, dict(s_bond=4, p_bond=4)),
                                 (3, 4, dict(s_bond=4, p_bond=4)), (4, 5, dict(s_bond=4, p_bond=4)),
                                 (5, 6, dict(s_bond=4, p_bond=4)), (6, 1, dict(s_bond=4, p_bond=4))])

        benzene0 = Graph()
        benzene0.add_edges_from([(1, 2, dict(s_bond=1, p_bond=1)), (2, 3, dict(s_bond=2, p_bond=2)),
                                 (3, 4, dict(s_bond=1, p_bond=1)), (4, 5, dict(s_bond=2, p_bond=2)),
                                 (5, 6, dict(s_bond=1, p_bond=1)), (6, 1, dict(s_bond=2, p_bond=2))])

        benzene1 = Graph()
        benzene1.add_edges_from([(1, 2, dict(s_bond=1, p_bond=1)), (2, 3, dict(s_bond=2, p_bond=2)),
                                 (3, 4, dict(s_bond=1, p_bond=1)), (4, 5, dict(s_bond=2, p_bond=2)),
                                 (5, 6, dict(s_bond=1, p_bond=1)), (6, 1, dict(s_bond=4, p_bond=4))])
        benzene2 = Graph()
        benzene2.add_edges_from([(1, 2, dict(s_bond=1, p_bond=1)), (2, 3, dict(s_bond=2, p_bond=2)),
                                 (3, 4, dict(s_bond=1, p_bond=1)), (4, 5, dict(s_bond=4, p_bond=4)),
                                 (5, 6, dict(s_bond=1, p_bond=1)), (6, 1, dict(s_bond=4, p_bond=4))])

        relabel_nodes(aromatic, relabel, copy=False)
        relabel_nodes(benzene0, relabel, copy=False)
        relabel_nodes(benzene1, relabel, copy=False)
        relabel_nodes(benzene2, relabel, copy=False)

        self.update_sp_marks(aromatic, copy=False)
        self.update_sp_marks(benzene0, copy=False)
        self.update_sp_marks(benzene1, copy=False)
        self.update_sp_marks(benzene2, copy=False)

        benzene0.__class__ = benzene1.__class__ = benzene2.__class__ = aromatic.__class__ = MoleculeContainer

        return [dict(substrats=benzene0, products=aromatic, meta=None),
                dict(substrats=benzene1, products=aromatic, meta=None),
                dict(substrats=benzene2, products=aromatic, meta=None)]
