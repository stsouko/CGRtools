# -*- coding: utf-8 -*-
#
# Copyright 2014-2016 Ramil Nugmanov <stsouko@live.ru>
# This file is part of cgrtools.
#
# cgrtools is free software; you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
import networkx as nx
from CGRtools.RDFrw import RDFread
from CGRtools.CGRreactor import CGRreactor, patcher
from CGRtools.CGRcore import CGRcore


class CGRcombo(CGRcore):
    def __init__(self, cgr_type='0', extralabels=False, b_templates=None, m_templates=None, speed=False,
                 isotop=False, element=True, deep=0, stereo=False):
        CGRcore.__init__(self, cgr_type=cgr_type, extralabels=extralabels)
        self.__bal = CGRbalanser(b_templates, balanse_groups=True, speed=speed, stereo=stereo, isotop=isotop,
                                 extralabels=extralabels, element=element,
                                 deep=deep) if b_templates is not None else None

        self.__map = CGRbalanser(m_templates, balanse_groups=False, speed=speed, stereo=stereo, isotop=isotop,
                                 extralabels=extralabels, element=element,
                                 deep=deep) if m_templates is not None else None

    def getCGR(self, data):
        g = super(CGRcombo, self).getCGR(data)
        g = self.__map.prepare(g) if self.__map is not None else g
        g = self.__bal.prepare(g) if self.__bal is not None else g
        return g


class CGRbalanser(CGRreactor):
    def __init__(self, templates, speed=False, balanse_groups=True, stereo=False, extralabels=False,
                 isotop=False, element=True, deep=0):
        CGRreactor.__init__(self, stereo=stereo, hyb=extralabels, neighbors=extralabels,
                            isotop=isotop, element=element, deep=deep)

        self.__searcher = self.searchtemplate(self.get_templates(templates), speed=speed, patch=True)
        self.__balanse_groups = balanse_groups

    def prepare(self, g):
        report = []
        if self.__balanse_groups:
            g = self.clonesubgraphs(g)

        while True:
            match = next(self.__searcher(g), None)
            if match:
                g = patcher(match)
                report.append(match['meta'].get('CGR_TEMPLATE'))
            else:
                g.graph.setdefault('CGR_REPORT', []).extend(report)
                return g

    @staticmethod
    def get_templates(templates):
        if templates:
            source = RDFread(templates)

            _templates = []
            for template in source.read():
                matrix = dict(meta=template['meta'])
                for i in ('products', 'substrats'):
                    x = nx.union_all(template[i])
                    matrix[i] = x

                for i in set(matrix['products']).intersection(matrix['substrats']):
                    for j in {'sp_charge', 'sp_hyb',  # todo: это надо сделать работать
                              'sp_neighbors', 'sp_stereo'}.intersection(matrix['products'].node[i]):
                        if isinstance(matrix['products'].node[i][j], list):
                            matrix['products'].node[i][j] = {x: y for x, y in zip(matrix['substrats'].node[i][j],
                                                                                  matrix['products'].node[i][j])}

                nx.relabel_nodes(matrix['substrats'], {x: x + 1000 for x in matrix['substrats']}, copy=False)
                nx.relabel_nodes(matrix['products'], {x: x + 1000 for x in matrix['products']}, copy=False)

                _templates.append(matrix)
            return _templates

        return None
