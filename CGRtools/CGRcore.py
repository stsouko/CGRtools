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
from CGRtools.CGRpreparer import CGRPreparer
from CGRtools.CGRreactor import CGRReactor, patcher
import networkx as nx
from networkx.algorithms import isomorphism as gis


class CGRcore(CGRPreparer, CGRReactor):
    def __init__(self, **kwargs):
        CGRPreparer.__init__(self, kwargs['type'], extralabels=kwargs.get('extralabels', False))
        CGRReactor.__init__(self, stereo=kwargs['stereo'], neighbors=True, hyb=True)

        ''' balancer rules.
        '''
        self.__searchpatch = self.__chkmap(self.searchtemplate(self.gettemplates(kwargs['b_templates']),
                                                               speed=kwargs.get('speed', False)),
                                           iswhitelist=False) if kwargs['b_templates'] else lambda x: x
        ''' mapper rules.
        '''
        self.chkmap = lambda x: x
        if kwargs['balance'] != 2:
            if kwargs['c_rules']:
                self.chkmap = self.__chkmap(self.searchtemplate(self.gettemplates(kwargs['c_rules'], isreaction=False),
                                                                patch=False, speed=kwargs.get('speed', False)))
            elif kwargs['e_rules']:
                self.chkmap = self.__chkmap(self.searchtemplate(self.gettemplates(kwargs['e_rules']),
                                                                speed=kwargs.get('speed', False)),
                                            iswhitelist=False)

        if kwargs['balance'] == 1 and self.acceptrepair():
            self.__step_first, self.__step_second = self.__compose, True
        elif kwargs['balance'] == 2 and self.acceptrepair():
            self.__step_first, self.__step_second = patcher, False
        else:
            self.__step_first, self.__step_second = self.__compose, False

        self.__disstemplate = self.__disscenter()

    __popdict = dict(products=dict(edge=('s_bond', 's_stereo'), node=('s_charge', 's_stereo', 's_neighbors', 's_hyb'),
                                   ext_node=('s_neighbors', 's_hyb')),
                     substrats=dict(edge=('p_bond', 'p_stereo'), node=('p_charge', 'p_stereo', 'p_neighbors', 'p_hyb'),
                                    ext_node=('p_neighbors', 'p_hyb')))

    __attrcompose = dict(edges=dict(substrats=dict(p_bond='s_bond', p_stereo='s_stereo'),
                                    products=dict(s_bond='p_bond', s_stereo='p_stereo')),
                         nodes=dict(substrats=dict(p_charge='s_charge', p_stereo='s_stereo',
                                                   p_neighbors='s_neighbors', p_hyb='s_hyb'),
                                    products=dict(s_charge='p_charge', s_stereo='p_stereo',
                                                  s_neighbors='p_neighbors', s_hyb='p_hyb')))

    def __compose(self, data):
        """ remove from union graphs of products or substrats data about substrats or products
        """
        common = set(data['substrats']).intersection(data['products'])
        matrix = dict(products=data['products'].copy(), substrats=data['substrats'].copy())
        extended_common = set()

        """ remove bond, stereo, neighbors and hybridization states for common atoms.
        """
        for i in ('substrats', 'products'):
            for n, m in matrix[i].edges(common):
                extended_common.update([n, m])
                for j in self.__popdict[i]['edge']:
                    matrix[i][n][m].pop(j, None)

            for n in common:
                for j in self.__popdict[i]['node']:
                    matrix[i].node[n].pop(j, None)

        """ remove neighbors and hybridization states for common frontier atoms.
        """
        bubble = extended_common.difference(common)
        for i in ('substrats', 'products'):
            for n in bubble.intersection(data[i]):
                for j in self.__popdict[i]['ext_node']:
                    matrix[i].node[n].pop(j, None)

        """ compose graphs. kostyl
        """
        g = matrix['substrats']
        g.add_nodes_from(matrix['products'].nodes(data=True))
        g.add_edges_from(matrix['products'].edges(data=True))

        """ update sp_* marks
        """
        for n in extended_common:
            label = g.node[n]
            for s, p, sp in (('s_%s' % mark, 'p_%s' % mark, 'sp_%s' % mark) for mark in
                             ('neighbors', 'hyb', 'stereo', 'charge')):
                ls = label.get(s)
                lr = label.get(p)
                if ls != lr:
                    label[sp] = (ls, lr)
                elif ls is not None:
                    label[sp] = ls
                else:
                    label.pop(sp, None)

        for n, m in g.edges(common):
            label = g[n][m]
            for s, p, sp in (('s_%s' % mark, 'p_%s' % mark, 'sp_%s' % mark) for mark in ('bond', 'stereo')):
                ls = label.get(s)
                lr = label.get(p)
                if ls != lr:
                    label[sp] = (ls, lr)
                elif ls is not None:
                    label[sp] = ls
                else:
                    label.pop(sp, None)
        return g

    @staticmethod
    def __disscenter():
        g1 = nx.Graph()
        g2 = nx.Graph()

        g1.add_edges_from([(1, 2, dict(s_bond=None))])
        g2.add_edges_from([(1, 2, dict(p_bond=None))])
        return dict(substrats=('s_bond', [g1]), products=('p_bond', [g2]))

    def __chkmap(self, rsearcher, iswhitelist=True):
        def searcher(g):
            if iswhitelist:
                if not rsearcher(g):
                    g.graph.setdefault('CGR_REPORT', []).append('MAPPING INCORRECT')
                return g

            while True:
                match = next(rsearcher(g), None)
                if match:
                    g = patcher(match)
                    g.graph.setdefault('CGR_REPORT', []).append(match['meta'].get('CGR_TEMPLATE'))
                else:
                    return g

        return searcher

    def getCGR(self, data):
        g = self.__step_first(self.prepare(data))

        ''' fear check. currently skip if white list fails
        '''
        g = self.chkmap(g)
        if self.__step_second:
            ''' reaction balancer.
            '''
            g = self.clonesubgraphs(g)
            g = self.__searchpatch(g)

        meta = data['meta'].copy()
        meta['CGR_REPORT'] = ';; '.join(g.graph.get('CGR_REPORT', []))
        g.graph['meta'] = meta
        return g

    def dissCGR(self, g):
        tmp = dict(substrats=[], products=[], meta=g.graph['meta'])
        for category, (edge, pattern) in self.__disstemplate.items():
            components, *_ = self.getbondbrokengraph(g, pattern, gis.categorical_edge_match(edge, None))
            for mol in components:
                for n, m, edge_attr in mol.edges(data=True):
                    for i, j in self.__attrcompose['edges'][category].items():
                        if j in edge_attr:
                            mol[n][m][i] = edge_attr[j]
                for n, node_attr in mol.nodes(data=True):
                    for i, j in self.__attrcompose['nodes'][category].items():
                        if j in node_attr:
                            mol.node[n][i] = node_attr[j]
                tmp[category].append(mol)
        return tmp
