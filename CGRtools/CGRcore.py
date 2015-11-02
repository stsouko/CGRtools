# -*- coding: utf-8 -*-
#
# Copyright 2014 Ramil Nugmanov <stsouko@live.ru>
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
import itertools
from CGRtools.CGRpreparer import CGRPreparer
from CGRtools.CGRreactor import CGRReactor
import networkx as nx
from networkx.algorithms import isomorphism as gis
from CGRtools.RDFread import RDFread
from CGRtools.SDFread import SDFread


class CGRcore(CGRPreparer, CGRReactor):
    def __init__(self, **kwargs):
        CGRPreparer.__init__(self, **kwargs)
        CGRReactor.__init__(self, self.__gettemplates(kwargs['b_templates']))

        self.chkmap = lambda x: x
        if kwargs['balance'] != 2:
            self.__maprepare = kwargs['map_repair']
            if kwargs['c_rules']:
                self.chkmap = self.__chkmap(self.searchtemplate(self.__gettemplates(kwargs['c_rules'])), True)
            elif kwargs['e_rules']:
                self.chkmap = self.__chkmap(self.searchtemplate(self.__gettemplates(kwargs['e_rules'])), False)

        if kwargs['balance'] == 1 and self.acceptrepair():
            self.__step_first, self.__step_second = self.__compose, True
        elif kwargs['balance'] == 2 and self.acceptrepair():
            self.__step_first, self.__step_second = self.__patch, False
        else:
            self.__step_first, self.__step_second = self.__compose, False

        self.__disstemplate = self.__disscenter()

    __popdict = dict(products=dict(edge=('s_bond', 's_stereo'), node=('s_charge', 's_stereo', 'x', 'y', 'z')),
                     substrats=dict(edge=('p_bond', 'p_stereo'), node=('p_charge', 'p_stereo')))

    __atrrset = {'p_bond', 'p_stereo', 's_bond', 's_stereo'}

    __attrcompose = dict(edges=dict(substrats=dict(p_bond='s_bond', p_stereo='s_stereo'),
                                    products=dict(s_bond='p_bond', s_stereo='p_stereo')),
                         nodes=dict(substrats=dict(p_charge='s_charge', p_stereo='s_stereo'),
                                    products=dict(s_charge='p_charge', s_stereo='p_stereo')))

    def __gettemplates(self, templates, isreaction=True):
        if templates:
            source = RDFread(templates) if isreaction else SDFread(templates)

            templates = []
            for template in source.readdata():
                if isreaction:
                    matrix = self.preparetemplate(template)
                    nx.relabel_nodes(matrix['substrats'], {x: x + 1000 for x in matrix['substrats']}, copy=False)
                    nx.relabel_nodes(matrix['products'], {x: x + 1000 for x in matrix['products']}, copy=False)
                else:
                    matrix = dict(substrats=template['structure'])
                templates.append(matrix)
            return templates

        return None

    def __compose(self, data):
        """ remove from union graphs of products or substrats data about substrats or products
        """
        common = set(data['substrats']).intersection(data['products'])
        matrix = dict(products=data['products'].copy(), substrats=data['substrats'].copy())

        """ remove bond and stereo states for common atoms.
        """
        for i in ('substrats', 'products'):
            for n, m in matrix[i].edges(common):
                for j in self.__popdict[i]['edge']:
                    matrix[i][n][m].pop(j)
            for n in common:
                for j in self.__popdict[i]['node']:
                    matrix[i].node[n].pop(j)

        """ compose graphs. kostyl
        """
        g = matrix['substrats']
        g.add_nodes_from(matrix['products'].nodes(data=True))
        g.add_edges_from(matrix['products'].edges(data=True))

        """ add lost data
        """
        for n, m, j in g.edges(common, data=True):
            for i in {'p_bond', 'p_stereo', 's_bond', 's_stereo'}.difference(j):
                g[n][m][i] = None
        return g

    @staticmethod
    def __patch(matrix):
        """ remove edges bw common nodes. add edges from template and replace nodes data
        """
        g = matrix['substrats'].copy()
        g.remove_edges_from(itertools.combinations(matrix['products'], 2))
        g = nx.compose(g, matrix['products'])
        return g

    def __disscenter(self):
        g1 = nx.Graph()
        g2 = nx.Graph()

        g1.add_edges_from([(1, 2, dict(s_bond=None))])
        g2.add_edges_from([(1, 2, dict(p_bond=None))])
        return dict(substrats=('s_bond', [g1]), products=('p_bond', [g2]))

    def __dissCGR(self, data):
        tmp = dict(substrats=[], products=[])
        for category, (edge, pattern) in self.__disstemplate.items():
            components, *_ = self.getbondbrokengraph(data, pattern, gis.categorical_edge_match(edge, None))
            for mol in components:
                for n, m, edge_attr in mol.edges(data=True):
                    for i, j in self.__attrcompose['edges'][category].items():
                        mol[n][m][i] = edge_attr[j]
                for n, node_attr in mol.nodes(data=True):
                    for i, j in self.__attrcompose['nodes'][category].items():
                        mol.node[n][i] = node_attr[j]
                tmp[category].append(self.getformattedcgr(mol))
        return tmp

    def __chkmap(self, rsearcher, rtype):
        def searcher(g):
            while True:
                match = rsearcher(g)
                if rtype:
                    if match:
                        return g
                else:
                    if not match:
                        return g
                    elif self.__maprepare:
                        g = self.__patch(match)
                        return g
                raise Exception

        return searcher

    def getCGR(self, data):
        g = self.__step_first(self.prepare(data))

        ''' fear check. try to remap
        '''
        g = self.chkmap(g)

        if self.__step_second:
            ''' reaction balancer.
            '''
            g = self.clonesubgraphs(g)
            while True:
                patch = self.searchpatch(g)
                if patch:
                    g = self.__patch(patch)
                else:
                    break

        return g

    def getFCGR(self, data):
        matrix = self.getformattedcgr(self.getCGR(data))
        matrix['meta'] = data['meta']
        return matrix

    def getFreaction(self, data):
        matrix = self.__dissCGR(self.getCGR(data))
        matrix['meta'] = data['meta']
        return matrix
