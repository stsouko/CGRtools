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
from collections import defaultdict
from networkx.algorithms import isomorphism as gis
import networkx as nx


class CGRReactor(object):
    def __init__(self, templates):
        self.__attrcompose = dict(edges=dict(substrats=dict(p_bond='s_bond', p_stereo='s_stereo'),
                                             products=dict(s_bond='p_bond', s_stereo='p_stereo')),
                                  nodes=dict(substrats=dict(p_charge='s_charge', p_stereo='s_stereo'),
                                             products=dict(s_charge='p_charge', s_stereo='p_stereo')))
        self.searchtemplate = self.__searchtemplate if templates else lambda _: None
        self.__templates = templates
        self.__rctemplate = self.__reactioncenter()
        self.__node_match = gis.categorical_node_match(['s_stereo', 'p_stereo', 'element', 'isotop', 's_charge',
                                                        'p_charge'], [None] * 6)
        self.__edge_match = gis.categorical_edge_match(['s_bond', 'p_bond', 's_stereo', 'p_stereo'], [None] * 4)

        self.__node_match_products = gis.categorical_node_match(['p_stereo', 'element', 'isotop', 'p_charge'],
                                                                [None] * 6)
        self.__edge_match_products = gis.categorical_edge_match(['p_bond', 'p_stereo'], [None] * 2)

        self.__edge_match_only_bond = gis.categorical_edge_match(['s_bond', 'p_bond'], [None] * 2)

    def __reactioncenter(self):
        g1 = nx.Graph()
        g2 = nx.Graph()
        g3 = nx.Graph()
        g1.add_edges_from([(1, 2, dict(s_bond=1, p_bond=None)), (2, 3, dict(s_bond=None, p_bond=1)),
                           (3, 4, dict(s_bond=1, p_bond=None))])
        g2.add_edges_from([(1, 2, dict(s_bond=1, p_bond=None)), (2, 3, dict(s_bond=None, p_bond=1))])

        g3.add_edges_from([(1, 2, dict(s_bond=None, p_bond=1))])
        return g1, g2, g3

    def __searchtemplate(self, data):
        for i in self.__templates:
            gm = gis.GraphMatcher(data, i['substrats'], node_match=self.__node_match, edge_match=self.__edge_match)
            if gm.subgraph_is_isomorphic():
                mapping = {y: x for x, y in gm.mapping.items()}
                mapping.update({x: y for x, y in zip(set(i['products']).difference(gm.mapping.values()),
                                                     set(range(1, 1000)).difference(data))})
                return dict(substrats=data, products=nx.relabel_nodes(i['products'], mapping, copy=True))

        return None

    @staticmethod
    def getbondbrokengraph(data, rc_templates, edge_match):
        lose_map = {}
        lose_bonds = defaultdict(dict)
        for i in rc_templates:
            gm = gis.GraphMatcher(data, i, edge_match=edge_match)
            for j in gm.subgraph_isomorphisms_iter():
                mapping = {y: x for x, y in j.items()}
                if mapping.get(3):
                    lose_map[mapping[2]] = mapping[1]
                    lose_bonds[mapping[2]][mapping[1]] = data.edge[mapping[1]][mapping[2]]
                    data.remove_edge(mapping[2], mapping[3])
                    if mapping.get(4):
                        lose_map[mapping[3]] = mapping[4]
                        lose_bonds[mapping[3]][mapping[4]] = data.edge[mapping[3]][mapping[4]]
                        data.remove_edge(mapping[3], mapping[4])

                data.remove_edge(mapping[1], mapping[2])
        components = list(nx.connected_component_subgraphs(data))
        return components, lose_bonds, lose_map

    def clonesubgraphs(self, data):
        r_group = {}
        x_group = {}
        r_group_clones = defaultdict(list)

        ''' search bond breaks and creations
        '''
        components, lose_bonds, lose_map = self.getbondbrokengraph(data.copy(), self.__rctemplate,
                                                                   self.__edge_match_only_bond)
        ''' extract subgraphs and sort by group type (R or X)
        '''
        setlose = set(lose_map.values())
        setlosekey = set(lose_map)
        newcomponents = []
        for i in components:
            x_terminal_atom = setlose.intersection(i)
            r_terminal_atom = setlosekey.intersection(i)

            if x_terminal_atom:
                x_group[x_terminal_atom.pop()] = i
            elif r_terminal_atom:
                r_group[r_terminal_atom.pop()] = i
            else:
                newcomponents.append(i)
        components = newcomponents
        ''' search similar R groups and patch.
        '''
        tmp = data.copy()
        for i in components:
            for k, j in r_group.items():
                gm = gis.GraphMatcher(j, i, node_match=self.__node_match_products,
                                      edge_match=self.__edge_match_products)
                ''' search for similar R-groups started from bond breaks.
                '''
                mapping = [x for x in gm.subgraph_isomorphisms_iter() if k in x]
                if mapping:
                    r_group_clones[k].append(gm.mapping[k])
                    break
        ''' add lose X groups to R groups
        '''
        for i, j in r_group_clones.items():
            for k in j:
                mapping = {x: y for x, y in zip(x_group[lose_map[i]], set(range(1, 1000)).difference(tmp))}
                X = nx.relabel_nodes(x_group[lose_map[i]], mapping, copy=True)
                tmp = nx.union(tmp, X)
                tmp.add_edge(k, mapping[lose_map[i]], attr_dict=lose_bonds[i][lose_map[i]])

        return tmp
