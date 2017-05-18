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
from networkx import Graph, compose, has_path, union, union_all, relabel_nodes, connected_component_subgraphs
from itertools import product, combinations
from networkx.algorithms import isomorphism as gis
from .files.RDFrw import RDFread


def patcher(matrix):
    """ remove edges bw common nodes. add edges from template and replace nodes data
    :param matrix: dict
    """
    s = matrix['substrats'].copy()
    p = matrix['products'].copy()

    common = set(p).intersection(s)
    for i in common:
        for j in {'s_charge', 's_hyb', 's_neighbors', 's_stereo',
                  'p_charge', 'p_hyb', 'p_neighbors', 'p_stereo'}.intersection(p.node[i]):
            if isinstance(p.node[i][j], dict):
                p.node[i][j] = p.node[i][j][s.node[i][j]]

    for m, n, a in p.edges(data=True):
        if m in common and n in common:
            for j in {'s_bond', 'p_bond', 's_stereo', 'p_stereo'}.intersection(a):
                if isinstance(a[j], dict):
                    a[j] = a[j][s.edge[m][n][j]]

    s.remove_edges_from(combinations(common, 2))
    composed = compose(s, p)
    composed.meta.update(s.meta)
    return composed


def list_eq(a, b):
    return True if b is None else a in b if isinstance(b, list) else a == b


def simple_eq(a, b):
    return True if b is None else a == b


class CGRreactor(object):
    def __init__(self, stereo=False, hyb=False, neighbors=False, isotope=False, element=True):
        self.__rc_template = self.__reaction_center()

        gnm_sp, gem_sp, pcem, pcnm = [], ['sp_bond'], ['p_bond'], ['element', 'p_charge']
        if isotope:
            gnm_sp.append('isotope')
            pcnm.append('isotope')
        if element:
            gnm_sp.extend(['sp_charge', 'element'])
        if neighbors:
            gnm_sp.append('sp_neighbors')
        if hyb:
            gnm_sp.append('sp_hyb')
        if stereo:
            gnm_sp.append('sp_stereo')
            gem_sp.append('sp_stereo')
            pcem.append('p_stereo')
            pcnm.append('p_stereo')

        self.__node_match = gis.generic_node_match(gnm_sp, [None] * len(gnm_sp), [list_eq] * len(gnm_sp))
        self.__edge_match = gis.generic_node_match(gem_sp, [None] * len(gem_sp), [list_eq] * len(gem_sp))

        self.__node_match_products = gis.categorical_node_match(pcnm, [None] * len(pcnm))
        self.__edge_match_products = gis.categorical_edge_match(pcem, [None] * len(pcem))
        self.__edge_match_only_bond = gis.categorical_edge_match(['s_bond', 'p_bond'], [None] * 2)

    @staticmethod
    def __reaction_center():
        g1 = Graph()
        g2 = Graph()
        g1.add_edges_from([(1, 2, dict(s_bond=1, p_bond=None)), (2, 3, dict(s_bond=None, p_bond=1))])
        g2.add_edges_from([(1, 2, dict(s_bond=None, p_bond=1))])
        return [g1, g2]

    def get_cgr_matcher(self, g, h):
        return gis.GraphMatcher(g, h, node_match=self.__node_match, edge_match=self.__edge_match)

    def get_template_searcher(self, templates):
        def searcher(g):
            for i in templates:
                gm = self.get_cgr_matcher(g, i['substrats'])
                for j in gm.subgraph_isomorphisms_iter():
                    res = dict(substrats=g, meta=i['meta'],
                               products=self.__remap_group(i['products'], g, {y: x for x, y in j.items()})[0])
                    yield res

        return searcher

    @staticmethod
    def get_bond_broken_graph(g, rc_templates, edge_match):
        g = g.copy()
        lose_bonds = {}
        for i in rc_templates:
            gm = gis.GraphMatcher(g, i, edge_match=edge_match)
            for j in gm.subgraph_isomorphisms_iter():
                mapping = {y: x for x, y in j.items()}
                if 3 in mapping:
                    lose_bonds[(mapping[2], mapping[1])] = g[mapping[1]][mapping[2]]
                    g.remove_edge(mapping[2], mapping[3])
                    g.remove_edge(mapping[1], mapping[2])
                elif not any(has_path(g, *x) for x in product((y for x in lose_bonds for y in x), mapping.values())):
                    # запилить проверку связности атомов 1 или 2 с lose_map атомами
                    g.remove_edge(mapping[1], mapping[2])
        components = list(connected_component_subgraphs(g))
        return components, lose_bonds

    def clone_subgraphs(self, g):
        r_group = []
        x_group = {}
        r_group_clones = []
        newcomponents = []

        ''' search bond breaks and creations
        '''
        components, lose_bonds = self.get_bond_broken_graph(g, self.__rc_template, self.__edge_match_only_bond)
        lose_map = {x: y for x, y in lose_bonds}
        ''' extract subgraphs and sort by group type (R or X)
        '''
        x_terminals = set(lose_map.values())
        r_terminals = set(lose_map)

        for i in components:
            x_terminal_atom = x_terminals.intersection(i)
            r_terminal_atom = r_terminals.intersection(i)

            if x_terminal_atom:
                x_group[x_terminal_atom.pop()] = i
            elif r_terminal_atom:
                r_group.append([r_terminal_atom, i])
            else:
                newcomponents.append(i)
        ''' search similar R groups and patch.
        '''
        tmp = g.copy()
        for i in newcomponents:
            for k, j in r_group:
                gm = gis.GraphMatcher(j, i, node_match=self.__node_match_products,
                                      edge_match=self.__edge_match_products)
                ''' search for similar R-groups started from bond breaks.
                '''
                mapping = next((x for x in gm.subgraph_isomorphisms_iter() if k.intersection(x)), None)
                if mapping:
                    r_group_clones.append([k, mapping])
                    tmp = compose(tmp, self.__remap_group(j, tmp, mapping)[0])
                    break
        ''' add lose X groups to R groups
        '''
        for i, j in r_group_clones:
            for k in i:
                remappedgroup, mapping = self.__remap_group(x_group[lose_map[k]], tmp, {})
                tmp = union(tmp, remappedgroup)
                tmp.add_edge(j[k], mapping[lose_map[k]], **lose_bonds[(k, lose_map[k])])

        return tmp

    @staticmethod
    def __remap_group(g, h, mapping):
        newmap = mapping.copy()
        newmap.update({x: y for x, y in zip(set(g).difference(newmap), set(range(1, 1000)).difference(h))})
        return relabel_nodes(g, newmap), newmap

    @staticmethod
    def get_templates(raw_templates):
        templates = []
        for template in raw_templates:
            matrix = dict(meta=template.meta.copy())
            for i in ('products', 'substrats'):
                matrix[i] = union_all([x.copy() for x in template[i]])

            common = set(matrix['products']).intersection(matrix['substrats'])
            for n in common:
                for j in {'s_charge', 's_hyb', 's_neighbors', 's_stereo',
                          'p_charge', 'p_hyb', 'p_neighbors', 'p_stereo'}.intersection(matrix['products'].node[n]):
                    if isinstance(matrix['products'].node[n][j], list):
                        matrix['products'].node[n][j] = {x: y for x, y in zip(matrix['substrats'].node[n][j],
                                                                              matrix['products'].node[n][j])}
                for j in ('x', 'y', 'z'):
                    matrix['products'].node[n].pop(j)

            for m, n, a in matrix['products'].edges(data=True):
                if m in common and n in common:
                    for j in {'s_bond', 'p_bond', 's_stereo', 'p_stereo'}.intersection(a):
                        if isinstance(a[j], list):
                            matrix['products'].edge[m][n][j] = {x: y for x, y in
                                                                zip(matrix['substrats'].edge[m][n][j], a[j])}

            relabel_nodes(matrix['substrats'], {x: x + 1000 for x in matrix['substrats']}, copy=False)
            relabel_nodes(matrix['products'], {x: x + 1000 for x in matrix['products']}, copy=False)

            templates.append(matrix)
        return templates
