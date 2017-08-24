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
from functools import reduce
from itertools import product, combinations
from networkx import Graph, compose, has_path
from networkx.algorithms import isomorphism as gis
from .containers import CGRTemplate, MatchContainer
from .core import CGRcore


def patcher(structure, patch):
    """ remove edges bw common nodes. add edges from template and replace nodes data
    :param structure: MoleculeContainer
    :param patch: MoleculeContainer with replacement data
    """
    s = structure.copy()
    p = patch.copy()

    common = set(p).intersection(s)
    for i in common:
        pni = p.nodes[i]
        for j in {'s_charge', 's_hyb', 's_neighbors', 'p_charge', 'p_hyb', 'p_neighbors'}.intersection(pni):
            if isinstance(pni[j], dict):
                pni[j] = pni[j][s.nodes[i][j]]

    for m, n, a in p.edges(data=True):
        if m in common and n in common:
            for j in {'s_bond', 'p_bond'}.intersection(a):
                if isinstance(a[j], dict):
                    a[j] = a[j][s[m][n][j]]

    s.remove_edges_from(combinations(common, 2))
    composed = compose(s, p)
    composed.meta.update(s.meta)

    for (a1, a2), x in s._stereo_dict.items():
        if s.has_edge(a1, a2):
            composed.add_stereo(a1, a2, x.get('s'), x.get('p'))
    for (a1, a2), x in p._stereo_dict.items():
        composed.add_stereo(a1, a2, x.get('s'), x.get('p'))
    return composed


class CGRreactor(object):
    """ CGR isomorphism based operations
    """
    def __init__(self, extralabels=False, isotope=False, element=True, stereo=False):
        self.__rc_template = self.__reaction_center()

        gnm_sp, gem_sp, pcem, pcnm = [], ['sp_bond'], ['p_bond'], ['element', 'p_charge']
        if isotope:
            gnm_sp.append('isotope')
            pcnm.append('isotope')
        if element:
            gnm_sp.extend(['sp_charge', 'element'])
        if extralabels:
            gnm_sp.append('sp_neighbors')
            gnm_sp.append('sp_hyb')

        self.__node_match = gis.generic_node_match(gnm_sp, [None] * len(gnm_sp), [self.__list_eq] * len(gnm_sp))
        self.__edge_match = gis.generic_node_match(gem_sp, [None] * len(gem_sp), [self.__list_eq] * len(gem_sp))

        self.__node_match_products = gis.categorical_node_match(pcnm, [None] * len(pcnm))
        self.__edge_match_products = gis.categorical_edge_match(pcem, [None] * len(pcem))
        self.__edge_match_only_bond = gis.categorical_edge_match(['s_bond', 'p_bond'], [None] * 2)

        self.__pickle = dict(stereo=stereo, extralabels=extralabels, isotope=isotope, element=element)

    def pickle(self):
        """ return config. for pickling
        """
        return self.__pickle.copy()

    @classmethod
    def unpickle(cls, config):
        """ return CGRreactor object instance
        """
        args = {'stereo', 'extralabels', 'isotope', 'element'}
        if args.difference(config):
            raise Exception('Invalid config')
        return cls(**{k: v for k, v in config.items() if k in args})

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
        def searcher(g, skip_intersection=True):
            if skip_intersection:
                found = set()

            for i in templates:
                gm = self.get_cgr_matcher(g, i.substrats)
                for j in gm.subgraph_isomorphisms_iter():
                    matched_atoms = set(j.keys())
                    if skip_intersection:
                        if found.intersection(matched_atoms):
                            continue
                        else:
                            found.update(matched_atoms)

                    yield MatchContainer(mapping=list(matched_atoms), meta=i.meta,
                                         patch=self.__remap_group(i.products, g, {y: x for x, y in j.items()})[0])

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
        return CGRcore.split(g), lose_bonds

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
                tmp = CGRcore.union(tmp, remappedgroup)
                tmp.add_edge(j[k], mapping[lose_map[k]], **lose_bonds[(k, lose_map[k])])

        return tmp

    @staticmethod
    def __list_eq(a, b):
        return True if b is None else a in b if isinstance(b, list) else a == b

    @staticmethod
    def __simple_eq(a, b):
        return True if b is None else a == b

    @staticmethod
    def __remap_group(g, h, mapping):
        newmap = mapping.copy()
        newmap.update({x: y for x, y in zip(set(g).difference(newmap), set(range(1, 1000)).difference(h))})
        return g.remap(newmap, copy=True), newmap

    @staticmethod
    def get_templates(raw_templates):
        templates = []
        for template in raw_templates:
            products = reduce(CGRcore.union, template.products).copy()
            substrats = reduce(CGRcore.union, template.substrats).copy()

            common = set(products).intersection(substrats)
            for n in common:
                for j in {'s_charge', 's_hyb', 's_neighbors',
                          'p_charge', 'p_hyb', 'p_neighbors'}.intersection(products.nodes[n]):
                    if isinstance(products.nodes[n][j], list):
                        products.nodes[n][j] = {x: y for x, y in zip(substrats.nodes[n][j], products.nodes[n][j])}
                for j in ('s_x', 's_y', 's_z', 'p_x', 'p_y', 'p_z'):
                    products.nodes[n].pop(j)

            for m, n, a in products.edges(data=True):
                if m in common and n in common:
                    for j in {'s_bond', 'p_bond'}.intersection(a):
                        if isinstance(a[j], list):
                            a[j] = {x: y for x, y in zip(substrats[m][n][j], a[j])}

            substrats.remap({x: x + 1000 for x in substrats})
            products.remap({x: x + 1000 for x in products})

            templates.append(CGRTemplate(substrats, products, template.meta.copy()))
        return templates
