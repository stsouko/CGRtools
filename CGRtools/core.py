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
from networkx import connected_components
from . import InvalidData
from .containers import MoleculeContainer, CGRContainer


class CGRcore(object):
    @staticmethod
    def split(m, meta=False):
        return [m.subgraph(c, meta=meta) for c in connected_components(m)]

    @classmethod
    def union(cls, m1, m2):
        if set(m1) & set(m2):
            raise InvalidData('The node sets of m1 and m2 are not disjoint.')

        u = CGRContainer() if isinstance(m1, CGRContainer) or isinstance(m2, CGRContainer) else MoleculeContainer()

        u.add_nodes_from(m1.nodes(data=True))
        u.add_nodes_from(m2.nodes(data=True))
        u.add_edges_from(m1.edges(data=True))
        u.add_edges_from(m2.edges(data=True))

        def fix(m):
            for n, m, attr in m.edges(data=True):
                u.add_edge(n, m, **cls.__fix_attr(attr, u._edge_marks))
            for n, attr in m.nodes(data=True):
                u.add_node(n, p_x=attr['s_x'], p_y=attr['s_y'], p_z=attr['s_z'], **cls.__fix_attr(attr, u._node_marks))

        if isinstance(u, CGRContainer) and not isinstance(m1, CGRContainer):
            fix(m1)
        elif isinstance(u, CGRContainer) and not isinstance(m2, CGRContainer):
            fix(m2)
        return u

    @staticmethod
    def __fix_attr(attr, marks):
        tmp = {}
        for s, p, sp in marks:
            if s in attr:
                tmp[p] = tmp[sp] = attr[s]
        return tmp

    @classmethod
    def compose(cls, m1, m2):
        """ remove from union graphs of products or reagents data about reagents or products
        """
        common = set(m1).intersection(m2)
        extended_common = set()
        h = CGRContainer()

        """ remove bond, neighbors and hybridization states for common atoms.
        """
        for i, g in (('reagents', m1), ('products' if isinstance(m2, CGRContainer) else 'non_cgr', m2)):
            pdi = cls.__popdict[i]
            ext_common = common.copy()
            e_pop, n_pop, x_pop = pdi['edge'], pdi['node'], pdi['ext_node']
            for n, m, attr in g.edges(common, data=True):
                ext_common.add(n)
                ext_common.add(m)
                bond = {e_pop[k]: v for k, v in attr.items() if k in e_pop}
                if bond:
                    h.add_edge(n, m, **bond)

            uniq = set(g).difference(ext_common)
            for n, m, attr in g.edges(uniq, data=True):
                h.add_edge(n, m, **attr)

            for n in common:
                h.add_node(n, **{n_pop[k]: v for k, v in g.nodes[n].items() if k in n_pop})

            for n in ext_common.difference(common):
                h.add_node(n, **{n_pop[k]: v for k, v in g.nodes[n].items() if k not in x_pop})

            for n in uniq:
                h.add_node(n, **g.nodes[n])

            extended_common.update(ext_common)

        """ update sp_* marks
        """
        h.fix_data(nodes_bunch=extended_common, edges_bunch=common)
        return h

    __popdict = dict(products=dict(edge=dict(p_bond='p_bond', p_stereo='p_stereo'),
                                   ext_node=('s_neighbors', 's_hyb', 'sp_neighbors', 'sp_hyb'),
                                   node=dict(p_charge='p_charge', p_neighbors='p_neighbors', p_hyb='p_hyb', p_x='p_x',
                                             p_y='p_y', p_z='p_z', p_stereo='p_stereo',
                                             mark='mark', element='element', map='map')),
                     reagents=dict(edge=dict(s_bond='s_bond', s_stereo='s_stereo'),
                                   ext_node=('p_neighbors', 'p_hyb', 'sp_neighbors', 'sp_hyb'),
                                   node=dict(s_charge='s_charge', s_neighbors='s_neighbors', s_hyb='s_hyb', s_x='s_x',
                                             s_y='s_y', s_z='s_z', s_stereo='s_stereo',
                                             mark='mark', element='element', map='map')),
                     non_cgr=dict(edge=dict(s_bond='p_bond', s_stereo='p_stereo'),
                                  ext_node=('sp_neighbors', 'sp_hyb'),
                                  node=dict(s_charge='p_charge', s_neighbors='p_neighbors', s_hyb='p_hyb', s_x='p_x',
                                            s_y='p_y', s_z='p_z', s_stereo='p_stereo',
                                            mark='mark', element='element', map='map')))
