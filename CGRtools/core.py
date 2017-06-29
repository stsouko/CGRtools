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
from .containers import MoleculeContainer


class CGRcore(object):
    @staticmethod
    def split(m, meta=False):
        return [m.subgraph(c, meta=meta) for c in connected_components(m)]

    @staticmethod
    def union(m1, m2):
        if set(m1) & set(m2):
            raise Exception('The node sets of m1 and m2 are not disjoint.')

        u = MoleculeContainer()
        u.add_nodes_from(m1.nodes(data=True))
        u.add_nodes_from(m2.nodes(data=True))
        u.add_edges_from(m1.edges(data=True))
        u.add_edges_from(m2.edges(data=True))

        for m in (m1, m2):
            for (a1, a2), x in m._stereo_dict.items():
                m.add_stereo(a1, a2, x.get('s'), x.get('p'))
        return u

    @classmethod
    def compose(cls, m1, m2):
        """ remove from union graphs of products or substrats data about substrats or products
        """
        common = set(m1).intersection(m2)
        extended_common = set()
        new_stereo = {}
        h = MoleculeContainer()

        """ remove bond, neighbors and hybridization states for common atoms.
        """
        for i, g in (('substrats', m1), ('products', m2)):
            ext_common = set()
            pop = cls.__popdict[i]['edge']
            s_pop = cls.__popdict[i]['stereo']
            for n, m, attr in g.edges(common, data=True):
                ext_common.update([n, m])
                bond = attr.get(pop)
                if bond:
                    h.add_edge(n, m, **{pop: bond})
                    if n in g.stereo and m in g.stereo[n]:
                        s = g.stereo[n][m]
                        stereo = s.get(s_pop)
                        if stereo:
                            new_stereo.setdefault(s['atoms'], {})[s_pop] = stereo

            uniq = set(g).difference(ext_common)
            for n, m, attr in g.edges(uniq, data=True):
                h.add_edge(n, m, **attr)
                if n in g.stereo and m in g.stereo[n]:
                    s = g.stereo[n][m]
                    new_stereo[s['atoms']] = s

            pop = cls.__popdict[i]['node']
            for n in common:
                h.add_node(n, **{k: v for k, v in g.node[n].items() if k in pop})

            pop = cls.__popdict[i]['ext_node']
            for n in ext_common.difference(common):
                h.add_node(n, **{k: v for k, v in g.node[n].items() if k not in pop})

            for n in uniq:
                h.add_node(n, **g.node[n])

            extended_common.update(ext_common)

        """ update sp_* marks
        """
        h.fix_sp_marks(nodes_bunch=extended_common, edges_bunch=common)

        for (a1, a2), x in new_stereo.items():
            h.add_stereo(a1, a2, x.get('s'), x.get('p'))
        return h

    __popdict = dict(products=dict(edge='p_bond', stereo='p',
                                   node=('p_charge', 'p_neighbors', 'p_hyb', 'p_x', 'p_y', 'p_z', 'mark', 'element'),
                                   ext_node=('s_neighbors', 's_hyb', 'sp_neighbors', 'sp_hyb')),
                     substrats=dict(edge='s_bond', stereo='s',
                                    node=('s_charge', 's_neighbors', 's_hyb', 's_x', 's_y', 's_z', 'mark', 'element'),
                                    ext_node=('p_neighbors', 'p_hyb', 'sp_neighbors', 'sp_hyb')))
