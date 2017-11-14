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
from itertools import cycle
from networkx import connected_components
from . import InvalidData
from .containers import MoleculeContainer, CGRContainer, ReactionContainer


class CGRcore(object):
    @staticmethod
    def split(m, meta=False):
        return [m.substructure(c, meta=meta) for c in connected_components(m)]

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

    @classmethod
    def compose(cls, m1, m2):
        """ remove from union graphs of products or reagents data about reagents or products
        """
        common = set(m1).intersection(m2)
        h = CGRContainer()
        unbalanced = {}

        """ remove bond, neighbors and hybridization states for common atoms.
        """
        for i, g in (('reagents', m1), ('products', m2)):
            is_cgr = isinstance(g, CGRContainer)
            pdi = cls.__popdict[i][is_cgr]
            ext_common = {}
            unbalanced[i] = ext_common
            e_pop, n_pop = pdi['edge'], pdi['node']

            for n, m, attr in g.edges(common, data=True):
                if n in common:
                    if m not in common:
                        ext_common.setdefault(m, []).append(n)
                elif m in common:
                    ext_common.setdefault(n, []).append(m)

                bond = {e_pop[k]: v for k, v in attr.items() if k in e_pop}
                if bond:
                    h.add_edge(n, m, **bond)

            unique = set(g).difference(common).difference(ext_common)

            for n, m, attr in g.edges(unique, data=True):
                if is_cgr:
                    tmp = attr
                else:
                    tmp = dict(p_bond=attr['s_bond'], **attr)
                    s = attr.get('s_stereo')
                    if s:
                        tmp['p_stereo'] = s
                h.add_edge(n, m, **tmp)

            for n in common:
                h.add_node(n, **{n_pop[k]: v for k, v in g.nodes[n].items() if k in n_pop})

            for n in unique:
                attr = g.nodes[n]
                h.add_node(n, **(attr if is_cgr else
                                 dict(p_charge=attr['s_charge'], p_x=attr['s_x'], p_y=attr['s_y'], p_z=attr['s_z'],
                                      p_neighbors=attr.get('s_neighbors'), p_hyb=attr.get('s_hyb'),
                                      p_radical=attr.get('s_radical'), p_stereo=attr.get('s_stereo'), **attr)))

            for n in ext_common:
                """ skip neighbors, hyb, stereo, charge, radical data from skin atoms
                """
                attr = g.nodes[n]
                if is_cgr:
                    tmp = {k: v for k, v in attr.items() if k not in pdi['ext_node']}
                elif i == 'reagents':
                    tmp = dict(p_x=attr['s_x'], p_y=attr['s_y'], p_z=attr['s_z'],
                               **{k: v for k, v in attr.items() if k in n_pop})
                else:
                    tmp = dict(s_x=attr['s_x'], s_y=attr['s_y'], s_z=attr['s_z'],
                               mark=attr['mark'], element=attr['element'], map=attr['map'],
                               **{n_pop[k]: v for k, v in attr.items() if k in n_pop})
                h.add_node(n, **tmp)

        """ calc unbalanced charges and radicals for skin atoms
        """
        reverse_ext = {}
        for i, e in unbalanced.items():
            for n, c in e.items():
                for j in c:
                    reverse_ext.setdefault(j, {}).setdefault(i, []).append(n)

        if reverse_ext:
            dnr = 0  # common atoms radicals or charges total changes
            for n in common:
                attr = h.nodes[n]
                dnr += cls.__radical_map[attr.get('p_radical')] - cls.__radical_map[attr.get('s_radical')]

        for n, sp in reverse_ext.items():
            atom = h.nodes[n]
            sh, ph = h.atom_implicit_h(n)
            dh = ph - sh
            dc = atom['p_charge'] - atom['s_charge']
            dr = cls.__radical_map[atom.get('p_radical')] - cls.__radical_map[atom.get('s_radical')]

            if not (dh or dc or dr):
                # common atom unchanged. Substitution, Elimination, Addition
                for m in sp.get('reagents', []):
                    attr = h.nodes[m]
                    attr.update(p_charge=attr['s_charge'], p_radical=attr.get('s_radical'))
                for m in sp.get('products', []):
                    attr = h.nodes[m]
                    attr.update(s_charge=attr['p_charge'], s_radical=attr.get('p_radical'))
            else:
                isp = dict(reagents=cycle(sp.get('reagents', [])), products=cycle(sp.get('products', [])))
                # radical balancing
                if dr > 0 and dnr > 0:
                    if dr <= dnr:
                        dnr -= dr
                    else:
                        dr = dnr
                        dnr = 0
                    # radical added or increased.
                    for x in range(dr):
                        m = next(isp['reagents'], None)  # homolysis
                        if m is not None:
                            attr = h.nodes[m]
                            attr['p_radical'] = \
                                cls.__radical_unmap[cls.__radical_map[attr.get('p_radical', attr.get('s_radical'))] + 1]
                        else:
                            m = next(isp['products'], None)  # radical addition
                            if m is not None:
                                attr = h.nodes[m]
                                attr['s_radical'] = \
                                    cls.__radical_unmap[cls.__radical_map[attr.get('s_radical',
                                                                                   attr.get('p_radical'))] + 1]

                elif dr < 0 and dnr < 0:
                    if dr >= dnr:
                        dnr -= dr
                    else:
                        dr = dnr
                        dnr = 0
                    # radical removed or decreased.
                    for x in range(-dr):
                        m = next(isp['products'], None)  # recombination
                        if m is not None:
                            attr = h.nodes[m]
                            attr['s_radical'] = \
                                cls.__radical_unmap[cls.__radical_map[attr.get('s_radical', attr.get('p_radical'))] + 1]
                        else:
                            m = next(isp['reagents'], None)  # radical elimination
                            if m is not None:
                                attr = h.nodes[m]
                                attr['p_radical'] = \
                                    cls.__radical_unmap[cls.__radical_map[attr.get('p_radical',
                                                                                   attr.get('s_radical'))] + 1]
                # protons and charge balancing
                if dh < 0 and dh < dc <= 0:  # deprotonation and charge decrease less.
                    for _, m in zip(range(dc - dh), isp['products']):  # electrophyle substitution
                        attr = h.nodes[m]
                        attr['s_charge'] = attr.get('s_charge', attr['p_charge']) + 1
                    for m in sp.get('reagents', []):
                        attr = h.nodes[m]
                        attr.update(s_charge=attr['p_charge'], s_radical=attr.get('p_radical'))

                elif dh > 0 and 0 <= dc < dh:  # charge increase and protonation less
                    for _, m in zip(range(dh - dc), isp['reagents']):  # cation elimination and anion protonation
                        attr = h.nodes[m]
                        attr['p_charge'] = attr.get('p_charge', attr['s_charge']) + 1
                    for m in sp.get('products', []):
                        attr = h.nodes[m]
                        attr.update(s_charge=attr['p_charge'], s_radical=attr.get('p_radical'))

                elif dc > 0 >= dh:  # charge increasing and deprotonation (if exists).
                    for x in range(dc - dh):
                        m = next(isp['reagents'], None)  # anion elimination
                        if m is not None:
                            attr = h.nodes[m]
                            attr['p_charge'] = attr.get('p_charge', attr['s_charge']) - 1
                        else:
                            m = next(isp['products'], None)  # cation addition
                            if m is not None:
                                attr = h.nodes[m]
                                attr['s_charge'] = attr.get('s_charge', attr['p_charge']) + 1

                elif dc < 0 <= dh:  # charge decreasing and protonation (if exists)
                    for _, m in zip(range(dh - dc), isp['reagents']):  # cation elimination
                        attr = h.nodes[m]
                        attr['p_charge'] = attr.get('p_charge', attr['s_charge']) + 1
                    for m in sp.get('products', []):
                        attr = h.nodes[m]
                        attr.update(s_charge=attr['p_charge'], s_radical=attr.get('p_radical'))

                else:
                    # restore charge and radical marks. we don't know what to do.
                    for m in sp.get('reagents', []):
                        attr = h.nodes[m]
                        attr.update(p_charge=attr['s_charge'], p_radical=attr.get('s_radical'))
                    for m in sp.get('products', []):
                        attr = h.nodes[m]
                        attr.update(s_charge=attr['p_charge'], s_radical=attr.get('p_radical'))

        """ update sp_* marks
        """
        h.fix_data()
        return h

    @classmethod
    def decompose(cls, g):
        tmp = ReactionContainer(meta=g.meta)

        x = g.copy()
        x.__class__ = MoleculeContainer
        for n, m in cls.__get_broken_paths(g, 's_bond'):
            x.remove_edge(n, m)

        for mol in cls.split(x):
            mol.fix_data()
            tmp.reagents.append(mol)

        x = g.copy()
        x.__class__ = MoleculeContainer
        for n, m in cls.__get_broken_paths(g, 'p_bond'):
            x.remove_edge(n, m)

        for mol in cls.split(x):
            for *_, edge_attr in mol.edges(data=True):
                edge_attr['s_bond'] = edge_attr.get('p_bond', None)
                s = edge_attr.get('p_stereo')
                if s:
                    edge_attr['s_stereo'] = s
                else:
                    edge_attr.pop('s_stereo', None)
            for _, node_attr in mol.nodes(data=True):
                for i, j in (('p_x', 's_x'), ('p_y', 's_y'), ('p_z', 's_z'), ('p_neighbors', 's_neighbors'),
                             ('p_hyb', 's_hyb'), ('p_charge', 's_charge'), ('p_radical', 's_radical')):
                    mark = node_attr.get(i)
                    if mark is not None:
                        node_attr[j] = mark
                    else:
                        node_attr.pop(j, None)

            mol.fix_data()
            tmp.products.append(mol)

        return tmp

    @staticmethod
    def __fix_attr(attr, marks):
        tmp = {}
        for s, p, sp in marks:
            if s in attr:
                tmp[p] = tmp[sp] = attr[s]
        return tmp

    @staticmethod
    def __get_broken_paths(g, edge):
        for m, n, attr in g.edges(data=True):
            if attr.get(edge) is None:
                yield n, m

    __tmp = dict(edge=dict(s_bond='s_bond', s_stereo='s_stereo'),
                 ext_node=('p_neighbors', 'p_hyb', 'p_stereo', 'p_charge', 'p_radical',
                           'sp_neighbors', 'sp_hyb', 'sp_stereo', 'sp_charge', 'sp_radical'),
                 node=dict(s_charge='s_charge', s_neighbors='s_neighbors', s_hyb='s_hyb',
                           s_x='s_x', s_y='s_y', s_z='s_z', s_stereo='s_stereo', s_radical='s_radical',
                           mark='mark', element='element', map='map'))

    __popdict = dict(reagents={True: __tmp, False: __tmp},
                     products={True: dict(edge=dict(p_bond='p_bond', p_stereo='p_stereo'),
                                          ext_node=('s_neighbors', 's_hyb', 's_stereo', 's_charge', 's_radical',
                                                    'sp_neighbors', 'sp_hyb', 'sp_stereo', 'sp_charge', 'sp_radical'),
                                          node=dict(p_charge='p_charge', p_neighbors='p_neighbors', p_hyb='p_hyb',
                                                    p_x='p_x', p_y='p_y', p_z='p_z', p_stereo='p_stereo',
                                                    p_radical='p_radical')),
                               False: dict(edge=dict(s_bond='p_bond', s_stereo='p_stereo'),
                                           node=dict(s_charge='p_charge', s_neighbors='p_neighbors', s_hyb='p_hyb',
                                                     s_x='p_x', s_y='p_y', s_z='p_z', s_stereo='p_stereo',
                                                     s_radical='p_radical'))})

    __radical_unmap = {None: None, 0: None, 1: 2, 2: 3}
    __radical_map = {1: 2, 2: 1, 3: 2, None: 0}
