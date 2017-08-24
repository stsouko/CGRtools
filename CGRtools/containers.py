# -*- coding: utf-8 -*-
#
#  Copyright 2017 Ramil Nugmanov <stsouko@live.ru>
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
from collections import namedtuple
from itertools import chain
from networkx import Graph, relabel_nodes
from networkx.readwrite.json_graph import node_link_graph, node_link_data
from .strings import get_morgan, get_cgr_string, hash_cgr_string


CGRTemplate = namedtuple('CGRTemplate', ['substrats', 'products', 'meta'])
MatchContainer = namedtuple('MatchContainer', ['mapping', 'meta', 'patch'])


class MergedReaction(namedtuple('MergedReaction', ['substrats', 'products'])):
    def copy(self):
        return MergedReaction(self.substrats.copy(), self.products.copy())


class MoleculeContainer(Graph):
    def __init__(self, meta=None):
        super(MoleculeContainer, self).__init__()
        if isinstance(meta, dict):
            self.__meta = meta

    def pickle(self, compress=True):
        """ return json serializable CGR
        :stereo: save stereo data
        :compres: return smaller json for molecules. don't use for queries structures!
        """
        g = Graph()
        if compress and not self.get_center_atoms():
            s_only = True
            node_marks = self.__node_s
            edge_marks = ('s_bond',)
        else:
            s_only = False
            node_marks = self.__node_sp
            edge_marks = ('s_bond', 'p_bond')

        g.add_nodes_from((n, {k: v for k, v in a.items() if k in node_marks}) for n, a in self.nodes(data=True))
        g.add_edges_from((n, m, {k: v for k, v in a.items() if k in edge_marks}) for n, m, a in self.edges(data=True))

        data = node_link_data(g, attrs=self.__attrs)
        data.update(meta=self.meta, s_only=s_only,
                    stereo=[[a1, a2, x.get('s'), x.get('p')] for (a1, a2), x in self._stereo_dict.items()])
        return data

    @classmethod
    def unpickle(cls, data):
        """ convert json serializable CGR into MoleculeContainer object instance 
        """
        g = node_link_graph(data, attrs=cls.__attrs)
        g.__class__ = cls
        g.meta.update(data['meta'])
        for s in data['stereo']:
            g.add_stereo(*s)

        if data['s_only']:
            for _, a in g.nodes(data=True):
                a.update(cls.__attr_sp_clone(a, cls.__node_marks))
                a.update(p_x=a['s_x'], p_y=a['s_y'], p_z=a['s_z'])
            for *_, a in g.edges(data=True):
                a.update(cls.__attr_sp_clone(a, (('s_bond', 'p_bond', 'sp_bond'),)))
        else:
            g.fix_sp_marks()
        return g

    def copy(self):
        copy = super(MoleculeContainer, self).copy()
        copy.__class__ = self.__class__
        copy.meta.update(self.meta)
        for (a1, a2), x in self._stereo_dict.items():
            copy.add_stereo(a1, a2, x.get('s'), x.get('p'))
        return copy

    def subgraph(self, nbunch, copy=True, meta=False):
        sub = super(MoleculeContainer, self).subgraph(nbunch)
        sub.__class__ = self.__class__
        if copy:
            sub = sub.copy()
        if meta:
            sub.meta.update(self.meta)
        for (a1, a2), x in self._stereo_dict.items():
            sub.add_stereo(a1, a2, x.get('s'), x.get('p'))
        return sub

    def remap(self, mapping, copy=False):
        g = relabel_nodes(self, mapping, copy=copy)
        old_stereo = self._stereo_dict
        if not copy:
            self.__stereo_dict = {}
            self.__stereo = {}

        for (a1, a2), x in old_stereo.items():
            g.add_stereo(mapping[a1], mapping[a2], x.get('s'), x.get('p'))
        return g

    @property
    def meta(self):
        if self.__meta is None:
            self.__meta = {}
        return self.__meta

    @property
    def stereo(self):
        if self.__stereo is None:
            self.__stereo = {}
        return self.__stereo

    @property
    def _stereo_dict(self):
        if self.__stereo_dict is None:
            self.__stereo_dict = {}
        return self.__stereo_dict

    def add_stereo(self, atom1, atom2, s=None, p=None):
        val = {}
        if s:
            val['s'] = s
        if p:
            val['p'] = p
        if val and self.has_edge(atom1, atom2):
            atoms = (atom1, atom2)
            val['atoms'] = atoms
            self._stereo_dict[atoms] = val
            self.stereo.setdefault(atom1, {})[atom2] = val
            self.stereo.setdefault(atom2, {})[atom1] = val

    def get_fear_hash(self, weights=None, isotope=False, stereo=False, hyb=False, element=True):
        return hash_cgr_string(self.get_fear(weights=weights, isotope=isotope, stereo=stereo, hyb=hyb, element=element))

    def get_fear(self, weights=None, isotope=False, stereo=False, hyb=False, element=True):
        """
        :param weights: dict of atoms in keys and orders in values 
        :param isotope: set isotope marks
        :param stereo: set stereo marks
        :param hyb: set hybridization mark of atom
        :param element: set elements marks
        :return: string representation of CGR
        """
        return get_cgr_string(self, weights or get_morgan(self, isotope=isotope, element=element),
                              isotope=isotope, stereo=stereo, hyb=hyb, element=element)

    def get_center_atoms(self, stereo=False):
        """ get atoms of reaction center (dynamic bonds, stereo or charges).
        """
        # todo: stereo
        nodes = set()
        for n, node_attr in self.nodes(data=True):
            if node_attr.get('s_charge') != node_attr.get('p_charge'):
                nodes.add(n)

        for *n, node_attr in self.edges(data=True):
            if node_attr.get('s_bond') != node_attr.get('p_bond'):
                nodes.update(n)

        return list(nodes)

    def get_environment(self, atoms, dante=False, deep=0):
        """ get subgraph with atoms and their neighbors
        :param atoms: list of core atoms in graph
        :param dante: if True return list of graphs containing atoms, atoms + first circle, atoms + 1st + 2nd, 
        etc up to deep
        :param deep: number of bonds between atoms and neighbors. 
        """
        nodes = [atoms]
        for i in range(deep):
            nodes.append(set(chain.from_iterable(self.edges(nodes[i]))) or nodes[i])

        if dante:
            centers = [self.subgraph(a) for a in nodes]
        else:
            centers = self.subgraph(nodes[-1])

        return centers

    def get_morgan(self, isotope=False, element=True):
        return get_morgan(self, isotope=isotope, element=element)

    def fix_sp_marks(self, copy=False, nodes_bunch=None, edges_bunch=None):
        g = self.copy() if copy else self
        for n, a in g.nodes(data=True) if nodes_bunch is None else ((x, g.nodes[x]) for x in nodes_bunch):
            a.update(self.__attr_renew(a, self.__node_marks))
        for *_, a in g.edges(nbunch=edges_bunch, data=True):
            a.update(self.__attr_renew(a, (('s_bond', 'p_bond', 'sp_bond'),)))
        return g if copy else None

    def reset_query_marks(self, copy=False):
        """
        set or reset hyb and neighbors marks to atoms. 
        :param copy: if True return copy of graph and keep existing without marks
        :return: graph if copy True else None
        """
        g = self.copy() if copy else self
        for i in g.nodes():
            label = dict(s_hyb=1, p_hyb=1, sp_hyb=1, s_neighbors=0, p_neighbors=0, sp_neighbors=0)
            #  hyb 1- sp3; 2- sp2; 3- sp1; 4- aromatic
            for b, h, n in (('s_bond', 's_hyb', 's_neighbors'), ('p_bond', 'p_hyb', 'p_neighbors')):
                for node, bond in g[i].items():
                    b_type = bond.get(b)
                    if b_type and g.nodes[node]['element'] != 'H':
                        label[n] += 1
                    if b_type in (1, None) or label[h] in (3, 4):
                        continue
                    elif b_type == 4:
                        label[h] = 4
                    elif b_type == 3 or (b_type == 2 and label[h] == 2):  # Если есть 3-я или две 2-х связи, то sp1
                        label[h] = 3
                    elif b_type == 2:  # Если есть 2-я связь, но до этого не было найдено другой 2-й, 3-й, или аром.
                        label[h] = 2

            for n, m, h in (('s_hyb', 'p_hyb', 'sp_hyb'), ('s_neighbors', 'p_neighbors', 'sp_neighbors')):
                label[h] = (label[n], label[m]) if label[n] != label[m] else label[n]

            g.nodes[i].update(label)
        return g if copy else None

    @staticmethod
    def __attr_sp_clone(attr, marks):
        new_attr = {}
        for s, p, sp in marks:
            ls = attr.get(s)
            if ls is not None:
                new_attr[p] = new_attr[sp] = ls
        return new_attr

    @staticmethod
    def __attr_renew(attr, marks):
        new_attr = {}
        for s, p, sp in marks:
            ls = attr.get(s)
            lp = attr.get(p)

            if isinstance(ls, list):
                if isinstance(lp, list):
                    if ls == lp:
                        new_attr[sp] = new_attr[s] = new_attr[p] = list(set(ls))
                        continue
                    new_attr[sp] = list(set((x, y) for x, y in zip(ls, lp) if x != y))
                else:
                    new_attr[sp] = list(set((x, lp) for x in ls if x != lp))

                new_attr[s] = [x for x, _ in new_attr[sp]]
                new_attr[p] = [x for _, x in new_attr[sp]]
            elif isinstance(lp, list):
                new_attr[sp] = list(set((ls, x) for x in lp if x != ls))
                new_attr[s] = [x for x, _ in new_attr[sp]]
                new_attr[p] = [x for _, x in new_attr[sp]]
            elif ls != lp:
                new_attr[sp] = (ls, lp)
            elif ls is not None:
                new_attr[sp] = ls
        return new_attr

    __node_marks = [('s_%s' % mark, 'p_%s' % mark, 'sp_%s' % mark) for mark in ('neighbors', 'hyb', 'charge')]
    __tmp_node = ['element', 'isotope', 'mark', 's_x', 's_y', 's_z']
    __node_s = [x for x, *_ in __node_marks] + __tmp_node
    __node_sp = [y for x in __node_marks for y in x[:2]] + __tmp_node + ['p_x', 'p_y', 'p_z']
    __attrs = dict(source='atom1', target='atom2', name='atom', link='bonds')
    __stereo = __meta = __stereo_dict = None


class ReactionContainer(dict):
    def __init__(self, substrats=None, products=None, meta=None):
        super(ReactionContainer, self).__init__(substrats=substrats or [], products=products or [], meta=meta or {})

    def pickle(self):
        """ return json serializable reaction
        """
        return dict(substrats=[x.pickle() for x in self.substrats], meta=self.meta,
                    products=[x.pickle() for x in self.products])

    @staticmethod
    def unpickle(data):
        """ convert json serializable reaction into ReactionContainer object instance 
        """
        return ReactionContainer(substrats=[MoleculeContainer.unpickle(x) for x in data['substrats']],
                                 products=[MoleculeContainer.unpickle(x) for x in data['products']], meta=data['meta'])

    @property
    def substrats(self):
        return self['substrats']

    @property
    def products(self):
        return self['products']

    @property
    def meta(self):
        return self['meta']

    def copy(self):
        return ReactionContainer(substrats=[x.copy() for x in self['substrats']], meta=self['meta'].copy(),
                                 products=[x.copy() for x in self['products']])
