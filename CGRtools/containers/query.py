# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Lesser General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Lesser General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from collections import namedtuple
from itertools import product
from .common import BaseContainer
from ..algorithms import CGRstring
from ..exceptions import InvalidData
from ..periodictable import elements


DynamicContainer = namedtuple('DynamicContainer', ['reagent', 'product'])


class QueryContainer(BaseContainer):
    def pickle(self):
        """return json serializable CGR"""
        data = super().pickle()
        data['is_query'] = True
        data['s_only'] = False
        return data

    @classmethod
    def unpickle(cls, data):
        """convert json serializable CGR or Molecule into MoleculeContainer or CGRcontainer object instance"""
        if data['s_only']:
            raise InvalidData('pickled data is invalid query. try MoleculeContainer.unpickle')
        elif not data.get('is_query'):
            raise InvalidData('pickled data is invalid query. try CGRContainer.unpickle')

        graph, meta = super().unpickle(data)

        g = cls(graph, meta)
        g.fix_data()
        return g

    def add_atom(self, *args, **kwargs):
        pass

    def add_bond(self, *args, **kwargs):
        pass

    def add_stereo(self, *args, **kwargs):
        pass

    def _prepare_stereo(self):
        return {}

    def _signature_generator(self, *args, **kwargs):
        return CGRstring(*args, is_cgr=True, **kwargs)

    @staticmethod
    def _attr_renew(attr, marks):
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
                if ls is not None:
                    new_attr[s] = ls
                if lp is not None:
                    new_attr[p] = lp
                new_attr[sp] = (ls, lp)
            elif ls is not None:
                new_attr[sp] = new_attr[s] = new_attr[p] = ls
        return new_attr

    @classmethod
    def _atom_container(cls, attrs):
        el = [elements[x] for x in attrs['element']] if isinstance(attrs['element'], list) else \
             [elements[attrs['element']]]
        if isinstance(attrs['s_charge'], list):
            ch = list(zip(attrs['s_charge'], attrs['p_charge']))
        else:
            ch = [(attrs['s_charge'], attrs['p_charge'])]

        sr = attrs.get('s_radical')
        if isinstance(sr, list):
            rd = list(zip(sr, attrs['p_radical']))
        elif sr or attrs.get('p_radical'):
            rd = [(sr, attrs.get('p_radical'))]
        else:
            rd = [(None, None)]

        it = attrs.get('isotope')
        if not isinstance(it, list):
            it = [it]

        res, rep = [], []
        for e, (sc, pc), (sr, pr), i in product(el, ch, rd, it):
            res.append(e(charge=sc, multiplicity=sr, isotope=i))
            rep.append(e(charge=pc, multiplicity=pr, isotope=i))

        if len(res) == 1:
            res = res[0]
            rep = rep[0]

        return DynamicContainer(res, rep)

    @classmethod
    def _stereo_container(cls, attrs):
        return DynamicContainer(attrs.get('s_stereo'), attrs.get('p_stereo'))

    @classmethod
    def _bond_container(cls, attrs):
        return DynamicContainer(attrs.get('s_bond'), attrs.get('p_bond'))

    _node_marks = tuple(('s_%s' % mark, 'p_%s' % mark, 'sp_%s' % mark)
                        for mark in ('charge', 'stereo', 'radical', 'neighbors', 'hyb'))
    _edge_marks = tuple(('s_%s' % mark, 'p_%s' % mark, 'sp_%s' % mark) for mark in ('bond', 'stereo'))

    __tmp1 = ('element', 'isotope', 'mark')
    __tmp2 = tuple(y for x in _node_marks for y in x[:2])
    __tmp3 = __tmp1 + __tmp2

    _node_base = __tmp1 + ('s_x', 's_y', 's_z', 'p_x', 'p_y', 'p_z')
    _node_save = __tmp2 + _node_base
    _edge_save = tuple(y for x in _edge_marks for y in x[:2])
    _atom_marks = {x: x for x in __tmp3}
    _bond_marks = {x: x for x in _edge_save}
