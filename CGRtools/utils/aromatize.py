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


def dearomatize(g):

    def trace_route(atom, trace):
        data = {}
        double = False
        for i, attr in g[atom].items():
            if not trace or i != trace[-1]:
                data[i] = attr['s_bond']
            if attr['s_bond'] == 2:
                double = True

        for i, bond in data.items():
            if i in trace and len(trace) > 5:
                if i == trace[-5]:
                    pass
                elif i == trace[-6]:
                    pass

            elif bond == 4:
                g[atom][i]['s_bond'] = g[atom][i]['p_bond'] = 1 if double else 2
                trace_route(i, trace + [atom])

    min_node = min(g)
    trace_route(min_node, [])
