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
from .CGRcore import CGRcore
from .CGRreactor import CGRreactor, patcher


class CGRcombo(CGRcore):
    def __init__(self, cgr_type='0', extralabels=False, b_templates=None, m_templates=None,
                 isotope=False, element=True, stereo=False):
        CGRcore.__init__(self, cgr_type=cgr_type, extralabels=extralabels)
        self.__bal = CGRbalancer(b_templates, balance_groups=True, stereo=stereo, isotope=isotope,
                                 extralabels=extralabels, element=element) if b_templates is not None else None

        self.__map = CGRbalancer(m_templates, balance_groups=False, stereo=stereo, isotope=isotope,
                                 extralabels=extralabels, element=element) if m_templates is not None else None

    def getCGR(self, data, is_merged=False):
        g = super(CGRcombo, self).getCGR(data, is_merged=is_merged)
        if self.__map is not None:
            g = self.__map.prepare(g)
        if self.__bal is not None:
            g = self.__bal.prepare(g)
        return g


class CGRbalancer(CGRreactor):
    def __init__(self, templates, balance_groups=True, stereo=False, extralabels=False, isotope=False, element=True):
        CGRreactor.__init__(self, stereo=stereo, hyb=extralabels, neighbors=extralabels,
                            isotope=isotope, element=element)

        self.__searcher = self.get_template_searcher(self.get_templates(templates))
        self.__balance_groups = balance_groups

    def prepare(self, g):
        report = []
        if self.__balance_groups:
            g = self.clone_subgraphs(g)

        while True:
            match = next(self.__searcher(g), None)
            if match:
                g = patcher(match)
                if 'CGR_TEMPLATE' in match['meta']:
                    report.append(match['meta']['CGR_TEMPLATE'])
            else:
                g.graph.setdefault('CGR_REPORT', []).extend(report)
                return g
