# -*- coding: utf-8 -*-
#
# Copyright 2014-2016 Ramil Nugmanov <stsouko@live.ru>
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
from .CGRcore import CGRcore
from .CGRreactor import CGRreactor, patcher


class CGRcombo(CGRcore):
    def __init__(self, cgr_type='0', extralabels=False, b_templates=None, m_templates=None, speed=False,
                 isotop=False, element=True, deep=0, stereo=False):
        CGRcore.__init__(self, cgr_type=cgr_type, extralabels=extralabels)
        self.__bal = CGRbalanser(b_templates, balanse_groups=True, speed=speed, stereo=stereo, isotop=isotop,
                                 extralabels=extralabels, element=element,
                                 deep=deep) if b_templates is not None else None

        self.__map = CGRbalanser(m_templates, balanse_groups=False, speed=speed, stereo=stereo, isotop=isotop,
                                 extralabels=extralabels, element=element,
                                 deep=deep) if m_templates is not None else None

    def getCGR(self, data):
        g = super(CGRcombo, self).getCGR(data)
        if self.__map is not None:
            g = self.__map.prepare(g)
        if self.__bal is not None:
            g = self.__bal.prepare(g)
        return g


class CGRbalanser(CGRreactor):
    def __init__(self, templates, speed=False, balanse_groups=True, stereo=False, extralabels=False,
                 isotop=False, element=True, deep=0):
        CGRreactor.__init__(self, stereo=stereo, hyb=extralabels, neighbors=extralabels,
                            isotop=isotop, element=element, deep=deep)

        self.__searcher = self.get_template_searcher(self.get_templates(templates), speed=speed, patch=True)
        self.__balanse_groups = balanse_groups

    def prepare(self, g):
        report = []
        if self.__balanse_groups:
            g = self.clonesubgraphs(g)

        while True:
            match = next(self.__searcher(g), None)
            if match:
                g = patcher(match)
                if 'CGR_TEMPLATE' in match['meta']:
                    report.append(match['meta']['CGR_TEMPLATE'])
            else:
                g.graph.setdefault('CGR_REPORT', []).extend(report)
                return g
