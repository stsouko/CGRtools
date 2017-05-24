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
from .containers import MoleculeContainer
from .core import CGRcore
from .reactor import CGRreactor, patcher


class CGRcombo(CGRcore):
    def __init__(self, cgr_type='0', extralabels=False, b_templates=None, m_templates=None,
                 isotope=False, element=True, stereo=False):
        CGRcore.__init__(self, cgr_type=cgr_type, extralabels=extralabels)
        if b_templates:
            self.__bal = CGRbalancer(b_templates, balance_groups=True, stereo=stereo, isotope=isotope,
                                     extralabels=extralabels, element=element)
        if m_templates:
            self.__map = CGRbalancer(m_templates, balance_groups=False, stereo=stereo, isotope=isotope,
                                     extralabels=extralabels, element=element)

        self.__pickle = dict(cgr_type=cgr_type, extralabels=extralabels, isotope=isotope,
                             element=element, stereo=stereo, b_templates=None, m_templates=None)
    __map = None
    __bal = None

    def pickle(self):
        """ remove attrs incorrectly dumped with dill
        """
        config = self.__pickle.copy()
        if self.__bal is not None:
            config['b_templates'] = self.__bal.pickle()['templates']
        if self.__map is not None:
            config['m_templates'] = self.__map.pickle()['templates']
        return config

    @staticmethod
    def unpickle(config):
        """ return CGRbalancer object instance
        """
        if {'cgr_type', 'stereo', 'extralabels', 'isotope', 'element', 'b_templates', 'm_templates'}.difference(config):
            raise Exception('Invalid config')
        config = config.copy()
        b_templates = [MoleculeContainer.unpickle(x) for x in config.pop('b_templates')]
        m_templates = [MoleculeContainer.unpickle(x) for x in config.pop('m_templates')]
        return CGRcombo(b_templates=b_templates, m_templates=m_templates, **config)

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

        self.__templates = templates
        self.__balance_groups = balance_groups

    __searcher = None

    def pickle(self):
        """ return config. for pickling
        """
        reactor = CGRreactor.pickle(self)
        return dict(templates=[x.pickle(compress=False) for x in self.__templates],
                    balance_groups=self.__balance_groups, stereo=reactor['stereo'],
                    extralabels=reactor['neighbors'], isotope=reactor['isotope'], element=reactor['element'])

    @staticmethod
    def unpickle(config):
        """ return CGRbalancer object instance
        """
        if {'templates', 'balance_groups', 'stereo', 'extralabels', 'isotope', 'element'}.difference(config):
            raise Exception('Invalid config')
        config = config.copy()
        templates = [MoleculeContainer.unpickle(x) for x in config.pop('templates')]
        return CGRbalancer(templates, **config)

    def prepare(self, g):
        if self.__searcher is None:
            self.__searcher = self.get_template_searcher(self.get_templates(self.__templates))

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
