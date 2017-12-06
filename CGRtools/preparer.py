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
from warnings import warn
from . import InvalidConfig, InvalidData
from .containers import MoleculeContainer, ReactionContainer, MergedReaction
from .core import CGRcore
from .reactor import CGRreactor


class CGRpreparer(CGRcore):
    def __init__(self, cgr_type='0', extralabels=False, balance=False, templates=None,
                 isotope=False, element=True, stereo=False):
        if templates:
            self.__bal = CGRstandardizer(templates, balance_groups=balance, stereo=stereo, isotope=isotope,
                                         extralabels=extralabels, element=element)

        self.__init_common(cgr_type, extralabels, isotope, element, stereo)

    def _init_unpickle(self, cgr_type, extralabels, templates, isotope, element, stereo):
        if templates:
            self.__bal = CGRbalancer.unpickle(dict(templates=templates, balance_groups=True, stereo=stereo,
                                                   isotope=isotope, extralabels=extralabels, element=element))

        self.__init_common(cgr_type, extralabels, isotope, element, stereo)

    def __init_common(self, cgr_type, extralabels, isotope, element, stereo):
        self.__cgr_type, self.__needed = self.__get_cgr_type(cgr_type)
        self.__extralabels = extralabels
        self.__pickle = dict(cgr_type=cgr_type, extralabels=extralabels, isotope=isotope,
                             element=element, stereo=stereo, templates=None)

    def pickle(self):
        """ remove attrs incorrectly dumped with dill
        """
        config = self.__pickle.copy()
        if self.__bal is not None:
            tmp = self.__bal.pickle()
            config.update(templates=tmp['templates'], balance=tmp['balance'])

        return config

    @classmethod
    def unpickle(cls, config):
        """ return CGRbalancer object instance
        """
        args = {'cgr_type', 'stereo', 'extralabels', 'isotope', 'element', 'templates'}
        if not args.issubset(config):
            raise InvalidConfig('Invalid config')
        obj = cls.__new__(cls)  # Does not call __init__
        obj._init_unpickle(**{k: v for k, v in config.items() if k in args})
        return obj

    def condense(self, data):
        """
        condense reaction container to cgr molecule container.
        :param data: reaction container or merge_mols structure.
        :return: molecule container.

        if cgr_type in (1, 2, 3, 4, 5, 6) return reagents or products union else return cgr. see CLI help.
        """
        is_merged = isinstance(data, MergedReaction)
        if self.__cgr_type in (1, 2, 3, 4, 5, 6):
            if is_merged:
                raise InvalidData('invalid data')
            g = self.__reaction_splitter(data)
            if self.__extralabels:
                g.reset_query_marks()
        else:
            res = data if is_merged else self.merge_mols(data)
            if self.__extralabels:
                res.reagents.reset_query_marks()
                res.products.reset_query_marks()

            g = self.compose(res.reagents, res.products)

        if not is_merged:
            g.meta.update(data.meta)

        if self.__bal is not None:
            g = self.__bal.prepare(g)

        return g

    @classmethod
    def dissociate(cls, g):
        return cls.decompose(g)

    def merge_mols(self, data):
        if self.__cgr_type == 0:
            reagents = self.__union_all(data.reagents)
            products = self.__union_all(data.products)

        elif self.__cgr_type == 7:
            reagents = self.__union_all(self.__get_mols(data.reagents, self.__needed['reagents']))
            products = self.__union_all(self.__get_mols(data.products, self.__needed['products']))
        elif self.__cgr_type == 8:
            reagents = self.__union_all(self.__exc_mols(data.reagents, self.__needed['reagents']))
            products = self.__union_all(self.__exc_mols(data.products, self.__needed['products']))
        elif self.__cgr_type == 9:
            reagents = self.__union_all(self.__exc_mols(data.reagents, self.__needed['reagents']))
            products = self.__union_all(self.__get_mols(data.products, self.__needed['products']))
        elif self.__cgr_type == 10:
            reagents = self.__union_all(self.__get_mols(data.reagents, self.__needed['reagents']))
            products = self.__union_all(self.__exc_mols(data.products, self.__needed['products']))
        else:
            raise InvalidData('Merging need reagents and products')

        res = MergedReaction(reagents=reagents, products=products)
        return res

    @staticmethod
    def __get_cgr_type(_type):
        needed = [int(x) for x in _type.split(',')]
        if needed[0] == 0:
            t = 0  # CGR
        elif needed[0] == 1:
            t = 1  # all reagents
        elif needed[0] == 2:
            t = 2  # all products
        elif not any(True for x in needed if -200 < x < -100) and not any(True for x in needed if -300 < x < -200) and \
                any(True for x in needed if 100 < x < 200) and any(True for x in needed if 200 < x < 300):
            t = 7  # CGR on included parts of reagents and products
        elif any(True for x in needed if -200 < x < -100) and any(True for x in needed if -300 < x < -200):
            t = 8  # CGR on excluded parts of reagents and products
        elif any(True for x in needed if -200 < x < -100) and not any(True for x in needed if -300 < x < -200) and \
                any(True for x in needed if 200 < x < 300):
            t = 9  # CGR on excluded part of reagents and included part of products
        elif not any(True for x in needed if -200 < x < -100) and any(True for x in needed if -300 < x < -200) and \
                any(True for x in needed if 100 < x < 200):
            t = 10  # CGR on excluded part of products and included part of reagents
        elif 100 < needed[0] < 200:
            t = 3  # only included part of reagents
        elif 200 < needed[0] < 300:
            t = 4  # only included part of products
        elif -200 < needed[0] < -100:
            t = 5  # only excluded part of reagents
        elif -300 < needed[0] < -200:
            t = 6  # only excluded part of products
        else:
            t = 0

        ind = dict(reagents=sorted([abs(x) - 101 for x in needed if 100 < abs(x) < 200], reverse=True),
                   products=sorted([abs(x) - 201 for x in needed if 200 < abs(x) < 300], reverse=True)) \
            if t > 2 else None

        return t, ind

    def __reaction_splitter(self, data):
        if self.__cgr_type == 1:
            g = self.__union_all(data.reagents)
        elif self.__cgr_type == 2:
            g = self.__union_all(data.products)
        elif self.__cgr_type == 3:
            g = self.__union_all(self.__get_mols(data.reagents, self.__needed['reagents']))
        elif self.__cgr_type == 4:
            g = self.__union_all(self.__get_mols(data.products, self.__needed['products']))
        elif self.__cgr_type == 5:
            g = self.__union_all(self.__exc_mols(data.reagents, self.__needed['reagents']))
        elif self.__cgr_type == 6:
            g = self.__union_all(self.__exc_mols(data.products, self.__needed['products']))
        else:
            raise InvalidData('Splitter Error')
        return g

    @staticmethod
    def __get_mols(data, needed):
        mols = []
        for x in needed:
            try:
                mols.append(data[x])
            except IndexError:
                pass
        return mols

    @staticmethod
    def __exc_mols(data, needed):
        mols = data.copy()
        for x in needed:
            try:
                mols.pop(x)
            except IndexError:
                pass
        return mols

    @classmethod
    def __union_all(cls, data):
        return reduce(cls.union, data) if data else MoleculeContainer()

    __bal = None

    def getCGR(self, data):  # Reverse compatibility
        warn('getCGR name is deprecated. use condense instead', DeprecationWarning)
        return self.condense(data)

    def dissCGR(self, data):  # Reverse compatibility
        warn('dissCGR name is deprecated. use dissociate instead', DeprecationWarning)
        return self.dissociate(data)


class CGRstandardizer(CGRreactor):
    def __init__(self, templates, balance_groups=False, **kwargs):
        CGRreactor.__init__(self, **kwargs)

        self.__templates = templates
        self.__balance_groups = balance_groups
        self.__searcher = self.get_template_searcher(self.prepare_templates(templates))

    def pickle(self):
        """ return config. for pickling
        """
        reactor = CGRreactor.pickle(self)
        return dict(templates=[x.pickle() for x in self.__templates], balance_groups=self.__balance_groups, **reactor)

    @classmethod
    def unpickle(cls, config):
        """ return CGRbalancer object instance
        """
        args = {'templates', 'balance_groups'}
        if not args.issubset(config):
            raise InvalidConfig('Invalid config')
        config = config.copy()
        templates = [ReactionContainer.unpickle(x) for x in config.pop('templates')]
        return cls(templates, **config)

    def prepare(self, g, copy=False):
        if copy:
            g = g.copy()

        report = []
        if self.__balance_groups:
            g = self.clone_subgraphs(g)

        while True:
            searcher = self.__searcher(g)
            first_match = next(searcher, None)
            if not first_match:
                if report:
                    g.graph.setdefault('CGR_REPORT', []).extend(report)
                return g

            g = CGRreactor.patcher(g, first_match.patch)
            if 'CGR_TEMPLATE' in first_match.meta:
                report.append(first_match.meta['CGR_TEMPLATE'])

            for match in searcher:
                g = CGRreactor.patcher(g, match.patch)
                if 'CGR_TEMPLATE' in match.meta:
                    report.append(match.meta['CGR_TEMPLATE'])


class CGRcombo:  # Reverse compatibility
    def __new__(cls, cgr_type='0', extralabels=False, isotope=False, element=True, stereo=False,
                b_templates=None, m_templates=None):
        warn('CGRcombo deprecated. use CGRpreparer instead')
        if b_templates or m_templates:
            warn('b_templates and m_templates now merged to single list of patterns with kwarg: templates. '
                 'for stoichemist balancing use balance kwarg in CGRpreparer')

        return CGRpreparer(cgr_type=cgr_type, extralabels=extralabels, balance=False, templates=None,
                           isotope=isotope, element=element, stereo=stereo)

    @classmethod
    def unpickle(cls, *args, **kwargs):
        raise Exception('CGRcombo unpickle incompatible with CGRpreparer')


class CGRbalancer:  # Reverse compatibility
    def __new__(cls, *args, **kwargs):
        warn('CGRbalancer deprecated. use CGRstandardizer instead')
        return CGRstandardizer(*args, **kwargs)

    @classmethod
    def unpickle(cls, *args, **kwargs):
        warn('CGRbalancer deprecated. use CGRstandardizer instead')
        return CGRstandardizer.unpickle(*args, **kwargs)
