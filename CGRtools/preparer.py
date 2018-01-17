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
from .containers import MoleculeContainer, ReactionContainer, MergedReaction
from .core import CGRcore
from .exceptions import InvalidConfig, InvalidData
from .reactor import CGRreactor


class CGRpreparer(CGRcore):
    """main object for CGR creation and manipulation"""
    def __init__(self, cgr_type='0', extralabels=False, templates=None, **kwargs):
        """
        :param cgr_type: control condensation procedure:

            * 0 - CGR
            * 1 - reagents only
            * 2 - products only
            * 101-199 - reagent 1,2 or later
            * 201+ - product 1,2… (e.g. 202 - second product)
            * -101/-199 or -201/-299 - exclude reagent or product

            comma-separated list of selected or excluded reagents/products also supported:

            * 101,102 - only first and second molecules of reagents
            * -101,[-]103 <second [-] no sense> - exclude 1st and 3rd molecules of reagents

            also supported CGR on parts of reagents or/and products molecules. e.g. 101,102,-201 - CGR on only first and
            second reagents molecules with all products molecules excluding first
        :param extralabels: [re]create query labels in Graphs before condensation. also see CGRreactor init and
            reset_query_marks in MoleculeContainer
        :param templates: see CGRstandardizer init. if not None extralabels automatically will be True
        :param kwargs: see CGRstandardizer init
        """
        if templates is not None:
            extralabels = True
            self.__bal = CGRstandardizer(templates, extralabels=extralabels, **kwargs)

        self.__init_common(cgr_type, extralabels)

    def _init_unpickle(self, cgr_type, extralabels, templates, **kwargs):
        if templates is not None:
            self.__bal = CGRstandardizer.unpickle(dict(templates=templates, extralabels=True, **kwargs))
        self.__init_common(cgr_type, extralabels)

    def __init_common(self, cgr_type, extralabels):
        self.__cgr_type, self.__needed = self.__get_cgr_type(cgr_type)
        self.__extralabels = extralabels
        self.__pickle = dict(cgr_type=cgr_type, extralabels=extralabels, templates=None)

    def pickle(self):
        """return json serializable config of object"""
        config = self.__pickle.copy()
        if self.__bal is not None:
            config.update(self.__bal.pickle())
        return config

    @classmethod
    def unpickle(cls, config):
        """restore CGRbalancer object instance from json-config"""
        args = {'cgr_type', 'extralabels', 'templates'}
        if not args.issubset(config):
            raise InvalidConfig('Invalid config')
        obj = cls.__new__(cls)
        obj._init_unpickle(**config)
        return obj

    def condense(self, data):
        """
        condense reaction container to CGR. see init for details about cgr_type

        :param data: ReactionContainer or MergedReaction object
        :return: CGRContainer
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

        res = MergedReaction(reagents=reagents, products=products, meta=data.meta)
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

    __bal = __visible = None

    def getCGR(self, data):  # Reverse compatibility
        """deprecated. see condense"""
        warn('getCGR name is deprecated. use condense instead', DeprecationWarning)
        return self.condense(data)

    def dissCGR(self, data):  # Reverse compatibility
        """deprecated. see dissociate"""
        warn('dissCGR name is deprecated. use dissociate instead', DeprecationWarning)
        return self.dissociate(data)


class CGRstandardizer(CGRreactor):
    """CGR standardization and reaction balancing"""
    def __init__(self, templates, balance_groups=False, **kwargs):
        """
        :param templates: CGRTemplates. rules for graph modifications. possible be False
        :param balance_groups: if True: for unbalanced reactions contains multiple attached functional groups in
            products and one of them described in reagents - will be restored information about all equal groups.
            for example:

                R + B1-X-> B'1-R'-B'2 + X'

            where B' is transformed B, R and X same.
            we know what B'1 and B'2 is equal and B'1 is transformed B1 =>
            this groups most likely appeared from a single reagent. we can add copy of B-X to reagents.
            results will be:

                R + B1-X1 + B2-X2 -> B'1-R'-B'2 + X'1 + X'2

        :param kwargs: see CGRreactor init
        """
        super().__init__(**kwargs)

        self.__templates = templates
        self.__balance_groups = balance_groups
        if templates:
            self.__searcher = self.get_template_searcher(self.prepare_templates(templates))
            self.__searching = True

    def pickle(self):
        """ return config. for pickling
        """
        reactor = super().pickle()
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

        while self.__searching:
            searcher = self.__searcher(g)
            first_match = next(searcher, None)
            if not first_match:
                if report:
                    g.graph.setdefault('CGR_REPORT', []).extend(report)
                break

            g = CGRreactor.patcher(g, first_match.patch)
            if 'CGR_TEMPLATE' in first_match.meta:
                report.append(first_match.meta['CGR_TEMPLATE'])

            for match in searcher:
                g = CGRreactor.patcher(g, match.patch)
                if 'CGR_TEMPLATE' in match.meta:
                    report.append(match.meta['CGR_TEMPLATE'])
        return g

    __searcher = None
    __searching = False


class CGRcombo:  # Reverse compatibility
    def __new__(cls, cgr_type='0', extralabels=False, isotope=False, element=True, stereo=False,
                b_templates=None, m_templates=None):
        warn('CGRcombo deprecated and can be work incorrectly. use CGRpreparer instead')
        if b_templates or m_templates:
            warn('b_templates and m_templates now merged to single list of patterns with kwarg: templates. '
                 'for stoichemist balancing use balance_groups kwarg in CGRpreparer')

        return CGRpreparer(cgr_type=cgr_type, extralabels=extralabels, balance_groups=False, templates=None,
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


__all__ = [CGRpreparer.__name__, CGRstandardizer.__name__]
