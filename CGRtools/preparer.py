# -*- coding: utf-8 -*-
#
#  Copyright 2014-2018 Ramil Nugmanov <stsouko@live.ru>
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
from functools import reduce
from warnings import warn
from .containers import MoleculeContainer, MergedReaction
from .core import CGRcore
from .exceptions import InvalidData


class CGRpreparer(CGRcore):
    def __init__(self, cgr_type='0', extralabels=False, balance=False):
        """
        main object for CGR creation and manipulation

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
        :param balance: do electron balancing for atom unbalanced reactions
        """
        self.__cgr_type_code = cgr_type
        self.__extralabels = extralabels
        self.__balance = balance
        self.__cgr_type, self.__needed = self.__get_cgr_type(cgr_type)

    def __getstate__(self):
        return dict(cgr_type=self.__cgr_type_code, extralabels=self.__extralabels)

    def __setstate__(self, state):
        self.__init__(**state)

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

            g = self.compose(res.reagents, res.products, self.__balance)

        g.meta.update(data.meta)
        return g

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


__all__ = [CGRpreparer.__name__]
