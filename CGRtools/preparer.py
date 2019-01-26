# -*- coding: utf-8 -*-
#
#  Copyright 2014-2019 Ramil Nugmanov <stsouko@live.ru>
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
from operator import or_
from .containers import MoleculeContainer, CGRContainer, ReactionContainer


class CGRpreparer:
    def __init__(self, cgr_type='0'):
        """
        CGR creation

        :param cgr_type: control condensation procedure:

            * 0 - CGR
            * 1 - reactants only
            * 2 - products only
            * 101-199 - reactant 1,2 or later
            * 201+ - product 1,2… (e.g. 202 - second product)
            * -101/-199 or -201/-299 - exclude reactant or product

            comma-separated list of selected or excluded reactants/products also supported:

            * 101,102 - only first and second molecules of reactants
            * -101,[-]103 <second [-] no sense> - exclude 1st and 3rd molecules of reactants

            also supported CGR on parts of reactants or/and products molecules. e.g. 101,102,-201 - CGR on only first and
            second reactants molecules with all products molecules excluding first
        """
        self.__cgr_type_code = cgr_type
        self.__cgr_type, self.__needed = self.__get_cgr_type(cgr_type)

    def compose(self, data):
        """
        condense reaction container to CGR. see init for details about cgr_type

        :param data: ReactionContainer
        :return: CGRContainer
        """
        g = self.__separate(data) if self.__cgr_type in (1, 2, 3, 4, 5, 6) else self.__condense(data)
        g.meta.update(data.meta)
        return g

    @staticmethod
    def decompose(data):
        if not isinstance(data, CGRContainer):
            raise TypeError('CGR only supported')
        r, p = ~data
        return ReactionContainer(r.split(), p.split(), meta=data.meta)

    @staticmethod
    def __get_cgr_type(_type):
        needed = [int(x) for x in _type.split(',')]
        if needed[0] == 0:
            t = 0  # CGR
        elif needed[0] == 1:
            t = 1  # all reactants
        elif needed[0] == 2:
            t = 2  # all products
        elif not any(True for x in needed if -200 < x < -100) and not any(True for x in needed if -300 < x < -200) and \
                any(True for x in needed if 100 < x < 200) and any(True for x in needed if 200 < x < 300):
            t = 7  # CGR on included parts of reactants and products
        elif any(True for x in needed if -200 < x < -100) and any(True for x in needed if -300 < x < -200):
            t = 8  # CGR on excluded parts of reactants and products
        elif any(True for x in needed if -200 < x < -100) and not any(True for x in needed if -300 < x < -200) and \
                any(True for x in needed if 200 < x < 300):
            t = 9  # CGR on excluded part of reactants and included part of products
        elif not any(True for x in needed if -200 < x < -100) and any(True for x in needed if -300 < x < -200) and \
                any(True for x in needed if 100 < x < 200):
            t = 10  # CGR on excluded part of products and included part of reactants
        elif 100 < needed[0] < 200:
            t = 3  # only included part of reactants
        elif 200 < needed[0] < 300:
            t = 4  # only included part of products
        elif -200 < needed[0] < -100:
            t = 5  # only excluded part of reactants
        elif -300 < needed[0] < -200:
            t = 6  # only excluded part of products
        else:
            t = 0

        ind = dict(reactants=sorted([abs(x) - 101 for x in needed if 100 < abs(x) < 200], reverse=True),
                   products=sorted([abs(x) - 201 for x in needed if 200 < abs(x) < 300], reverse=True)) \
            if t > 2 else None

        return t, ind

    def __condense(self, data):
        if self.__cgr_type == 0:
            reactants = self.__unite(data.reactants)
            products = self.__unite(data.products)
        elif self.__cgr_type == 7:
            reactants = self.__unite(self.__include(data.reactants, self.__needed['reactants']))
            products = self.__unite(self.__include(data.products, self.__needed['products']))
        elif self.__cgr_type == 8:
            reactants = self.__unite(self.__exclude(data.reactants, self.__needed['reactants']))
            products = self.__unite(self.__exclude(data.products, self.__needed['products']))
        elif self.__cgr_type == 9:
            reactants = self.__unite(self.__exclude(data.reactants, self.__needed['reactants']))
            products = self.__unite(self.__include(data.products, self.__needed['products']))
        else:  # 10
            reactants = self.__unite(self.__include(data.reactants, self.__needed['reactants']))
            products = self.__unite(self.__exclude(data.products, self.__needed['products']))

        return reactants ^ products

    def __separate(self, data):
        if self.__cgr_type == 1:
            g = self.__unite(data.reactants)
        elif self.__cgr_type == 2:
            g = self.__unite(data.products)
        elif self.__cgr_type == 3:
            g = self.__unite(self.__include(data.reactants, self.__needed['reactants']))
        elif self.__cgr_type == 4:
            g = self.__unite(self.__include(data.products, self.__needed['products']))
        elif self.__cgr_type == 5:
            g = self.__unite(self.__exclude(data.reactants, self.__needed['reactants']))
        else:  # 6
            g = self.__unite(self.__exclude(data.products, self.__needed['products']))
        return g

    @staticmethod
    def __include(data, needed):
        mols = []
        for x in needed:
            try:
                mols.append(data[x])
            except IndexError:
                pass
        return mols

    @staticmethod
    def __exclude(data, needed):
        mols = data.copy()
        for x in needed:
            try:
                mols.pop(x)
            except IndexError:
                pass
        return mols

    @staticmethod
    def __unite(data):
        return reduce(or_, data) if data else MoleculeContainer()


__all__ = ['CGRpreparer']
