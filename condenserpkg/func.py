# -*- coding: utf-8 -*-
#
#  func.py
#
#  Copyright 2013 nougmanoff <stsouko@live.ru>
#  This file is part of condenser.
#
#  condenser is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
from collections import defaultdict
import copy

import numpy


def main():
    print
    "модуль генерации конденсированных графов реакций."
    return 0


class Condenser(object):
    def __init__(self, charge, type):
        self.__matrix = {}
        self.__charge = charge
        self.__cgrtype = type

    def calc(self, data):
        self.__data = data
        self.__creatematrix()  #подготовка матриц
        self.__matrix['substrats']['meta'] = data['meta']
        self.__scan()
        return copy.deepcopy(self.__matrix)

    def __creatematrix(self):
        if self.__cgrtype == 1:
            self.__getmatrix(range(0, self.__data['substrats']), 'substrats')
            self.__matrix['products'] = self.__matrix['substrats']
        elif self.__cgrtype == 2:
            self.__getmatrix(range(self.__data['substrats'], self.__data['substrats'] + self.__data['products']),
                             'products')
            self.__matrix['substrats'] = self.__matrix['products']
        elif self.__cgrtype == 0:
            self.__getmatrix(range(0, self.__data['substrats']), 'substrats')
            self.__getmatrix(range(self.__data['substrats'], self.__data['substrats'] + self.__data['products']),
                             'products')
            self.__sortmatrix()
        elif 10 < self.__cgrtype < 20:
            self.__getmatrix(range(self.__cgrtype - 11, self.__cgrtype - 10), 'substrats')
            self.__matrix['products'] = self.__matrix['substrats']
        elif 20 < self.__cgrtype < 30:
            self.__getmatrix(range(self.__data['substrats'] + self.__cgrtype - 21,
                                   self.__data['substrats'] + self.__cgrtype - 20), 'products')
            self.__matrix['substrats'] = self.__matrix['products']
        elif -20 < self.__cgrtype < -10:
            self.__getmatrix(set(range(0, self.__data['substrats'])) - {-self.__cgrtype - 11}, 'substrats')
            self.__matrix['products'] = self.__matrix['substrats']
        elif -30 < self.__cgrtype < -20:
            self.__getmatrix(
                set(range(self.__data['substrats'], self.__data['substrats'] + self.__data['products'])) - {
                    self.__data['substrats'] - self.__cgrtype - 21}, 'products')
            self.__matrix['substrats'] = self.__matrix['products']

    def __getmatrix(self, mollist, category):
        atomlist = []
        #объединим атомы в 1 список. и создадим пустую матрицу.
        for i in mollist:
            atomlist += self.__data['molecules'][i]['atomlist']
        length = len(atomlist)
        self.__length = length
        tempMatrix = {'bondmatrix': numpy.zeros((length, length), dtype=int), 'maps': atomlist, 'diff': []}
        step = 0
        for i in mollist:
            matrixlength = len(self.__data['molecules'][i]['atomlist'])
            bonds = self.__data['molecules'][i]['bondmatrix']
            tempMatrix['bondmatrix'][step:step + matrixlength, step:step + matrixlength] = bonds
            step += matrixlength
        self.__matrix[category] = tempMatrix

    def __sortmatrix(self):

        maps = {'substrats': [int(x['map']) for x in self.__matrix['substrats']['maps']],
                'products': [int(x['map']) for x in self.__matrix['products']['maps']]}

        length = max((max(maps['products']), max(maps['substrats'])))
        self.__length = length

        for j in ('substrats', 'products'):
            tmp = numpy.zeros((length, length), dtype=int)

            shape = self.__matrix[j]['bondmatrix'].shape[0]
            for i in xrange(len(maps[j])):
                tmp[maps[j][i] - 1, 0:shape] = self.__matrix[j]['bondmatrix'][i, :]
            self.__matrix[j]['bondmatrix'] = numpy.zeros((length, length), dtype=int)
            for i in xrange(len(maps[j])):
                self.__matrix[j]['bondmatrix'][:, maps[j][i] - 1] = tmp[:, i]

            atomlist = [0] * length
            for i in self.__matrix[j]['maps']:
                atomlist[i['map']-1] = i
            self.__matrix[j]['maps'] = atomlist

        if 0 in self.__matrix['substrats']['maps']:
            lostlist = []
            for i in range(length):
                if self.__matrix['substrats']['maps'][i] == 0:
                    self.__matrix['substrats']['maps'][i] = self.__matrix['products']['maps'][i]
                    lostlist.append(i)
            self.__repmatrix(lostlist, 'products', 'substrats')
        if 0 in self.__matrix['products']['maps']:
            lostlist = []
            for i in range(length):
                if self.__matrix['products']['maps'][i] == 0:
                    self.__matrix['products']['maps'][i] = self.__matrix['substrats']['maps'][i]
                    lostlist.append(i)
            self.__repmatrix(lostlist, 'substrats', 'products')

    def __repmatrix(self, lost, s, t):
        for i in lost:
            for j in lost:
                self.__matrix[t]['bondmatrix'][i, j] = self.__matrix[s]['bondmatrix'][i, j]

    __cgr = {0: {0: 0, 1: 81, 2: 82, 3: 83, 4: 84, 5: 85, 6: 86, 7: 87},
             1: {0: 18, 1: 1, 2: 12, 3: 13, 4: 14, 5: 15, 6: 16, 7: 17},
             2: {0: 28, 1: 21, 2: 2, 3: 23, 4: 24, 5: 25, 6: 26, 7: 27},
             3: {0: 38, 1: 31, 2: 32, 3: 3, 4: 34, 5: 35, 6: 36, 7: 37},
             4: {0: 48, 1: 41, 2: 42, 3: 43, 4: 4, 5: 45, 6: 46, 7: 47},
             5: {0: 58, 1: 51, 2: 52, 3: 53, 4: 54, 5: 5, 6: 56, 7: 57},
             6: {0: 68, 1: 61, 2: 62, 3: 63, 4: 64, 5: 65, 6: 6, 7: 67},
             7: {0: 78, 1: 71, 2: 72, 3: 73, 4: 74, 5: 75, 6: 76, 7: 7}}

    def __chdiffpseudo(self, length):
        for i in xrange(length):
            s = self.__matrix['substrats']['maps'][i]['charge']
            p = self.__matrix['products']['maps'][i]['charge']
            self.__matrix['substrats']['maps'][i]['charge'] = str(self.__cgr[s][p])

    def __chdiffreal(self, length):
        for i in xrange(length):
            s = int(self.__matrix['substrats']['maps'][i]['element'])
            p = int(self.__matrix['products']['maps'][i]['element'])
            self.__matrix['substrats']['maps'][i]['element'] = str(p - s)

    def __chreplacereal(self, length):
        for i in xrange(length):
            self.__matrix['substrats']['maps'][i]['element'] = self.__matrix['products']['maps'][i]['element']

    def __chreplacepseudo(self, length):
        for i in xrange(length):
            self.__matrix['substrats']['maps'][i]['charge'] = self.__matrix['products']['maps'][i]['charge']

    def __scan(self):
        length = self.__length
        if self.__charge in (4, 6, 7):
            self.__chreplacereal(length)
        if self.__charge in (5, 6, 8):
            self.__chreplacepseudo(length)
        if self.__charge in (2, 3, 8):
            self.__chdiffreal(length)
        if self.__charge in (1, 3, 7):
            self.__chdiffpseudo(length)

        for i in xrange(length - 1):
            for j in xrange(i + 1, length):
                if self.__matrix['substrats']['bondmatrix'][i][j] != 0 or \
                                self.__matrix['products']['bondmatrix'][i][j] != 0:
                    self.__matrix['substrats']['diff'] += [(i + 1, j + 1,
                                                            self.__cgr[self.__matrix['substrats']['bondmatrix'][i][j]][
                                                                self.__matrix['products']['bondmatrix'][i][j]])]


if __name__ == '__main__':
    main()
