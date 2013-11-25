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
    def __init__(self):
        self.__matrix = {}

    def calc(self, data):
        self.__creatematrix(data) #подготовка матриц
        self.__matrix['substrats']['meta'] = data['meta']
        self.__scan()
        return copy.deepcopy(self.__matrix)

    def __creatematrix(self, data):
        for j in [(0, data['substrats'], 'substrats'),
                  (data['substrats'], data['substrats'] + data['products'], 'products')]:
            atomlist = []
            for i in range(j[0], j[1]): #объединим атомы в 1 список. и создадим пустую матрицу.
                atomlist += data['molecules'][i]['atomlist']
            length = len(atomlist)
            tempMatrix = {'bondmatrix': numpy.zeros((length, length), dtype=int), 'maps': atomlist, 'diff': []}
            step = 0
            for i in range(j[0], j[1]):
                matrixlength = len(data['molecules'][i]['atomlist'])
                bonds = data['molecules'][i]['bondmatrix']
                tempMatrix['bondmatrix'][step:step + matrixlength, step:step + matrixlength] = bonds
                step += matrixlength
            self.__matrix[j[2]] = tempMatrix
        self.__sortmatrix()

    def __sortmatrix(self):
        length = len(self.__matrix['substrats']['maps'])
        for j in ['substrats', 'products']:
            tmp = numpy.empty((length, length), dtype=int)
            for i in range(length):
                tmp[int(self.__matrix[j]['maps'][i]['map']) - 1, :] = self.__matrix[j]['bondmatrix'][i, :]
            for i in range(length):
                self.__matrix[j]['bondmatrix'][:, int(self.__matrix[j]['maps'][i]['map']) - 1] = tmp[:, i]
            self.__matrix[j]['maps'].sort(key=lambda a: int(a['map']))

    __cgr = {0: {0: 0, 1: 81, 2: 82, 3: 83, 4: 84},
             1: {0: 18, 1: 1, 2: 12, 3: 13, 4: 14},
             2: {0: 28, 1: 21, 2: 2, 3: 23, 4: 24},
             3: {0: 38, 1: 31, 2: 32, 3: 3, 4: 34},
             4: {0: 48, 1: 41, 2: 42, 3: 43, 4: 4}}

    def __scan(self):
        length = len(self.__matrix['substrats']['maps'])
        for i in xrange(length - 1):
            for j in xrange(i + 1, length):
                if self.__matrix['substrats']['bondmatrix'][i][j] != 0 or \
                                self.__matrix['products']['bondmatrix'][i][j] != 0:
                    self.__matrix['substrats']['diff'] += [(i + 1, j + 1,
                                                            self.__cgr[self.__matrix['substrats']['bondmatrix'][i][j]][
                                                                self.__matrix['products']['bondmatrix'][i][j]])]


if __name__ == '__main__':
    main()
