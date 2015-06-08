#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2014, 2015 Ramil Nugmanov <stsouko@live.ru>
# This file is part of FEAR (Find Errors in Automapped Reactions).
#
# FEAR is free software; you can redistribute it and/or modify
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
import copy
#import numpy

#from .mendel import replace
#from .maprules import rules


class FEAR:
    def __init__(self, **kwargs):
        self.__rulesC = kwargs['c_rules']
        self.__rulesE = kwargs['e_rules']
        self.__repare = kwargs['map_repair']
        
        self.chkmap = self.__chkmap if self.__rulesC or self.__rulesE else lambda x: x

    def __chkmap(self, data):
        print('#$#$', data)
        return None

    def __diff(self):
        diff = []
        for x, y in list(self.__data['diff'].items()):
            for w, z in list(y.items()):
                if (w + 1, x + 1, z) not in diff:
                    diff += [(x + 1, w + 1, z)]
        return diff

    def __reactcenter(self):
        atoms = []
        report = []
        newstruct = 0

        def getsmarts(sub, level):
            resout = []
            while sub:
                self.__nextnumb = self.__numb()
                res = self.__getSmiles(sub[:1], sub[0], level=level)
                resout.append(''.join(res[1]))
                sub = list(set(sub) - res[0])
            return resout

        rv = numpy.absolute(self.__data['products']['bondmatrix'] - self.__data['substrats']['bondmatrix'])
        for i in range(self.__length):
            if i not in atoms:
                sub = list(self.__reactcenteratoms([i], i))
                atoms.extend(sub)

                if len(sub) > 1:
                    srv = rv[numpy.ix_(sub, sub)].sum()
                    hashes = []
                    for j in [2]:#range(2, -1, -1):
                        mhash = self.__getMorgan(sub, level=j)
                        mhashkey = tuple(sorted(mhash.values()))
                        ruledataC = self.__rulesC.get(mhashkey, [0])
                        ruledataE = self.__rulesE.get(mhashkey, [0])
                        if ruledataC[0]:
                            hashes = [ruledataC[1]]
                            newhash = 0
                            break
                        elif ruledataE[0]:
                            hashes = [ruledataC[1]]
                            newhash = 1
                            break
                        else:
                            newhash = 1
                            self.__recsub = sub
                            self.__smartstarget = 'substrats'
                            self.__nextmap = self.__numb()
                            smirks = '.'.join(getsmarts(sub, j))
                            self.__smartstarget = 'products'
                            smirks += '>>' + '.'.join(getsmarts(sub, j))
                            hashes.append("new: %s'%.1f + %s" % (','.join(str(x) for x in mhashkey), srv, smirks))
                    report.extend(hashes)
                    newstruct += newhash
        return (True, report) if not newstruct else (False, report)

    __tosmiles = {1: '-', 2: '=', 3: '#', 1.5: ':'}

    def __numb(self):
        i = 1
        while True:
            yield i
            i += 1

    def __getSmiles(self, trace, inter, level=None):
        strace = set(trace)
        rule = {'substrats': lambda x: x < 80, 'products': lambda x: x % 10 != 8}
        iterlist = set([x for x, y in list(self.__data['diff'][inter].items()) if rule[self.__smartstarget](y)]).intersection(self.__recsub).difference(trace[-2:])
        #print iterlist, inter, self.__data['diff'][inter]
        if level == 1:
            symbol = '#%d' % self.__replace[self.__data[self.__smartstarget]['atomlist'][inter]['element']]
        elif level == 2: # с типами атомов
            symbol = '#%d%s' % (self.__replace[self.__data[self.__smartstarget]['atomlist'][inter]['element']], self.__types[self.__data[self.__smartstarget]['atomlist'][inter]['type']])
        else:
            symbol = '*'
        if self.__smartstarget == 'substrats':
            self.__mapsrepl[inter + 1] = next(self.__nextmap)
        smi = ['[%s:%d]' % (symbol, self.__mapsrepl[inter + 1])]
        #print smi
        concat = []
        stoplist = []
        for i in iterlist:
            if i in strace:
                if i not in stoplist:
                    # костыль для циклов. чтоб не было 2х проходов.
                    cyc = next(self.__nextnumb)
                    concat += [(i, cyc, inter)]
                    #print concat
                    smi[0] += '%s%d' % (self.__tosmiles[self.__data[self.__smartstarget]['bondmatrix'][inter][i]], cyc)
                    #print smi, '!'
                continue
            deep = self.__getSmiles(copy.copy(trace + [i]), i, level)
            strace.update(deep[0])
            #print strace
            #print deep
            if deep[2]:
                #print inter, i, '1'
                concat += deep[2]
                #print deep[2]
                for j in deep[2]:
                    if j[0] == inter:
                        #print '2', inter
                        stoplist += [j[2]]
                        #print stoplist
                        smi[0] += '%s%d' % (self.__tosmiles[self.__data[self.__smartstarget]['bondmatrix'][inter][i]], j[1])
            smi += ['(%s' % self.__tosmiles[self.__data[self.__smartstarget]['bondmatrix'][inter][i]]] + deep[1] + [')']

            #strace.update(self.__reactcenteratoms(copy.copy(trace + [i]), i))
        return strace, smi, concat

    def __reactcenteratoms(self, trace, inter):
        '''
        поиск всех атомов подструктур-реакционных центров
        '''
        strace = {inter}
        iterlist = set(self.__data['diff'][inter]).difference(trace)

        for i in iterlist:
            if i in strace:
                # костыль для циклов. чтоб не было 2х проходов.
                continue
            if self.__data['diff'][inter][i] > 10:
                strace.update(self.__reactcenteratoms(copy.copy(trace + [i]), i))
        return strace

    __types = {4: 'A', 3: 'T', 2: 'D', 1: 'S'}

    def __getMorgan(self, sub, level=None):
        '''
        модифицированный алгоритм моргана используется для хэширования подструктур.
        изначально все атомы имеют веса согласно сорту атомов.
        итеративно значения весов атомов умножатся на кратности образуемых связей с последующим суммированием полученных
        значений к текущему весу.

        после возвращается хэш подструктуры.
        '''
        if level == 1:
            weights = {x: self.__replace[self.__data['products']['atomlist'][x]['element']] * 1000 for x in sub}
        elif level == 2:
            weights = {x: self.__replace[self.__data['products']['atomlist'][x]['element']] * 1000 + self.__data['products']['atomlist'][x]['type'] * 1000000 for x in sub}
        else:
            weights = dict.fromkeys(sub, 0)
        for i in sub:
            l = sum([x for x in list(self.__data['diff'][i].values()) if x > 10]) + self.__data['diff'][i][i]
            weights[i] += l
        return weights


    def __accept(self, numb, length):
        '''
        функция остановки для алгоритма моргана. останавливает поиск после первого простоя улучшения.
        '''
        numl = len(set(numb.values()))
        if numl == length:
            return False
        if self.__buffer == numl:
            return False
        self.__buffer = numl
        return True
