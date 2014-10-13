# -*- coding: utf-8 -*-
#
#  func.py
#
#  Copyright 2014 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of condenser.
#
#  condenser is free software; you can redistribute it and/or modify
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

from collections import defaultdict
import copy
from weightable import replace
import numpy


def main():
    print "This file is part of condenser."
    return 0


class Condenser(object):
    def __init__(self, charge, type, repare):
        self.__matrix = {}
        self.__charge = charge
        self.__cgrtype = type
        self.__repare = repare
        self.__replace = replace()

    def calc(self, data):
        self.__moltorepare = False
        self.__data = data
        self.__creatematrix()  #подготовка матриц
        self.__matrix['substrats']['meta'] = data['meta']
        self.__scan()
        return copy.deepcopy(self.__matrix)

    def __creatematrix(self):
        if self.__cgrtype == 1:
            self.__getmatrix(self.__data['substrats'], 'substrats')
            self.__matrix['products'] = self.__matrix['substrats']
        elif self.__cgrtype == 2:
            self.__getmatrix(self.__data['products'], 'products')
            self.__matrix['substrats'] = self.__matrix['products']
        elif self.__cgrtype == 0:
            self.__getmatrix(self.__data['substrats'], 'substrats')
            self.__getmatrix(self.__data['products'], 'products')
            self.__sortmatrix()
        elif 10 < self.__cgrtype < 20:
            try:
                mols = self.__data['substrats'][self.__cgrtype - 11]
            except:
                mols = self.__data['substrats'][-1]
            self.__getmatrix(mols, 'substrats')
            self.__matrix['products'] = self.__matrix['substrats']
        elif 20 < self.__cgrtype < 30:
            try:
                mols = self.__data['products'][self.__cgrtype - 21]
            except:
                mols = self.__data['products'][-1]
            self.__getmatrix(mols, 'products')
            self.__matrix['substrats'] = self.__matrix['products']
        elif -20 < self.__cgrtype < -10:
            try:
                mols = self.__data['substrats'][:-self.__cgrtype - 11] + self.__data['substrats'][-self.__cgrtype - 10:]
            except:
                mols = self.__data['substrats'][:-1]
            self.__getmatrix(mols, 'substrats')
            self.__matrix['products'] = self.__matrix['substrats']
        elif -30 < self.__cgrtype < -20:
            try:
                mols = self.__data['products'][:-self.__cgrtype - 21] + self.__data['products'][-self.__cgrtype - 20:]
            except:
                mols = self.__data['products'][:-1]
            self.__getmatrix(mols, 'products')
            self.__matrix['substrats'] = self.__matrix['products']

    def __getmatrix(self, data, category):
        '''
        постройка матриц связности
        '''
        atomlist = []
        #объединим атомы в 1 список. и создадим пустую матрицу.
        for i in data:
            #self.__mols[role].append(range(len(atomlist), len(atomlist) + len(i['atomlist'])))
            atomlist += i['atomlist']
        length = len(atomlist)
        self.__length = length
        tempMatrix = {'bondmatrix': numpy.zeros((length, length), dtype=float), 'maps': atomlist, 'diff': []}
        step = 0
        for i in data:
            matrixlength = len(i['atomlist'])
            tempMatrix['bondmatrix'][step:step + matrixlength, step:step + matrixlength] = i['bondmatrix']
            step += matrixlength

        self.__matrix[category] = tempMatrix

    def __sortmatrix(self):
        '''
        атомы отсортировываются по порядку. для атомов представленных только в продуктах или реагентах достраивается
        симметричная матрица с разрывом в месте стыка с основной молекулой.
        '''
        maps = {'substrats': [int(x['map']) for x in self.__matrix['substrats']['maps']],
                'products': [int(x['map']) for x in self.__matrix['products']['maps']]}

        length = max((max(maps['products']), max(maps['substrats'])))
        lose = sorted(list(set(range(1, length + 1)).difference(maps['products']).difference(maps['substrats'])))[::-1]
        if lose:
            for k in maps.keys():
                for i in lose:
                    newmaps = []
                    for j in maps[k]:
                        newmaps.append(j if j < i else j - 1)
                    maps[k] = newmaps
            self.__length = length - len(lose)
        else:
            self.__length = length

        for j in ('substrats', 'products'):
            tmp = numpy.zeros((self.__length, self.__length), dtype=float)

            shape = self.__matrix[j]['bondmatrix'].shape[0]
            for i in xrange(len(maps[j])):
                tmp[maps[j][i] - 1, 0:shape] = self.__matrix[j]['bondmatrix'][i, :]
            self.__matrix[j]['bondmatrix'] = numpy.zeros((self.__length, self.__length), dtype=float)
            for i in xrange(len(maps[j])):
                self.__matrix[j]['bondmatrix'][:, maps[j][i] - 1] = tmp[:, i]

            atomlist = [0] * self.__length
            for i, k in zip(self.__matrix[j]['maps'], maps[j]):
                i['map'] = k
                atomlist[k - 1] = i
            self.__matrix[j]['maps'] = atomlist

        if 0 in self.__matrix['substrats']['maps']:
            '''
            метод восстановления используется только для недостающих атомов реагентов.
            информация о факте отвалившегося фрагмента в продуктах самодостаточна для корректного ААО или полного
            восстановления структуры уходящей группы.
            '''
            self.__moltorepare = True
            self.__searchlost('products', 'substrats')
        if 0 in self.__matrix['products']['maps']:
            self.__moltocharge = True
            self.__searchlost('substrats', 'products')

    def __searchlost(self, s, t):
        '''
        внедрение в матрицу продуктов/субстратов атомов из субстратов/продуктов
        '''
        lostlist = []
        for i in xrange(self.__length):
            if self.__matrix[t]['maps'][i] == 0:
                self.__matrix[t]['maps'][i] = self.__matrix[s]['maps'][i]
                #print s, self.__matrix[s]['maps'][i]
                lostlist.append(i)
        for i in lostlist:
            for j in lostlist:
                self.__matrix[t]['bondmatrix'][i, j] = self.__matrix[s]['bondmatrix'][i, j]

    __cgr = {0: {0: 0, 1: 81, 2: 82, 3: 83, 4: 84, 1.5: 84, 5: 85, 6: 86, 7: 87},
             1: {0: 18, 1: 1, 2: 12, 3: 13, 4: 14, 1.5: 14, 5: 15, 6: 16, 7: 17},
             2: {0: 28, 1: 21, 2: 2, 3: 23, 4: 24, 1.5: 24, 5: 25, 6: 26, 7: 27},
             3: {0: 38, 1: 31, 2: 32, 3: 3, 4: 34, 1.5: 34, 5: 35, 6: 36, 7: 37},
             4: {0: 48, 1: 41, 2: 42, 3: 43, 4: 4, 1.5: 4, 5: 45, 6: 46, 7: 47},
             1.5: {0: 48, 1: 41, 2: 42, 3: 43, 4: 4, 1.5: 4, 5: 45, 6: 46, 7: 47},
             5: {0: 58, 1: 51, 2: 52, 3: 53, 4: 54, 1.5: 54, 5: 5, 6: 56, 7: 57},
             6: {0: 68, 1: 61, 2: 62, 3: 63, 4: 64, 1.5: 64, 5: 65, 6: 6, 7: 67},
             7: {0: 78, 1: 71, 2: 72, 3: 73, 4: 74, 1.5: 74, 5: 75, 6: 76, 7: 7}}

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

        connections = {x: {} for x in range(length)}

        for i in xrange(length - 1):
            for j in xrange(i + 1, length):
                if self.__matrix['substrats']['bondmatrix'][i][j] != 0 or \
                                self.__matrix['products']['bondmatrix'][i][j] != 0:
                    diff = self.__cgr[self.__matrix['substrats']['bondmatrix'][i][j]][
                        self.__matrix['products']['bondmatrix'][i][j]]
                    self.__matrix['substrats']['diff'] += [(i + 1, j + 1, diff)]
                    connections[i][j] = connections[j][i] = diff

        self.__connections = copy.deepcopy(connections)

        if self.__repare and self.__moltorepare:
            groups = defaultdict(list)
            backup = copy.deepcopy(self.__matrix)
            atomlist = set(range(length))
            self.__smiless = {}
            while atomlist:
                atom = list(atomlist)[0]
                block = self.__subsearch([atom], atom)
                atomlist.difference_update(block)
                smiles = self.__getSmiles(block)
                #print block, smiles
                groups[tuple(sorted(smiles.values()))].append(block)
                self.__smiless[tuple(block)] = smiles
            #print groups
            try:
                self.__balancer(groups)
            except:
                self.__connections = connections
                self.__matrix = backup
            diff = []
            for x, y in self.__connections.items():
                for w, z in y.items():
                    if (w + 1, x + 1, z) not in diff:
                        diff += [(x + 1, w + 1, z)]
            self.__matrix['substrats']['diff'] = diff

        #место для зарядового балансировщика!!!
        if self.__repare and self.__moltocharge:
            print self.__connections
            self.__chargebalancer()

        if self.__charge in (4, 6, 7):
            self.__chreplacereal(length)
        if self.__charge in (5, 6, 8):
            self.__chreplacepseudo(length)
        if self.__charge in (2, 3, 8):
            self.__chdiffreal(length)
        if self.__charge in (1, 3, 7):
            self.__chdiffpseudo(length)

    def __chargebalancer(self):
        #print self.__matrix, self.__connections
        atoms = []
        for i in xrange(self.__length):
            if i not in atoms:
                sub = list(self.__reactcenteratoms([i], i))
                atoms.extend(sub)
                if len(sub) > 1:
                    #локализовали реакционный центр. осталось раскидать электроны.
                    print '%', sub

                    for x in sub:
                        print x, self.__matrix['substrats']['maps'][x]['charge'], self.__matrix['products']['maps'][x]['charge']
                    pass

    def __reactcenteratoms(self, trace, inter):
        '''
        поиск всех атомов подструктур-реакционных центров
        '''
        strace = {inter}
        iterlist = set(self.__connections[inter]).difference(trace)

        for i in iterlist:
            if i in strace:
                # костыль для циклов. чтоб не было 2х проходов.
                continue
            if self.__connections[inter][i] > 10:
                strace.update(self.__reactcenteratoms(copy.copy(trace + [i]), i))
        return strace

    def __balancer(self, groups):
        '''
        для всех подструктур с идентичными хэшами находится фрагмент, информация о котором имеется в реагентах.
        после выбирается атом входа для обхода по подструктурам для выявления динамических связей.
        '''
        for i, j in groups.items():
            target = None
            if len(j) > 1:
                for fragment in j:
                    for substrat in self.__data['substrats']:
                        #print fragment, set([int(x['map']) - 1 for x in self.__data['molecules'][substrat]['atomlist']]).intersection(fragment)
                        if set([int(x['map']) - 1 for x in substrat['atomlist']]).intersection(fragment):
                            target = fragment
                            #print '!', target
                            self.__target = self.__smiless[tuple(target)]
                            break
                #target not found problem
                if not target:
                    #print "@"
                    continue
                #print target
                for fragment in j:
                    if fragment == target:
                        continue
                    fstart = min(self.__smiless[tuple(fragment)], key=self.__smiless[tuple(fragment)].get)
                    tstart = min(self.__smiless[tuple(target)], key=self.__smiless[tuple(target)].get)
                    #print fstart, tstart
                    self.__fragment = self.__smiless[tuple(fragment)]
                    #print self.__fragment
                    self.__builder(tstart, fstart, [])

    def __builder(self, target, fragment, trace):
        '''
        рекурсивная функция для прохода по фрагментам шаблона и восстанавливаемого для копирования связей.
        но есть нюанс. связь концов циклических структур не копируется.
        но поскольку это далекая часть от реакционного центра, то монопениссуально.
        '''
        #print '_+', target, fragment
        #self.__fragment.pop(fragment)
        titerlist = set(self.__connections[target]) #.difference(trace)
        # if len(titerlist) == 1 and self.__connections[target][list(titerlist)[0]] % 10 == 8:
        #     print '*'
        #     return 0
        #print titerlist
        for i in titerlist:
            #print "!", i
            if i in trace:
                continue
            j = self.__connections[target][i]
            #print '^', j, i
            if j not in (81, 82, 83, 84, 85, 86, 87):
                if self.__target.get(i):
                    #print self.__target[i], next(key for key, val in self.__fragment.items() if val == self.__target[i])
                    fi = next((key for key in self.__connections[fragment] if self.__fragment.get(key) == self.__target[i] and key not in trace), next((key for key, val in self.__fragment.items() if val == self.__target[i]), None))
                    #print "#", target, i, fragment, fi, j
                    self.__connections[fragment][fi] = self.__connections[fi][fragment] = j
                    trace += self.__builder(i, fi, trace + [target, fragment])
                else:
                    #запилятор уходящей группы
                    #print "============================="
                    self.__nucleofug(target, fragment, i, {})
        return trace

    def __nucleofug(self, target, fragment, vector, trace):
        '''
        полностью копирует информацию об уходящей группе из известной. копируется состояние в реагентах.
        '''
        mapn = len(self.__matrix['substrats']['maps'])
        #print "^", target, vector, fragment, self.__connections[target][vector], mapn

        data = copy.copy(self.__matrix['substrats']['maps'][vector])
        data['map'] = mapn + 1
        self.__matrix['substrats']['maps'] += [data]
        self.__matrix['products']['maps'] += [data]
        self.__connections[mapn] = {}
        self.__connections[mapn][fragment] = self.__connections[fragment][mapn] = self.__connections[target][vector]

        trace.update({vector: mapn})

        for i in set(self.__connections[vector]) - {target}:
            #print i
            if i in trace:
                #print trace
                # это затычка для замыкания циклов. будет создавать несколько лишних циклов. но это не критично.
                old = trace[i]
                self.__connections[mapn][old] = self.__connections[old][mapn] = self.__connections[i][vector]
                continue
            trace.update(self.__nucleofug(vector, mapn, i, trace))
        return trace

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

    def __getSmiles(self, sub):
        '''
        модифицированный алгоритм моргана используется для хэширования подструктур.
        изначально все атомы имеют веса согласно сорту атомов.
        итеративно значения весов атомов умножатся на кратности образуемых связей с последующим суммированием полученных
        значений к текущему весу.

        после возвращается хэш подструктуры.
        '''
        newweights = {x: self.__replace[self.__matrix['products']['maps'][x]['element']] for x in xrange(self.__length)}
        length = len(sub)
        self.__buffer = 0

        while self.__accept(newweights, length):
            weights = copy.copy(newweights)
            for k, i in enumerate(sub):
                l = sum([weights[x] * self.__matrix['products']['bondmatrix'][i][x] for x in self.__connections[i].keys()])
                newweights[i] += l
        return {x: y for x, y in newweights.items() if x in sub}

    def __subsearch(self, trace, inter):
        '''
        поиск всех подструктур граничащих на разрывах и соединениях
        '''
        strace = {inter}
        iterlist = set(self.__connections[inter]).difference(trace)

        for i in iterlist:
            if i in strace:
                # костыль для циклов. чтоб не было 2х проходов.
                continue
            if self.__connections[inter][i] not in (81, 82, 83, 84, 85, 86, 87, 18, 28, 38, 48, 58, 68, 78):
                strace.update(self.__subsearch(copy.copy(trace + [i]), i))
        return strace


if __name__ == '__main__':
    main()
