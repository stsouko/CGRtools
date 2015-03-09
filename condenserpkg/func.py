# -*- coding: utf-8 -*-
#
# func.py
#
# Copyright 2014 Ramil Nugmanov <stsouko@live.ru>
# This file is part of condenser.
#
# condenser is free software; you can redistribute it and/or modify
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

from collections import defaultdict
import copy
import itertools
from weightable import replace
import numpy


def main():
    print "This file is part of condenser."
    return 0


class Condenser(object):
    def __init__(self, charge, type, repare, stereo):
        self.__matrix = {}
        self.__charge = charge
        self.__cgrtype = type
        self.__stereo = stereo
        self.__repare = repare
        self.__cgr = self.__tocgr()
        self.__replace = replace()

    def calc(self, data):
        self.__moltorepare = False
        self.__data = data
        self.__creatematrix()  #подготовка матриц
        #        print(self.__matrix)
        self.__scan()
        #print(self.__matrix)
        out = self.__matrix['substrats']
        return dict(meta=data['meta'], stereo=out['stereo'], dynbonds=out['dynbonds'], dyncharges=out['dyncharge'],
                    maps=out['maps'], diff=out['diff'], products_maps=self.__matrix['products']['maps'])

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
            self.__getmatrix([mols], 'substrats')
            self.__matrix['products'] = self.__matrix['substrats']
        elif 20 < self.__cgrtype < 30:
            try:
                mols = self.__data['products'][self.__cgrtype - 21]
            except:
                mols = self.__data['products'][-1]
            self.__getmatrix([mols], 'products')
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

    __toMDL = {-3: 7, -2: 6, -1: 5, 0: 0, 1: 3, 2: 2, 3: 1}

    __fromMDL = {0: 0, 1: 3, 2: 2, 3: 1, 4: 0, 5: -1, 6: -2, 7: -3}

    __bondlabels = {'0': 0, '1': 1, '2': 2, '3': 3, 'a': 4, 'c': 5, 'n': 0, 'h': 6}

    def __getmatrix(self, data, category):
        '''
        постройка матриц связности
        '''
        atomlist = []
        datF = {'dynatom': [], 'dynbond': [], 'dynatomstereo': [], 'dynbondstereo': [],
                'bondstereo': defaultdict(dict), 'atomstereo': {}, 'extrabond': []}

        #объединим атомы в 1 список. и создадим пустую матрицу.
        for i in data:
            dat = {'dynatom': [], 'dynbond': [], 'dynatomstereo': [], 'dynbondstereo': [],
                   'bondstereo': [], 'atomstereo': [], 'extrabond': []}

            shift = len(atomlist)
            for j in i['DAT'].values():
                dat[j['type']].append([sorted([x + shift for x in j['atoms']]), j['value']])

            atomlist.extend(i['atomlist'])

            for j in dat['dynbondstereo']:
                # если на вход дали кгр, то надо выдрать состояние в соответствии с положением (субстрат или продукт)
                val = j[1].lower().split('>')
                datF['bondstereo'][j[0][0]][j[0][1]] = datF['bondstereo'][j[0][1]][j[0][0]] = val[
                    0] if category == 'substrats' else val[-1]

            for j in dat['dynatomstereo']:
                val = j[1].lower().split('>')
                datF['atomstereo'][j[0][0]] = val[0] if category == 'substrats' else val[-1]

            for j in dat['bondstereo']:
                datF['bondstereo'][j[0][0]][j[0][1]] = datF['bondstereo'][j[0][1]][j[0][0]] = j[1].lower()

            for j in dat['atomstereo']:
                datF['atomstereo'][j[0][0]] = j[1]

            for j in dat['dynatom']:
                x = j[0][0] - 1  # atom index
                key = j[1][0].lower()
                diff = int(j[1][1:])
                if key == 'c':  #update atom charges from CGR
                    charge = self.__fromMDL.get(i['atomlist'][x]['charge'], 0)
                    v = 10 if category == 'substrats' else charge + diff
                    if v != 10:
                        i['atomlist'][x]['charge'] = self.__toMDL.get(v, 0)
                elif key == 'h':
                    pass
                elif key == 'r':
                    pass

            for j in dat['dynbond']:  #update bondmatrix from CGR reactant bonds
                x, y = j[0][0] - 1 - shift, j[0][1] - 1 - shift
                raw = j[1].lower().split('>')
                v = self.__bondlabels.get(raw[0] if category == 'substrats' else raw[-1], 0)
                i['bondmatrix'][x, y] = i['bondmatrix'][y, x] = v

            for j in dat['extrabond']:
                x, y = j[0][0] - 1 - shift, j[0][1] - 1 - shift
                raw = j[1].lower()
                i['bondmatrix'][x, y] = i['bondmatrix'][y, x] = self.__bondlabels.get(raw)

        length = len(atomlist)
        self.__length = length

        tempMatrix = dict(bondmatrix=numpy.zeros((length, length), dtype=int), maps=atomlist, diff=[],
                          atomstereo=datF['atomstereo'],
                          bondstereo=datF['bondstereo'],
                          stereo=[], dynbonds=[], dyncharge=[])
        step = 0
        for i in data:
            matrixlength = len(i['atomlist'])
            tempMatrix['bondmatrix'][step:step + matrixlength, step:step + matrixlength] = i['bondmatrix']
            step += matrixlength

        self.__matrix[category] = tempMatrix

    def __sortmatrix(self):
        """
        атомы отсортировываются по порядку. для атомов представленных только в продуктах или реагентах достраивается
        симметричная матрица с разрывом в месте стыка с основной молекулой.
        """
        maps = {'substrats': [int(x['map']) for x in self.__matrix['substrats']['maps']],
                'products': [int(x['map']) for x in self.__matrix['products']['maps']]}

        length = max((max(maps['products']), max(maps['substrats'])))

        for j in ('substrats', 'products'):
            if 0 in maps[j]:
                nmaps = []
                for x in maps[j]:
                    if x:
                        nmaps.append(x)
                    else:
                        length += 1
                        nmaps.append(length)
                maps[j] = nmaps

        lose = sorted(list(set(range(1, length + 1)).difference(maps['products']).difference(maps['substrats'])),
                      reverse=True)
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
            dat = {}
            for i, k in enumerate(zip(self.__matrix[j]['maps'], maps[j])):
                k[0]['map'] = k[1]
                atomlist[k[1] - 1] = k[0]
                dat[i + 1] = k[1] - 1
            self.__matrix[j]['maps'] = atomlist

            newbondstereoblock = defaultdict(dict)
            for x, y in self.__matrix[j]['bondstereo'].items():
                x = dat[x]
                for z, w in y.items():
                    z = dat[z]
                    newbondstereoblock[z][x] = newbondstereoblock[x][z] = w
            self.__matrix[j]['bondstereo'] = newbondstereoblock

            newatomstereoblock = {}
            for x, y in self.__matrix[j]['atomstereo'].items():
                newatomstereoblock[dat[x]] = y
            self.__matrix[j]['atomstereo'] = newatomstereoblock

        self.__lostlist = dict.fromkeys(('products', 'substrats'))

        if 0 in self.__matrix['substrats']['maps']:
            '''
            метод восстановления используется только для недостающих атомов реагентов.
            информация о факте отвалившегося фрагмента в продуктах самодостаточна для корректного ААО или полного
            восстановления структуры уходящей группы.
            '''
            self.__moltorepare = True
            self.__searchlost('products', 'substrats')
        if 0 in self.__matrix['products']['maps']:
            self.__moltocharge = False
            self.__searchlost('substrats', 'products')

    def __searchlost(self, s, t):
        """
        внедрение в матрицу продуктов/субстратов атомов из субстратов/продуктов
        """
        lostlist = []
        for i in xrange(self.__length):
            if self.__matrix[t]['maps'][i] == 0:
                self.__matrix[t]['maps'][i] = self.__matrix[s]['maps'][i]
                #print s, self.__matrix[s]['maps'][i]
                if i in self.__matrix[s]['atomstereo']:
                    self.__matrix[t]['atomstereo'][i] = self.__matrix[s]['atomstereo'][i]
                lostlist.append(i)
        # надо сохранить список несбалансированных атомов. и когда их нельзя сбалансировать за счет списывания у многостадийных реакций, то надо внести фиктивные атомы итд.
        self.__lostlist[t] = lostlist
        for i in lostlist:
            first = self.__matrix[s]['bondstereo'].get(i)
            for j in lostlist:
                if first:
                    second = first.get(j)
                    if second:
                        self.__matrix[t]['bondstereo'][i][j] = second

                self.__matrix[t]['bondmatrix'][i, j] = self.__matrix[s]['bondmatrix'][i, j]

    def __tocgr(self):
        ways = ['0', '1', '2', '3', 'a', 'c', 'h']
        bondcodes = {0: '0', 1: '1', 2: '2', 3: '3', 4: '4', 5: 'c', 6: 'h'}
        cgrdict = defaultdict(dict)
        for x, y in itertools.permutations(ways, 2):
            cgrdict[self.__bondlabels[x]][self.__bondlabels[y]] = ('%s>%s' % (x, y))
        for x, y in cgrdict.items():
            y[x] = bondcodes[x]
        return cgrdict

    def __chdiffpseudo(self, length):
        for i in xrange(length):
            s = self.__fromMDL.get(self.__matrix['substrats']['maps'][i]['charge'], 0)
            p = self.__fromMDL.get(self.__matrix['products']['maps'][i]['charge'], 0)
            diff = p - s
            if diff:
                self.__matrix['substrats']['dyncharge'].append([i + 1, 'c%+d' % diff, 'dyncharge'])

    def __chreplacepseudo(self, length):
        for i in xrange(length):
            self.__matrix['substrats']['maps'][i]['charge'] = self.__matrix['products']['maps'][i]['charge']

    def __getStereo(self, target, i, stype):
        if stype == 'atomstereo':
            x = self.__matrix[target][stype].get(i, 'n')
        elif stype == 'bondstereo':
            x = self.__matrix[target][stype].get(i[0])
            if x:
                x = x.get(i[1], 'n')
            else:
                x = 'n'
        else:
            x = 'n'
        return x

    def __scan(self):
        length = self.__length
        connections = {x: {} for x in range(length)}

        for i in xrange(length - 1):
            if self.__stereo:
                a, b = self.__getStereo('substrats', i, 'atomstereo'), self.__getStereo('products', i, 'atomstereo')
                if a == b:
                    if a != 'n':
                        stereo = a
                        stype = 'atomstereo'
                    else:
                        stereo = False
                else:
                    stereo = '%s>%s' % (a, b)
                    stype = 'dynatomstereo'
                if stereo:
                    self.__matrix['substrats']['stereo'] += [(i + 1, stereo, stype)]
            for j in xrange(i + 1, length):
                if self.__matrix['substrats']['bondmatrix'][i][j] != 0 or \
                                self.__matrix['products']['bondmatrix'][i][j] != 0:
                    diff = self.__cgr[self.__matrix['substrats']['bondmatrix'][i][j]][
                        self.__matrix['products']['bondmatrix'][i][j]]

                    if not diff.isdigit():
                        if len(diff) == 1:
                            btype = 'extrabond'
                        else:
                            btype = 'dynbond'
                        self.__matrix['substrats']['dynbonds'].append((i + 1, j + 1, diff, btype))
                        diff = '8'

                    if self.__stereo:
                        a = self.__getStereo('substrats', [i, j], 'bondstereo')
                        b = self.__getStereo('products', [i, j], 'bondstereo')
                        if a == b:
                            if a != 'n':
                                stereo = a
                                stype = 'bondstereo'
                            else:
                                stereo = None
                        else:
                            stereo = '%s>%s' % (a, b)
                            stype = 'dynbondstereo'
                        if stereo:
                            self.__matrix['substrats']['stereo'] += [(i + 1, j + 1, stereo, stype)]
                    self.__matrix['substrats']['diff'] += [(i + 1, j + 1, diff)]
                    connections[i][j] = connections[j][i] = diff

        self.__connections = copy.deepcopy(connections)

        if self.__repare and (self.__moltorepare or self.__moltocharge):
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
            #сначала попробуем восстановить инфу о продуктах. это позволит для многостадийных реакций лучше восстановить баланс.
            if self.__moltocharge:
                try:
                    self.__chargebalancer()
                except:
                    self.__connections = connections
                    self.__matrix = backup
            if self.__moltorepare:
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

        if self.__charge == 2:
            self.__chreplacepseudo(length)
        elif self.__charge == 1:
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
                        print x, self.__matrix['substrats']['maps'][x]['charge'], self.__matrix['products']['maps'][x][
                            'charge']
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
        titerlist = set(self.__connections[target])  #.difference(trace)
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
                    fi = next((key for key in self.__connections[fragment] if
                               self.__fragment.get(key) == self.__target[i] and key not in trace),
                              next((key for key, val in self.__fragment.items() if val == self.__target[i]), None))
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
                l = sum(
                    [weights[x] * self.__matrix['products']['bondmatrix'][i][x] for x in self.__connections[i].keys()])
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
