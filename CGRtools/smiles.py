#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2015 Ramil Nugmanov <stsouko@live.ru>
# This file is part of fragger.
#
#  fragger is free software; you can redistribute it and/or modify
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


class smiles():
    tosmileskeys = {0: '', 1: '-', 2: '=', 3: '#', 1.5: ':', 4: ':', 8: '~', 9: '.'}

    def __tosmiles(self, bond):
        if bond > 10:
            b1 = bond // 10 % 10
            b2 = bond % 10
            s = [self.tosmileskeys[b1 if b1 != 8 else 9], self.tosmileskeys[b2 if b2 != 8 else 9]]
            if bond > 100:
                s.append(self.tosmileskeys[bond // 100])
            s = ''.join(s)
        else:
            s = self.tosmileskeys[bond]
        return s

    def get_name(self, sub):
        self.__nextnumb = self.__numb()
        self.__sub = set(sub)
        self.__levels = self.__getMorgan(sub)
        inter = min(self.__levels, key=self.__levels.get)
        res = self.__getSmiles([inter], inter)
        return ''.join(res[1])

    def __numb(self):
        i = 1
        while True:
            yield i
            i += 1

    def __getSmiles(self, trace, inter):
        strace = set(trace)
        iterlist = set(self.connections[inter]).intersection(self.__sub).difference(trace[-2:])

        # solution for unbounded fragments
        if iterlist and strace != self.__sub:
            remained = {k: self.__levels[k] for k in self.__sub.difference(strace)}
            iterlist = {min(remained, key=remained.get)}

        # get atom labels
        smi = [self.table[inter]]
        concat = []
        stoplist = []
        iterlen = len(iterlist) - 1
        for b, i in enumerate(sorted(list(iterlist), key=self.__levels.get)):
            if i in strace:
                if i not in stoplist:
                    # костыль для циклов. чтоб не было 2х проходов.
                    cyc = next(self.__nextnumb)
                    concat += [(i, cyc, inter)]
                    smi[0] += '%s%d' % (self.__tosmiles(self.connections[inter][i]), cyc)
                if b == iterlen and len(smi) > 3:
                    smi[-1] = smi[-3] = ''
                continue
            deep = self.__getSmiles(copy.copy(trace + [i]), i)
            strace.update(deep[0])

            for j in deep[2]:
                if j[0] == inter:
                    stoplist += [j[2]]
                    smi[0] += '%s%d' % (self.__tosmiles(self.connections[inter][j[2]]), j[1])
                else:
                    concat.append(j)

            smi += ['(' if iterlen - b else '', '%s' % self.__tosmiles(self.connections[inter][i]) + deep[1],
                    ')' if iterlen - b else '']

        return strace, ''.join(smi), concat

    def __getMorgan(self):
        """
        модифицированный алгоритм моргана используется для хэширования подструктур.
        итеративно значения весов атомов умножатся на кратности образуемых связей с последующим суммированием полученных
        значений к текущему весу.
        после возвращается хэш подструктуры.
        """
        newweights = {x: self.data['atoms'][x]['weight'] for x in self.__sub}
        self.__buffer = 0
        self.__pososhok = False
        while self.__accept(newweights):
            weights = copy.copy(newweights)
            for i in self.__sub:
                l = sum([weights[x] * self.connections[i][x] / 10 for x in
                         set(self.connections[i]).intersection(self.__sub)])
                newweights[i] += l
        return newweights

    def __accept(self, numb):
        """
        функция остановки для алгоритма моргана. останавливает поиск после первого простоя улучшения.
        """
        numl = len(set(numb.values()))

        if self.__buffer == numl:
            return False
        self.__buffer = numl
        return True


if __name__ == '__main__':
    print('This file is part of fragger')
