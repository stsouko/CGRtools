#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2014 Ramil Nugmanov <stsouko@live.ru>
# This file is part of FEAR (Fix Errors in Automapped Reactions).
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
from collections import defaultdict

import numpy


def main():
    print("RDF parser module")
    tst = RDFread('test.rdf')
    print tst.readdata()
    return 0


class RDFread(object):
    def __init__(self, file):
        self.__RDFfile = file

    def __postProc(self, atomlist, matrix):  #обработчик сырых данных. ставит НЭП на диагональ
        table = dict(H=1, D=1, He=2, Li=1, Be=2, B=3, C=4, N=5, O=6, F=7, Ne=8, Na=1, Mg=2, Al=3, Si=4, P=5, S=6,
                     Cl=7, Ar=8, K=1, Ca=2, Sc=3, Ti=4, V=5, Cr=6, Mn=7, Fe=8, Co=9, Ni=10, Cu=11, Zn=12,
                     Ga=3, Ge=4, As=5, Se=6, Br=7, Kr=8, Rb=1, Sr=2, Y=39, Zr=40, Nb=41, Mo=42, Tc=43, Ru=44,
                     Rh=45, Pd=46, Ag=47, Cd=48, In=49, Sn=50, Sb=5, Te=6, I=7, Xe=8, Cs=1, Ba=2, La=3, Ce=58,
                     Pr=59, Nd=60, Pm=61, Sm=62, Eu=63, Gd=64, Tb=65, Dy=66, Ho=67, Er=68, Tm=69, Yb=70, Lu=71, Hf=72,
                     Ta=73, W=74, Re=75, Os=76, Ir=77, Pt=78, Au=79, Hg=80, Tl=81, Pb=82, Bi=83, Po=84, At=85, Rn=86,
                     Fr=1, Ra=88, Ac=89, Th=90, Pa=91, U=92, Np=93, Pu=94, Am=95, Cm=96, Bk=97)
        charge = {0: 0, 1: 3, 2: 2, 3: 1, 4: 0, 5: -1, 6: -2, 7: -3}  # формат записи зарядов в MDL
        bondSumm = matrix.sum(axis=0)  #посчитаем число связей у каждого атома.
        for i in range(len(atomlist)):  #найдем число несвязанных электронов для каждого атома.
            lonepaire = table[atomlist[i]['element']] - bondSumm[i] + matrix[i, i] / 2 - charge[atomlist[i]['charge']]
            matrix[i][i] += lonepaire if lonepaire > 0 else 0

    def chkRDF(self):
        try:
            f = open(self.__RDFfile)
            if "$RDFILE 1" in f.next():
                return True
            else:
                return False
        except IOError:
            return False

    __aromatize = {1: 1, 2: 2, 3: 3, 4: 1.5}

    def readdata(self):
        '''парсер RDF файлов. возвращает пакет данных вида
        {'substrats':substrats, 'products':products, 'molecules':{'номер':{'atomlist':atomlist, 'bondmatrix':matrix}}}
        '''
        with open(self.__RDFfile) as f:
            ir = -1
            im = -1
            atomcount = -1
            bondcount = -1
            failkey = True
            reaction = None
            meta = None
            for n, line in enumerate(f):
                if not failkey and "$RXN" not in line[0:4]:
                    continue
                elif "$RXN" in line[0:4]:
                    if reaction:
                        yield reaction
                    reaction = {'substrats': [], 'products': [], 'meta': defaultdict(str)}
                    meta = None
                    ir = n + 4
                    failkey = True
                elif n == ir:
                    try:
                        substrats, products = int(line[0:3]), int(line[3:6])
                    except:
                        failkey = False
                        reaction = None
                elif "$MOL" in line[0:4]:
                    molecule = {'atomlist': {}, 'bondmatrix': None, 'bondlist': []}
                    im = n + 4
                    pnum = 0
                elif n == im:
                    try:
                        atomcount = int(line[0:3]) + im
                        bondcount = int(line[3:6]) + atomcount
                        if atomcount == bondcount:
                            molecule['bondmatrix'] = numpy.zeros((1, 1), dtype=float)
                    except:
                        failkey = False
                        reaction = None
                elif n <= atomcount:
                    if line[31:33].strip() != "H":
                        molecule['atomlist'][pnum] = dict(element=line[31:34].strip(), izotop=line[34:36].strip(),
                                                          charge=int(line[38:39]), map=int(line[60:63]),
                                                          mark=line[51:54].strip(),
                                                          x=float(line[0:10]),
                                                          y=float(line[10:20]),
                                                          z=float(line[20:30]))
                    pnum += 1
                elif n <= bondcount:
                    try:
                        if molecule['bondmatrix'] is None:
                            mlen = len(molecule['atomlist'])
                            molecule['bondmatrix'] = numpy.zeros((mlen, mlen), dtype=float)
                        a1, a2 = int(line[0:3]) - 1, int(line[3:6]) - 1
                        if a1 in molecule['atomlist'] and a2 in molecule['atomlist']:
                            molecule['bondmatrix'][a1, a2] = molecule['bondmatrix'][a2, a1] = self.__aromatize[int(line[6:9])]
                            molecule['bondlist'].append((a1 + 1, a2 + 1, int(line[6:9]), int(line[9:12])))
                        elif a1 in molecule['atomlist']:
                            molecule['bondmatrix'][a1, a1] += 2
                        else:
                            molecule['bondmatrix'][a2, a2] += 2
                    except:
                        failkey = False
                        reaction = None
                elif "M  END" in line:
                    molecule['atomlist'] = molecule['atomlist'].values()
                    try:
                        self.__postProc(molecule['atomlist'], molecule['bondmatrix'])
                        if len(reaction['substrats']) < substrats:
                            reaction['substrats'].append(molecule)
                        else:
                            reaction['products'].append(molecule)
                    except:
                        failkey = False
                        reaction = None
                elif '$DTYPE' in line:
                    meta = line[7:].strip()
                elif '$RFMT' not in line and meta:
                    reaction['meta'][meta] += line.strip("$DATUM").strip() + ' '
            else:
                if reaction:
                    yield reaction

if __name__ == '__main__':
    main()
