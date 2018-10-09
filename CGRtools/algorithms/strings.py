# -*- coding: utf-8 -*-
#
#  Copyright 2017, 2018 Ramil Nugmanov <stsouko@live.ru>
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
from hashlib import md5, sha256
from itertools import chain, count
from warnings import warn


def hash_cgr_string(string):
    """
    concatenated md5 and sha256 hashes of cgr string 
    :param string: 
    :return: 48 bytes length string  
    """
    bs = string.encode()
    return md5(bs).digest() + sha256(bs).digest()


class CGRstring:
    def __init__(self, element=True, isotope=False, stereo=False, hybridization=False, neighbors=False, is_cgr=False):
        self.__isotope = element and isotope
        self.__stereo = stereo
        self.__hyb = hybridization
        self.__neighbors = neighbors
        self.__element = element
        self.__get_smi = self.__cgr_smi if is_cgr else self.__mol_smi
        self.__is_cgr = is_cgr

    def __call__(self, g, weights):
        self.__weights = weights
        self.__countmap, self.__countcyc = count(1), count(1)
        self.__g = g

        has_next = set(g)
        ssmiles, psmiles = [], []
        while has_next:
            self.__visited = set()
            firstatom = self.__get_next_atom(has_next)
            smirks = self.__do_cgr_smarts({firstatom}, firstatom, firstatom)
            ssmiles.append(''.join(smirks[1]))
            if self.__is_cgr:
                psmiles.append(''.join(smirks[2]))
            has_next.difference_update(self.__visited)

        jssmiles = '.'.join(ssmiles)
        if self.__is_cgr:
            jpsmiles = '.'.join(psmiles)
            return '%s>>%s' % (jssmiles, jpsmiles) if jssmiles != jpsmiles else jssmiles

        return jssmiles

    def __get_next_atom(self, atoms):
        if len(atoms) == 1:
            nextatom = list(atoms)[0]
        else:
            weights = self.__weights
            nextatom = sorted((weights[i], i) for i in atoms)[-1][1]

        self.__visited.add(nextatom)
        return nextatom

    def __mol_smi(self, gni):
        return self.__get_sp_smi(gni)

    def __cgr_smi(self, gni):
        return self.__get_sp_smi(gni, True)

    def __get_sp_smi(self, gni, cgr=False):
        smi = []
        pmi = [] if cgr else None

        if self.__stereo:
            if gni.get('s_stereo'):
                smi.append(self.__stereo_types[gni['s_stereo']])
            if cgr and gni.get('s_stereo'):
                pmi.append(self.__stereo_types[gni['p_stereo']])

        if self.__hyb:
            s = gni.get('s_hyb')
            if cgr:
                p = gni.get('p_hyb')
                if isinstance(s, list):
                    tmp = sorted((self.__hyb_types[x], self.__hyb_types[y]) for x, y in zip(s, p))
                    smi.append('<%s>' % ''.join(x for x, _ in tmp))
                    pmi.append('<%s>' % ''.join(x for _, x in tmp))
                else:
                    if s:
                        smi.append(self.__hyb_types[s])
                    if p:
                        pmi.append(self.__hyb_types[p])
            elif s:
                smi.append(self.__hyb_types[s])

        if self.__neighbors:
            s = gni.get('s_neighbors')
            if cgr:
                p = gni.get('p_neighbors')
                if isinstance(s, list):
                    tmp = sorted((x is None and 'n' or str(x), y is None and 'n' or str(y)) for x, y in zip(s, p))
                    smi.append('<%s>' % ''.join(x for x, _ in tmp))
                    pmi.append('<%s>' % ''.join(x for _, x in tmp))
                else:
                    if s is not None:
                        smi.append(str(s))
                    if p is not None:
                        pmi.append(str(p))
            elif s is not None:
                smi.append(str(s))

        if smi:
            smi.append(';')
            smi.insert(0, ';')
        if pmi:
            pmi.append(';')
            pmi.insert(0, ';')

        if self.__element:
            ge = gni.get('element')
            if cgr:
                if isinstance(ge, list):
                    ge = list(','.join(sorted(ge)))
                    smi = ge + smi
                    pmi = ge + pmi
                else:
                    if not ge:
                        ge = '*'
                    smi.insert(0, ge)
                    pmi.insert(0, ge)
            else:
                smi.insert(0, ge)

            s = gni.get('s_charge')
            if cgr:
                p = gni.get('p_charge')
                if isinstance(s, list):
                    tmp = [(self.__charge_to_string(x), self.__charge_to_string(y)) for x, y in sorted(zip(s, p))]
                    smi.append('<%s>' % ''.join(x for x, _ in tmp))
                    pmi.append('<%s>' % ''.join(x for _, x in tmp))
                else:
                    if s:
                        smi.append(self.__charge_to_string(s))
                    if p:
                        pmi.append(self.__charge_to_string(p))
            elif s:
                smi.append(self.__charge_to_string(s))

            s = gni.get('s_radical')
            if cgr:
                p = gni.get('p_radical')
                if isinstance(s, list):
                    tmp = [(self.__radical_map[x], self.__radical_map[y]) for x, y in sorted(zip(s, p))]
                    smi.append('<%s>' % ''.join(x for x, _ in tmp))
                    pmi.append('<%s>' % ''.join(x for _, x in tmp))
                else:
                    if s:
                        smi.append(self.__radical_map[s])
                    if p:
                        pmi.append(self.__radical_map[p])
            elif s:
                smi.append(self.__radical_map[s])

        else:
            smi.insert(0, '*')
            if cgr:
                pmi.insert(0, '*')

        if self.__isotope:
            s = gni.get('isotope')
            if cgr:
                if isinstance(s, list):
                    s = '<%s>' % ','.join(str(x) for x in s)
                    smi.insert(0, s)
                    pmi.insert(0, s)
                elif s:
                    s = str(s)
                    smi.insert(0, s)
                    pmi.insert(0, s)
            elif s:
                smi.insert(0, str(s))

        if len(smi) != 1 or smi[0] == '*':
            smi.insert(0, '[')
            smi.append(']')

        if cgr:
            if len(pmi) != 1 or pmi[0] == '*':
                pmi.insert(0, '[')
                pmi.append(']')

        return smi, pmi

    def __do_cgr_smarts(self, trace, inter, prev):
        g = self.__g
        countcyc = self.__countcyc
        stereo = self.__stereo

        smis, smip = self.__get_smi(self.__g.nodes[inter])
        concat = []
        stoplist = []
        iterlist = set(g.neighbors(inter)).difference([prev])
        while iterlist:
            i = self.__get_next_atom(iterlist)
            iterlist.discard(i)
            gii = g[inter][i]
            if i in trace:
                if i not in stoplist:  # костыль для циклов. чтоб не было 2х проходов.
                    cyc = next(countcyc)
                    concat.append((i, cyc, inter))
                    if smip:
                        sb, pb = self.__bond_sp(gii.get('s_bond'), gii.get('p_bond'))
                        smip.append('%s%s%d' % (pb, stereo and gii.get('p_stereo') or '', cyc))
                    else:
                        sb = self.__bond_s(gii.get('s_bond'))
                    smis.append('%s%s%d' % (sb, stereo and gii.get('s_stereo') or '', cyc))
                continue

            deep0, deep1, deep2, deep3 = self.__do_cgr_smarts(set(chain(trace, [i])), i, inter)
            trace.update(deep0)
            if deep3:
                for j0, j1, j2 in deep3:
                    if j0 == inter:
                        gij = g[inter][j2]
                        stoplist.append(j2)
                        if smip:
                            sb, pb = self.__bond_sp(gij.get('s_bond'), gij.get('p_bond'))
                            smip.append('%s%s%d' % (pb, stereo and gij.get('p_stereo') or '', j1))
                        else:
                            sb = self.__bond_s(gij.get('s_bond'))
                        smis.append('%s%s%d' % (sb, stereo and gij.get('s_stereo') or '', j1))
                    else:
                        concat.append((j0, j1, j2))

            if smip:
                sb, pb = self.__bond_sp(gii.get('s_bond'), gii.get('p_bond'))
                if iterlist:
                    smip.append('(')
                smip.append(pb + (stereo and gii.get('p_stereo') or ''))
                smip.extend(deep2)
                if iterlist:
                    smip.append(')')
            else:
                sb = self.__bond_s(gii.get('s_bond'))

            if iterlist:
                smis.append('(')
            smis.append(sb + (stereo and gii.get('s_stereo') or ''))
            smis.extend(deep1)
            if iterlist:
                smis.append(')')

        return trace, smis, smip, concat

    @classmethod
    def __bond_s(cls, s):
        if isinstance(s, list):
            return '<%s>' % ''.join(sorted(cls.__to_smiles_map[x] for x in s))
        return cls.__to_smiles_map[s]

    @classmethod
    def __bond_sp(cls, s, p):
        if isinstance(s, list):
            tmp = sorted((cls.__to_smiles_map[x], cls.__to_smiles_map[y]) for x, y in zip(s, p))
            return '<%s>' % ''.join(x for x, _ in tmp), '<%s>' % ''.join(x for _, x in tmp)
        return cls.__to_smiles_map[s], cls.__to_smiles_map[p]

    @classmethod
    def __charge_to_string(cls, c):
        if c == 1:
            c = '+'
        elif c == -1:
            c = '-'
        elif c:
            c = '%+d' % c
        else:
            c = '0'
        return c

    __to_smiles_map = {1: '-', 2: '=', 3: '#', 4: ':', None: '.', 9: '~'}
    __hyb_types = {4: 'a', 3: 't', 2: 'd', 1: 's', None: 'n'}
    __stereo_types = {1: '@', -1: '@@'}
    __radical_map = {1: '*', 2: '*2', 3: '*3', None: 'n'}


def get_cgr_string(g, weights, isotope=False, stereo=False, hyb=False, element=True, is_cgr=False):
    warn('deprecated function. use CGRstring class instead', DeprecationWarning)
    s = CGRstring(isotope, stereo, hyb, element, is_cgr)
    return s(g, weights)
