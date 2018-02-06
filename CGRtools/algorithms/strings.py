# -*- coding: utf-8 -*-
#
#  Copyright 2017 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
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
        self.__visited = visited = set()
        self.__countmap, self.__countcyc = count(1), count(1)
        self.__g = g

        has_next = g
        ssmiles, psmiles = [], []
        while has_next:
            firstatom = self.__get_next_atom(has_next)
            smirks = self.__do_cgr_smarts({firstatom}, firstatom, firstatom)
            ssmiles.append(''.join(smirks[1]))
            if self.__is_cgr:
                psmiles.append(''.join(smirks[2]))
            has_next = set(g).difference(visited)

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
        return self.__get_sp_smi(gni, 's_stereo', 's_hyb', 's_neighbors', 's_charge'), None

    def __cgr_smi(self, gni):
        return self.__get_sp_smi(gni, 's_stereo', 's_hyb', 's_neighbors', 's_charge'), \
               self.__get_sp_smi(gni, 'p_stereo', 'p_hyb', 'p_neighbors', 'p_charge')

    def __get_sp_smi(self, gni, stereo, hyb, neighbors, charge):
        smi = []
        if self.__stereo and gni.get(stereo):
            smi.append(self.__stereo_types[gni[stereo]])
        if self.__hyb and gni.get(hyb):
            smi.append(self.__hyb_types[gni[hyb]])
        if self.__neighbors and gni.get(neighbors):
            smi.append(str(gni[neighbors]))
        if smi:
            smi.append(';')
            smi.insert(0, ';')

        if self.__element:
            ge = gni.get('element') or '*'
            gcs = gni.get(charge)
            if gcs:
                if abs(gcs) == 1:
                    gcs = '+' if gcs > 0 else '-'
                else:
                    gcs = '%+d' % gcs
                smi.append(gcs)

            smi.insert(0, ge)
        else:
            smi.insert(0, '*')

        if self.__isotope and gni.get('isotope'):
            smi.insert(0, str(gni['isotope']))

        if len(smi) != 1 or smi[0] == '*':
            smi.insert(0, '[')
            smi.append(']')

        return ''.join(smi)

    def __do_cgr_smarts(self, trace, inter, prev):
        g = self.__g
        countcyc = self.__countcyc
        to_smiles = self.__to_smiles
        stereo = self.__stereo

        s_atom, p_atom = self.__get_smi(self.__g.nodes[inter])
        smis, smip = [s_atom], [p_atom]
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
                    smis.append('%s%s%d' % (to_smiles[gii.get('s_bond')], stereo and gii.get('s_stereo') or '', cyc))
                    if smip:
                        smip.append('%s%s%d' % (to_smiles[gii.get('p_bond')],
                                                stereo and gii.get('p_stereo') or '', cyc))
                continue

            deep0, deep1, deep2, deep3 = self.__do_cgr_smarts(set(chain(trace, [i])), i, inter)
            trace.update(deep0)
            if deep3:
                for j0, j1, j2 in deep3:
                    if j0 == inter:
                        gij = g[inter][j2]
                        stoplist.append(j2)
                        smis.append('%s%s%d' % (to_smiles[gij.get('s_bond')],
                                                stereo and gij.get('s_stereo') or '', j1))
                        if smip:
                            smip.append('%s%s%d' % (to_smiles[gij.get('p_bond')],
                                                    stereo and gij.get('p_stereo') or '', j1))
                    else:
                        concat.append((j0, j1, j2))
            smis.extend(['(' if iterlist else '',
                         '%s%s' % (to_smiles[gii.get('s_bond')], stereo and gii.get('s_stereo') or '')] + deep1 +
                        [')' if iterlist else ''])
            if smip:
                smip.extend(['(' if iterlist else '',
                             '%s%s' % (to_smiles[gii.get('p_bond')], stereo and gii.get('p_stereo') or '')] + deep2 +
                            [')' if iterlist else ''])
        return trace, smis, smip, concat

    __to_smiles = {1: '-', 2: '=', 3: '#', 4: ':', None: '.', 9: '~'}
    __hyb_types = {4: 'a', 3: 't', 2: 'd', 1: 's', None: ''}
    __stereo_types = {1: '@', -1: '@@'}


def get_cgr_string(g, weights, isotope=False, stereo=False, hyb=False, element=True, is_cgr=False):
    warn('deprecated function. use CGRstring class instead', DeprecationWarning)
    s = CGRstring(isotope, stereo, hyb, element, is_cgr)
    return s(g, weights)
