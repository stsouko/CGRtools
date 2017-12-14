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
    def __init__(self, isotope=False, stereo=False, hyb=False, element=True, is_cgr=False):
        self.__isotope = element and isotope
        self.__stereo = stereo
        self.__hyb = hyb
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
        if self.__stereo:
            ss = gni.get('s_stereo')
            if ss:
                s_sh = (';%%s%s;' % self.__hyb_types[gni['s_hyb']] if self.__hyb and gni.get('s_hyb') else ';%s;') % ss
            else:
                s_sh = ''
        elif self.__hyb:
            s_sh = ';%s;' % self.__hyb_types[gni['s_hyb']] if gni.get('s_hyb') else ''
        else:
            s_sh = ''

        if self.__isotope and gni.get('isotope'):
            s_ish = '[%s%%s%s%%s]' % (gni['isotope'], s_sh)
        else:
            s_ish = '[%%s%s%%s]' % s_sh if s_sh else ''

        if self.__element:
            ge = gni.get('element') or '*'
            gcs = gni.get('s_charge') and '%+d' % gni['s_charge'] or ''
            smis = [s_ish and s_ish % (ge, gcs) or gcs and '[%s%s]' % (ge, gcs) or ge != '*' and ge or '[*]']
        else:
            smis = [s_ish % ('*', '') if s_ish else '[*]']

        return smis, None

    def __cgr_smi(self, gni):
        if self.__stereo:
            ss = gni.get('s_stereo')
            ps = gni.get('p_stereo')
            if ss:
                s_sh = (';%%s%s;' % self.__hyb_types[gni['s_hyb']] if self.__hyb and gni.get('s_hyb') else ';%s;') % ss
            else:
                s_sh = ''
            if ps:
                p_sh = (';%%s%s;' % self.__hyb_types[gni['p_hyb']] if self.__hyb and gni.get('p_hyb') else ';%s;') % ps
            else:
                p_sh = ''
        elif self.__hyb:
            s_sh = ';%s;' % self.__hyb_types[gni['s_hyb']] if gni.get('s_hyb') else ''
            p_sh = ';%s;' % self.__hyb_types[gni['p_hyb']] if gni.get('p_hyb') else ''
        else:
            s_sh = p_sh = ''

        if self.__isotope and gni.get('isotope'):
            s_ish = '[%s%%s%s%%s]' % (gni['isotope'], s_sh)
            p_ish = '[%s%%s%s%%s]' % (gni['isotope'], p_sh)
        else:
            s_ish = '[%%s%s%%s]' % s_sh if s_sh else ''
            p_ish = '[%%s%s%%s]' % p_sh if p_sh else ''

        if self.__element:
            ge = gni.get('element') or '*'
            gcs = gni.get('s_charge') and '%+d' % gni['s_charge'] or ''
            gcp = gni.get('p_charge') and '%+d' % gni['p_charge'] or ''
            smis = [s_ish and s_ish % (ge, gcs) or gcs and '[%s%s]' % (ge, gcs) or ge != '*' and ge or '[*]']
            smip = [p_ish and p_ish % (ge, gcp) or gcp and '[%s%s]' % (ge, gcp) or ge != '*' and ge or '[*]']
        else:
            smis = [s_ish % ('*', '') if s_ish else '[*]']
            smip = [p_ish % ('*', '') if p_ish else '[*]']

        return smis, smip

    def __do_cgr_smarts(self, trace, inter, prev):
        g = self.__g
        countcyc = self.__countcyc
        to_smiles = self.__to_smiles
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


def get_cgr_string(g, weights, isotope=False, stereo=False, hyb=False, element=True, is_cgr=False):
    warn('deprecated function. use CGRstring class instead', DeprecationWarning)
    s = CGRstring(isotope, stereo, hyb, element, is_cgr)
    return s(g, weights)
