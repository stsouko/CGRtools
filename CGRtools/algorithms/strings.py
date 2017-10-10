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


def hash_cgr_string(string):
    """
    concatenated md5 and sha256 hashes of cgr string 
    :param string: 
    :return: 48 bytes length string  
    """
    bs = string.encode()
    return md5(bs).digest() + sha256(bs).digest()


def get_cgr_string(g, weights, isotope=False, stereo=False, hyb=False, element=True):
    visited = set()
    countmap, countcyc = count(1), count(1)
    stereo = False

    def getnextatom(atoms):
        if len(atoms) == 1:
            nextatom = list(atoms)[0]
        else:
            nextatom = sorted((weights[i], i) for i in atoms)[-1][1]

        visited.add(nextatom)
        return nextatom

    def dosmarts(trace, inter, prev):
        gni = g.nodes[inter]
        iso = isotope and element and gni.get('isotope')

        if stereo:
            ss = gni.get('s_stereo')
            ps = gni.get('p_stereo')
            if ss:
                s_sh = (';%%s%s;' % hyb_types[gni['s_hyb']] if hyb and gni.get('s_hyb') else ';%s;') % ss
            else:
                s_sh = ''
            if ps:
                p_sh = (';%%s%s;' % hyb_types[gni['p_hyb']] if hyb and gni.get('p_hyb') else ';%s;') % ps
            else:
                p_sh = ''
        elif hyb:
            s_sh = ';%s;' % hyb_types[gni['s_hyb']] if gni.get('s_hyb') else ''
            p_sh = ';%s;' % hyb_types[gni['p_hyb']] if gni.get('p_hyb') else ''
        else:
            s_sh = p_sh = ''

        if iso:
            s_ish = '[%s%%s%s%%s]' % (gni['isotope'], s_sh)
            p_ish = '[%s%%s%s%%s]' % (gni['isotope'], p_sh)
        else:
            s_ish = '[%%s%s%%s]' % s_sh if s_sh else ''
            p_ish = '[%%s%s%%s]' % p_sh if p_sh else ''

        if element:
            ge = gni.get('element') or '*'
            gcs = gni.get('s_charge') and '%+d' % gni['s_charge'] or ''
            gcp = gni.get('p_charge') and '%+d' % gni['p_charge'] or ''
            smis = [s_ish and s_ish % (ge, gcs) or gcs and '[%s%s]' % (ge, gcs) or ge != '*' and ge or '[*]']
            smip = [p_ish and p_ish % (ge, gcp) or gcp and '[%s%s]' % (ge, gcp) or ge != '*' and ge or '[*]']
        else:
            smis = [s_ish % ('*', '') if s_ish else '[*]']
            smip = [p_ish % ('*', '') if p_ish else '[*]']

        concat = []
        stoplist = []
        iterlist = set(g.neighbors(inter)).difference([prev])
        while iterlist:
            i = getnextatom(iterlist)
            iterlist.discard(i)
            gii = g[inter][i]
            if i in trace:
                if i not in stoplist:  # костыль для циклов. чтоб не было 2х проходов.
                    cyc = next(countcyc)
                    concat.append((i, cyc, inter))
                    smis.append('%s%s%d' % (to_smiles[gii.get('s_bond')], stereo and gii.get('s_stereo') or '', cyc))
                    smip.append('%s%s%d' % (to_smiles[gii.get('p_bond')], stereo and gii.get('p_stereo') or '', cyc))
                continue

            deep0, deep1, deep2, deep3 = dosmarts(set(chain(trace, [i])), i, inter)
            trace.update(deep0)
            if deep3:
                concat.extend(deep3)
                for j0, j1, j2 in deep3:
                    if j0 == inter:
                        gij = g[inter][j2]
                        stoplist.append(j2)
                        smis.append('%s%s%d' % (to_smiles[gij.get('s_bond')],
                                                stereo and gij.get('s_stereo') or '', j1))
                        smip.append('%s%s%d' % (to_smiles[gij.get('p_bond')],
                                                stereo and gij.get('p_stereo') or '', j1))
            smis.extend(['(' if iterlist else '',
                         '%s%s' % (to_smiles[gii.get('s_bond')], stereo and gii.get('s_stereo') or '')] + deep1 +
                        [')' if iterlist else ''])
            smip.extend(['(' if iterlist else '',
                         '%s%s' % (to_smiles[gii.get('p_bond')], stereo and gii.get('p_stereo') or '')] + deep2 +
                        [')' if iterlist else ''])
        return trace, smis, smip, concat

    has_next = g
    ssmiles, psmiles = [], []
    while has_next:
        firstatom = getnextatom(has_next)
        smirks = dosmarts({firstatom}, firstatom, firstatom)
        ssmiles.append(''.join(smirks[1]))
        psmiles.append(''.join(smirks[2]))
        has_next = set(g).difference(visited)
    jssmiles = '.'.join(ssmiles)
    jpsmiles = '.'.join(psmiles)
    return '%s>>%s' % (jssmiles, jpsmiles) if jssmiles != jpsmiles else jssmiles


to_smiles = {1: '-', 2: '=', 3: '#', 4: ':', None: '.', 9: '~'}
hyb_types = {4: 'a', 3: 't', 2: 'd', 1: 's', None: ''}
