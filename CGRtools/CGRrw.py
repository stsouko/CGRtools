#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
# Copyright 2014-2016 Ramil Nugmanov <stsouko@live.ru>
# This file is part of FEAR (C) Ramil Nugmanov <stsouko@live.ru>.
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
from itertools import count, product
from collections import defaultdict
import periodictable as pt

toMDL = {-3: 7, -2: 6, -1: 5, 0: 0, 1: 3, 2: 2, 3: 1}
fromMDL = {0: 0, 1: 3, 2: 2, 3: 1, 4: 0, 5: -1, 6: -2, 7: -3}
bondlabels = {'0': None, '1': 1, '2': 2, '3': 3, '4': 4, '9': 9}


class CGRRead:
    def __init__(self):
        self.__prop = {}
        self.__mendeleyset = set(pt.elements)

    def collect(self, line):
        if 'M  ALS' in line:
            self.__prop[line[3:10]] = dict(atoms=[int(line[7:10])],
                                           type='atomlist' if line[14] == 'F' else 'atomnotlist',
                                           value=[line[16 + x*4: 20 + x*4].strip() for x in range(int(line[10:13]))])
        elif 'M  ISO' in line:
            self.__prop[line[3:10]] = dict(atoms=[int(line[10:13])], type='isotop', value=int(line[14:17]))

        elif 'M  STY' in line:
            for i in range(int(line[8])):
                if 'DAT' in line[10 + 8 * i:17 + 8 * i]:
                    self.__prop[int(line[10 + 8 * i:13 + 8 * i])] = {}
        elif 'M  SAL' in line:
            if int(line[7:10]) in self.__prop:
                key = []
                for i in range(int(line[10:13])):
                    key.append(int(line[14 + 4 * i:17 + 4 * i]))
                self.__prop[int(line[7:10])]['atoms'] = sorted(key)
        elif 'M  SDT' in line and int(line[7:10]) in self.__prop:
            key = line.split()[-1].lower()
            if key not in self.__cgrkeys:
                self.__prop.pop(int(line[7:10]))
            else:
                self.__prop[int(line[7:10])]['type'] = key
        elif 'M  SED' in line and int(line[7:10]) in self.__prop:
            self.__prop[int(line[7:10])]['value'] = line[10:].strip().replace('/', '').lower()

    __cgrkeys = dict(dynatom=1, atomstereo=1, dynatomstereo=1,
                     bondstereo=2, dynbondstereo=2, extrabond=2, dynbond=2,
                     atomhyb=1, atomneighbors=1, dynatomhyb=1, dynatomneighbors=1,
                     atomnotlist=1, atomlist=1, isotop=1)

    def getdata(self):
        prop = []
        for i in self.__prop.values():
            if len(i['atoms']) == self.__cgrkeys[i['type']]:
                prop.append(i)

        self.__prop = {}
        return prop

    def cgr_dat(self, g, k, atom1, atom2):

        def _parsedyn(target, name):
            tmp = ([], [])
            for x in k['value'].split(','):
                s, *_, p = x.split('>')
                tmp[0].append((int(s) if s != 'n' else None) if name != 'bond' else bondlabels[s])
                tmp[1].append((int(p) if p != 'n' else None) if name != 'bond' else bondlabels[p])

            if len(tmp[0]) == 1:
                target['sp_%s' % name] = tmp[0][0], tmp[1][0]
                for sp_key, sp_ind in (('s', 0), ('p', 1)):
                    if tmp[sp_ind][0] is not None:
                        target['%s_%s' % (sp_key, name)] = tmp[sp_ind][0]
                    else:
                        target.pop('%s_%s' % (sp_key, name), None)
            else:
                target['sp_%s' % name] = list(zip(*tmp))
                target.pop('s_%s' % name, None)
                target.pop('p_%s' % name, None)

        def _parselist(target, name):
            tmp = [bondlabels[x] if name == 'bond' else (int(x) if x != 'n' else None) for x in k['value'].split(',')]
            if len(tmp) == 1:
                target['s_%s' % name] = target['p_%s' % name] = target['sp_%s' % name] = tmp[0]
            else:
                target['sp_%s' % name] = tmp
                target.pop('s_%s' % name, None)
                target.pop('p_%s' % name, None)

        if k['type'] == 'dynatomstereo':
            _parsedyn(g.node[atom1], 'stereo')

        elif k['type'] == 'atomstereo':
            _parselist(g.node[atom1], 'stereo')

        elif k['type'] == 'dynatom':
            key = k['value'][0]
            diff = k['value'][1:].split(',')
            if key == 'c':  # update atom charges from CGR
                base = fromMDL.get(g.node[atom1]['s_charge'], 0)

                if len(diff) > 1:
                    g.node[atom1].pop('s_charge')
                    g.node[atom1].pop('p_charge')
                    if int(diff[0]) == 0:  # for list of charges c0,1,-2,+3...
                        g.node[atom1]['sp_charge'] = [toMDL.get(base + int(x), 0) for x in diff]
                    elif (len(diff) - 1) % 2 == 0:  # for dyn charges c1,-1,0... -1,0 is dyn group relatively to base
                        s = [toMDL.get(base + int(x), 0) for x in [0] + diff[1::2]]
                        p = [toMDL.get(fromMDL.get(x, 0) + int(y), 0) for x, y in zip(s, diff[::2])]
                        g.node[atom1]['sp_charge'] = list(zip(s, p))

                else:
                    g.node[atom1]['p_charge'] = toMDL.get(base + int(diff[0]), 0)
                    g.node[atom1]['sp_charge'] = (g.node[atom1]['s_charge'], g.node[atom1]['p_charge'])

            elif key == '*':
                pass  # not implemented

        elif k['type'] == 'dynbond':
            _parsedyn(g.edge[atom1][atom2], 'bond')

        elif k['type'] == 'bondstereo':
            _parselist(g.edge[atom1][atom2], 'stereo')

        elif k['type'] == 'dynbondstereo':
            _parsedyn(g.edge[atom1][atom2], 'stereo')

        elif k['type'] == 'extrabond':
            _parselist(g.edge[atom1][atom2], 'bond')

        elif k['type'] == 'atomlist':
            g.node[atom1]['element'] = k['value']

        elif k['type'] == 'atomnotlist':
            g.node[atom1]['element'] = list(self.__mendeleyset.difference(k['value']))

        elif k['type'] == 'atomhyb':
            _parselist(g.node[atom1], 'hyb')

        elif k['type'] == 'atomneighbors':
            _parselist(g.node[atom1], 'neighbors')

        elif k['type'] == 'dynatomhyb':
            _parsedyn(g.node[atom1], 'hyb')

        elif k['type'] == 'dynatomneighbors':
            _parsedyn(g.node[atom1], 'neighbors')

        elif k['type'] == 'isotop':
            g.node[atom1]['isotop'] = k['value']

    @staticmethod
    def parsecolors(g, key, colors, remapped):
        adhoc, before = [], []
        for x in colors:
            if (len(x) == 81 or len(x) == 75 and not before) and x[-1] == '+':
                before.append(x[:-1])
            else:
                adhoc.append(''.join(before + [x]))
                before = []

        for population, *keys in (x.split() for x in adhoc):
            for atom, val in (x.split(':') for x in keys):
                atom = remapped[int(atom)]
                if 'dyn' not in key[:3]:
                    g.node[atom].setdefault('s_%s' % key, {})[int(population)] = val
                    g.node[atom].setdefault('p_%s' % key, {})[int(population)] = val
                else:
                    v1, v2 = val.split('>')
                    g.node[atom].setdefault('s_%s' % key[3:], {})[int(population)] = v1
                    g.node[atom].setdefault('p_%s' % key[3:], {})[int(population)] = v2


class CGRWrite:
    def __init__(self):
        self.__cgr = self.__tocgr()

    def getformattedcgr(self, g, flushmap=False):
        data = dict(atoms=[], bonds=[], CGR_DAT=[], extended=[], meta=g.graph.get('meta', {}).copy(), colors={})
        renum = {}
        colors = {}
        for n, (i, j) in enumerate(g.nodes(data=True), start=1):
            renum[i] = n
            meta = self.__charge(j['s_charge'], j['p_charge'])
            if meta:
                meta['atoms'] = (n,)
                data['CGR_DAT'].append(meta)

            meta = self.__getstate(j.get('s_stereo'), j.get('p_stereo'), 'atomstereo', 'dynatomstereo')
            if meta:
                meta['atoms'] = (n,)
                data['CGR_DAT'].append(meta)

            meta = self.__getstate(j.get('s_hyb'), j.get('p_hyb'), 'atomhyb', 'dynatomhyb')
            if meta:
                meta['atoms'] = (n,)
                data['CGR_DAT'].append(meta)

            meta = self.__getstate(j.get('s_neighbors'), j.get('p_neighbors'), 'atomneighbors', 'dynatomneighbors')
            if meta:
                meta['atoms'] = (n,)
                data['CGR_DAT'].append(meta)

            if 'isotop' in j:
                data['extended'].append(dict(atoms=(n,), value=j.get('isotop'), type='isotop'))

            for k in ('PHTYP', 'FFTYP', 'PCTYP', 'EPTYP', 'HBONDCHG', 'CNECHG'):
                for part, s_val in j.get('s_%s' % k, {}).items():
                    p_val = j['p_%s' % k][part]
                    meta = self.__getstate(s_val, p_val, k, 'dyn%s' % k)
                    if meta:
                        colors.setdefault(meta['type'], {}).setdefault(part, []).append('%d:%s' % (n, meta['value']))

            data['atoms'].append(dict(map=j['mark'] if flushmap else i, charge=j['s_charge'],
                                      element=j.get('element', 'A'), x=j['x'], y=j['y'], z=j['z'], mark=j['mark']))

        for i, j in colors.items():
            data['colors'][i] = '\n'.join('%s %s' % (x, ' '.join(y)) for x, y in j.items())

        for i, l, j in g.edges(data=True):
            bond, cbond, btype = self.__cgr[j.get('s_bond')][j.get('p_bond')]
            if btype:
                data['CGR_DAT'].append({'value': cbond, 'type': btype, 'atoms': (renum[i], renum[l])})

            data['bonds'].append((renum[i], renum[l], bond))

            meta = self.__getstate(j.get('s_stereo'), j.get('p_stereo'), 'bondstereo', 'dynbondstereo')
            if meta:
                meta['atoms'] = (renum[i], renum[l])
                data['CGR_DAT'].append(meta)

        return data

    def getformattedtext(self, data):
        text = []
        for i in data['extended']:
            if i['type'] == 'isotop':
                text.append('M  ISO  1 %3d %3d' % (i['atoms'][0], i['value']))
            elif i['type'] == 'atomlist':
                pass

        for j in count():
            sty = data['CGR_DAT'][j * 8:j * 8 + 8]
            if sty:
                stydat = ' '.join(['%3d DAT' % (x + 1 + j * 8) for x in range(len(sty))])
                text.append('M  STY  %d %s\n' % (len(sty), stydat))
            else:
                break
        for i, j in enumerate(data['CGR_DAT'], start=1):
            cx, cy = self.__getposition(j['atoms'], data['atoms'])
            text.append('M  SAL %3d%3d %s\n' % (i, len(j['atoms']), ' '.join(['%3d' % x for x in j['atoms']])))
            text.append('M  SDT %3d %s\n' % (i, j['type']))
            text.append('M  SDD %3d %10.4f%10.4f    DAU   ALL  0       0\n' % (i, cx, cy))
            text.append('M  SED %3d %s\n' % (i, j['value']))
        return ''.join(text)

    @staticmethod
    def __getstate(s, p, t1, t2):
        rep = {None: 'n'}
        s = rep.get(s, s)
        p = rep.get(p, p)
        if s == p:
            stereo = s
            stype = t1
        else:
            stereo = '%s>%s' % (s, p)
            stype = t2
        if stereo != 'n':
            return {'value': stereo, 'type': stype}
        else:
            return None

    @staticmethod
    def __charge(s, p):
        ss = fromMDL.get(s)
        pp = fromMDL.get(p)
        diff = pp - ss
        if diff:
            meta = {'value': 'c%+d' % diff, 'type': 'dynatom'}
        else:
            meta = None
        return meta

    @staticmethod
    def __tocgr():
        ways = [None, 0, 1, 2, 3, 4, 9]
        rep = {None: 0}
        cgrdict = defaultdict(dict)
        for x, y in product(ways, repeat=2):
            cgrdict[x][y] = ('8', '%s>%s' % (rep.get(x, x), rep.get(y, y)), 'dynbond') if x != y else \
                ('8', 's', 'extrabond') if x == 9 else (str(rep.get(x, x)), None, False)
        return cgrdict

    @staticmethod
    def __getposition(inp, atoms):
        cord = []
        for i in inp:
            cord.append(atoms[i - 1])
        if len(cord) > 1:
            x = (cord[-1]['x'] + cord[0]['x']) / 2 + .2
            y = (cord[-1]['y'] + cord[0]['y']) / 2
            dy = cord[-1]['y'] - cord[0]['y']
            dx = cord[-1]['x'] - cord[0]['x']
            if dx > 0:
                if dy > 0:
                    y -= .2
                else:
                    y += .2
            elif dx < 0:
                if dy < 0:
                    y -= .2
                else:
                    y += .2
        else:
            x, y = cord[0]['x'] + .25, cord[0]['y']

        return x, y
