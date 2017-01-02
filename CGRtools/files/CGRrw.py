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
from itertools import count, product, chain
from collections import defaultdict
import periodictable as pt

toMDL = {-3: 7, -2: 6, -1: 5, 0: 0, 1: 3, 2: 2, 3: 1}
fromMDL = {0: 0, 1: 3, 2: 2, 3: 1, 4: 0, 5: -1, 6: -2, 7: -3}
bondlabels = {'0': None, '1': 1, '2': 2, '3': 3, '4': 4, '9': 9, 'n': None, 's': 9}
mendeleyset = set(x.symbol for x in pt.elements)


class CGRread:
    __prop = {}

    def collect(self, line):
        if line.startswith('M  ALS'):
            self.__prop[line[3:10]] = dict(atoms=[int(line[7:10])],
                                           type='atomlist' if line[14] == 'F' else 'atomnotlist',
                                           value=[line[16 + x*4: 20 + x*4].strip() for x in range(int(line[10:13]))])
        elif line.startswith('M  ISO'):
            self.__prop[line[3:10]] = dict(atoms=[int(line[10:13])], type='isotop', value=line[14:17].strip())

        elif line.startswith('M  STY'):
            for i in range(int(line[8])):
                if 'DAT' in line[10 + 8 * i:17 + 8 * i]:
                    self.__prop[int(line[10 + 8 * i:13 + 8 * i])] = {}
        elif line.startswith('M  SAL') and int(line[7:10]) in self.__prop:
            key = []
            for i in range(int(line[10:13])):
                key.append(int(line[14 + 4 * i:17 + 4 * i]))
            self.__prop[int(line[7:10])]['atoms'] = sorted(key)
        elif line.startswith('M  SDT') and int(line[7:10]) in self.__prop:
            key = line.split()[-1].lower()
            if key not in self.__cgrkeys:
                self.__prop.pop(int(line[7:10]))
            else:
                self.__prop[int(line[7:10])]['type'] = key
        elif line.startswith('M  SED') and int(line[7:10]) in self.__prop:
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

        self.__prop.clear()
        return prop

    @staticmethod
    def cgr_dat(g, k, atom1, atom2):

        def _parsedyn(target, name):
            tmp = ([], [])
            for x in k['value'].split(','):
                s, *_, p = x.split('>')
                if s != p:
                    tmp[0].append((s != 'n' and int(s) or None) if name != 'bond' else bondlabels[s])
                    tmp[1].append((p != 'n' and int(p) or None) if name != 'bond' else bondlabels[p])

            if len(tmp[0]) == 1:
                target['sp_%s' % name] = tmp[0][0], tmp[1][0]
                for sp_key, sp_ind in (('s', 0), ('p', 1)):
                    if tmp[sp_ind][0] is not None:
                        target['%s_%s' % (sp_key, name)] = tmp[sp_ind][0]
                    else:
                        target.pop('%s_%s' % (sp_key, name), None)
            elif len(tmp[0]) > 1:
                tmp2 = list(zip(*tmp))
                if len(set(tmp2)) == len(tmp2):
                    target['sp_%s' % name] = tmp2
                    target['s_%s' % name] = tmp[0]
                    target['p_%s' % name] = tmp[1]

        def _parselist(target, name):
            tmp = [(x != 'n' and int(x) or None) if name != 'bond' else bondlabels[x] for x in k['value'].split(',')]
            if len(tmp) == 1:
                if tmp[0] is not None:
                    target['s_%s' % name] = target['p_%s' % name] = target['sp_%s' % name] = tmp[0]
            elif len(set(tmp)) == len(tmp):
                target['s_%s' % name] = target['p_%s' % name] = target['sp_%s' % name] = tmp

        if k['type'] == 'dynatomstereo':
            _parsedyn(g.node[atom1], 'stereo')

        elif k['type'] == 'atomstereo':
            _parselist(g.node[atom1], 'stereo')

        elif k['type'] == 'dynatom':
            key = k['value'][0]
            diff = [int(x) for x in k['value'][1:].split(',')]
            if key == 'c':  # update atom charges from CGR
                base = g.node[atom1]['s_charge']

                if len(diff) > 1:
                    if diff[0] == 0:  # for list of charges c0,1,-2,+3...
                        tmp = [base + x for x in diff]
                        if len(set(tmp)) == len(tmp):
                            g.node[atom1]['s_charge'] = g.node[atom1]['p_charge'] = g.node[atom1]['sp_charge'] = tmp
                    elif len(diff) % 2 == 1:  # for dyn charges c1,-1,0... -1,0 is dyn group relatively to base
                        s = [base] + [base + x for x in diff[1::2]]
                        p = [base + x for x in diff[::2]]
                        tmp = list(zip(s, p))
                        if len(set(tmp)) == len(tmp):
                            g.node[atom1]['sp_charge'] = tmp
                            g.node[atom1]['s_charge'] = s
                            g.node[atom1]['p_charge'] = p
                else:
                    tmp = diff[0]
                    if tmp:
                        g.node[atom1]['p_charge'] = base + tmp
                        g.node[atom1]['sp_charge'] = (g.node[atom1]['s_charge'], g.node[atom1]['p_charge'])

            else:
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
            g.node[atom1]['element'] = list(mendeleyset.difference(k['value']))

        elif k['type'] == 'atomhyb':
            _parselist(g.node[atom1], 'hyb')

        elif k['type'] == 'atomneighbors':
            _parselist(g.node[atom1], 'neighbors')

        elif k['type'] == 'dynatomhyb':
            _parsedyn(g.node[atom1], 'hyb')

        elif k['type'] == 'dynatomneighbors':
            _parsedyn(g.node[atom1], 'neighbors')

        elif k['type'] == 'isotop':
            tmp = k['value'].split(',')
            g.node[atom1]['isotop'] = [int(x) for x in tmp] if len(tmp) > 1 else int(tmp[0])

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


class CGRwrite:
    def __init__(self, extralabels=False, mark_to_map=False):
        self.__cgr = self.__tocgr()
        self.__mark_to_map = mark_to_map
        tmp = ['stereo'] + (['hyb', 'neighbors'] if extralabels else [])
        self.__atomprop = [('s_%s' % x, 'p_%s' % x, 'sp_%s' % x, 'atom%s' % x, 'dynatom%s' % x) for x in tmp]

    def getformattedcgr(self, g):
        data = dict(meta=g.meta.copy())
        cgr_dat, extended, atoms, bonds = [], [], [], []
        renum, colors = {}, {}
        for n, (i, j) in enumerate(g.nodes(data=True), start=1):
            renum[i] = n
            charge, meta = self.__charge(j.get('s_charge'), j.get('p_charge'))
            if meta:
                meta['atoms'] = (n,)
                cgr_dat.append(meta)

            for s_key, p_key, _, a_key, d_key in self.__atomprop:
                meta = self.__getstate(j.get(s_key), j.get(p_key),  a_key, d_key)
                if meta:
                    meta['atoms'] = (n,)
                    cgr_dat.append(meta)

            if 'isotop' in j:
                extended.append(dict(atoms=(n,), value=j['isotop'], type='isotop'))

            for k in ('PHTYP', 'FFTYP', 'PCTYP', 'EPTYP', 'HBONDCHG', 'CNECHG'):
                for part, s_val in j.get('s_%s' % k, {}).items():
                    p_val = j['p_%s' % k][part]
                    meta = self.__getstate(s_val, p_val, k, 'dyn%s' % k)
                    if meta:
                        colors.setdefault(meta['type'], {}).setdefault(part, []).append('%d:%s' % (n, meta['value']))

            if isinstance(j.get('element'), list):
                extended.append(dict(atoms=(n,), value=j['element'], type='atomlist'))
                j['element'] = 'L'

            atoms.append(dict(map=j['mark'] if self.__mark_to_map else i, charge=charge,
                              element=j.get('element', 'A'), x=j['x'], y=j['y'], z=j['z'], mark=j['mark']))

        data['colors'] = {i: '\n'.join('%s %s' % (x, ' '.join(y)) for x, y in j.items()) for i, j in colors.items()}

        for i, l, j in g.edges(data=True):
            bond, cbond, btype = self.__cgr[j.get('s_bond')][j.get('p_bond')]
            if btype:
                cgr_dat.append({'value': cbond, 'type': btype, 'atoms': (renum[i], renum[l])})

            bonds.append((renum[i], renum[l], bond))

            meta = self.__getstate(j.get('s_stereo'), j.get('p_stereo'), 'bondstereo', 'dynbondstereo')
            if meta:
                meta['atoms'] = (renum[i], renum[l])
                cgr_dat.append(meta)

        data['CGR'] = ''.join(chain(("\n  CGRtools. (c) Dr. Ramil I. Nugmanov\n"
                                     "\n%3s%3s  0  0  0  0            999 V2000\n" % (len(atoms), len(bonds)),),
                                    ("%(x)10.4f%(y)10.4f%(z)10.4f %(element)-3s 0%(charge)3s  0  0  0  0  0%(mark)3s  "
                                     "0%(map)3s  0  0\n" % i for i in atoms),
                                    ("%3d%3d%3s  0  0  0  0\n" % i for i in bonds),
                                    self.__formatprop(extended, cgr_dat, atoms)))
        return data

    def __formatprop(self, extended, cgr_dat, atoms):
        text = []
        for i in extended:
            if i['type'] == 'isotop':
                if isinstance(i['value'], list):
                    i['value'] = ','.join(str(x) for x in i['value'])
                    cgr_dat.insert(0, i)
                else:
                    text.append('M  ISO  1 %3d %3d\n' % (i['atoms'][0], i['value']))
            elif i['type'] == 'atomlist':
                atomslist, _type = (mendeleyset.difference(i['value']), 'T') \
                    if len(i['value']) > len(mendeleyset) / 2 else (i['value'], 'F')

                text.append('M  ALS %3d%3d %s %s\n' % (i['atoms'][0], len(atomslist), _type,
                                                       ''.join('%-4s' % x for x in atomslist)))

        for j in count():
            sty = cgr_dat[j * 8:j * 8 + 8]
            if sty:
                stydat = ' '.join(['%3d DAT' % (x + 1 + j * 8) for x in range(len(sty))])
                text.append('M  STY  %d %s\n' % (len(sty), stydat))
            else:
                break

        for i, j in enumerate(cgr_dat, start=1):
            cx, cy = self.__getposition([atoms[i - 1] for i in j['atoms']])
            text.append('M  SAL %3d%3d %s\n' % (i, len(j['atoms']), ' '.join(['%3d' % x for x in j['atoms']])))
            text.append('M  SDT %3d %s\n' % (i, j['type']))
            text.append('M  SDD %3d %10.4f%10.4f    DAU   ALL  0       0\n' % (i, cx, cy))
            text.append('M  SED %3d %s\n' % (i, j['value']))
        return text

    @staticmethod
    def __getstate(s, p, t1, t2):
        if isinstance(s, list):
            if s == p:
                _val = ','.join(x and str(x) or 'n' for x in s)
                _type = t1
            else:
                _val = ','.join('%s>%s' % (x or 'n', y or 'n') for x, y in zip(s, p) if x != y)
                _type = t2

        elif s != p:
            _val = '%s>%s' % (s or 'n', p or 'n')
            _type = t2
        else:
            _val = s
            _type = t1

        return _val and {'value': _val, 'type': _type}

    @staticmethod
    def __charge(s, p):
        charge = s
        if isinstance(s, list):
            charge = s[0]
            if s == p:
                meta = dict(type='dynatom', value='c0%s' % ','.join(str(x) for x in s[1:])) if len(s) > 1 else None
            elif s[0] != p[0]:
                meta = dict(type='dynatom', value='c%s' % ','.join(chain(['%+d' % (p[0] - charge)],
                                                                         ('%+d,%+d' % (x - charge, y - charge)
                                                                          for x, y in zip(s[1:], p[1:]) if x != y))))
            else:
                meta = None
        elif p != s:
            meta = dict(value='c%+d' % (p - s), type='dynatom')
        else:
            meta = None
        return toMDL.get(charge, 0), meta

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
    def __getposition(cord):
        if len(cord) > 1:
            x = (cord[-1]['x'] + cord[0]['x']) / 2 + .2
            y = (cord[-1]['y'] + cord[0]['y']) / 2
            dy = cord[-1]['y'] - cord[0]['y']
            dx = cord[-1]['x'] - cord[0]['x']
            if dx > 0:
                y += -.2 if dy > 0 else .2
            elif dx < 0:
                y += -.2 if dy < 0 else .2
        else:
            x, y = cord[0]['x'] + .25, cord[0]['y']

        return x, y
