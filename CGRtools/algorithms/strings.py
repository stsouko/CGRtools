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
        self.__get_atom = self.__cgr_atom if is_cgr else self.__mol_atom
        self.__add_bond = self.__cgr_bond if is_cgr else self.__mol_bond
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

    def __mol_atom(self, atom):
        smi = []

        if self.__stereo and atom.stereo:
            smi.append(self.__stereo_types[atom.stereo])

        if self.__hyb and atom.hybridization:
            smi.append(self.__hyb_types[atom.hybridization])

        if self.__neighbors and atom.neighbors is not None:
            smi.append(str(atom.neighbors))

        if smi:
            smi.append(';')
            smi.insert(0, ';')

        if self.__element:
            smi.insert(0, atom.symbol)

            if atom.charge:
                smi.append(self.__charge_to_string(atom.charge))

            if atom.multiplicity:
                smi.append(self.__radical_map[atom.multiplicity])
        else:
            smi.insert(0, '*')

        if self.__isotope and atom.isotope:
            smi.insert(0, str(atom.isotope))

        if len(smi) != 1 or smi[0] == '*':
            smi.insert(0, '[')
            smi.append(']')

        return smi

    def __cgr_atom(self, atom):
        smi = []
        pmi = []

        if self.__stereo:
            if atom.get('s_stereo'):
                smi.append(self.__stereo_types[atom['s_stereo']])
            if atom.get('s_stereo'):
                pmi.append(self.__stereo_types[atom['p_stereo']])

        if self.__hyb:
            s = atom.get('s_hyb')
            p = atom.get('p_hyb')
            if isinstance(s, list):
                tmp = sorted((self.__hyb_types[x], self.__hyb_types[y]) for x, y in zip(s, p))
                smi.append('<%s>' % ''.join(x for x, _ in tmp))
                pmi.append('<%s>' % ''.join(x for _, x in tmp))
            else:
                if s:
                    smi.append(self.__hyb_types[s])
                if p:
                    pmi.append(self.__hyb_types[p])

        if self.__neighbors:
            s = atom.get('s_neighbors')
            p = atom.get('p_neighbors')
            if isinstance(s, list):
                tmp = sorted((x is None and 'n' or str(x), y is None and 'n' or str(y)) for x, y in zip(s, p))
                smi.append('<%s>' % ''.join(x for x, _ in tmp))
                pmi.append('<%s>' % ''.join(x for _, x in tmp))
            else:
                if s is not None:
                    smi.append(str(s))
                if p is not None:
                    pmi.append(str(p))

        if smi:
            smi.append(';')
            smi.insert(0, ';')
        if pmi:
            pmi.append(';')
            pmi.insert(0, ';')

        if self.__element:
            ge = atom.get('element')
            if isinstance(ge, list):
                ge = list(','.join(sorted(ge)))
                smi = ge + smi
                pmi = ge + pmi
            else:
                if not ge:
                    ge = '*'
                smi.insert(0, ge)
                pmi.insert(0, ge)

            s = atom.get('s_charge')
            p = atom.get('p_charge')
            if isinstance(s, list):
                tmp = [(self.__charge_to_string(x), self.__charge_to_string(y)) for x, y in sorted(zip(s, p))]
                smi.append('<%s>' % ''.join(x for x, _ in tmp))
                pmi.append('<%s>' % ''.join(x for _, x in tmp))
            else:
                if s:
                    smi.append(self.__charge_to_string(s))
                if p:
                    pmi.append(self.__charge_to_string(p))

            s = atom.get('s_radical')
            p = atom.get('p_radical')
            if isinstance(s, list):
                tmp = [(self.__radical_map[x], self.__radical_map[y]) for x, y in sorted(zip(s, p))]
                smi.append('<%s>' % ''.join(x for x, _ in tmp))
                pmi.append('<%s>' % ''.join(x for _, x in tmp))
            else:
                if s:
                    smi.append(self.__radical_map[s])
                if p:
                    pmi.append(self.__radical_map[p])
        else:
            smi.insert(0, '*')
            pmi.insert(0, '*')

        if self.__isotope:
            s = atom.get('isotope')
            if isinstance(s, list):
                s = '<%s>' % ','.join(str(x) for x in s)
                smi.insert(0, s)
                pmi.insert(0, s)
            elif s:
                s = str(s)
                smi.insert(0, s)
                pmi.insert(0, s)

        if len(smi) != 1 or smi[0] == '*':
            smi.insert(0, '[')
            smi.append(']')
        if len(pmi) != 1 or pmi[0] == '*':
            pmi.insert(0, '[')
            pmi.append(']')

        return smi, pmi

    def __mol_bond(self, smiles, bond, closure=''):
        if isinstance(bond.order, list):
            sb = '<%s>' % ''.join(sorted(self.__to_smiles_map[x] for x in bond.order))
        else:
            sb = self.__to_smiles_map[bond.order]
        smiles.append('%s%s%s' % (sb, self.__stereo and bond.stereo or '', closure))

    def __cgr_bond(self, smiles, bond, closure=''):
        if isinstance(bond.order, list):
            tmp = sorted((self.__to_smiles_map[x], self.__to_smiles_map[y]) for x, y in zip(bond.order, bond.p_order))
            sb, pb = '<%s>' % ''.join(x for x, _ in tmp), '<%s>' % ''.join(x for _, x in tmp)
        else:
            sb, pb = self.__to_smiles_map[bond.order], self.__to_smiles_map[bond.p_order]

        smiles[0].append('%s%s%s' % (pb, self.__stereo and bond.stereo or '', closure))
        smiles[1].append('%s%s%s' % (pb, self.__stereo and bond.p_stereo or '', closure))

    def __do_cgr_smarts(self, trace, inter, prev):
        g = self.__g
        countcyc = self.__countcyc

        smiles = self.__get_atom(g._node[inter])
        concat = []
        stoplist = []
        iterlist = set(g._adj[inter]).difference([prev])
        while iterlist:
            i = self.__get_next_atom(iterlist)
            iterlist.discard(i)
            gii = g._adj[inter][i]
            if i in trace:
                if i not in stoplist:  # костыль для циклов. чтоб не было 2х проходов.
                    cyc = next(countcyc)
                    concat.append((i, cyc, inter))
                    self.__add_bond(smiles, gii, cyc)
                continue

            deep0, deep1, deep2 = self.__do_cgr_smarts(set(chain(trace, [i])), i, inter)
            trace.update(deep0)
            if deep2:
                for j0, j1, j2 in deep2:
                    if j0 == inter:
                        gij = g[inter][j2]
                        stoplist.append(j2)
                        self.__add_bond(smiles, gij, j1)
                    else:
                        concat.append((j0, j1, j2))

            if self.__is_cgr:
                if iterlist:
                    smiles[0].append('(')
                    smiles[1].append('(')
                self.__add_bond(smiles, gii)
                smiles[0].extend(deep1[0])
                smiles[1].extend(deep1[1])
                if iterlist:
                    smiles[0].append(')')
                    smiles[1].append(')')
            else:
                if iterlist:
                    smiles.append('(')
                self.__add_bond(smiles, gii)
                smiles.extend(deep1)
                if iterlist:
                    smiles.append(')')

        return trace, smiles, concat

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
