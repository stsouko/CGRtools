# -*- coding: utf-8 -*-
#
#  Copyright 2017-2019 Ramil Nugmanov <stsouko@live.ru>
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
from collections import defaultdict
from hashlib import sha512
from itertools import count
from ..attributes import Atom, DynAtom
from ..cache import cached_method


hybridization_str = {4: 'a', 3: 't', 2: 'd', 1: 's', None: 'n'}
multiplicity_str = {1: '*', 2: '*2', 3: '*3', None: 'n'}
charge_str = {-3: '-3', -2: '-2', -1: '-', 0: '0', 1: '+', 2: '+2', 3: '+3'}
order_str = {1: '-', 2: '=', 3: '#', 4: ':', 5: '~', None: '.'}
stereo_str = {1: '@', -1: '@@'}


class StringCommon:
    @cached_method
    def __bytes__(self):
        sha512(str(self).encode()).digest()

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(self) == str(other)

    def _flatten(self, weights):
        node = self._node
        adj = self._adj
        cycles = BroodMother()
        atoms_set = set(self._node)
        string = []

        while True:
            start = min(atoms_set, key=weights)
            visited, edges, disconnected = self.__dfs(start, weights, len(atoms_set))
            atoms_set.difference_update(visited)
            reverse = defaultdict(list)
            atoms = defaultdict(list)
            for parent, children in disconnected.items():
                atoms[parent].append(node[parent])
                for child in children:
                    bond = (adj[parent][child], next(cycles))
                    atoms[parent].append(bond)
                    reverse[child].append(bond)

            for parent, children in reverse.items():
                if parent not in atoms:
                    atoms[parent].append(node[parent])
                atoms[parent].extend(children)

            for i in visited - atoms.keys():
                atoms[i] = node[i]

            stack = [[start, [atoms[start]], 0]]
            while True:
                tail, smiles, closure = stack[-1]
                if tail in edges:
                    children = edges[tail]
                    if len(children) > 1:
                        child = children[-1]
                        stack_len = len(stack)
                        stack.append([child, [adj[tail][child], atoms[child]], 0])
                        for child in children[:-1]:
                            stack.append([child, ['(', adj[tail][child], atoms[child]], stack_len])
                    else:
                        child = children[-1]
                        stack[-1][0] = child
                        smiles.append(adj[tail][child])
                        smiles.append(atoms[child])
                elif closure:
                    stack.pop()
                    smiles.append(')')
                    stack[closure - 1][1].extend(smiles)
                elif len(stack) > 2:
                    stack.pop()
                    stack[-1][0] = tail
                    stack[-1][1].extend(smiles)
                elif len(stack) == 2:
                    stack[0][1].extend(smiles)
                    smiles = stack[0][1]
                    break
                else:
                    break
            string.extend(smiles)
            if atoms_set:
                string.append('.')
            else:
                break
        return string

    def __dfs(self, start, weights, depth_limit):
        """
        modified NX dfs
        """
        adj = self._adj

        stack = [(start, depth_limit, iter(sorted(adj[start], key=weights)))]
        visited = {start}
        disconnected = defaultdict(list)
        edges = defaultdict(list)

        while stack:
            parent, depth_now, children = stack[-1]
            try:
                child = next(children)
            except StopIteration:
                stack.pop()
            else:
                if child not in visited:
                    edges[parent].append(child)
                    visited.add(child)
                    if depth_now > 1:
                        front = adj[child].keys() - {parent}
                        if front:
                            stack.append((child, depth_now - 1, iter(sorted(front, key=weights))))
                elif child not in disconnected:
                    disconnected[parent].append(child)

        return visited, edges, disconnected


class SMILES(StringCommon):
    @cached_method
    def __str__(self):
        smiles = []
        for x in self._flatten(self.atoms_order.__getitem__):
            if isinstance(x, str):
                smiles.append(x)
            elif isinstance(x, list):
                smiles.append(self.__format_atom(x[0]))
                for b, c in sorted(x[1:], key=lambda e: int(e[1])):
                    smiles.append(order_str[b.order])
                    smiles.append(str(c))
            elif isinstance(x, Atom):
                smiles.append(self.__format_atom(x))
            else:
                smiles.append(order_str[x.order])
        return ''.join(smiles)

    @staticmethod
    def __format_atom(atom):
        if atom.isotope != atom.common_isotope:
            smi = [str(atom.isotope), atom.element]
        else:
            smi = [atom.element]

        if atom.charge:
            smi.append(charge_str[atom.charge])
        if atom.multiplicity:
            smi.append(multiplicity_str[atom.multiplicity])

        if len(smi) != 1 or atom.element not in ('C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B'):
            smi.insert(0, '[')
            smi.append(']')

        return ''.join(smi)


class SMILES_CGR(StringCommon):
    @cached_method
    def __str__(self):
        smiles = []
        p_smiles = []
        for x in self._flatten(self.atoms_order.__getitem__):
            if isinstance(x, str):
                smiles.append(x)
                p_smiles.append(x)
            elif isinstance(x, list):
                a, p_a = self.__format_atom(x[0])
                smiles.append(a)
                p_smiles.append(p_a)
                for b, c in sorted(x[1:], key=lambda e: int(e[1])):
                    smiles.append(order_str[b.order])
                    smiles.append(str(c))
                    p_smiles.append(order_str[b.p_order])
                    p_smiles.append(str(c))
            elif isinstance(x, DynAtom):
                a, p_a = self.__format_atom(x)
                smiles.append(a)
                p_smiles.append(p_a)
            else:
                smiles.append(order_str[x.order])
                p_smiles.append(order_str[x.p_order])
        return f'{"".join(smiles)}>>{"".join(p_smiles)}'

    @staticmethod
    def __format_atom(atom):
        if atom.isotope != atom.common_isotope:
            smi = [str(atom.isotope), atom.element]
            p_smi = [str(atom.isotope), atom.element]
        else:
            smi = [atom.element]
            p_smi = [atom.element]

        if atom.charge:
            smi.append(charge_str[atom.charge])
        if atom.p_charge:
            p_smi.append(charge_str[atom.p_charge])
        if atom.multiplicity:
            smi.append(multiplicity_str[atom.multiplicity])
        if atom.p_multiplicity:
            p_smi.append(multiplicity_str[atom.p_multiplicity])

        if atom.element not in {'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B'}:
            smi.insert(0, '[')
            smi.append(']')
            p_smi.insert(0, '[')
            p_smi.append(']')
        else:
            if len(smi) != 1:
                smi.insert(0, '[')
                smi.append(']')
            if len(p_smi) != 1:
                p_smi.insert(0, '[')
                p_smi.append(']')

        return ''.join(smi), ''.join(p_smi)


class BroodMother:
    __slots__ = ('__count', '__map')

    def __init__(self, start=1, step=1):
        self.__count = count(start, step)
        self.__map = {}

    def __next__(self):
        return Brood(self)

    def __getitem__(self, brood):
        try:
            return self.__map[brood]
        except KeyError:
            return self.__map.setdefault(brood, next(self.__count))


class Brood:
    __slots__ = '__mother'

    def __init__(self, mother):
        self.__mother = mother

    def __str__(self):
        return str(self.__mother[self])

    def __int__(self):
        return self.__mother[self]


'''

        if atom.hybridization:
            smi.append(';')
            smi.append(hybridization_str[atom.hybridization])
            if atom.neighbors is not None:
                smi.append(str(atom.neighbors))
            smi.append(';')
        elif atom.neighbors is not None:
            smi.append(';')
            smi.append(str(atom.neighbors))
            smi.append(';')


    def stringify(self, atom=True, isotope=True, stereo=True, hybridization=True, neighbors=True):
        smi = []
        if stereo and self.stereo:
            smi.append(self._stereo_str[self.stereo])
        if hybridization:
            if len(self.hybridization) > 1:
                smi.append('<%s>' % ''.join(self._hybridization_str[x] for x in self.hybridization))
            elif self.hybridization:
                smi.append(self._hybridization_str[self.hybridization[0]])
        if neighbors:
            if len(self.neighbors) > 1:
                smi.append('<%s>' % ''.join(str(x) for x in self.neighbors))
            elif self.neighbors:
                smi.append(str(self.neighbors[0]))
        if smi:
            smi.append(';')
            smi.insert(0, ';')

        if atom:
            if self.element == ('A',):
                atom = False
                smi.insert(0, '*')
            elif len(self.element) > 1:
                atom = False
                smi.insert(0, ','.join(self.element))
            else:
                if self.element[0] not in ('C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B'):
                    atom = False
                smi.insert(0, self.element[0])

            if len(self.charge) > 1:
                smi.append('<%s>' % ''.join(self._charge_str[x] for x in self.charge))
            elif self.charge != (0,):
                smi.append(self._charge_str[self.charge[0]])

            if len(self.multiplicity) > 1:
                smi.append('<%s>' % ''.join(self._multiplicity_str[x] for x in self.multiplicity))
            elif self.multiplicity:
                smi.append(self._multiplicity_str[self.multiplicity[0]])

            if isotope:
                if len(self.isotope) > 1:
                    smi.insert(0, '<%s>' % ','.join(str(x) for x in self.isotope))
                elif self.isotope:
                    smi.insert(0, str(self.isotope[0]))
        else:
            smi.insert(0, '*')

        if len(smi) != 1 or not atom:
            smi.insert(0, '[')
            smi.append(']')

        return ''.join(smi)
'''

__all__ = ['SMILES_CGR', 'SMILES']
