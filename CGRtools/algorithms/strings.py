# -*- coding: utf-8 -*-
#
#  Copyright 2017-2019 Ramil Nugmanov <stsouko@live.ru>
#  Copyright 2019 Timur Gimadiev <timur.gimadiev@gmail.com>
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
from ..attributes import Atom, DynAtom, QueryAtom, DynQueryAtom, DynBond
from ..cache import cached_method, cached_args_method


hybridization_str = {4: 'a', 3: 't', 2: 'd', 1: 's', None: 'n'}
multiplicity_str = {1: '*', 2: '*2', 3: '*3', None: 'n'}
charge_str = {-3: '-3', -2: '-2', -1: '-', 0: '0', 1: '+', 2: '+2', 3: '+3'}
order_str = {1: '-', 2: '=', 3: '#', 4: ':', 5: '~', None: '.'}
stereo_str = {1: '@', -1: '@@'}
dyn_order_str = {(None, 1): "[.>-]", (None, 2): "[.>=]", (None, 3): "[.>#]", (None, 4): "[.>:]", (None, 5): "[.>~]",
                 (1, None): "[->.]", (1, 1): "", (1, 2): "[->=]", (1, 3): "[->#]", (1, 4): "[->:]", (1, 5): "[->~]",
                 (2, None): "[=>.]", (2, 1): "[=>-]", (2, 2): "=", (2, 3): "[=>#]", (2, 4): "[=>:]", (2, 5): "[=>~]",
                 (3, None): "[#>.]", (3, 1): "[#>-]", (3, 2): "[#>=]", (3, 3): "#", (3, 4): "[#>:]", (3, 5): "[#>~]",
                 (4, None): "[:>.]", (4, 1): "[:>-]", (4, 2): "[:>=]", (4, 3): "[:>#]", (4, 4): ":", (4, 5): "[:>~]",
                 (5, None): "[~>.]", (5, 1): "[~>-]", (5, 2): "[~>=]", (5, 3): "[~>#]", (5, 4): "[~>:]", (5, 5): "[~]"}
dyn_charge_str = {(-3, -3): "-3", (-3, -2): "-3>-2", (-3, -1): "-3>-", (-3, 0): "-3>0", (-3, 1): "-3>+",
                  (-3, 2): "-3>+2", (-3, 3): "-3>3", (-2, -3): "-2>-3", (-2, -2): "-2", (-2, -1): "-2>-",
                  (-2, 0): "-2>0", (-2, 1): "-2>+", (-2, 2): "-2>+2", (-2, 3): "-2>+3", (-1, -3): "->-3",
                  (-1, -2): "->-2", (-1, -1): "-", (-1, 0): "->0", (-1, 1): "->+", (-1, 2): "->+2", (-1, 3): "->+3",
                  (0, -3): "0>-3", (0, -2): "0>-2", (0, -1): "0>-", (0, 0): "", (0, 1): "0>+", (0, 2): "0>+2",
                  (0, 3): "0>+3", (1, -3): "+>-3", (1, -2): "+>-2", (1, -1): "+>-", (1, 0): "+>0", (1, 1): "+",
                  (1, 2): "+>+2", (1, 3): "+>+3", (2, -3): "+2>-3", (2, -2): "+2>-2", (2, -1): "+2>-", (2, 0): "+2>0",
                  (2, 1): "+2>+", (2, 2): "+2", (2, 3): "+2>+3", (3, -3): "+3>-3", (3, -2): "+3>-2", (3, -1): "+3>-",
                  (3, 0): "+3>0", (3, 1): "+3>+", (3, 2): "+3>+2", (3, 3): "+3"}
dyn_multiplicity_str = {(1, 1): "*", (1, 2): "*>*2", (1, 3): "*>*3", (1, None): "*>n", (2, 1): "*2>*", (2, 2): "*2",
                        (2, 3): "*2>*3", (2, None): "*2>n", (3, 1): "*3>*", (3, 2): "*3>*2", (3, 3): "*3",
                        (3, None): "*3>n", (None, 1): "n>1", (None, 2): "n>2", (None, 3): "n>3", (None, None): ""}


class HashableSmiles:
    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(self) == str(other)


class StringCommon:
    @cached_method
    def __bytes__(self):
        return sha512(str(self).encode()).digest()

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


class Smiles(StringCommon, HashableSmiles):
    @cached_method
    def __str__(self):
        return format(self)

    @cached_args_method
    def __format__(self, format_spec):
        """
        format molecule as SMILES string

        :param format_spec: if == 'n' add neighbors count of atoms. don't forget to call reset query marks before.
        if == 'h' add hybridizations of atoms. if 'nh' or 'hn' add both.
        """
        if not format_spec:
            neighbors = False
            hybridization = False
        elif format_spec == 'n':
            neighbors = True
            hybridization = False
        elif format_spec == 'h':
            neighbors = False
            hybridization = True
        elif format_spec in ('hn', 'nh'):
            neighbors = True
            hybridization = True
        else:
            raise ValueError('invalid format_spec')
        return self._format_string(self.atoms_order.__getitem__, neighbors, hybridization)

    def _format_string(self, order, neighbors, hybridization):
        smiles = []
        for x in self._flatten(order):
            if isinstance(x, str):
                smiles.append(x)
            elif isinstance(x, list):
                smiles.append(self.__format_atom(x[0], neighbors, hybridization))
                for b, c in sorted(x[1:], key=lambda e: int(e[1])):
                    smiles.append(order_str[b.order])
                    smiles.append(str(c))
            elif isinstance(x, Atom):
                smiles.append(self.__format_atom(x, neighbors, hybridization))
            else:
                smiles.append(order_str[x.order])
        return ''.join(smiles)

    @staticmethod
    def __format_atom(atom, neighbors, hybridization):
        if atom.isotope != atom.common_isotope:
            smi = [str(atom.isotope), atom.element]
        else:
            smi = [atom.element]

        if hybridization:
            smi.append(';')
            smi.append(hybridization_str[atom.hybridization])
            if neighbors:
                smi.append(str(atom.neighbors))
            smi.append(';')
        elif neighbors:
            smi.append(';')
            smi.append(str(atom.neighbors))
            smi.append(';')

        if atom.charge:
            smi.append(charge_str[atom.charge])
        if atom.multiplicity:
            smi.append(multiplicity_str[atom.multiplicity])

        if len(smi) != 1 or atom.element not in {'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B'}:
            smi.insert(0, '[')
            smi.append(']')

        return ''.join(smi)


class SmilesCGR(StringCommon, HashableSmiles):
    @cached_method
    def __str__(self):
        return format(self)

    @cached_args_method
    def __format__(self, format_spec):
        """
        format CGR as SMIRKS string

        :param format_spec: if 's' in fromat_spec only representation of CGR as smiles will be shown.
        No hybridization and neighbors count will be used. High priority option.
        if == 'n' add neighbors count of atoms. don't forget to call reset query marks before.
        if == 'h' add hybridization of atoms. if 'nh' or 'hn' add both.
        """
        if format_spec and 's' in format_spec:
            return self._format_string_cgr(self.atoms_order.__getitem__)
        if not format_spec:
            neighbors = False
            hybridization = False
        elif format_spec == 'n':
            neighbors = True
            hybridization = False
        elif format_spec == 'h':
            neighbors = False
            hybridization = True
        elif format_spec in ('hn', 'nh'):
            neighbors = True
            hybridization = True
        else:
            raise ValueError('invalid format_spec')
        return self._format_string(self.atoms_order.__getitem__, neighbors, hybridization)

    def _format_string(self, order, neighbors, hybridization):
        smiles = []
        p_smiles = []
        for x in self._flatten(order):
            if isinstance(x, str):
                smiles.append(x)
                p_smiles.append(x)
            elif isinstance(x, list):
                a, p_a = self.__format_atom(x[0], neighbors, hybridization)
                smiles.append(a)
                p_smiles.append(p_a)
                for b, c in sorted(x[1:], key=lambda e: int(e[1])):
                    smiles.append(order_str[b.order])
                    smiles.append(str(c))
                    p_smiles.append(order_str[b.p_order])
                    p_smiles.append(str(c))
            elif isinstance(x, DynAtom):
                a, p_a = self.__format_atom(x, neighbors, hybridization)
                smiles.append(a)
                p_smiles.append(p_a)
            else:
                smiles.append(order_str[x.order])
                p_smiles.append(order_str[x.p_order])
        return f'{"".join(smiles)}>>{"".join(p_smiles)}'

    @staticmethod
    def __format_atom(atom, neighbors, hybridization):
        if atom.isotope != atom.common_isotope:
            smi = [str(atom.isotope), atom.element]
            p_smi = [str(atom.isotope), atom.element]
        else:
            smi = [atom.element]
            p_smi = [atom.element]

        if hybridization:
            smi.append(';')
            smi.append(hybridization_str[atom.hybridization])
            p_smi.append(';')
            p_smi.append(hybridization_str[atom.p_hybridization])
            if neighbors:
                smi.append(str(atom.neighbors))
                p_smi.append(str(atom.p_neighbors))
            smi.append(';')
            p_smi.append(';')
        elif neighbors:
            smi.append(';')
            p_smi.append(';')
            smi.append(str(atom.neighbors))
            p_smi.append(str(atom.p_neighbors))
            smi.append(';')
            p_smi.append(';')

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

    def _format_string_cgr(self, order):
        smiles = []
        for x in self._flatten(order):
            if isinstance(x, str):
                smiles.append(x)
            elif isinstance(x, list):
                smiles.append(self.__format_atom_cgr(x[0]))
                for b, c in sorted(x[1:], key=lambda e: int(e[1])):
                    smiles.append(dyn_order_str[(b.order, b.p_order)])
                    smiles.append(str(c))
            elif isinstance(x, DynAtom):
                smiles.append(self.__format_atom_cgr(x))
            else:
                smiles.append(dyn_order_str[(x.order, x.p_order)])
        return "".join(smiles)

    @staticmethod
    def __format_atom_cgr(atom):
        if atom.isotope != atom.common_isotope:
            smi = [str(atom.isotope), atom.element]
        else:
            smi = [atom.element]
        if atom.charge or atom.p_charge:
            smi.append(dyn_charge_str[(atom.charge, atom.p_charge)])
        if atom.multiplicity or atom.p_multiplicity:
            smi.append(dyn_multiplicity_str[(atom.multiplicity, atom.p_multiplicity)])
        if len(smi) != 1 or atom.element not in {'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B'}:
            smi.insert(0, '[')
            smi.append(']')
        return ''.join(smi)


class SmilesQuery(StringCommon):
    @cached_method
    def __str__(self):
        smiles = []
        for x in self._flatten(lambda x: x):
            if isinstance(x, str):
                smiles.append(x)
            elif isinstance(x, list):
                smiles.append(self.__format_atom(x[0]))
                for b, c in sorted(x[1:], key=lambda e: int(e[1])):
                    smiles.append(order_str[b.order])
                    smiles.append(str(c))
            elif isinstance(x, QueryAtom):
                smiles.append(self.__format_atom(x))
            else:
                smiles.append(order_str[x.order])
        return ''.join(smiles)

    @staticmethod
    def __format_atom(atom):
        if atom.isotope:
            smi = ['[', str(atom.isotope)]
        else:
            smi = ['[']

        if not atom.element:
            smi.append('*')
        elif len(atom.element) > 1:
            smi.append(','.join(atom.element))
        else:
            smi.extend(atom.element)

        if len(atom.hybridization) > 1:
            smi.append(';<%s>' % ''.join(hybridization_str[x] for x in atom.hybridization))
            if len(atom.neighbors) > 1:
                smi.append('<%s>' % ''.join(str(x) for x in atom.neighbors))
            elif atom.neighbors:
                smi.append(str(atom.neighbors[0]))
            smi.append(';')
        elif atom.hybridization:
            smi.append(';')
            smi.append(hybridization_str[atom.hybridization[0]])
            if len(atom.neighbors) > 1:
                smi.append('<%s>' % ''.join(str(x) for x in atom.neighbors))
            elif atom.neighbors:
                smi.append(str(atom.neighbors[0]))
            smi.append(';')
        elif len(atom.neighbors) > 1:
            smi.append(';<%s>;' % ''.join(str(x) for x in atom.neighbors))
        elif atom.neighbors:
            smi.append(';')
            smi.append(str(atom.neighbors[0]))
            smi.append(';')

        if atom.charge:
            smi.append(charge_str[atom.charge])
        if atom.multiplicity:
            smi.append(multiplicity_str[atom.multiplicity])

        smi.append(']')
        return ''.join(smi)


class SmilesQueryCGR(StringCommon):
    @cached_method
    def __str__(self):
        smiles = []
        p_smiles = []
        for x in self._flatten(lambda x: x):
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
            elif isinstance(x, DynQueryAtom):
                a, p_a = self.__format_atom(x)
                smiles.append(a)
                p_smiles.append(p_a)
            else:
                smiles.append(order_str[x.order])
                p_smiles.append(order_str[x.p_order])
        return f'{"".join(smiles)}>>{"".join(p_smiles)}'

    @staticmethod
    def __format_atom(atom):
        if atom.isotope:
            smi = ['[', str(atom.isotope)]
            p_smi = ['[', str(atom.isotope)]
        else:
            smi = ['[']
            p_smi = ['[']

        if not atom.element:
            smi.append('*')
            p_smi.append('*')
        elif len(atom.element) > 1:
            e = ','.join(atom.element)
            smi.append(e)
            p_smi.append(e)
        else:
            smi.extend(atom.element)
            p_smi.extend(atom.element)

        if len(atom.hybridization) > 1:
            smi.append(';<%s>' % ''.join(hybridization_str[x] for x in atom.hybridization))
            if len(atom.neighbors) > 1:
                smi.append('<%s>' % ''.join(str(x) for x in atom.neighbors))
            elif atom.neighbors:
                smi.append(str(atom.neighbors[0]))
            smi.append(';')
        elif atom.hybridization:
            smi.append(';')
            smi.append(hybridization_str[atom.hybridization[0]])
            if len(atom.neighbors) > 1:
                smi.append('<%s>' % ''.join(str(x) for x in atom.neighbors))
            elif atom.neighbors:
                smi.append(str(atom.neighbors[0]))
            smi.append(';')
        elif len(atom.neighbors) > 1:
            smi.append(';<%s>;' % ''.join(str(x) for x in atom.neighbors))
        elif atom.neighbors:
            smi.append(';')
            smi.append(str(atom.neighbors[0]))
            smi.append(';')

        if len(atom.p_hybridization) > 1:
            p_smi.append(';<%s>' % ''.join(hybridization_str[x] for x in atom.p_hybridization))
            if len(atom.p_neighbors) > 1:
                p_smi.append('<%s>' % ''.join(str(x) for x in atom.p_neighbors))
            elif atom.p_neighbors:
                p_smi.append(str(atom.p_neighbors[0]))
            p_smi.append(';')
        elif atom.p_hybridization:
            p_smi.append(';')
            p_smi.append(hybridization_str[atom.p_hybridization[0]])
            if len(atom.p_neighbors) > 1:
                p_smi.append('<%s>' % ''.join(str(x) for x in atom.p_neighbors))
            elif atom.p_neighbors:
                p_smi.append(str(atom.p_neighbors[0]))
            p_smi.append(';')
        elif len(atom.p_neighbors) > 1:
            p_smi.append(';<%s>;' % ''.join(str(x) for x in atom.p_neighbors))
        elif atom.p_neighbors:
            p_smi.append(';')
            p_smi.append(str(atom.p_neighbors[0]))
            p_smi.append(';')

        if atom.charge:
            smi.append(charge_str[atom.charge])
        if atom.p_charge:
            p_smi.append(charge_str[atom.p_charge])
        if atom.multiplicity:
            smi.append(multiplicity_str[atom.multiplicity])
        if atom.p_multiplicity:
            p_smi.append(multiplicity_str[atom.p_multiplicity])

        smi.append(']')
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


__all__ = ['Smiles', 'SmilesCGR', 'SmilesQuery', 'SmilesQueryCGR', 'HashableSmiles']
