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
from collections import defaultdict
from hashlib import sha512
from itertools import count
from ..cache import cached_method


class StringCommon:
    @cached_method
    def __str__(self):
        self._stringify()

    @cached_method
    def __bytes__(self):
        sha512(str(self).encode()).digest()

    def __hash__(self):
        return hash(str(self))

    def __eq__(self, other):
        return str(self) == str(other)

    def _dfs(self, start, weights, depth_limit):
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

    def stringify(self, atom=True, isotope=True, stereo=False, hybridization=False, neighbors=False):
        smi = []

        if stereo and self.__stereo:
            smi.append(self._stereo_str[self.__stereo])
        if hybridization and self.__hybridization:
            smi.append(self._hybridization_str[self.__hybridization])
        if neighbors and self.__neighbors is not None:
            smi.append(str(self.__neighbors))
        if smi:
            smi.append(';')
            smi.insert(0, ';')

        if atom:
            if self.element not in ('C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B'):
                atom = False
            smi.insert(0, self.element)
            if self.charge:
                smi.append(self._charge_str[self.charge])
            if self.multiplicity:
                smi.append(self._multiplicity_str[self.multiplicity])
            if isotope and self.isotope != self.common_isotope:
                smi.insert(0, str(self.isotope))
        else:
            smi.insert(0, '*')

        if len(smi) != 1 or not atom:
            smi.insert(0, '[')
            smi.append(']')

        return ''.join(smi)

    def _stringify(self):
        weights = self.atoms_order.__getitem__
        cycles = BroodMother()
        atoms_set = set(self._node)
        smiles = []

        while atoms_set:
            start = min(atoms_set, key=weights)
            visited, edges, disconnected = self._dfs(start, weights, len(atoms_set))
            atoms_set.difference_update(visited)
            smiles.append(self.__stringify(start, visited, edges, disconnected, cycles))
        return '.'.join(smiles)

    def __stringify(self, start, visited, edges, disconnected, cycles):
        node = self._node
        adj = self._adj

        reverse = defaultdict(list)
        atoms = defaultdict(Atom)
        for parent, children in disconnected.items():
            atoms[parent].atom = str(node[parent])
            for child in children:
                bond = (adj[parent][child].stringify(f_stereo), next(cycles))
                atoms[parent].connections.append(bond)
                reverse[child].append(bond)

        for parent, children in reverse.items():
            if parent not in atoms:
                atoms[parent].atom = node[parent].stringify(f_atom, f_isotope, f_stereo, f_hybridization, f_neighbors)
            atoms[parent].connections.extend(children)

        for i in visited - atoms.keys():
            atoms[i].atom = node[i].stringify(f_atom, f_isotope, f_stereo, f_hybridization, f_neighbors)

        stack = [[start, [atoms[start]], 0]]
        while True:
            tail, smiles, closure = stack[-1]
            if tail in edges:
                children = edges[tail]
                if len(children) > 1:
                    child = children[-1]
                    stack_len = len(stack)
                    stack.append([child, [adj[tail][child].stringify(f_stereo), atoms[child]], 0])
                    for child in children[:-1]:
                        stack.append([child, ['(', adj[tail][child].stringify(f_stereo), atoms[child]], stack_len])
                else:
                    child = children[-1]
                    stack[-1][0] = child
                    smiles.append(adj[tail][child].stringify(f_stereo))
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
        return ''.join(str(x) for x in smiles)


class SMILES_CGR(StringCommon):

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

    def _stringify_full(self, weights, atom=True, isotope=True, stereo=False, hybridization=False, neighbors=False):
        cycles = BroodMother()
        atoms_set = set(self._node)
        r_smiles, p_smiles = [], []

        while atoms_set:
            start = min(atoms_set, key=weights)
            visited, edges, disconnected = self._dfs(start, weights, len(atoms_set))
            atoms_set.difference_update(visited)
            tmp0, tmp1 = self.__stringify(start, visited, edges, disconnected, cycles, atom, isotope, stereo,
                                          hybridization, neighbors)
            r_smiles.append(tmp0)
            p_smiles.append(tmp1)
        return f"{'.'.join(r_smiles)}>>{'.'.join(p_smiles)}"

    def _stringify_augmented(self, start, depth_limit, weights, atom=True, isotope=True, stereo=False,
                             hybridization=False, neighbors=False):
        visited, edges, disconnected = self._dfs(start, weights, depth_limit)
        return '>>'.join(self.__stringify(start, visited, edges, disconnected, BroodMother(), atom, isotope, stereo,
                                          hybridization, neighbors))

    def _stringify_chain(self, start, stop, atom=True, isotope=True, stereo=False,
                         hybridization=False, neighbors=False):
        """
        container specific signature generation
        """
        pass

    def __stringify(self, start, visited, edges, disconnected, cycles,
                    f_atom, f_isotope, f_stereo, f_hybridization, f_neighbors):
        node = self._node
        adj = self._adj

        reverse = defaultdict(list)
        atoms = defaultdict(Atom)
        p_reverse = defaultdict(list)
        p_atoms = defaultdict(Atom)
        for parent, children in disconnected.items():
            atoms[parent].atom, p_atoms[parent].atom = node[parent].stringify(f_atom, f_isotope, f_stereo,
                                                                              f_hybridization, f_neighbors)
            for child in children:
                nc = next(cycles)
                bond, p_bond = adj[parent][child].stringify(f_stereo)
                bond = (bond, nc)
                atoms[parent].connections.append(bond)
                reverse[child].append(bond)
                p_bond = (p_bond, nc)
                p_atoms[parent].connections.append(p_bond)
                p_reverse[child].append(p_bond)

        for (parent, children), p_children in zip(reverse.items(), p_reverse.values()):
            if parent not in atoms:
                atoms[parent].atom, p_atoms[parent].atom = node[parent].stringify(f_atom, f_isotope, f_stereo,
                                                                                  f_hybridization, f_neighbors)
            atoms[parent].connections.extend(children)
            p_atoms[parent].connections.extend(p_children)

        for i in visited - atoms.keys():
            atoms[i].atom, p_atoms[i].atom = node[i].stringify(f_atom, f_isotope, f_stereo, f_hybridization,
                                                               f_neighbors)

        stack = [[start, [atoms[start]], [p_atoms[start]], 0]]
        while True:
            tail, smiles, p_smiles, closure = stack[-1]
            if tail in edges:
                children = edges[tail]
                if len(children) > 1:
                    child = children[-1]
                    stack_len = len(stack)
                    bond, p_bond = adj[tail][child].stringify(f_stereo)
                    stack.append([child, [bond, atoms[child]], [p_bond, p_atoms[child]], 0])
                    for child in children[:-1]:
                        bond, p_bond = adj[tail][child].stringify(f_stereo)
                        stack.append([child, ['(', bond, atoms[child]], ['(', p_bond, p_atoms[child]], stack_len])
                else:
                    child = children[-1]
                    stack[-1][0] = child
                    bond, p_bond = adj[tail][child].stringify(f_stereo)
                    smiles.append(bond)
                    smiles.append(atoms[child])
                    p_smiles.append(p_bond)
                    p_smiles.append(p_atoms[child])
            elif closure:
                stack.pop()
                smiles.append(')')
                p_smiles.append(')')
                stack[closure - 1][1].extend(smiles)
                stack[closure - 1][2].extend(p_smiles)
            elif len(stack) > 2:
                stack.pop()
                stack[-1][0] = tail
                stack[-1][1].extend(smiles)
                stack[-1][2].extend(p_smiles)
            elif len(stack) == 2:
                stack[0][1].extend(smiles)
                stack[0][2].extend(p_smiles)
                smiles = stack[0][1]
                p_smiles = stack[0][2]
                break
            else:
                break

        return ''.join(str(x) for x in smiles), ''.join(str(x) for x in p_smiles)


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


class Atom:
    __slots__ = ('atom', 'connections')

    def __init__(self):
        self.connections = []

    def __str__(self):
        if self.connections:
            return self.atom + ''.join(f'{o}{c}' for o, c in sorted(self.connections, key=lambda x: int(x[1])))
        return self.atom


__all__ = ['SMILES_CGR', 'SMILES']
