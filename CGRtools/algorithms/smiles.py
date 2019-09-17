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
from CachedMethods import cached_method, cached_args_method, cached_property
from collections import defaultdict
from hashlib import sha512
from itertools import count


charge_str = {-3: '-3', -2: '-2', -1: '-', 0: '0', 1: '+', 2: '++', 3: '+3'}
order_str = {1: '', 2: '=', 3: '#', 4: ':', 8: '~', None: '.'}
organic_set = {'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B'}
hybridization_str = {4: 'a', 3: 't', 2: 'd', 1: 's', None: 'n'}
stereo_str = {1: '@', -1: '@@'}
dyn_order_str = {(None, 1): "[.>-]", (None, 2): "[.>=]", (None, 3): "[.>#]", (None, 4): "[.>:]", (None, 8): "[.>~]",
                 (1, None): "[->.]", (1, 1): "", (1, 2): "[->=]", (1, 3): "[->#]", (1, 4): "[->:]", (1, 8): "[->~]",
                 (2, None): "[=>.]", (2, 1): "[=>-]", (2, 2): "=", (2, 3): "[=>#]", (2, 4): "[=>:]", (2, 8): "[=>~]",
                 (3, None): "[#>.]", (3, 1): "[#>-]", (3, 2): "[#>=]", (3, 3): "#", (3, 4): "[#>:]", (3, 8): "[#>~]",
                 (4, None): "[:>.]", (4, 1): "[:>-]", (4, 2): "[:>=]", (4, 3): "[:>#]", (4, 4): ":", (4, 8): "[:>~]",
                 (8, None): "[~>.]", (8, 1): "[~>-]", (8, 2): "[~>=]", (8, 3): "[~>#]", (8, 4): "[~>:]", (8, 8): "~"}
dyn_charge_str = {(-3, -3): "-3", (-3, -2): "-3>-2", (-3, -1): "-3>-", (-3, 0): "-3>0", (-3, 1): "-3>+",
                  (-3, 2): "-3>+2", (-3, 3): "-3>3", (-2, -3): "-2>-3", (-2, -2): "-2", (-2, -1): "-2>-",
                  (-2, 0): "-2>0", (-2, 1): "-2>+", (-2, 2): "-2>+2", (-2, 3): "-2>+3", (-1, -3): "->-3",
                  (-1, -2): "->-2", (-1, -1): "-", (-1, 0): "->0", (-1, 1): "->+", (-1, 2): "->+2", (-1, 3): "->+3",
                  (0, -3): "0>-3", (0, -2): "0>-2", (0, -1): "0>-", (0, 0): "", (0, 1): "0>+", (0, 2): "0>+2",
                  (0, 3): "0>+3", (1, -3): "+>-3", (1, -2): "+>-2", (1, -1): "+>-", (1, 0): "+>0", (1, 1): "+",
                  (1, 2): "+>+2", (1, 3): "+>+3", (2, -3): "+2>-3", (2, -2): "+2>-2", (2, -1): "+2>-", (2, 0): "+2>0",
                  (2, 1): "+2>+", (2, 2): "+2", (2, 3): "+2>+3", (3, -3): "+3>-3", (3, -2): "+3>-2", (3, -1): "+3>-",
                  (3, 0): "+3>0", (3, 1): "+3>+", (3, 2): "+3>+2", (3, 3): "+3"}
dyn_radical_str = {(True, True): "*", (True, False): "*>n", (False, True): "n>*"}


class Smiles:
    __slots__ = ()

    @cached_method
    def __str__(self):
        return self._smiles(self.atoms_order.get)

    def __format__(self, format_spec):
        if format_spec == "ac":
            return self._smiles(self.atoms_order.get, asymmetric_closures=True)
        return self._smiles(self.atoms_order.get)

    def __eq__(self, other):
        return str(self) == str(other)

    @cached_method
    def __hash__(self):
        return hash(str(self))

    @cached_method
    def __bytes__(self):
        return sha512(str(self).encode()).digest()

    def _smiles(self, weights, *, asymmetric_closures=False):
        bonds = self._bonds
        atoms_set = set(self._atoms)
        cycles = count()
        casted_cycles = {}
        string = []
        if asymmetric_closures:
            visited_bond = set()

        while True:
            start = min(atoms_set, key=weights)

            # modified NX dfs with cycle detection
            stack = [(start, len(atoms_set), iter(sorted(bonds[start], key=weights)))]
            visited = {start: [None]}  # predecessors for stereo. atom: (visited[atom], *edges[atom])
            disconnected = set()
            edges = defaultdict(list)
            tokens = defaultdict(list)
            while stack:
                parent, depth_now, children = stack[-1]
                try:
                    child = next(children)
                except StopIteration:
                    stack.pop()
                else:
                    if child not in visited:
                        edges[parent].append(child)
                        visited[child] = [parent]
                        if depth_now > 1:
                            front = bonds[child].keys() - {parent}
                            if front:
                                stack.append((child, depth_now - 1, iter(sorted(front, key=weights))))
                    elif child not in disconnected:
                        disconnected.add(parent)
                        cycle = next(cycles)
                        tokens[parent].append((child, cycle))
                        tokens[child].append((parent, cycle))

            # flatten directed graph: edges
            stack = [[start, 0, [start]]]
            while True:
                tail, closure, smiles = stack[-1]
                if tail in edges:
                    children = edges[tail]
                    if len(children) > 1:  # has side chain
                        child = children[-1]
                        stack_len = len(stack)
                        stack.append([child, 0, [(tail, child), child]])  # end of current chain
                        for child in children[-2::-1]:  # start side chains
                            stack.append([child, stack_len, ['(', (tail, child), child]])
                    else:  # chain grow
                        child = children[-1]
                        stack[-1][0] = child
                        smiles.append((tail, child))
                        smiles.append(child)
                elif closure:  # end of side chain
                    stack.pop()
                    smiles.append(')')
                    stack[closure - 1][2].extend(smiles)
                elif len(stack) > 2:
                    stack.pop()
                    stack[-1][0] = tail
                    stack[-1][2].extend(smiles)
                elif len(stack) == 2:
                    stack[0][2].extend(smiles)
                    smiles = stack[0][2]
                    break
                else:
                    break

            # prepare new neighbors order for stereo sign calculation
            for token in smiles:
                if token in tokens:
                    tokens[token].sort(key=lambda x: casted_cycles.get(x[1]) or
                                                     casted_cycles.setdefault(x[1], len(casted_cycles) + 1))
                    visited[token].extend(n for n, _ in tokens[token])
                if token in edges:
                    visited[token].extend(edges[token])
            for token in smiles:
                if isinstance(token, int):  # atoms
                    string.append(self._format_atom(token, visited))
                    if token in tokens:
                        for m, c in tokens[token]:
                            if asymmetric_closures:
                                if (token, m) not in visited_bond:
                                    string.append(self._format_bond(token, m, visited))
                                    visited_bond.add((m, token))
                            else:
                                string.append(self._format_bond(token, m, visited))
                            string.append(str(casted_cycles[c]))
                elif token in ('(', ')'):
                    string.append(token)
                else:  # bonds
                    string.append(self._format_bond(*token, visited))

            atoms_set.difference_update(visited)
            if atoms_set:
                string.append('.')
            else:
                break
        return ''.join(string)


class MoleculeSmiles(Smiles):
    __slots__ = ()

    @cached_property
    def __aromatic_atoms(self):
        aromatics = set()
        for ring in self.aromatic_rings:
            aromatics.update(ring)
        return aromatics

    def _format_atom(self, n, adjacency):
        atom = self._atoms[n]
        if atom.isotope:
            smi = [str(atom.isotope), atom.atomic_symbol]
        else:
            smi = [atom.atomic_symbol]

        # todo: stereo
        # smi.append()
        if atom.charge:
            h = atom.implicit_hydrogens
            if h == 1:
                smi.append('H')
            elif h:
                smi.append(f'H{h}')
            smi.append(charge_str[atom.charge])
        elif atom.is_radical:
            h = atom.implicit_hydrogens
            if h == 1:
                smi.append('H')
            elif h:
                smi.append(f'H{h}')
        elif n in self.__aromatic_atoms and atom.atomic_symbol in ('N', 'P') and atom.implicit_hydrogens:
            smi.append('H')

        if len(smi) != 1 or atom.atomic_symbol not in organic_set:
            smi.insert(0, '[')
            smi.append(']')
        return ''.join(smi)

    def _format_bond(self, n, m, adjacency):
        return order_str[self._bonds[n][m].order]


class CGRSmiles(Smiles):
    __slots__ = ()

    def _format_atom(self, n, adjacency):
        atom = self._atoms[n]
        if atom.isotope:
            smi = [str(atom.isotope), atom.atomic_symbol]
        else:
            smi = [atom.atomic_symbol]

        if atom.charge or atom.p_charge:
            smi.append(dyn_charge_str[(atom.charge, atom.p_charge)])
        if atom.is_radical or atom.p_is_radical:
            smi.append(dyn_radical_str[(atom.is_radical, atom.p_is_radical)])

        if len(smi) != 1 or atom.atomic_symbol not in organic_set:
            smi.insert(0, '[')
            smi.append(']')
        return ''.join(smi)

    def _format_bond(self, n, m, adjacency):
        bond = self._bonds[n][m]
        return dyn_order_str[(bond.order, bond.p_order)]


class QuerySmiles(Smiles):
    __slots__ = ()

    def _format_atom(self, n, adjacency):
        atom = self._atoms[n]
        if atom.isotope:
            smi = ['[', str(atom.isotope), atom.atomic_symbol]
        else:
            smi = ['[', atom.atomic_symbol]

        if len(atom.hybridization) > 1:
            smi.append(';')
            smi.append(''.join(hybridization_str[x] for x in atom.hybridization))
            if len(atom.neighbors) > 1:
                smi.append(''.join(str(x) for x in atom.neighbors))
            elif atom.neighbors:
                smi.append(str(atom.neighbors[0]))
            smi.append(';')
        elif atom.hybridization:
            smi.append(';')
            smi.append(hybridization_str[atom.hybridization[0]])
            if len(atom.neighbors) > 1:
                smi.append(''.join(str(x) for x in atom.neighbors))
            elif atom.neighbors:
                smi.append(str(atom.neighbors[0]))
            smi.append(';')
        elif len(atom.neighbors) > 1:
            smi.append(';')
            smi.append(''.join(str(x) for x in atom.neighbors))
            smi.append(';')
        elif atom.neighbors:
            smi.append(';')
            smi.append(str(atom.neighbors[0]))
            smi.append(';')

        if atom.charge:
            smi.append(charge_str[atom.charge])
        if atom.is_radical:
            smi.append('*')

        smi.append(']')
        return ''.join(smi)

    def _format_bond(self, n, m, adjacency):
        return order_str[self._bonds[n][m].order]


class QueryCGRSmiles(Smiles):
    __slots__ = ()

    def _format_atom(self, n, adjacency):
        atom = self._atoms[n]
        if atom.isotope:
            smi = ['[', str(atom.isotope), atom.atomic_symbol]
        else:
            smi = ['[', atom.atomic_symbol]

        if len(atom.hybridization) > 1:
            smi.append(';')
            smi.append(''.join(hybridization_str[x] for x in atom.hybridization))
            smi.append('>')
            smi.append(''.join(hybridization_str[x] for x in atom.p_hybridization))
            if len(atom.neighbors) > 1:
                smi.append(''.join(str(x) for x in atom.neighbors))
                smi.append('>')
                smi.append(''.join(str(x) for x in atom.p_neighbors))
            elif atom.neighbors:
                smi.append(str(atom.neighbors[0]))
                smi.append('>')
                smi.append(str(atom.p_neighbors[0]))
            smi.append(';')
        elif atom.hybridization:
            smi.append(';')
            smi.append(hybridization_str[atom.hybridization[0]])
            smi.append('>')
            smi.append(hybridization_str[atom.p_hybridization[0]])
            if len(atom.neighbors) > 1:
                smi.append(''.join(str(x) for x in atom.neighbors))
                smi.append('>')
                smi.append(''.join(str(x) for x in atom.p_neighbors))
            elif atom.neighbors:
                smi.append(str(atom.neighbors[0]))
                smi.append('>')
                smi.append(str(atom.p_neighbors[0]))
            smi.append(';')
        elif len(atom.neighbors) > 1:
            smi.append(';')
            smi.append(''.join(str(x) for x in atom.neighbors))
            smi.append('>')
            smi.append(''.join(str(x) for x in atom.p_neighbors))
            smi.append(';')
        elif atom.neighbors:
            smi.append(';')
            smi.append(str(atom.neighbors[0]))
            smi.append('>')
            smi.append(str(atom.p_neighbors[0]))
            smi.append(';')

        if atom.charge or atom.p_charge:
            smi.append(dyn_charge_str[(atom.charge, atom.p_charge)])
        if atom.is_radical or atom.p_is_radical:
            smi.append(dyn_radical_str[(atom.is_radical, atom.p_is_radical)])

        smi.append(']')
        return ''.join(smi)

    def _format_bond(self, n, m, adjacency):
        bond = self._bonds[n][m]
        return dyn_order_str[(bond.order, bond.p_order)]


__all__ = ['MoleculeSmiles', 'CGRSmiles', 'QuerySmiles', 'QueryCGRSmiles']
