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
from CachedMethods import cached_method, cached_property
from collections import defaultdict
from hashlib import sha512
from itertools import count, product


charge_str = {-4: '-4', -3: '-3', -2: '-2', -1: '-', 0: '0', 1: '+', 2: '++', 3: '+3', 4: '+4'}
order_str = {1: '', 2: '=', 3: '#', 4: ':', 8: '~', None: '.'}
organic_set = {'C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B'}
hybridization_str = {4: 'a', 3: 't', 2: 'd', 1: 's', None: 'n'}
dyn_order_str = {(None, 1): "[.>-]", (None, 2): "[.>=]", (None, 3): "[.>#]", (None, 4): "[.>:]", (None, 8): "[.>~]",
                 (1, None): "[->.]", (1, 1): "", (1, 2): "[->=]", (1, 3): "[->#]", (1, 4): "[->:]", (1, 8): "[->~]",
                 (2, None): "[=>.]", (2, 1): "[=>-]", (2, 2): "=", (2, 3): "[=>#]", (2, 4): "[=>:]", (2, 8): "[=>~]",
                 (3, None): "[#>.]", (3, 1): "[#>-]", (3, 2): "[#>=]", (3, 3): "#", (3, 4): "[#>:]", (3, 8): "[#>~]",
                 (4, None): "[:>.]", (4, 1): "[:>-]", (4, 2): "[:>=]", (4, 3): "[:>#]", (4, 4): ":", (4, 8): "[:>~]",
                 (8, None): "[~>.]", (8, 1): "[~>-]", (8, 2): "[~>=]", (8, 3): "[~>#]", (8, 4): "[~>:]", (8, 8): "~"}

dyn_charge_str = {(i, j): f'{charge_str[i]}>{charge_str[j]}' if i != j else charge_str[i]
                  for i, j in product(range(-4, 5), repeat=2)}
dyn_charge_str[(0, 0)] = ''

dyn_radical_str = {(True, True): "*", (True, False): "*>n", (False, True): "n>*"}


class Smiles:
    __slots__ = ()

    @cached_method
    def __str__(self):
        return self._smiles(self.atoms_order.get)

    def __format__(self, format_spec):
        """
        signature generation options

        :param format_spec: string with keys:
            a - generate asymmetric closures
            !s - disable stereo marks
            !h - disable hybridization marks in queries
            !n - disable neighbors marks in queries

            combining possible. order independent. another key ignored
        """
        if format_spec:
            kwargs = {}
            if 'a' in format_spec:
                kwargs['asymmetric_closures'] = True
            if '!s' in format_spec:
                kwargs['stereo'] = False
            if '!h' in format_spec:
                kwargs['hybridization'] = False
            if '!n' in format_spec:
                kwargs['neighbors'] = False
            return self._smiles(self.atoms_order.get, **kwargs)
        return str(self)

    def __eq__(self, other):
        return isinstance(other, Smiles) and str(self) == str(other)

    @cached_method
    def __hash__(self):
        return hash(str(self))

    @cached_method
    def __bytes__(self):
        return sha512(str(self).encode()).digest()

    def _smiles(self, weights, *, asymmetric_closures=False, **kwargs):
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
            visited = {start: []}  # predecessors for stereo. atom: (visited[atom], *edges[atom])
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
                    if smiles[-2] == '(':
                        smiles.pop(-2)
                    else:
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
                    string.append(self._format_atom(token, adjacency=visited, **kwargs))
                    if token in tokens:
                        for m, c in tokens[token]:
                            if asymmetric_closures:
                                if (token, m) not in visited_bond:
                                    string.append(self._format_bond(token, m, adjacency=visited, **kwargs))
                                    visited_bond.add((m, token))
                            else:
                                string.append(self._format_bond(token, m, adjacency=visited, **kwargs))
                            c = casted_cycles[c]
                            string.append(str(c) if c < 10 else f'%{c}')
                elif token in ('(', ')'):
                    string.append(token)
                else:  # bonds
                    string.append(self._format_bond(*token, adjacency=visited, **kwargs))

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

    def _format_atom(self, n, **kwargs):
        atom = self._atoms[n]
        charge = self._charges[n]
        ih = self._hydrogens[n]
        if atom.isotope:
            smi = [str(atom.isotope), atom.atomic_symbol]
        else:
            smi = [atom.atomic_symbol]

        if kwargs.get('stereo', True) and n in self._atoms_stereo:  # carbon only
            smi.append('@' if self._translate_tetrahedron_stereo(n, kwargs['adjacency'][n]) else '@@')
            if ih:
                smi.append('H')
            smi.insert(0, '[')
            smi.append(']')
        elif charge:
            if ih == 1:
                smi.append('H')
            elif ih:
                smi.append(f'H{ih}')
            smi.append(charge_str[charge])
            smi.insert(0, '[')
            smi.append(']')
        elif self._radicals[n]:
            if ih == 1:
                smi.append('H')
            elif ih:
                smi.append(f'H{ih}')
            smi.insert(0, '[')
            smi.append(']')
        elif n in self.__aromatic_atoms and atom.atomic_symbol in ('N', 'P'):  # heterocycles
            if ih == 1:
                smi.append('H]')
                smi.insert(0, '[')
            elif ih:
                smi.append(f'H{ih}]')
                smi.insert(0, '[')
        elif atom.atomic_symbol not in organic_set:
            if ih == 1:
                smi.append('H')
            elif ih:
                smi.append(f'H{ih}')
            smi.insert(0, '[')
            smi.append(']')
        return ''.join(smi)

    def _format_bond(self, n, m, **kwargs):
        return order_str[self._bonds[n][m].order]


class CGRSmiles(Smiles):
    __slots__ = ()

    def _format_atom(self, n, **kwargs):
        atom = self._atoms[n]
        charge = self._charges[n]
        is_radical = self._radicals[n]
        p_charge = self._p_charges[n]
        p_is_radical = self._p_radicals[n]
        if atom.isotope:
            smi = [str(atom.isotope), atom.atomic_symbol]
        else:
            smi = [atom.atomic_symbol]

        if charge or p_charge:
            smi.append(dyn_charge_str[(charge, p_charge)])
        if is_radical or p_is_radical:
            smi.append(dyn_radical_str[(is_radical, p_is_radical)])

        if len(smi) != 1 or atom.atomic_symbol not in organic_set:
            smi.insert(0, '[')
            smi.append(']')
        return ''.join(smi)

    def _format_bond(self, n, m, **kwargs):
        bond = self._bonds[n][m]
        return dyn_order_str[(bond.order, bond.p_order)]


class QuerySmiles(Smiles):
    __slots__ = ()

    def _format_atom(self, n, **kwargs):
        atom = self._atoms[n]
        charge = self._charges[n]
        hybridization = self._hybridizations[n]
        neighbors = self._neighbors[n]

        if atom.isotope:
            smi = ['[', str(atom.isotope), atom.atomic_symbol]
        else:
            smi = ['[', atom.atomic_symbol]

        if kwargs.get('stereo', True) and n in self._atoms_stereo:  # carbon only
            smi.append('@' if self._translate_tetrahedron_stereo(n, kwargs['adjacency'][n]) else '@@')

        if kwargs.get('hybridization', True) and hybridization:
            smi.append(';')
            smi.append(''.join(hybridization_str[x] for x in hybridization))
            if kwargs.get('neighbors', True) and neighbors:
                smi.append(''.join(str(x) for x in neighbors))
            smi.append(';')
        elif kwargs.get('neighbors', True) and neighbors:
            smi.append(';')
            smi.append(''.join(str(x) for x in neighbors))
            smi.append(';')

        if charge:
            smi.append(charge_str[charge])
        if self._radicals[n]:
            smi.append('*')

        smi.append(']')
        return ''.join(smi)

    def _format_bond(self, n, m, **kwargs):
        return order_str[self._bonds[n][m].order]


class QueryCGRSmiles(Smiles):
    __slots__ = ()

    def _format_atom(self, n, **kwargs):
        atom = self._atoms[n]
        charge = self._charges[n]
        hybridization = self._hybridizations[n]
        neighbors = self._neighbors[n]
        is_radical = self._radicals[n]
        p_charge = self._p_charges[n]
        p_hybridization = self._p_hybridizations[n]
        p_neighbors = self._p_neighbors[n]
        p_is_radical = self._p_radicals[n]

        if atom.isotope:
            smi = ['[', str(atom.isotope), atom.atomic_symbol]
        else:
            smi = ['[', atom.atomic_symbol]

        if kwargs.get('hybridization', True) and hybridization:
            smi.append(';')
            smi.append(''.join(hybridization_str[x] for x in hybridization))
            smi.append('>')
            smi.append(''.join(hybridization_str[x] for x in p_hybridization))
            if kwargs.get('neighbors', True) and neighbors:
                smi.append(''.join(str(x) for x in neighbors))
                smi.append('>')
                smi.append(''.join(str(x) for x in p_neighbors))
            smi.append(';')
        elif kwargs.get('neighbors', True) and neighbors:
            smi.append(';')
            smi.append(''.join(str(x) for x in neighbors))
            smi.append('>')
            smi.append(''.join(str(x) for x in p_neighbors))
            smi.append(';')

        if charge or p_charge:
            smi.append(dyn_charge_str[(charge, p_charge)])
        if is_radical or p_is_radical:
            smi.append(dyn_radical_str[(is_radical, p_is_radical)])

        smi.append(']')
        return ''.join(smi)

    def _format_bond(self, n, m, **kwargs):
        bond = self._bonds[n][m]
        return dyn_order_str[(bond.order, bond.p_order)]


__all__ = ['MoleculeSmiles', 'CGRSmiles', 'QuerySmiles', 'QueryCGRSmiles']
