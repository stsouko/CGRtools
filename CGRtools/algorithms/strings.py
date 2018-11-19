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


class StringCommon:
    def _dfs(self, start, weights=None, depth_limit=None):
        """
        modified NX dfs
        """
        adj = self._adj
        if weights is None:
            def weights(x):
                return x
        else:
            weights = weights.get

        if depth_limit is None:
            depth_limit = len(adj)
        stack = [(start, depth_limit, iter(sorted(adj[start], key=weights)))]
        visited = {start}
        cycles = defaultdict(list)
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
                elif child not in cycles:
                    cycles[parent].append(child)

        return visited, edges, cycles


class StringMolecule:
    def _stringify(self, start=None, depth_limit=None, weights=None, atom=True, isotope=False, stereo=False,
                   hybridization=False, neighbors=False):
        if start is None:
            depth_limit = None
            g_set = set(self._node)
            smiles = []
            while g_set:
                if weights is None:
                    start = min(g_set)
                else:
                    start = min(g_set, key=weights.get)

                visited, edges, cycles = self._dfs(start, weights, depth_limit)
                smiles.append(self.__stringify(start, visited, edges, cycles, atom, isotope, stereo,
                              hybridization, neighbors))
                g_set.difference_update(visited)
            return '.'.join(smiles)

        visited, edges, cycles = self._dfs(start, weights, depth_limit)
        return self.__stringify(start, visited, edges, cycles, atom, isotope, stereo, hybridization, neighbors)

    def __stringify(self, start, visited, edges, cycles, f_atom, f_isotope, f_stereo, f_hybridization, f_neighbors):
        node = self._node
        adj = self._adj
        c = 1
        atoms = {}
        reverse = defaultdict(list)
        for parent, children in cycles.items():
            tmp = [node[parent].stringify(f_atom, f_isotope, f_stereo, f_hybridization, f_neighbors)]
            for child in children:
                bond = f'{adj[parent][child].stringify(f_stereo)}{c}'
                tmp.append(bond)
                reverse[child].append(bond)
                c += 1
            atoms[parent] = ''.join(tmp)

        for parent, children in reverse.items():
            tmp = [node[parent].stringify(f_atom, f_isotope, f_stereo, f_hybridization, f_neighbors)]
            tmp.extend(children)
            atoms[parent] = ''.join(tmp)

        for i in visited - atoms.keys():
            atoms[i] = node[i].stringify(f_atom, f_isotope, f_stereo, f_hybridization, f_neighbors)

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
                        stack.append([child, ['(', adj[tail][child].stringify(f_stereo), atoms[child]],
                                      stack_len])
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

        return ''.join(smiles)


class StringCGR:
    def _stringify(self, start=None, depth_limit=None, weights=None, atom=True, isotope=False, stereo=False,
                   hybridization=False, neighbors=False):
        if start is None:
            depth_limit = None
            g_set = set(self._node)
            r_smiles, p_smiles = [], []
            while g_set:
                if weights is None:
                    start = min(g_set)
                else:
                    start = min(g_set, key=weights.get)

                visited, edges, cycles = self._dfs(start, weights, depth_limit)
                tmp = self.__stringify(start, visited, edges, cycles, atom, isotope, stereo,
                                       hybridization, neighbors)
                r_smiles.append(tmp[0])
                p_smiles.append(tmp[1])
                g_set.difference_update(visited)
            return f"{'.'.join(r_smiles)}>>{'.'.join(p_smiles)}"

        visited, edges, cycles = self._dfs(start, weights, depth_limit)
        return '>>'.join(self.__stringify(start, visited, edges, cycles, atom, isotope, stereo, hybridization,
                                          neighbors))

    def __stringify(self, start, visited, edges, cycles, f_atom, f_isotope, f_stereo, f_hybridization, f_neighbors):
        node = self._node
        adj = self._adj
        c = 1
        atoms = {}
        r_reverse = defaultdict(list)
        p_reverse = defaultdict(list)
        for parent, children in cycles.items():
            r_tmp, p_tmp = node[parent].stringify(f_atom, f_isotope, f_stereo, f_hybridization, f_neighbors)
            r_tmp = [r_tmp]
            p_tmp = [p_tmp]
            for child in children:
                r_bond, p_bond = adj[parent][child].stringify(f_stereo)
                sc = str(c)
                r_tmp.append(r_bond)
                p_tmp.append(p_bond)
                r_tmp.append(sc)
                p_tmp.append(sc)
                r_reverse[child].append(r_bond)
                r_reverse[child].append(sc)
                p_reverse[child].append(p_bond)
                p_reverse[child].append(sc)
                c += 1
            atoms[parent] = (''.join(r_tmp), ''.join(p_tmp))

        for (parent, r_children), p_children in zip(r_reverse.items(), p_reverse.values()):
            r_tmp, p_tmp = node[parent].stringify(f_atom, f_isotope, f_stereo, f_hybridization, f_neighbors)
            r_tmp = [r_tmp]
            p_tmp = [p_tmp]
            r_tmp.extend(r_children)
            p_tmp.extend(p_children)
            atoms[parent] = (''.join(r_tmp), ''.join(p_tmp))

        for i in visited - atoms.keys():
            atoms[i] = node[i].stringify(f_atom, f_isotope, f_stereo, f_hybridization, f_neighbors)

        r_atom, p_atom = atoms[start]
        stack = [[start, [r_atom], [p_atom], 0]]
        while True:
            tail, r_smiles, p_smiles, closure = stack[-1]
            if tail in edges:
                children = edges[tail]
                if len(children) > 1:
                    child = children[-1]
                    stack_len = len(stack)
                    r_bond, p_bond = adj[tail][child].stringify(f_stereo)
                    r_atom, p_atom = atoms[child]
                    stack.append([child, [r_bond, r_atom], [p_bond, p_atom], 0])
                    for child in children[:-1]:
                        r_bond, p_bond = adj[tail][child].stringify(f_stereo)
                        r_atom, p_atom = atoms[child]
                        stack.append([child, ['(', r_bond, r_atom], ['(', p_bond, p_atom], stack_len])
                else:
                    child = children[-1]
                    stack[-1][0] = child
                    r_bond, p_bond = adj[tail][child].stringify(f_stereo)
                    r_atom, p_atom = atoms[child]
                    r_smiles.append(r_bond)
                    r_smiles.append(r_atom)
                    p_smiles.append(p_bond)
                    p_smiles.append(p_atom)
            elif closure:
                stack.pop()
                r_smiles.append(')')
                p_smiles.append(')')
                stack[closure - 1][1].extend(r_smiles)
                stack[closure - 1][2].extend(p_smiles)
            elif len(stack) > 2:
                stack.pop()
                stack[-1][0] = tail
                stack[-1][1].extend(r_smiles)
                stack[-1][2].extend(p_smiles)
            elif len(stack) == 2:
                stack[0][1].extend(r_smiles)
                stack[0][2].extend(p_smiles)
                r_smiles = stack[0][1]
                p_smiles = stack[0][2]
                break
            else:
                break

        return ''.join(r_smiles), ''.join(p_smiles)


__all__ = ['StringCommon', 'StringCGR', 'StringMolecule']
