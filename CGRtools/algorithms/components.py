# -*- coding: utf-8 -*-
#
#  Copyright 2019, 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from CachedMethods import cached_property, FrozenDict
from collections import defaultdict, deque
from itertools import chain
from typing import Tuple, Dict, Set, Any, Union
from ..exceptions import ValenceError


class GraphComponents:
    __slots__ = ()

    @cached_property
    def connected_components(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Isolated components of single graph. E.g. salts as ion pair.
        """
        if not self._atoms:
            return ()
        return self._connected_components(self._bonds)

    @staticmethod
    def _connected_components(bonds: Dict[int, Union[Set[int], Dict[int, Any]]]) -> Tuple[Tuple[int, ...], ...]:
        atoms = set(bonds)
        components = []
        while atoms:
            start = atoms.pop()
            seen = {start}
            queue = deque([start])
            while queue:
                current = queue.popleft()
                for i in bonds[current]:
                    if i not in seen:
                        queue.append(i)
                        seen.add(i)
            components.append(tuple(seen))
            atoms.difference_update(seen)
        return tuple(components)

    @property
    def connected_components_count(self) -> int:
        """
        Number of components in graph
        """
        return len(self.connected_components)

    @cached_property
    def skin_atoms(self) -> Tuple[int, ...]:
        """
        Atoms of rings and rings linkers [without terminal atoms]
        """
        return tuple(self._skin_graph(self._bonds))

    @cached_property
    def skin_graph(self):
        """
        Graph without terminal atoms. Only rings and linkers
        """
        return FrozenDict((n, frozenset(ms)) for n, ms in self._skin_graph(self._bonds).items())

    @staticmethod
    def _skin_graph(bonds: Dict[int, Union[Set[int], Dict[int, Any]]]) -> Dict[int, Set[int]]:
        """
        Graph without terminal nodes. Only rings and linkers
        """
        bonds = {n: set(ms) for n, ms in bonds.items() if ms}
        while True:  # skip not-cycle chains
            try:
                n = next(n for n, ms in bonds.items() if len(ms) <= 1)
            except StopIteration:
                break
            for m in bonds.pop(n):
                bonds[m].discard(n)
        return bonds

    @cached_property
    def connected_rings(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Rings groups with common atoms. E.g. naphthalene has two connected rings. Rings not atom ordered like sssr.
        """
        rings = self.sssr
        if len(rings) <= 1:
            return rings

        rings = [set(r) for r in rings]
        out = []
        for i in range(len(rings)):
            r = rings[i]
            for x in rings[i + 1:]:
                if not r.isdisjoint(x):
                    x.update(r)
                    break
            else:  # isolated ring[s] found
                out.append(tuple(r))
        return tuple(out)

    @cached_property
    def ring_atoms(self):
        """
        Atoms in rings
        """
        return tuple({x for x in self.sssr for x in x})

    @cached_property
    def rings_count(self):
        """
        SSSR rings count.
        """
        bonds = self._bonds
        return sum(len(x) for x in bonds.values()) // 2 - len(bonds) + self.connected_components_count


class StructureComponents:
    __slots__ = ()

    @cached_property
    def aromatic_rings(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Aromatic rings atoms numbers
        """
        bonds = self._bonds
        return tuple(ring for ring in self.sssr if bonds[ring[0]][ring[-1]].order == 4
                     and all(bonds[n][m].order == 4 for n, m in zip(ring, ring[1:])))

    @cached_property
    def cumulenes(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Alkenes, allenes and cumulenes atoms numbers
        """
        return self._cumulenes()

    @cached_property
    def connected_rings_cumulenes(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Connected ring with attached cumulenes.
        """
        components = self.connected_rings
        if not components:
            return ()

        components = [set(r) for r in components]
        components.extend(set(c) for c in self.cumulenes)

        out = []
        for i in range(len(components)):
            c = components[i]
            for x in components[i + 1:]:
                if not c.isdisjoint(x):
                    x.update(c)
                    break
            else:  # isolated ring[s] found
                out.append(tuple(c))
        return tuple(out)

    @cached_property
    def tetrahedrons(self) -> Tuple[int, ...]:
        """
        Carbon sp3 atoms numbers
        """
        atoms = self._atoms
        bonds = self._bonds
        tetra = []
        for n, atom in atoms.items():
            if atom.atomic_number == 6 and not self._charges[n]:
                env = bonds[n]
                if all(x.order == 1 for x in env.values()):
                    b_sum = sum(x.order for x in env.values())
                    if b_sum > 4:
                        raise ValenceError(f'carbon atom: {n} has invalid valence = {b_sum}')
                    tetra.append(n)
        return tetra

    def _cumulenes(self, heteroatoms=False):
        atoms = self._atoms
        bonds = self._bonds

        if heteroatoms:
            atoms_numbers = {5, 6, 7, 8, 14, 15, 16, 33, 34, 52}
        else:
            atoms_numbers = {6}

        adj = defaultdict(set)  # carbon double bonds adjacency matrix
        for n, atom in atoms.items():
            if atom.atomic_number in atoms_numbers:
                adj_n = adj[n].add
                for m, bond in bonds[n].items():
                    order = bond.order
                    if order == 2 and atoms[m].atomic_number in atoms_numbers:
                        adj_n(m)
        if not adj:
            return ()

        for n, ms in adj.items():
            if len(ms) > 2:
                raise ValenceError(f'atom: {n} has invalid valence')

        terminals = [x for x, y in adj.items() if len(y) == 1]
        cumulenes = []
        while terminals:
            m = terminals.pop(0)
            path = [m]
            while m not in terminals:
                n, m = m, adj[m].pop()
                adj[m].discard(n)
                path.append(m)
            terminals.remove(m)
            if sum(1 for b in chain(bonds[path[0]].values(), bonds[path[-1]].values()) if b.order == 2) == 2:
                # check for carbon only double-bonded chains.
                cumulenes.append(tuple(path))
        return cumulenes


__all__ = ['GraphComponents', 'StructureComponents']
