# -*- coding: utf-8 -*-
#
#  Copyright 2019, 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2020 Ravil Mukhametgaleev <sonic-mc@mail.ru>
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
from collections import defaultdict, deque, ChainMap
from itertools import chain, product
from typing import List, Tuple, Dict, Set, Any, Union, TYPE_CHECKING, Iterator
from ..containers import molecule  # cyclic imports resolve
from ..exceptions import ValenceError


if TYPE_CHECKING:
    from CGRtools import ReactionContainer


class GraphComponents:
    __slots__ = ()

    @cached_property
    def connected_components(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Isolated components of single graph. E.g. salts as ion pair.
        """
        if not self._atoms:
            return ()
        return tuple(tuple(x) for x in self._connected_components(self._bonds))

    @staticmethod
    def _connected_components(bonds: Dict[int, Union[Set[int], Dict[int, Any]]]) -> List[Set[int]]:
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
            components.append(seen)
            atoms.difference_update(seen)
        return components

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

    @cached_property
    def atoms_rings(self) -> Dict[int, Tuple[Tuple[int, ...]]]:
        """
        Dict of atoms rings which contains it.
        """
        rings = defaultdict(list)
        for r in self.sssr:
            for n in r:
                rings[n].append(r)
        return {n: tuple(rs) for n, rs in rings.items()}

    def _augmented_substructure(self, atoms, deep):
        atoms = set(atoms)
        bonds = self._bonds
        if atoms - self._atoms.keys():
            raise ValueError('invalid atom numbers')
        nodes = [atoms]
        for i in range(deep):
            n = {y for x in nodes[-1] for y in bonds[x]} | nodes[-1]
            if n in nodes:
                break
            nodes.append(n)
        return nodes


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
        charges = self._charges
        radicals = self._radicals

        tetra = []
        for n, atom in atoms.items():
            if atom.atomic_number == 6 and not charges[n] and not radicals[n]:
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

        adj = defaultdict(set)  # double bonds adjacency matrix
        if heteroatoms:
            atoms_numbers = {5, 6, 7, 8, 14, 15, 16, 33, 34, 52}
            for n, atom in atoms.items():
                if atom.atomic_number in atoms_numbers:
                    adj_n = adj[n].add
                    for m, bond in bonds[n].items():
                        if bond.order == 2 and atoms[m].atomic_number in atoms_numbers:
                            adj_n(m)
        else:
            for n, atom in atoms.items():
                if atom.atomic_number == 6:
                    adj_n = adj[n].add
                    for m, bond in bonds[n].items():
                        if bond.order == 2 and atoms[m].atomic_number == 6:
                            adj_n(m)
        if not adj:
            return ()

        terminals = [x for x, y in adj.items() if len(y) == 1]
        cumulenes = []
        while terminals:
            n = terminals.pop(0)
            m = adj[n].pop()
            path = [n, m]
            while m not in terminals:
                adj_m = adj[m]
                if len(adj_m) > 2:  # not cumulene. SO3 etc.
                    cumulenes.extend(zip(path, path[1:]))  # keep single double bonds.
                    break
                adj_m.discard(n)
                n, m = m, adj_m.pop()
                path.append(m)
            else:
                terminals.remove(m)
                adj[m].pop()
                cumulenes.append(tuple(path))
        return cumulenes


class CGRComponents:
    __slots__ = ()

    @cached_property
    def centers_list(self) -> Tuple[Tuple[int, ...], ...]:
        """ Get a list of lists of atoms of reaction centers
        """
        radicals = self._radicals
        p_charges = self._p_charges
        p_radicals = self._p_radicals
        center = set()
        adj = defaultdict(set)

        for n, c in self._charges.items():
            if c != p_charges[n] or radicals[n] != p_radicals[n]:
                center.add(n)

        for n, m_bond in self._bonds.items():
            for m, bond in m_bond.items():
                if bond.order != bond.p_order:
                    adj[n].add(m)
        center.update(adj)

        # changes in condensed aromatic rings.
        if self.aromatic_rings:
            adj_update = defaultdict(set)
            for r in self.aromatic_rings:
                if not center.isdisjoint(r):
                    n = r[0]
                    m = r[-1]
                    if n in adj and m in adj[n]:
                        for n, m in zip(r, r[1:]):
                            if m not in adj[n]:
                                adj_update[n].add(m)
                                adj_update[m].add(n)
                    elif any(n in adj and m in adj[n] for n, m in zip(r, r[1:])):
                        adj_update[n].add(m)
                        adj_update[m].add(n)
                        for n, m in zip(r, r[1:]):
                            if m not in adj[n]:
                                adj_update[n].add(m)
                                adj_update[m].add(n)
            for n, ms in adj_update.items():
                adj[n].update(ms)

        out = []
        while center:
            n = center.pop()
            if n in adj:
                seen = set()
                nextlevel = {n}
                while nextlevel:
                    thislevel = nextlevel
                    nextlevel = set()
                    for v in thislevel:
                        if v not in seen:
                            seen.add(v)
                            nextlevel.update(adj[v])
                out.append(tuple(seen))
                center.difference_update(seen)
            else:
                out.append((n,))
        return tuple(out)

    @cached_property
    def center_atoms(self) -> Tuple[int, ...]:
        """ Get list of atoms of reaction center (atoms with dynamic: bonds, charges, radicals).
        """
        radicals = self._radicals
        p_charges = self._p_charges
        p_radicals = self._p_radicals

        center = set()
        for n, c in self._charges.items():
            if c != p_charges[n] or radicals[n] != p_radicals[n]:
                center.add(n)

        for n, m_bond in self._bonds.items():
            if any(bond.order != bond.p_order for bond in m_bond.values()):
                center.add(n)

        return tuple(center)

    @cached_property
    def center_bonds(self) -> Tuple[Tuple[int, int], ...]:
        """ Get list of bonds of reaction center (bonds with dynamic orders).
        """
        return tuple((n, m) for n, m, bond in self.bonds() if bond.order != bond.p_order)

    @cached_property
    def aromatic_rings(self) -> Tuple[Tuple[int, ...], ...]:
        """
        existed or formed aromatic rings atoms numbers
        """
        adj = self._bonds
        return tuple(ring for ring in self.sssr if
                     adj[ring[0]][ring[-1]].order == 4 and all(adj[n][m].order == 4 for n, m in zip(ring, ring[1:])) or
                     adj[ring[0]][ring[-1]].p_order == 4 and all(adj[n][m].p_order == 4 for n, m in zip(ring, ring[1:]))
                     )


class ReactionComponents:
    __slots__ = ()

    @cached_property
    def centers_list(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Reaction centers with leaving and coming groups.
        """
        if not self.reactants or not self.products:
            return ()  # no rc
        elif not isinstance(self.reactants[0], molecule.MoleculeContainer):
            raise TypeError('Only Molecules supported')

        cgr = self.compose()
        bonds = cgr._bonds
        all_groups = {x for x in self.reactants for x in x} ^ {x for x in self.products for x in x}
        all_groups = cgr._connected_components({n: bonds[n].keys() & all_groups for n in all_groups})
        centers_list = list(cgr.centers_list)

        for x in all_groups:
            intersection = []
            for i, y in enumerate(centers_list):
                if not x.isdisjoint(y):
                    intersection.append(i)
            if intersection:
                for i in reversed(intersection):
                    x.update(centers_list.pop(i))
                centers_list.append(x)
        return tuple(tuple(x) for x in centers_list)

    @cached_property
    def extended_centers_list(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Additionally to `centers_list` include:
        * First environment of dynamic atoms.
        * Whole formed cycles. For condensed cycles smallest is taken.
        * Whole aromatic cycle with at least one dynamic atom.
        * Whole small (3, 4) cycle with at least one dynamic atom.
        * Double or triple bonds connected to previous atoms.

        Note for multiple RCs intersection possible. Use `enumerate_centers` to prevent unobvious RCs.
        """
        cgr = self.compose()
        bonds = cgr._bonds
        center_atoms = set(cgr.center_atoms)

        formed_rings = {}
        small_aromatic_rings = set()
        for r in cgr.sssr:
            if len(r) < 5 and not center_atoms.isdisjoint(r):
                small_aromatic_rings.add(frozenset(r))

            n, m = r[0], r[-1]
            if bonds[n][m].order is None:
                nm = frozenset((n, m))
                if nm not in formed_rings or len(formed_rings[nm]) > len(r):
                    formed_rings[nm] = r
            for n, m in zip(r, r[1:]):
                if bonds[n][m].order is None:
                    nm = frozenset((n, m))
                    if nm not in formed_rings or len(formed_rings[nm]) > len(r):
                        formed_rings[nm] = r

        for r in cgr.aromatic_rings:
            if not center_atoms.isdisjoint(r):
                small_aromatic_rings.add(frozenset(r))

        out = []
        for rc in self.centers_list:
            c = center_atoms.intersection(rc)
            fe = {m for n in c for m in bonds[n]}  # add first environment
            # add double or triple bond to first env
            fe |= {m for n in fe for m, b in bonds[n].items() if m not in fe and b.order in (2, 3)}
            fe.update(rc)  # add leaving|coming groups and alone center atoms

            for rk, r in formed_rings.items():  # add formed rings to RC
                if c.issuperset(rk):
                    fe.update(r)
            for r in small_aromatic_rings:  # add small or aromatic rings with dyn atoms
                if not c.isdisjoint(r):
                    fe.update(r)
            out.append(tuple(fe))
        return tuple(out)

    def enumerate_centers(self) -> Iterator['ReactionContainer']:
        """
        Get all possible single stage reactions from multistage.
        Note multicomponent molecules (salts etc) can be treated incorrectly.
        """
        if len(self.centers_list) > 1:
            centers_list = self.centers_list

            charges = ChainMap(*(x._charges for x in self.reactants))
            radicals = ChainMap(*(x._radicals for x in self.reactants))
            bonds = ChainMap(*(x._bonds for x in self.reactants))
            atoms = ChainMap(*(x._atoms for x in self.reactants))
            p_charges = ChainMap(*(x._charges for x in self.products))
            p_radicals = ChainMap(*(x._radicals for x in self.products))
            p_bonds = ChainMap(*(x._bonds for x in self.products))
            p_atoms = ChainMap(*(x._atoms for x in self.products))

            centers = {x for x in centers_list for x in x}
            common = {x for x in chain(self.reactants, self.products) for x in x if x not in centers}
            reactants = {x for x in self.reactants for x in x}
            products = {x for x in self.products for x in x}

            common_molecule = molecule.MoleculeContainer()
            for n in common:
                common_molecule.add_atom(atoms[n].copy(), n, charge=charges[n], is_radical=radicals[n])
            seen = set()
            for n in common:
                seen.add(n)
                for m, b in bonds[n].items():
                    if m not in seen and m in common:
                        common_molecule.add_bond(n, m, b.copy())

            products_bonds = {}
            reactants_bonds = {}
            common_bonds = []
            seen = set()
            p_seen = set()
            for c in centers_list:
                not_rc = centers.difference(c)
                reactants_bonds[c] = (c_bonds, c_atoms) = [], reactants.intersection(c)
                for n in c_atoms:
                    seen.add(n)
                    for m, b in bonds[n].items():
                        if m not in seen and m in reactants:
                            if m in not_rc:
                                common_bonds.append((n, m, b))
                            else:
                                c_bonds.append((n, m, b))
                products_bonds[c] = (c_bonds, c_atoms) = [], products.intersection(c)
                for n in c_atoms:
                    p_seen.add(n)
                    for m, b in p_bonds[n].items():
                        if m not in p_seen and m in products and m not in not_rc:
                            c_bonds.append((n, m, b))

            for rc in range(len(centers_list)):
                not_rc = centers_list[:rc] + centers_list[rc + 1:]
                rc = centers_list[rc]
                for combo in list(product((False, True), repeat=len(not_rc))):
                    r = common_molecule.copy()
                    p = common_molecule.copy()

                    for n in reactants_bonds[rc][1]:
                        r.add_atom(atoms[n].copy(), n, charge=charges[n], is_radical=radicals[n])
                    for n in products_bonds[rc][1]:
                        p.add_atom(p_atoms[n].copy(), n, charge=p_charges[n], is_radical=p_radicals[n])

                    for is_p, center in zip(combo, not_rc):
                        if is_p:
                            c_bonds, c_atoms = products_bonds[center]
                            for n in c_atoms:
                                r.add_atom(p_atoms[n].copy(), n, charge=p_charges[n], is_radical=p_radicals[n])
                                p.add_atom(p_atoms[n].copy(), n, charge=p_charges[n], is_radical=p_radicals[n])
                        else:
                            c_bonds, c_atoms = reactants_bonds[center]
                            for n in c_atoms:
                                r.add_atom(atoms[n].copy(), n, charge=charges[n], is_radical=radicals[n])
                                p.add_atom(atoms[n].copy(), n, charge=charges[n], is_radical=radicals[n])
                        for n, m, b in c_bonds:
                            r.add_bond(n, m, b.copy())
                            p.add_bond(n, m, b.copy())

                    for n, m, b in products_bonds[rc][0]:
                        p.add_bond(n, m, b.copy())
                    for n, m, b in reactants_bonds[rc][0]:
                        r.add_bond(n, m, b.copy())
                    for n, m, b in common_bonds:
                        r.add_bond(n, m, b.copy())
                        p.add_bond(n, m, b.copy())
                    yield self.__class__(r.split(), p.split(), [x.copy() for x in self.reagents])
        else:
            cp = self.copy()
            cp.meta.clear()
            yield cp


__all__ = ['GraphComponents', 'StructureComponents', 'CGRComponents', 'ReactionComponents']
