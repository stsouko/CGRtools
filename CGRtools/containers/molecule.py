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
from networkx.algorithms.isomorphism import GraphMatcher
from networkx.classes.function import frozen
from typing import List
from .common import BaseContainer
from ..algorithms import Aromatize, Calculate2D, Compose, DepictMolecule, Morgan, Smiles, Standardize
from ..attributes import Atom, Bond
from ..cache import cached_args_method, cached_property
from ..periodictable import H


class MoleculeContainer(Aromatize, Calculate2D, Compose, Morgan, Smiles, Standardize, DepictMolecule, BaseContainer):
    """
    storage for Molecules

    new molecule object creation or copy of data object
    """
    node_attr_dict_factory = Atom
    edge_attr_dict_factory = Bond

    def reset_query_marks(self):
        """
        set or reset hyb and neighbors marks to atoms.
        """
        for i, atom in self.atoms():
            neighbors = 0
            hybridization = 1
            # hybridization 1- sp3; 2- sp2; 3- sp1; 4- aromatic
            for j, bond in self._adj[i].items():
                if self._node[j].element != 'H':
                    neighbors += 1

                if hybridization in (3, 4):
                    continue
                order = bond.order
                if order == 4:
                    hybridization = 4
                elif order == 3:
                    hybridization = 3
                elif order == 2:
                    if hybridization == 2:
                        hybridization = 3
                    else:
                        hybridization = 2

            atom._neighbors = neighbors
            atom._hybridization = hybridization
        self.flush_cache()

    def implicify_hydrogens(self):
        """
        remove explicit hydrogen if possible

        :return: number of removed hydrogens
        """
        explicit = defaultdict(list)
        c = 0
        for n, atom in self.atoms():
            if atom.element == 'H':
                for m in self.neighbors(n):
                    if self._node[m].element != 'H':
                        explicit[m].append(n)

        for n, h in explicit.items():
            atom = self._node[n]
            len_h = len(h)
            for i in range(len_h, 0, -1):
                hi = h[:i]
                if atom.get_implicit_h([y.order for x, y in self._adj[n].items() if x not in hi]) == i:
                    for x in hi:
                        self.remove_node(x)
                        c += 1
                    break

        self.flush_cache()
        return c

    def explicify_hydrogens(self):
        """
        add explicit hydrogens to atoms

        :return: number of added atoms
        """
        tmp = []
        for n, atom in self.atoms():
            if atom.element != 'H':
                for _ in range(atom.get_implicit_h([x.order for x in self._adj[n].values()])):
                    tmp.append(n)
        for n in tmp:
            self.add_bond(n, self.add_atom(H), Bond())

        self.flush_cache()
        return len(tmp)

    def substructure(self, atoms, meta=False, as_view=True):
        """
        create substructure containing atoms from nbunch list

        :param atoms: list of atoms numbers of substructure
        :param meta: if True metadata will be copied to substructure
        :param as_view: If True, the returned graph-view provides a read-only view
            of the original structure scaffold without actually copying any data.
        """
        s = super().substructure(atoms, meta, as_view)
        if as_view:
            s.check_valence = s.explicify_hydrogens = s.implicify_hydrogens = s.reset_query_marks = frozen
            s.standardize = s.aromatize = frozen
        return s

    @cached_args_method
    def atom_implicit_h(self, atom):
        return self._node[atom].get_implicit_h([x.order for x in self._adj[atom].values()])

    @cached_args_method
    def atom_explicit_h(self, atom):
        return sum(self._node[x].element == 'H' for x in self.neighbors(atom))

    @cached_args_method
    def atom_total_h(self, atom):
        return self.atom_explicit_h(atom) + self.atom_implicit_h(atom)

    def check_valence(self):
        """
        check valences of all atoms

        :return: list of invalid atoms
        """
        return [x for x, atom in self.atoms() if not atom.check_valence(self.environment(x))]

    @cached_property
    def aromatic_rings(self) -> List[List[int]]:
        """
        aromatic rings atoms numbers
        """
        adj = self._bonds
        return [ring for ring in self.sssr if len(ring) in (5, 6, 7) and adj[ring[0]][ring[-1]].order == 4
                and all(adj[n][m].order == 4 for n, m in zip(ring, ring[1:]))]

    @cached_property
    def cumulenes(self) -> List[List[int]]:
        """
        alkenes, allenes and cumulenes atoms numbers
        """
        atoms = self._atoms
        adj = defaultdict(set)  # carbon double bonds adjacency matrix
        for n, m_bond in self._bonds.items():
            if atoms[n].element == 'C':
                adj_n = adj[n].add
                for m, bond in m_bond.items():
                    if bond.order == 2 and atoms[m].element == 'C':
                        adj_n(m)
        if not adj:
            return []

        terminals = {x for x, y in adj.items() if len(y) == 1}
        cumulenes = []
        while terminals:
            m = terminals.pop()
            path = [m]
            cumulenes.append(path)
            while m not in terminals:
                n, m = m, adj[m].pop()
                adj[m].discard(n)
                path.append(m)
            terminals.discard(m)
        return cumulenes

    @cached_property
    def tetrahedrons(self) -> List[int]:
        """
        carbon sp3 atoms numbers
        """
        atoms = self._atoms
        return [n for n, env in self._bonds.items()
                if atoms[n].element == 'C' and all(x.order == 1 for x in env.values())]

    def _matcher(self, other):
        """
        return VF2 GraphMatcher

        MoleculeContainer < MoleculeContainer
        MoleculeContainer < CGRContainer
        """
        if isinstance(other, (self._get_subclass('CGRContainer'), MoleculeContainer)):
            return GraphMatcher(other, self, lambda x, y: x == y, lambda x, y: x == y)
        raise TypeError('only cgr-cgr possible')

    def __setstate__(self, state):
        if '_BaseContainer__meta' in state:  # 2.8 reverse compatibility
            node = {}
            for n, atom in state['_node'].items():
                a = node[n] = Atom()
                a.element = atom['element']
                a.mapping = atom['map']
                a.x = atom['s_x']
                a.y = atom['s_y']
                a.z = atom['s_z']
                if atom['s_charge']:
                    a.charge = atom['s_charge']
                if 's_radical' in atom:
                    a.multiplicity = atom['s_radical']
                if 'isotope' in atom:
                    a.isotope = atom['isotope']

            adj = defaultdict(dict)
            seen = set()
            for n, m_bond in state['_adj'].items():
                seen.add(n)
                for m, bond in m_bond.items():
                    if m not in seen:
                        b = adj[n][m] = adj[m][n] = Bond()
                        b.order = bond['s_bond'] if bond['s_bond'] != 9 else 5
            state = {'meta': state['_BaseContainer__meta'], 'node': node, 'adj': dict(adj)}
        super().__setstate__(state)


__all__ = ['MoleculeContainer']
