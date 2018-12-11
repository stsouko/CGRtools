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
from networkx.algorithms.isomorphism import GraphMatcher
from .common import BaseContainer
from ..algorithms import Aromatize, StringMolecule, CGRCompose
from ..attributes import Atom, Bond
from ..periodictable import H


class MoleculeContainer(StringMolecule, Aromatize, CGRCompose, BaseContainer):
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
        for i, atom in self._node.items():
            neighbors = 0
            hybridization = 1
            # hybridization 1- sp3; 2- sp2; 3- sp1; 4- aromatic
            for j, bond in self._adj[i].items():
                if self._node[j] != 'H':
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
        for n, atom in self._node.items():
            if atom == 'H':
                m = next(self.neighbors(n))
                if self._node[m] != 'H':
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
        for n, atom in self._node.items():
            if atom != 'H':
                for _ in range(atom.get_implicit_h([x.order for x in self._adj[n].values()])):
                    tmp.append(n)
        for n in tmp:
            self.add_bond(n, self.add_atom(H), Bond())

        self.flush_cache()
        return len(tmp)

    def atom_implicit_h(self, atom):
        return self._node[atom].get_implicit_h([x.order for x in self._adj[atom].values()])

    def atom_explicit_h(self, atom):
        return sum(self._node[x] == 'H' for x in self.neighbors(atom))

    def atom_total_h(self, atom):
        return self.atom_explicit_h(atom) + self.atom_implicit_h(atom)

    def check_valence(self):
        """
        check valences of all atoms

        :return: list of invalid atoms
        """
        return [f'atom {x} has invalid valence' for x, atom in self._node.items()
                if not atom.check_valence(self.environment(x))]

    def _matcher(self, other):
        """
        MoleculeContainer < MoleculeContainer
        MoleculeContainer < CGRContainer
        """
        if isinstance(other, (self._get_subclass('CGRContainer'), MoleculeContainer)):
            return GraphMatcher(other, self, lambda x, y: x == y, lambda x, y: x == y)
        raise TypeError('only cgr-cgr possible')

    _visible = ()


__all__ = ['MoleculeContainer']
