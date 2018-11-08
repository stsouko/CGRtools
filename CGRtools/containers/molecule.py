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
from .common import BaseContainer
from ..algorithms import CGRstring, pyramid_volume, aromatize
from ..algorithms.morgan import initial_weights
from ..attributes import Atom, Bond
from ..exceptions import InvalidStereo
from ..periodictable import H


class MoleculeContainer(BaseContainer):
    """
    storage for Molecules

    new molecule object creation or copy of data object
    """
    node_attr_dict_factory = Atom
    edge_attr_dict_factory = Bond

    def __dir__(self):
        if self.__visible is None:
            self.__visible = super().__dir__() + [self.explicify_hydrogens.__name__, self.implicify_hydrogens.__name__,
                                                  self.atom_implicit_h.__name__, self.reset_query_marks.__name__,
                                                  self.aromatize.__name__, self.check_valence.__name__]
        return self.__visible

    def substructure(self, *args, **kwargs):
        sub = super().substructure(*args, **kwargs)
        sub._fix_stereo_stage_2(self._fix_stereo_stage_1())
        return sub

    def add_bond(self, atom1, atom2, bond):
        stereo = self._fix_stereo_stage_1()
        super().add_bond(atom1, atom2, bond)
        self._fix_stereo_stage_2(stereo)

    def delete_atom(self, n):
        """
        implementation of atom removing
        """
        stereo = self._fix_stereo_stage_1()
        self.remove_node(n)
        self._fix_stereo_stage_2(stereo)

    def delete_bond(self, n, m):
        """
        implementation of bond removing
        """
        stereo = self._fix_stereo_stage_1()
        self.remove_edge(n, m)
        self._fix_stereo_stage_2(stereo)

    def add_stereo(self, atom1, atom2, mark):
        if mark not in (1, -1):
            raise ValueError('stereo mark invalid')
        if not self.has_edge(atom1, atom2):
            raise KeyError('atom or bond not found')

        if self._node[atom1].stereo:
            raise self._stereo_exception3

        if self._node[atom1].z or any(self._node[x].z for x in self._adj[atom1]):
            raise self._stereo_exception1

        implicit = self.atom_implicit_h(atom1)
        if implicit + self.atom_explicit_h(atom1) > 1:
            raise self._stereo_exception4

        total = implicit + len(self._adj[atom1])
        if total == 4 and all(x.order == 1 for x in self._adj[atom1].values()):  # tetrahedron
            self.__tetrahedron_parse(atom1, atom2, mark, implicit)
        else:
            raise self._stereo_exception2

    def reset_query_marks(self, copy=False):
        """
        set or reset hyb and neighbors marks to atoms.

        :param copy: if True return copy of graph and keep existing as is
        :return: graph if copy True else None
        """
        g = self.copy() if copy else self
        for i, atom in g._node.items():
            neighbors = 0
            hybridization = 1
            # hybridization 1- sp3; 2- sp2; 3- sp1; 4- aromatic
            for j, bond in g._adj[i].items():
                if g._node[j] != 'H':
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

            atom.neighbors = neighbors
            atom.hybridization = hybridization

        if copy:
            return g
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

    def aromatize(self) -> int:
        """
        convert structure to aromatic form

        :return: number of processed rings
        """
        res = aromatize(self)
        if res:
            self.flush_cache()
        return res

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

    def __tetrahedron_parse(self, atom1, atom2, mark, implicit):
        weights = self.get_morgan(stereo=True)

        neighbors = list(self._adj[atom1])
        if len(neighbors) != len(set(weights[x] for x in neighbors)):
            raise InvalidStereo('stereo impossible. neighbors equivalent')

        order = sorted(((x, self._node[x]) for x in neighbors), key=lambda x: weights[x[0]])
        if implicit:
            central = self._node[atom1]
            vol = pyramid_volume((central.x, central.y, 0), *((atom.x, atom.y, 0 if x != atom2 else mark)
                                                              for x, atom in order))
        else:
            vol = pyramid_volume(*((atom.x, atom.y, 0 if x != atom2 else mark) for x, atom in order))

        if not vol:
            raise InvalidStereo('unknown')

        self._node[atom1].stereo = vol > 0 and 1 or -1
        self.flush_cache()

    def _prepare_stereo(self):
        _stereo_cache = {}
        nodes = list(self._node.items())
        while True:
            failed = []
            for i, tmp in enumerate(nodes, start=1):
                n, atom = tmp
                s = atom.stereo
                if not s:
                    continue
                neighbors = list(self._adj[n])
                len_n = len(neighbors)
                if len_n in (3, 4):  # tetrahedron
                    if atom.z or any(self._node[x].z for x in neighbors):
                        continue  # 3d molecules ignored

                    weights = self.get_morgan(stereo=True)
                    order = sorted(neighbors, key=weights.get)
                    for _ in range(len_n):
                        if (order[0], n) in _stereo_cache:
                            order.append(order.pop(0))
                        else:
                            failed.append(tmp)
                            break
                    else:
                        failed.insert(0, tmp)
                        failed.extend(nodes[i:])
                        _stereo_cache = None
                        nodes = failed
                        break

                    if len_n == 4:
                        zero = self._node[order[0]]
                        zero = (zero.x, zero.y, 1)
                        first = self._node[order[1]]
                        first = (first.x, first.y, 0)
                    else:
                        zero = (atom.x, atom.y, 0)
                        first = self._node[order[0]]
                        first = (first.x, first.y, 1)

                    second = self._node[order[-2]]
                    third = self._node[order[-1]]
                    vol = pyramid_volume(zero, first, (second.x, second.y, 0), (third.x, third.y, 0))

                    _stereo_cache[(n, order[0])] = 1 if vol > 0 and s == 1 or vol < 0 and s == -1 else -1
            else:
                return _stereo_cache

    def _signature_generator(self, *args, **kwargs):
        return CGRstring(*args, **kwargs)

    @property
    def _morgan_init(self):
        return initial_weights

    def _fix_stereo_stage_1(self):
        pass

    def _fix_stereo_stage_2(self, stereo):
        pass

    __visible = None
    _stereo_exception1 = InvalidStereo('molecule have 3d coordinates. bond up/down stereo unusable')
    _stereo_exception2 = InvalidStereo('unsupported stereo or stereo impossible. '
                                       'single bonded tetrahedron only supported')
    _stereo_exception3 = InvalidStereo('atom has stereo. change impossible')
    _stereo_exception4 = InvalidStereo('stereo impossible. too many H atoms')
