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
from ..exceptions import InvalidStereo


class StereoMolecule:
    def add_stereo(self, atom1, atom2, mark):
        if mark not in (1, -1):
            raise ValueError('stereo mark invalid')
        if not self.has_edge(atom1, atom2):
            raise KeyError('atom or bond not found')

        if any(x.z for x in self._node.values()):
            raise InvalidStereo('molecule have 3d coordinates. bond up/down stereo unusable')

        implicit = self.atom_implicit_h(atom1)
        total = implicit + len(self._adj[atom1])

        if total == 4:  # tetrahedron
            if implicit:
                self._node[atom1].stereo = self.__pyramid_calc(atom1, atom2, mark)
            else:
                self._node[atom1].stereo = self.__tetrahedron_calc(atom1, atom2, mark)
        else:
            raise InvalidStereo('unsupported stereo or stereo impossible. tetrahedron only supported')

    def _wedge_map(self):
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
                    for _ in range(len_n):
                        if (neighbors[0], n) in _stereo_cache:
                            neighbors.append(neighbors.pop(0))
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
                        zero = self._node[neighbors[0]]
                        zero = (zero.x, zero.y, 1)
                        first = self._node[neighbors[1]]
                        first = (first.x, first.y, 0)
                    else:
                        zero = (atom.x, atom.y, 0)
                        first = self._node[neighbors[0]]
                        first = (first.x, first.y, 1)

                    second = self._node[neighbors[-2]]
                    third = self._node[neighbors[-1]]
                    vol = self._pyramid_volume(zero, first, (second.x, second.y, 0), (third.x, third.y, 0))

                    _stereo_cache[(n, neighbors[0])] = 1 if vol > 0 and s == 1 or vol < 0 and s == -1 else -1
            else:
                return _stereo_cache

    def __pyramid_calc(self, atom1, atom2, mark):
        center = self._node[atom1]
        order = ((self._node[x], mark if x == atom2 else 0) for x in self._adj[atom1])
        vol = self._pyramid_volume((center.x, center.y, 0), *((atom.x, atom.y, z) for atom, z in order))
        if vol > 0:
            return 1
        elif vol < 0:
            return -1
        raise InvalidStereo('unknown')

    def __tetrahedron_calc(self, atom1, atom2, mark):
        order = ((self._node[x], mark if x == atom2 else 0) for x in self._adj[atom1])
        vol = self._pyramid_volume(*((atom.x, atom.y, z) for atom, z in order))
        if vol > 0:
            return 1
        elif vol < 0:
            return -1
        raise InvalidStereo('unknown')

    @staticmethod
    def _pyramid_volume(n, u, v, w):
        nx, ny, nz = n
        ux, uy, uz = u
        vx, vy, vz = v
        wx, wy, wz = w
        ux -= nx
        uy -= ny
        uz -= nz
        vx -= nx
        vy -= ny
        vz -= nz
        wx -= nx
        wy -= ny
        wz -= nz
        return ux * (vy * wz - vz * wy) + uy * (vz * wx - vx * wz) + uz * (vx * wy - vy * wx)


class StereoCGR(StereoMolecule):
    def add_stereo(self, atom1, atom2, mark, p_mark=None):
        return  # disabled. need fix.
        if mark not in (1, -1) and p_mark not in (1, -1):
            raise InvalidData('stereo marks invalid')
        if not self.has_edge(atom1, atom2):
            raise InvalidAtom('atom or bond not found')

        n_atom1 = self.nodes[atom1]
        if n_atom1.get('s_stereo') or n_atom1.get('p_stereo'):
            raise self._stereo_exception3

        tmp_s = [(x, y['s_bond']) for x, y in self[atom1].items() if y.get('s_bond')]
        tmp_p = [(x, y['p_bond']) for x, y in self[atom1].items() if y.get('p_bond')]
        neighbors = [x for x, _ in tmp_s]

        if mark and (n_atom1['s_z'] or any(self.nodes[x]['s_z'] for x in neighbors)):
            raise self._stereo_exception1
        elif p_mark and (n_atom1['p_z'] or any(self.nodes[x]['p_z'] for x in neighbors)):
            raise self._stereo_exception1

        neighbors_e = [self.nodes[x]['element'] for x in neighbors]
        implicit_s, implicit_p = self.atom_implicit_h(atom1)
        if mark and (implicit_s > 1 or implicit_s == 1 and 'H' in neighbors_e or neighbors_e.count('H') > 1):
            raise self._stereo_exception4
        elif p_mark and (implicit_p > 1 or implicit_p == 1 and 'H' in neighbors_e or neighbors_e.count('H') > 1):
            raise self._stereo_exception4

        if mark:
            bonds = [x for _, x in tmp_s]
            total = implicit_s + len(neighbors)
            if total == 4:  # tetrahedron
                self._tetrahedron_parse(atom1, atom2, mark, neighbors, bonds, implicit_s)
            else:
                raise self._stereo_exception2
        if p_mark:
            bonds = [x for _, x in tmp_p]
            total = implicit_p + len(neighbors)
            if total == 4:  # tetrahedron
                self._tetrahedron_parse(atom1, atom2, p_mark, neighbors, bonds, implicit_p, label='p')
            else:
                raise self._stereo_exception2

    def _wedge_map(self):
        return {}
