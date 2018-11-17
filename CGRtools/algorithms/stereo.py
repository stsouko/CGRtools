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
        return {}
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
