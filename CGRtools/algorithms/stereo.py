# -*- coding: utf-8 -*-
#
#  Copyright 2019 Ramil Nugmanov <stsouko@live.ru>
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
from CachedMethods import cached_property
from collections import defaultdict
from itertools import combinations, product
from logging import info
from ..exceptions import AtomNotFound, NotChiral, IsChiral, ValenceError


def _pyramid_sign(n, u, v, w):
    #
    #  |   n /
    #  |   |\
    #  |   | \
    #  |  /|  \
    #  | / u---v
    #  |/___\_/___
    #        w
    #
    nx, ny, nz = n
    ux, uy, uz = u
    vx, vy, vz = v
    wx, wy, wz = w

    q1x = ux - nx
    q1y = uy - ny
    q1z = uz - nz
    q2x = vx - nx
    q2y = vy - ny
    q2z = vz - nz
    q3x = wx - nx
    q3y = wy - ny
    q3z = wz - nz

    vol = q1x * (q2y * q3z - q2z * q3y) + q1y * (q2z * q3x - q2x * q3z) + q1z * (q2x * q3y - q2y * q3x)
    if vol > 0:
        return 1
    elif vol < 0:
        return -1
    return 0


def _dihedral_sign(n, u, v, w):
    # n    w
    # |   /
    # u--v
    nx, ny, nz = n
    ux, uy, uz = u
    vx, vy, vz = v
    wx, wy, wz = w

    q1x = ux - nx
    q1y = uy - ny
    q1z = uz - nz
    q2x = vx - ux
    q2y = vy - uy
    q2z = vz - uz
    q3x = wx - vx
    q3y = wy - vy
    q3z = wz - vz

    # cross vectors
    q1q2x = q1y * q2z - q1z * q2y
    q1q2y = q1z * q2x - q1x * q2z
    q1q2z = q1x * q2y - q1y * q2x
    q2q3x = q2y * q3z - q2z * q3y
    q2q3y = q2z * q3x - q2x * q3z
    q2q3z = q2x * q3y - q2y * q3x

    # angle calculation
    # len_q1q2 = sqrt(q1q2x ** 2 + q1q2y ** 2 + q1q2z ** 2)
    # n1x = q1q2x / len_q1q2
    # n1y = q1q2y / len_q1q2
    # n1z = q1q2z / len_q1q2
    # len_q2q3 = sqrt(q2q3x ** 2 + q2q3y ** 2 + q2q3z ** 2)
    # u1x = q2q3x / len_q2q3
    # u1y = q2q3y / len_q2q3
    # u1z = q2q3z / len_q2q3
    # len_q2 = sqrt(q2x ** 2 + q2y ** 2 + q2z ** 2)
    # u3x = q2x / len_q2
    # u3y = q2y / len_q2
    # u3z = q2z / len_q2
    # u2x = u3y * u1z - u3z * u1y
    # u2y = u3z * u1x - u3x * u1z
    # u2z = u3x * u1y - u3y * u1x
    # cos_theta = n1x * u1x + n1y * u1y + n1z * u1z
    # sin_theta = n1x * u2x + n1y * u2y + n1z * u2z
    # return -atan2(sin_theta, cos_theta)

    dot = q1q2x * q2q3x + q1q2y * q2q3y + q1q2z * q2q3z
    if dot > 0:
        return 1
    elif dot < 0:
        return -1
    return 0


class Stereo:
    __slots__ = ()

    def get_mapping(self, other, **kwargs):
        stereo = self._atoms_stereo
        if stereo:
            tetrahedrons = self._tetrahedrons
            for mapping in super().get_mapping(other, **kwargs):
                for n in stereo.keys() & mapping.keys():
                    m = mapping[n]
                    if m not in other._atoms_stereo:  # self stereo atom not stereo in other
                        break
                    # translate stereo mark in other in order of self tetrahedron
                    if stereo[n] != other._translate_tetrahedron_stereo(m, [mapping[x] for x in tetrahedrons[n]]):
                        break
                else:
                    yield mapping
        else:
            yield from super().get_mapping(other, **kwargs)

    @cached_property
    def _wedge_map(self):
        plane = self._plane
        wedge = []
        for n, s in self._atoms_stereo.items():
            order = sorted(self._tetrahedrons[n], key=self.atoms_order.get)
            s = self._translate_tetrahedron_stereo(n, order)
            # need recalculation if XY changed
            if len(order) == 3:
                v = _pyramid_sign((*plane[n], 0),
                                  (*plane[order[0]], 1), (*plane[order[1]], 0), (*plane[order[2]], 0))
            else:
                v = _pyramid_sign((*plane[order[3]], 0),
                                  (*plane[order[0]], 1), (*plane[order[1]], 0), (*plane[order[2]], 0))
            if not v:
                info(f'need 2d clean. wedge stereo ambiguous for atom {{{n}}}')
            v = v > 0
            if s:
                wedge.append((n, order[0], 1 if v else -1))
            else:
                wedge.append((n, order[0], -1 if v else 1))
        return tuple(wedge)

    def _translate_tetrahedron_stereo(self, n, env):
        order = self._tetrahedrons[n]
        if len(order) == 3:
            if len(env) == 4:
                atoms = self._atoms
                order = (*order, next(x for x in env if atoms[x].atomic_number == 1))  # see translate scheme
            elif len(env) != 3:
                raise ValueError('invalid atoms list')
        elif len(env) not in (3, 4):
            raise ValueError('invalid atoms list')
        s = self._atoms_stereo[n]

        translate = tuple(order.index(x) + 1 for x in env[:3])
        if _tetrahedron_translate[translate]:
            return not s
        return s


class MoleculeStereo(Stereo):
    __slots__ = ()

    def add_wedge(self, n: int, m: int, mark: bool):
        if n not in self._atoms:
            raise AtomNotFound
        if n in self._atoms_stereo:
            raise IsChiral

        plane = self._plane
        if n in self._chiral_atoms:
            if m not in self._bonds[n]:
                raise AtomNotFound

            if self._atoms[m].atomic_number == 1:
                s = _pyramid_sign((*plane[m], mark), *((*plane[x], 0) for x in self._tetrahedrons[n]))
            else:
                order = [(*plane[x], mark if x == m else 0) for x in self._tetrahedrons[n]]
                if len(order) == 3:
                    s = _pyramid_sign((*plane[n], 0), *order)
                else:
                    s = _pyramid_sign(order[-1], *order[:3])
            if s:
                self._atoms_stereo[n] = s > 0
                self.flush_cache()
        else:  # only tetrahedrons supported
            raise NotChiral

    def add_atom_stereo(self, n, env, mark: bool):
        if n not in self._atoms:
            raise AtomNotFound
        if n in self._atoms_stereo:
            raise IsChiral
        if len(env) != len(set(env)):
            raise ValueError('invalid environment')
        if not isinstance(mark, bool):
            raise TypeError('stereo mark should be bool')

        if n in self._chiral_atoms:
            if set(env) != set(self._bonds[n]):
                raise AtomNotFound

            translate = tuple(env.index(x) + 1 for x in self._tetrahedrons[n][:3])
            if _tetrahedron_translate[translate]:
                mark = not mark

            self._atoms_stereo[n] = mark
            self.flush_cache()
        else:  # only tetrahedrons supported
            raise NotChiral

    @cached_property
    def _tetrahedrons(self):
        #    2
        #    |
        # 1--K--3
        #    |
        #    4?
        atoms = self._atoms
        bonds = self._bonds
        tetrahedrons = {}
        for n in self.tetrahedrons:
            env = tuple(x for x in bonds[n] if atoms[x].atomic_number != 1)
            if len(env) in (3, 4):
                tetrahedrons[n] = env
        return tetrahedrons

    def _fix_stereo(self):
        if self._atoms_stereo:
            atoms = self._atoms
            stereo = {k: v for k, v in self._atoms_stereo.items() if k in atoms}
            self._atoms_stereo = new_stereo = {}

            old_stereo = 0
            while stereo and len(stereo) != old_stereo:
                old_stereo = len(stereo)
                failed_stereo = {}
                chiral = self._chiral_atoms
                self.__dict__.clear()
                for n, s in stereo.items():
                    if n in chiral:
                        new_stereo[n] = s
                    else:
                        failed_stereo[n] = s
                stereo = failed_stereo

    @cached_property
    def _chiral_atoms(self):
        morgan = self.atoms_order
        atoms_stereo = self._atoms_stereo
        tetrahedrons = self._tetrahedrons
        if atoms_stereo:
            grouped_stereo = defaultdict(list)
            morgan_update = {}

            for n in atoms_stereo:
                grouped_stereo[morgan[n]].append(n)  # collect equal stereo atoms
            for group in grouped_stereo.values():
                if len(group) % 2 == 0:  # only even number of equal stereo atoms give new stereo center
                    s = [n for n in group
                         if self._translate_tetrahedron_stereo(n, sorted(tetrahedrons[n], key=morgan.get))]
                    if 0 < len(s) < len(group):  # RS pair required
                        for n in s:
                            morgan_update[n] = -morgan[n]
            if morgan_update:
                morgan = self._morgan(self._sorted_primed({**morgan, **morgan_update}))

        return {n for n, env in self._tetrahedrons.items()
                if n not in atoms_stereo and len(set(morgan[x] for x in env)) == len(env)}


class QueryStereo(Stereo):
    __slots__ = ()

    @cached_property
    def _tetrahedrons(self):
        #    2
        #    |
        # 1--K--3
        #    |
        #    4?
        atoms = self._atoms
        bonds = self._bonds
        tetrahedrons = {}
        for n, atom in atoms.items():
            if atom.atomic_number == 6 and not self._charges[n]:
                env = bonds[n]
                b_sum = sum(x.order for x in env.values())
                if b_sum > 4:
                    raise ValenceError(f'carbon atom: {n} has invalid valence = {b_sum}')
                if all(x.order == 1 for x in env.values()):
                    env = tuple(x for x in env if atoms[x].atomic_number != 1)
                    if len(env) in (3, 4):
                        tetrahedrons[n] = env
        return tetrahedrons


class NotUsed:
    @cached_property
    def __cumulenes(self):
        # 5       4
        #  \     /
        #   2---3
        #  /     \
        # 1       6
        bonds = self._bonds
        atoms = self._atoms
        cumulenes = {}
        for path in self.cumulenes:
            n, m = path[0], path[-1]
            n1, m1 = path[1], path[-2]

            nn = [x for x in bonds[n] if x != n1 and atoms[x].atomic_number != 1]
            mn = [x for x in bonds[m] if x != m1 and atoms[x].atomic_number != 1]
            if nn and mn:
                sn = nn[1] if len(nn) == 2 else None
                sm = mn[1] if len(mn) == 2 else None
                cumulenes[(n, m)] = (nn[0], mn[0], sn, sm)
        return cumulenes

    def __atropoisomers(self):
        #     ___
        #    |   |
        # 7--3   5--9
        #     \ /
        #      1[n]
        #      |
        #      2[m]
        #     / \
        # 8--4   6--10
        #    |___|
        #
        bonds = self._bonds
        if len(self.sssr) < 2:
            return {}
        aromatic = self.aromatic_rings
        if len(aromatic) < 2:
            return {}

        reduced_aromatic = [[x for x in x if len(bonds[x]) == 3] for x in aromatic]  # remove :[CH]: atoms
        connections = {}
        for rings in combinations(range(len(aromatic)), 2):
            ring1, ring2 = rings
            for n, m in product(reduced_aromatic[ring1], reduced_aromatic[ring2]):
                if n in bonds[m]:
                    if rings in connections:  # remove condensed rings or twice-bonded rings
                        del connections[rings]
                        break  # skip rings
                    connections[rings] = (n, m)

        atropos = {}
        for (ring1, ring2), (n, m) in connections.items():
            # neighbors of connection atoms in rings
            a3, a5 = (x for x in bonds[n] if x != m)
            a4, a6 = (x for x in bonds[m] if x != n)
            # substituents of neighbors
            a7 = next((x for x in bonds[a3] if x not in aromatic[ring1]), None)
            a9 = next((x for x in bonds[a5] if x not in aromatic[ring1]), None)
            a8 = next((x for x in bonds[a4] if x not in aromatic[ring2]), None)
            a10 = next((x for x in bonds[a6] if x not in aromatic[ring2]), None)

            # skip rings without substituents
            # todo: rings bounded with chain
            if not a7 and (not a9 or not a8 or not a10):
                continue
            elif not a9 and (not a8 or not a10):
                continue
            elif not a8 and not a10:
                continue
            atropos[(n, m)] = (a3, a4, a5, a6, a7, a8, a9, a10)
        return atropos


# 2  3
#  \ |
#   \|
#    0---4
#   /
#  /
# 1
_tetrahedron_translate = {(1, 2, 3): False, (2, 3, 1): False, (3, 1, 2): False,
                          (1, 3, 2): True, (2, 1, 3): True, (3, 2, 1): True,
                          (1, 4, 2): False, (4, 2, 1): False, (2, 1, 4): False,
                          (1, 2, 4): True, (2, 4, 1): True, (4, 1, 2): True,
                          (1, 3, 4): False, (3, 4, 1): False, (4, 1, 3): False,
                          (1, 4, 3): True, (4, 3, 1): True, (3, 1, 4): True,
                          (2, 4, 3): False, (4, 3, 2): False, (3, 2, 4): False,
                          (2, 3, 4): True, (3, 4, 2): True, (4, 2, 3): True}
# 5       4
#  \     /
#   2---3
#  /     \
# 1       6
_alkene_translate = {(1, 2, 3, 4): False, (4, 3, 2, 1): False, (1, 2, 3, 6): True, (6, 3, 2, 1): True,
                     (5, 2, 3, 6): False, (6, 3, 2, 5): False, (5, 2, 3, 4): True, (4, 3, 2, 5): True}


__all__ = ['MoleculeStereo', 'QueryStereo']
