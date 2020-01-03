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
from CachedMethods import cached_property
from collections import defaultdict
from itertools import combinations, product
from logging import info
from typing import Set, Tuple, Dict, Optional, Union
from ..exceptions import AtomNotFound, NotChiral, IsChiral


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


def _cis_trans_sign(n, u, v, w):
    # n      w
    #  \    /
    #   u--v
    #  /    \
    # x      x
    nx, ny = n
    ux, uy = u
    vx, vy = v
    wx, wy = w

    q1x = ux - nx
    q1y = uy - ny
    q2x = vx - ux
    q2y = vy - uy
    q3x = wx - vx
    q3y = wy - vy

    # cross vectors
    q1q2z = q1x * q2y - q1y * q2x
    q2q3z = q2x * q3y - q2y * q3x

    dot = q1q2z * q2q3z
    if dot > 0:
        return 1
    elif dot < 0:
        return -1
    return 0


def _allene_sign(n, u, v, w):
    # n    w
    # |   /
    # u--v
    nx, ny, nz = n
    ux, uy = u
    vx, vy = v
    wx, wy, wz = w

    q1x = ux - nx
    q1y = uy - ny
    q1z = -nz
    q2x = vx - ux
    q2y = vy - uy
    q3x = wx - vx
    q3y = wy - vy
    q3z = wz

    # cross vectors
    q1q2x = -q1z * q2y
    q1q2y = q1z * q2x
    q1q2z = q1x * q2y - q1y * q2x
    q2q3x = q2y * q3z
    q2q3y = -q2x * q3z
    q2q3z = q2x * q3y - q2y * q3x

    q1q2q3x = q1q2y * q2q3z - q1q2z * q2q3y
    q1q2q3y = q1q2z * q2q3x - q1q2x * q2q3z

    dot = q1q2q3x * q2x + q1q2q3y * q2y
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
                    if stereo[n] != other._translate_tetrahedron_sign(m, [mapping[x] for x in tetrahedrons[n]]):
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
            s = self._translate_tetrahedron_sign(n, order)
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

    def _translate_tetrahedron_sign(self, n, env):
        """
        get sign of chiral tetrahedron atom for specified neighbors order

        :param n: stereo atom
        :param env: neighbors order
        """
        s = self._atoms_stereo[n]
        order = self._tetrahedrons[n]
        if len(order) == 3:
            if len(env) == 4:  # hydrogen atom passed to env
                atoms = self._atoms
                # hydrogen always last in order
                order = (*order, next(x for x in env if atoms[x].atomic_number == 1))  # see translate scheme
            elif len(env) != 3:  # pyramid or tetrahedron expected
                raise ValueError('invalid atoms list')
        elif len(env) not in (3, 4):  # pyramid or tetrahedron expected
            raise ValueError('invalid atoms list')

        translate = tuple(order.index(x) for x in env[:3])
        if _tetrahedron_translate[translate]:
            return not s
        return s

    def _translate_cis_trans_sign(self, n, m, nn, nm):
        """
        get sign for specified opposite neighbors

        :param n: first double bonded atom
        :param m: last double bonded atom
        :param nn: neighbor of first atom
        :param nm: neighbor of last atom
        """
        cis_trans_stereo = self._cis_trans_stereo
        try:
            k = (n, m)
            s = cis_trans_stereo[k]
        except KeyError:
            k = (m, n)
            s = cis_trans_stereo[k]

        order = self._cis_trans[k]
        translate = (order.index(nn), order.index(nm))
        if _alkene_translate[translate]:
            return not s
        return s

    def _translate_allene_sign(self, c, nn, nm):
        """
        get sign for specified opposite neighbors

        :param c: central double bonded atom
        :param nn: neighbor of first double bonded atom
        :param nm: neighbor of last double bonded atom
        """
        s = self._allenes_stereo[c]
        order = self._allenes[c]
        translate = (order.index(nn), order.index(nm))
        if _alkene_translate[translate]:
            return not s
        return s

    @cached_property
    def _cumulenes(self) -> Dict[Tuple[int, ...], Tuple[int, int, Optional[int], Optional[int]]]:
        # 5       4
        #  \     /
        #   2---3
        #  /     \
        # 1       6
        bonds = self._bonds
        atoms = self._atoms
        cumulenes = {}
        for path in self.cumulenes:
            n1, m1 = path[1], path[-2]
            nn = [x for x in bonds[path[0]] if x != n1 and atoms[x].atomic_number != 1]
            mn = [x for x in bonds[path[-1]] if x != m1 and atoms[x].atomic_number != 1]
            if nn and mn:
                sn = nn[1] if len(nn) == 2 else None
                sm = mn[1] if len(mn) == 2 else None
                cumulenes[path] = (nn[0], mn[0], sn, sm)
        return cumulenes

    @cached_property
    def _tetrahedrons(self) -> Dict[int, Union[Tuple[int, int, int], Tuple[int, int, int, int]]]:
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

    @cached_property
    def _cis_trans(self) -> Dict[Tuple[int, int], Tuple[int, int, Optional[int], Optional[int]]]:
        return {(n, m): env for (n, *mid, m), env in self._cumulenes.items() if not len(mid) % 2}

    @cached_property
    def _allenes(self) -> Dict[int, Tuple[int, int, Optional[int], Optional[int]]]:
        """
        dict of allenes centers valued with neighbors of terminal atoms
        """
        return {path[len(path) // 2]: env for path, env in self._cumulenes.items() if len(path) % 2}

    @cached_property
    def _allenes_centers(self) -> Dict[int, int]:
        """
        allene terminal atom to center mapping
        """
        terminals = {}
        for path, env in self._cumulenes.items():
            if len(path) % 2:
                c = path[len(path) // 2]
                terminals[path[0]] = terminals[path[-1]] = c
        return terminals

    @cached_property
    def _allenes_terminals(self) -> Dict[int, Tuple[int, int]]:
        """
        allene center atom to terminals mapping
        """
        terminals = {}
        for path, env in self._cumulenes.items():
            if len(path) % 2:
                c = path[len(path) // 2]
                terminals[c] = (path[0], path[-1])
        return terminals


class MoleculeStereo(Stereo):
    __slots__ = ()

    def add_wedge(self, n: int, m: int, mark: bool):
        if n not in self._atoms:
            raise AtomNotFound
        if n in self._atoms_stereo:
            raise IsChiral

        plane = self._plane
        if n in self._chiral_tetrahedrons:
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
                del self.__dict__['_MoleculeStereo__chiral_centers']
        else:
            c = self._allenes_centers.get(n)
            if c:
                if c in self._allenes_stereo:
                    raise IsChiral

                order = self._allenes[c]
                t1, t2 = self._allenes_terminals[c]
                w = order.index(m)
                if w == 0:
                    m1 = order[1]
                    r = False
                elif w == 1:
                    m1 = order[0]
                    r = False
                elif w == 2:
                    m1 = order[1]
                    r = True
                else:
                    m1 = order[0]
                    r = True
                s = _dihedral_sign((*plane[m], mark), (*plane[t1], 0), (*plane[t2], 0), (*plane[m1], 0))
                if s:
                    self._allenes_stereo[c] = s < 0 if r else s > 0
                    del self.__dict__['_MoleculeStereo__chiral_centers']
            else:
                # only tetrahedrons and allenes supported
                raise NotChiral

    def calculate_cis_trans_from_2d(self):
        cis_trans_stereo = self._cis_trans_stereo
        plane = self._plane
        while True:
            if not self._chiral_cis_trans:
                break
            stereo = {}
            for nm in self._chiral_cis_trans:
                n, m = nm
                n1, m1 = self._cis_trans[nm][:2]
                s = _cis_trans_sign(plane[n1], plane[n], plane[m], plane[m1])
                if s:
                    stereo[nm] = s > 0
            if stereo:
                cis_trans_stereo.update(stereo)
                del self.__dict__['_MoleculeStereo__chiral_centers']
            else:
                break

    def add_atom_stereo(self, n, env, mark: bool):
        if n not in self._atoms:
            raise AtomNotFound
        if n in self._atoms_stereo:
            raise IsChiral
        if len(env) != len(set(env)):
            raise ValueError('invalid environment')
        if not isinstance(mark, bool):
            raise TypeError('stereo mark should be bool')

        if n in self._chiral_tetrahedrons:
            if set(env) != set(self._bonds[n]):
                raise AtomNotFound

            translate = tuple(env.index(x) for x in self._tetrahedrons[n][:3])
            if _tetrahedron_translate[translate]:
                mark = not mark

            self._atoms_stereo[n] = mark
            del self.__dict__['_MoleculeStereo__chiral_centers']
        else:  # only tetrahedrons supported
            raise NotChiral

    def _fix_stereo(self):
        if self._atoms_stereo:
            atoms = self._atoms
            stereo = {k: v for k, v in self._atoms_stereo.items() if k in atoms}
            self._atoms_stereo = new_stereo = {}

            old_stereo = 0
            while stereo and len(stereo) != old_stereo:
                old_stereo = len(stereo)
                failed_stereo = {}
                chiral = self._chiral_tetrahedrons
                del self.__dict__['_MoleculeStereo__chiral_centers']
                for n, s in stereo.items():
                    if n in chiral:
                        new_stereo[n] = s
                    else:
                        failed_stereo[n] = s
                stereo = failed_stereo

    @property
    def _chiral_tetrahedrons(self) -> Set[int]:
        return self.__chiral_centers[0]

    @property
    def _chiral_cis_trans(self) -> Set[Tuple[int, int]]:
        return self.__chiral_centers[1]

    @property
    def _chiral_allenes(self) -> Set[int]:
        return self.__chiral_centers[2]

    @cached_property
    def __chiral_centers(self):
        atoms_stereo = self._atoms_stereo
        cis_trans_stereo = self._cis_trans_stereo
        allenes_stereo = self._allenes_stereo

        morgan = self.atoms_order
        tetrahedrons = self._tetrahedrons.copy()
        cumulenes = self._cumulenes.copy()

        morgan_update = {}
        while True:
            # get obviously chiral centers
            chiral_t = {n for n, env in tetrahedrons.items() if len(set(morgan[x] for x in env)) == len(env)}
            chiral_c = set()
            chiral_a = set()
            for path, (n1, m1, n2, m2) in cumulenes.items():
                if morgan[n1] != morgan.get(n2, 0) and morgan[m1] != morgan.get(m2, 0):
                    if len(path) % 2:
                        chiral_a.add(path)
                    else:
                        chiral_c.add(path)

            if atoms_stereo:
                grouped_stereo = defaultdict(list)
                for n in chiral_t:
                    if n in atoms_stereo:
                        grouped_stereo[morgan[n]].append(n)  # collect equal stereo atoms
                for group in grouped_stereo.values():
                    if not len(group) % 2:  # only even number of equal stereo atoms give new stereo center
                        s = [n for n in group
                             if self._translate_tetrahedron_sign(n, sorted(tetrahedrons[n], key=morgan.get))]
                        if 0 < len(s) < len(group):  # RS pair required
                            for n in s:
                                morgan_update[n] = -morgan[n]
                    for n in group:  # remove seen stereo atoms
                        del tetrahedrons[n]
                        chiral_t.discard(n)

            if cis_trans_stereo:
                grouped_stereo = defaultdict(list)
                for path in chiral_c:
                    n, *_, m = path
                    mn, mm = morgan[n], morgan[m]
                    if (n, m) in cis_trans_stereo or (m, n) in cis_trans_stereo:
                        if mn <= mm:
                            grouped_stereo[mn].append((n, m, cumulenes[path], path))
                        else:
                            grouped_stereo[mm].append((m, n, cumulenes[path], path))
                for group in grouped_stereo.values():
                    if not len(group) % 2:  # only even number of equal stereo bonds give new stereo center
                        s = []
                        for n, m, (n1, m1, n2, m2), _ in group:
                            if n2 is None:
                                a = n1
                            else:
                                a = min(n1, n2, key=morgan.get)
                            if m2 is None:
                                b = m1
                            else:
                                b = min(m1, m2, key=morgan.get)
                            if self._translate_cis_trans_sign(n, m, a, b):
                                s.append(n)
                        if 0 < len(s) < len(group):  # RS pair required
                            for n in s:
                                morgan_update[n] = -morgan[n]
                    for *_, path in group:  # remove seen stereo atoms
                        del cumulenes[path]
                        chiral_c.discard(path)

            if allenes_stereo:
                grouped_stereo = defaultdict(list)
                for path in chiral_a:
                    c = path[len(path) // 2]
                    grouped_stereo[morgan[c]].append((c, cumulenes[path], path))
                for group in grouped_stereo.values():
                    if not len(group) % 2:  # only even number of equal stereo bonds give new stereo center
                        s = []
                        for c, (n1, m1, n2, m2), _ in group:
                            if n2 is None:
                                a = n1
                            else:
                                a = min(n1, n2, key=morgan.get)
                            if m2 is None:
                                b = m1
                            else:
                                b = min(m1, m2, key=morgan.get)
                            if self._translate_allene_sign(c, a, b):
                                s.append(c)
                        if 0 < len(s) < len(group):  # RS pair required
                            for c in s:
                                morgan_update[c] = -morgan[c]
                    for *_, path in group:  # remove seen stereo atoms
                        del cumulenes[path]
                        chiral_a.discard(path)

            if morgan_update:
                morgan = self._morgan(self._sorted_primed({**morgan, **morgan_update}))
                morgan_update = {}
            else:
                return chiral_t, {(n, m) for n, *_, m in chiral_c}, {path[len(path) // 2] for path in chiral_a}


class QueryStereo(Stereo):  # todo: implement add_wedge, calculate_cis_trans_from_2d
    __slots__ = ()


class NotUsed:
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


# 1  2
#  \ |
#   \|
#    n---3
#   /
#  /
# 0
_tetrahedron_translate = {(0, 1, 2): False, (1, 2, 0): False, (2, 0, 1): False,
                          (0, 2, 1): True, (1, 0, 2): True, (2, 1, 0): True,
                          (0, 3, 1): False, (3, 1, 0): False, (1, 0, 3): False,
                          (0, 1, 3): True, (1, 3, 0): True, (3, 0, 1): True,
                          (0, 2, 3): False, (2, 3, 0): False, (3, 0, 2): False,
                          (0, 3, 2): True, (3, 2, 0): True, (2, 0, 3): True,
                          (1, 3, 2): False, (3, 2, 1): False, (2, 1, 3): False,
                          (1, 2, 3): True, (2, 3, 1): True, (3, 1, 2): True}
# 2       1
#  \     /
#   n---m
#  /     \
# 0       3
_alkene_translate = {(0, 1): False, (1, 0): False, (0, 3): True, (3, 0): True,
                     (2, 3): False, (3, 2): False, (2, 1): True, (1, 2): True}


__all__ = ['MoleculeStereo', 'QueryStereo']
