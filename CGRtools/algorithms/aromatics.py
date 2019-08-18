# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
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
from itertools import product
from typing import List, Tuple, Optional
from ..exceptions import InvalidAromaticRing


_pyrole_atoms = ('N', 'O', 'S', 'Se', 'P')


class Aromatize:
    __slots__ = ()

    def dummy_thiele(self):
        """
        convert structure to aromatic form (dummy algorithm. don't detect quinones)

        :return: number of processed rings
        """
        bonds = self._bonds
        atoms = self._atoms
        unsaturated = {n for n, m_bond in bonds.items() if any(bond.order in (2, 4) for bond in m_bond.values())}
        aromatics = []

        for ring in self.sssr:
            lr = len(ring)
            if lr in (5, 6, 7):
                if unsaturated.issuperset(ring):  # benzene, azulene, pyridine and quinones
                    aromatics.append(ring)
                else:  # pyrole, cyclopentapyridine and quinones
                    sr = set(ring)
                    if len(unsaturated & sr) == lr - 1:
                        atom = atoms[(sr - unsaturated).pop()]
                        if atom.atomic_symbol in _pyrole_atoms or atom.charge == -1 and atom.atomic_symbol == 'C':
                            aromatics.append(ring)

        if aromatics:
            for ring in aromatics:
                for n, m in zip(ring, ring[1:]):
                    b = bonds[n][m]
                    if b.order != 4:
                        b._Bond__order = 4
                b = bonds[ring[0]][ring[-1]]
                if b.order != 4:
                    b._Bond__order = 4
            self.flush_cache()
            return len(aromatics)
        return 0

    def thiele(self) -> int:
        """
        convert structure to aromatic form

        :return: number of processed rings
        """
        total = self.dummy_aromatize()
        bonds = self._bonds
        atoms = self._atoms
        patch = set()
        double_bonded = {n for n, m_bond in bonds.items() if any(bond.order == 2 for bond in m_bond.values())}

        pyroles = set()
        quinones = []
        condensed_rings = defaultdict(lambda: defaultdict(list))
        for ring in self.aromatic_rings:
            ring = tuple(ring)
            lr = len(ring)
            if lr in (5, 6):
                pyroles.update(n for n in ring if atoms[n].atomic_symbol in _pyrole_atoms or
                               atoms[n].charge == -1 and atoms[n].atomic_symbol == 'C')
            if not double_bonded.isdisjoint(ring):
                quinones.append(ring)

            for n, m in zip(ring, ring[1:]):  # fill condensed rings graph
                condensed_rings[n][m].append(ring)
                condensed_rings[m][n].append(ring)
            n, *_, m = ring
            condensed_rings[n][m].append(ring)
            condensed_rings[m][n].append(ring)

        while quinones:
            ring = quinones.pop()
            for n, m in zip(ring, ring[1:]):  # remove from condensed rings graph
                condensed_rings[n][m].remove(ring)
                condensed_rings[m][n].remove(ring)
            n, *_, m = ring
            condensed_rings[n][m].remove(ring)
            condensed_rings[m][n].remove(ring)

            doubles = [n for n, m in enumerate(ring) if m in double_bonded]
            start = doubles[0]
            if start:  # reorder double bonded to starting position
                ordered_ring = ring[start:] + ring[:start]
            else:
                ordered_ring = ring

            lr = len(ring)
            if lr == 6:
                unbalanced_ring = len(doubles) in (1, 3)
            else:  # bis- or tetra- azulene 7-ring and pyrole or bis azulenes 5-rings quinones
                unbalanced_ring = len(doubles) in (2, 4)

            pyrole_like = not pyroles.isdisjoint(ring)

            bond = 1
            n = ordered_ring[0]
            for m in ordered_ring[1:]:
                if bond == 1:
                    if m not in double_bonded:
                        if unbalanced_ring:
                            if pyrole_like:
                                if m in pyroles:
                                    unbalanced_ring = False
                                else:
                                    bond = 2
                            elif condensed_rings[n][m]:
                                unbalanced_ring = False
                            else:
                                bond = 2
                        else:
                            bond = 2
                    if not condensed_rings[n][m]:
                        patch.add((n, m, 1))
                    elif n in double_bonded:  # found new quinone ring (Y)
                        q = condensed_rings[n][m][0]
                        if q not in quinones:
                            quinones.append(q)
                else:
                    if m in double_bonded:
                        raise InvalidAromaticRing(ring)
                    if unbalanced_ring and not pyrole_like and condensed_rings[n][m]:
                        unbalanced_ring = False
                    else:
                        bond = 1
                    if not condensed_rings[n][m]:
                        patch.add((n, m, 2))
                        double_bonded.add(n)
                        double_bonded.add(m)
                        if condensed_rings[n][p]:
                            q = condensed_rings[n][p][0]
                            if q not in quinones:
                                quinones.append(q)
                p, n = n, m
            else:
                m = ordered_ring[0]
                if bond != 1:
                    raise InvalidAromaticRing(ring)
                patch.add((n, m, 1))

        if patch:
            for n, m, b in patch:
                bonds[n][m]._Bond__order = b
            self.flush_cache()
        return total

    def kekule(self):
        """
        convert structure to kekule form.

        only one of possible double/single bonds positions will be set.
        for enumerate bonds positions use `enumerate_kekule`
        """
        kekule = next(self.__kekule_full(), None)
        if kekule:
            self.__kekule_patch(kekule)
            self.flush_cache()

    def enumerate_kekule(self):
        """
        enumerate all possible kekule forms of molecule
        """
        for form in self.__kekule_full():
            copy = self.copy()
            copy._Aromatize__kekule_patch(form)
            yield copy

    def __kekule_patch(self, patch):
        bonds = self._bonds
        for n, m, b in patch:
            bonds[n][m]._Bond__order = b

    def __kekule_full(self):
        atoms = self._atoms
        charges = self._charges
        radicals = self._radicals
        bonds = self._bonds

        rings = defaultdict(set)  # aromatic skeleton
        double_bonded = set()
        for n, m_bond in bonds.items():
            for m, bond in m_bond.items():
                if bond.order == 4:
                    rings[n].add(m)
                elif bond.order == 2:
                    double_bonded.add(n)
        if not rings:
            return
        elif any(len(ms) not in (2, 3) for ms in rings.values()):
            raise InvalidAromaticRing('not in ring aromatic bond or hypercondensed rings')

        double_bonded &= rings.keys()
        if any(len(rings[n]) != 2 for n in double_bonded):  # double bonded never condensed
            raise InvalidAromaticRing('quinone valence error')
        if any(atoms[n].atomic_number not in (6, 15, 16, 24) or charges[n] for n in double_bonded):
            raise InvalidAromaticRing('quinone should be neutral S, Se, C, P atom')

        pyroles = set()
        for n, ms in rings.items():
            an = atoms[n].atomic_number
            ac = charges[n]
            if an == 6:  # carbon
                if ac == 0:
                    continue
                elif ac == -1:
                    if radicals[n]:
                        if len(bonds[n]) == 2:  # anion-radical
                            double_bonded.add(n)
                        else:
                            raise InvalidAromaticRing
                    else:
                        pyroles.add(n)
                elif ac != 1 or radicals[n] or len(bonds[n]) != 2:  # not benzene cation
                    raise InvalidAromaticRing
            elif an in (7, 15):
                if ac == 0:  # pyrole or pyridine. include radical pyrole
                    if radicals[n] and len(bonds[n]) != 2:
                        raise InvalidAromaticRing
                    elif len(bonds[n]) == 3:  # pyrole or P-oxyde only possible
                        double_bonded.add(n)
                    else:
                        pyroles.add(n)
                elif ac == -1:  # pyrole only
                    if radicals[n] or len(bonds[n]) != 2:
                        raise InvalidAromaticRing
                    double_bonded.add(n)
                elif ac != 1:
                    raise InvalidAromaticRing
                elif radicals[n]:
                    if len(bonds[n]) != 2:  # not cation-radical pyridine
                        raise InvalidAromaticRing
                elif len(bonds[n]) == 2:  # pyrole cation
                    double_bonded.add(n)
                elif len(bonds[n]) != 3:  # not pyridine oxyde
                    raise InvalidAromaticRing
            elif an == 8:  # furan
                if len(bonds[n]) == 2 and (ac == 0 and not radicals[n] or ac == 1 and radicals[n]):
                    double_bonded.add(n)
                else:
                    raise InvalidAromaticRing
            elif an in (16, 24):  # thiophene or sulphoxyde or sulphone
                if n not in double_bonded:
                    if len(bonds[n]) == 2 and (ac == 0 and not radicals[n] or ac == 1 and radicals[n]):
                        double_bonded.add(n)
                    else:
                        raise InvalidAromaticRing
            elif an == 5:  # boron
                if ac == 0:
                    if len(bonds[n]) == 3 and not radicals[n] or len(bonds[n]) == 2:
                        double_bonded.add(n)
                    else:
                        raise InvalidAromaticRing
                elif ac in (-1, 1) and len(bonds[n]) == 2 and not radicals[n]:
                    double_bonded.add(n)
                else:
                    raise InvalidAromaticRing
            else:
                raise InvalidAromaticRing(f'only B, C, N, P, O, S, Se possible, not: {atoms[n].atomic_symbol}')

        atoms = set(rings)
        components = []
        while atoms:
            start = atoms.pop()
            component = {n: rings[n] for n in self.__component(rings, start)}
            components.append(component)
            atoms.difference_update(component)

        for keks in product(*(self.__kekule_component(c, double_bonded & c.keys(), pyroles & c.keys())
                              for c in components)):
            yield list(*keks)

    @staticmethod
    def __kekule_component(rings, double_bonded, pyroles):
        # (current atom, previous atom, bond between cp atoms, path deep for cutting [None if cut impossible])
        stack: List[List[Tuple[int, int, int, Optional[int]]]]
        if double_bonded:  # start from double bonded if exists
            start = next(iter(double_bonded))
            stack = [[(next(iter(rings[start])), start, 1, 0)]]
        else:  # select not pyrole not condensed atom
            start = next(n for n, ms in rings.items() if len(ms) == 2 and n not in pyroles)
            stack = [[(next_atom, start, 1, 0)] for next_atom in rings[start]]

        size = sum(len(x) for x in rings.values()) // 2
        path = []
        hashed_path = set()
        while stack:
            atom, prev_atom, bond, _ = stack[-1].pop()
            path.append((atom, prev_atom, bond))
            hashed_path.add(atom)

            if len(path) == size:
                yield path
                del stack[-1]
                if stack:
                    path = path[:stack[-1][-1][-1]]
                    hashed_path = {x for x, *_ in path}
            elif atom != start:  # we finished. next step is final closure
                for_stack = []
                closures = []
                loop = 0
                for next_atom in rings[atom]:
                    if next_atom == prev_atom:  # only forward. behind us is the homeland
                        continue
                    elif next_atom == start:
                        loop = next_atom
                    elif next_atom in hashed_path:  # closure found
                        closures.append(next_atom)
                    else:
                        for_stack.append(next_atom)

                if loop:  # we found starting point.
                    if bond == 2:  # finish should be single bonded
                        if double_bonded:  # ok
                            stack[-1].insert(0, (loop, atom, 1, None))
                            if not for_stack:
                                continue
                        else:
                            del stack[-1]
                            if stack:
                                path = path[:stack[-1][-1][-1]]
                                hashed_path = {x for x, *_ in path}
                            continue
                    elif atom in double_bonded:  # we on quinone atom. finish should be single bonded
                        stack[-1].insert(0, (loop, atom, 1, None))
                    elif double_bonded:  # we in quinone ring. finish should be single bonded
                        if for_stack:  # need path for storing double bond
                            stack[-1].insert(0, (loop, atom, 1, None))
                        else:
                            del stack[-1]
                            if stack:
                                path = path[:stack[-1][-1][-1]]
                                hashed_path = {x for x, *_ in path}
                            continue
                    else:  # finish should be double bonded
                        stack[-1].insert(0, (loop, atom, 2, None))
                        if not for_stack:
                            continue
                        bond = 2  # grow should be single bonded

                if bond == 2 or atom in double_bonded:  # double in - single out. quinone has two single bonds
                    for next_atom in closures:
                        path.append((next_atom, atom, 1))  # closures always single-bonded
                        stack[-1].remove((atom, next_atom, 1, None))  # remove fork from stack
                    for next_atom in for_stack:
                        stack[-1].append((next_atom, atom, 1, None))
                elif len(for_stack) == 1:  # easy path grow. next bond double
                    next_atom = for_stack[0]
                    if next_atom in double_bonded:  # need double bond, but next atom quinone
                        del stack[-1]
                        if stack:
                            path = path[:stack[-1][-1][-1]]
                            hashed_path = {x for x, *_ in path}
                    else:
                        stack[-1].append((next_atom, atom, 2, None))
                        for next_atom in closures:
                            path.append((next_atom, atom, 1))  # closures always single-bonded
                            stack[-1].remove((atom, next_atom, 1, None))  # remove fork from stack
                elif for_stack:  # fork
                    next_atom1, next_atom2 = for_stack
                    if next_atom1 in double_bonded:  # quinone next from fork
                        if next_atom2 in double_bonded:  # bad path
                            del stack[-1]
                            if stack:
                                path = path[:stack[-1][-1][-1]]
                                hashed_path = {x for x, *_ in path}
                        else:
                            stack[-1].append((next_atom1, atom, 1, None))
                            stack[-1].append((next_atom2, atom, 2, None))
                    elif next_atom2 in double_bonded:  # quinone next from fork
                        stack[-1].append((next_atom2, atom, 1, None))
                        stack[-1].append((next_atom1, atom, 2, None))
                    else:  # new path
                        opposite = stack[-1].copy()
                        stack[-1].append((next_atom1, atom, 1, None))
                        stack[-1].append((next_atom2, atom, 2, len(path)))
                        opposite.append((next_atom2, atom, 1, None))
                        opposite.append((next_atom1, atom, 2, None))
                        stack.append(opposite)
                elif closures:  # need double bond, but closure should be single bonded
                    del stack[-1]
                    if stack:
                        path = path[:stack[-1][-1][-1]]
                        hashed_path = {x for x, *_ in path}

    @staticmethod
    def __component(bonds, start):
        seen = {start}
        queue = [start]
        while queue:
            start = queue.pop(0)
            yield start
            for i in bonds[start] - seen:
                queue.append(i)
                seen.add(i)


__all__ = ['Aromatize']
