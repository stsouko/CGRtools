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
from ..exceptions import InvalidAromaticRing


_pyrole_atoms = ('N', 'O', 'S', 'Se', 'P')


class Aromatize:
    __slots__ = ()

    def dummy_aromatize(self):
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

    def aromatize(self) -> int:
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

    def dearomatize(self):
        rings = defaultdict(set)  # aromatic skeleton
        double_bonded = set()
        for n, m_bond in self._bonds.items():
            for m, bond in m_bond.items():
                if bond.order == 4:
                    rings[n].add(m)
                elif bond.order == 2 and n not in double_bonded:
                    double_bonded.add(n)
        atoms = set(rings)
        double_bonded &= rings.keys()
        while atoms:
            if double_bonded:
                enter = double_bonded.pop()
                atoms.discard(enter)
            else:
                enter = atoms.pop()


__all__ = ['Aromatize']
