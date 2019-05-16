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
from ..periodictable import C


_pyrole_atoms = ('N', 'O', 'S', 'Se', 'P', C(-1))


class Aromatize:
    def dummy_aromatize(self):
        """
        convert structure to aromatic form (dummy algorithm. don't detect quinones)

        :return: number of processed rings
        """
        adj = self._bonds
        atom = self._atoms
        total = 0
        unsaturated = {n for n, m_bond in adj.items() if any(bond.order in (2, 4) for bond in m_bond.values())}

        for ring in self.sssr:
            lr = len(ring)
            if lr in (5, 6, 7) and unsaturated.issuperset(ring):
                for n, m in zip(ring, ring[1:]):
                    b = adj[n][m]
                    if b.order != 4:
                        b.order = 4
                b = adj[ring[0]][ring[-1]]
                if b.order != 4:
                    b.order = 4
                total += 1
            elif lr == 5:
                sr = set(ring)
                if len(unsaturated & sr) == 4 and atom[(sr - unsaturated).pop()]._atom in _pyrole_atoms:
                    for n, m in zip(ring, ring[1:]):
                        b = adj[n][m]
                        if b.order != 4:
                            b.order = 4
                    b = adj[ring[0]][ring[-1]]
                    if b.order != 4:
                        b.order = 4
                    total += 1
        if total:
            self.flush_cache()
        return total

    def aromatize(self) -> int:
        """
        convert structure to aromatic form

        :return: number of processed rings
        """
        total = self.dummy_aromatize()
        adj = self._bonds
        atom = self._atoms
        patch = set()
        double_bonded = {n for n, m_bond in adj.items() if any(bond.order == 2 for bond in m_bond.values())}

        pyroles = set()
        quinones = []
        azulenes = set()
        condensed_rings = defaultdict(lambda: defaultdict(list))
        for ring in self.aromatic_rings:
            ring = tuple(ring)
            if not double_bonded.isdisjoint(ring):  # search quinones
                quinones.append(ring)

            lr = len(ring)
            if lr == 5:
                pyroles.update(n for n in ring if atom[n]._atom in _pyrole_atoms)
            elif lr == 7:
                azulenes.update(ring)

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
            if lr == 7:
                unbalansed_ring = len(doubles) in (2, 4)  # bis- or tetra- azulene 7-ring quinones
            elif lr == 6:
                unbalansed_ring = len(doubles) % 2
            elif not azulenes.isdisjoint(ring):
                unbalansed_ring = len(doubles) == 2  # bis- azulene 5-ring quinones
            else:  # pyroles
                unbalansed_ring = len(doubles) == 1

            bond = 1
            n = ordered_ring[0]
            for m in ordered_ring[1:]:
                if bond == 1:
                    if not (m in double_bonded or m in pyroles or
                            unbalansed_ring and condensed_rings[n][m] and not condensed_rings[n][p]):
                        bond = 2
                        # single bond followed by double if common atom:
                        # not already has double bond [quinone] or
                        # not pyrole atom with LP or
                        # unbalansed_ring has 2 condensed rings in a row
                    if not condensed_rings[n][m]:
                        patch.add((n, m, 1))
                    elif n in double_bonded:  # found new quinone ring (Y)
                        q = condensed_rings[n][m][0]
                        if q not in quinones:
                            quinones.insert(0, q)  # low priority
                else:
                    if m in double_bonded:
                        raise InvalidAromaticRing(ring)
                    if not (unbalansed_ring and condensed_rings[n][m] and not condensed_rings[n][p]):
                        bond = 1
                    if not condensed_rings[n][m]:
                        patch.add((n, m, 2))
                        double_bonded.add(n)
                        double_bonded.add(m)
                        if condensed_rings[n][p]:
                            q = condensed_rings[n][p][0]
                            if q in quinones:  # up priority
                                quinones.remove(q)
                                quinones.append(q)
                            else:
                                quinones.insert(0, q)
                p, n = n, m
            else:
                m = ordered_ring[0]
                if bond != 1:
                    raise InvalidAromaticRing(ring)
                patch.add((n, m, 1))

        if patch:
            for n, m, b in patch:
                adj[n][m].order = b
            self.flush_cache()
        return total

    def dearomatize(self):
        raise NotImplementedError
        adj = defaultdict(set)  # aromatic skeleton
        for n, m_bond in self._bonds.items():
            for m, bond in m_bond.items():
                if bond.order == 4:
                    adj[n].add(m)


__all__ = ['Aromatize']
