# -*- coding: utf-8 -*-
#
#  Copyright 2021 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from lazy_object_proxy import Proxy
from typing import TYPE_CHECKING


if TYPE_CHECKING:
    from CGRtools import MoleculeContainer


def _freaks():
    from ...containers.query import QueryContainer
    rules = []

    q = QueryContainer()
    q.add_atom('N', neighbors=2)
    q.add_atom('A')
    q.add_atom('A')
    q.add_atom('A')
    q.add_atom('A')
    q.add_bond(1, 2, 1)
    q.add_bond(2, 3, (2, 4))
    q.add_bond(3, 4, 1)
    q.add_bond(4, 5, 4)
    q.add_bond(1, 5, 1)
    rules.append(q)
    return rules


freak_rules = Proxy(_freaks)


class Thiele:
    __slots__ = ()

    def thiele(self: 'MoleculeContainer', *, fix_tautomers=True, fix_metal_organics=True) -> bool:
        """
        Convert structure to aromatic form (Huckel rule ignored). Return True if found any kekule ring.
        Also marks atoms as aromatic.

        :param fix_tautomers: try to fix condensed rings with pyroles.
            N1C=CC2=NC=CC2=C1>>N1C=CC2=CN=CC=C12
        :param fix_metal_organics: create neutral form of ferrocenes and imidazolium complexes.
        """
        atoms = self._atoms
        bonds = self._bonds
        sh = self.hybridization
        charges = self._charges
        hydrogens = self._hydrogens

        rings = defaultdict(set)  # aromatic? skeleton. include quinones
        tetracycles = []
        pyroles = set()
        acceptors = set()
        donors = []
        freaks = []
        fixed_charges = {}
        for ring in self.sssr:
            lr = len(ring)
            if not 3 < lr < 8:  # skip 3-membered and big rings
                continue
            sp2 = sum(sh(n) == 2 and atoms[n].atomic_number in (5, 6, 7, 8, 15, 16) for n in ring)
            if sp2 == lr:  # benzene like
                if lr == 4:  # two bonds condensed aromatic rings
                    tetracycles.append(ring)
                else:
                    if fix_tautomers and lr % 2:  # find potential pyroles
                        try:
                            n = next(n for n in ring if atoms[n].atomic_number == 7 and not charges[n])
                        except StopIteration:
                            pass
                        else:
                            acceptors.add(n)
                    n, *_, m = ring
                    rings[n].add(m)
                    rings[m].add(n)
                    for n, m in zip(ring, ring[1:]):
                        rings[n].add(m)
                        rings[m].add(n)
            elif 4 < lr == sp2 + 1:  # pyroles, furanes, etc
                try:
                    n = next(n for n in ring if sh(n) == 1)
                except StopIteration:  # exotic, just skip
                    continue
                an = atoms[n].atomic_number
                if lr == 7 and an != 5:  # skip electron-rich 7-membered rings
                    continue
                elif an in (5, 7, 8, 15, 16, 34) and not charges[n]:
                    if fix_tautomers and lr == 6 and an == 7 and len(bonds[n]) == 2:
                        donors.append(n)
                    elif fix_metal_organics and lr == 5 and an == 7:
                        try:  # check for imidazolium. CN1C=C[N+](C)=C1[Cu,Ag,Au-]X
                            m = next(m for m in bonds[n] if atoms[m].atomic_number == 6 and len(bonds[m]) == 3 and
                                     all(atoms[x].atomic_number in (29, 47, 79, 46) and charges[x] < 0 or
                                         atoms[x].atomic_number == 7 and charges[x] == 1
                                         for x in bonds[m] if x != n))
                        except StopIteration:
                            pass
                        else:
                            for x in bonds[m]:
                                if charges[x] < 0:
                                    if x in fixed_charges:
                                        fixed_charges[x] += 1
                                    else:
                                        fixed_charges[x] = charges[x] + 1
                                else:
                                    fixed_charges[x] = 0
                    pyroles.add(n)
                    n, *_, m = ring
                    rings[n].add(m)
                    rings[m].add(n)
                    for n, m in zip(ring, ring[1:]):
                        rings[n].add(m)
                        rings[m].add(n)
                elif an == 6 and lr == 5 and charges[n] == -1:  # ferrocene, etc.
                    if fix_metal_organics:
                        try:
                            m = next(m for m, b in bonds[n].items() if b == 8 and charges[m] > 0 and
                                     atoms[m].atomic_number in
                                     (22, 23, 24, 25, 26, 27, 28, 40, 41, 42, 44, 72, 74, 75, 77))
                        except StopIteration:
                            pass
                        else:
                            fixed_charges[n] = 0  # remove charges in thiele form
                            if m in fixed_charges:
                                fixed_charges[m] -= 1
                            else:
                                fixed_charges[m] = charges[m] - 1
                    pyroles.add(n)
                    n, *_, m = ring
                    rings[n].add(m)
                    rings[m].add(n)
                    for n, m in zip(ring, ring[1:]):
                        rings[n].add(m)
                        rings[m].add(n)
            # like N1C=Cn2cccc12
            elif lr == 5 and sum(atoms[x].atomic_number == 7 and not charges[x] for x in ring) > 1:
                freaks.append(ring)
        if not rings:
            return False
        double_bonded = {n for n in rings if any(m not in rings and b.order == 2 for m, b in bonds[n].items())}

        # fix_tautomers
        if fix_tautomers and acceptors and donors:
            for start in donors:
                stack = [(start, n, 0, 2) for n in rings[start] if n not in double_bonded]
                path = []
                seen = {start}
                while stack:
                    last, current, depth, order = stack.pop()
                    if len(path) > depth:
                        seen.difference_update(x for _, x, _ in path[depth:])
                        path = path[:depth]
                    path.append((last, current, order))
                    if current in acceptors:  # we found
                        if order == 1:
                            acceptors.discard(current)
                            pyroles.discard(start)
                            pyroles.add(current)
                            hydrogens[current] = 1
                            hydrogens[start] = 0
                            break
                        else:
                            continue

                    depth += 1
                    seen.add(current)
                    new_order = 1 if order == 2 else 2
                    stack.extend((current, n, depth, new_order) for n in rings[current] if
                                 n not in seen and n not in double_bonded and bonds[current][n].order == order)
                for n, m, o in path:
                    bonds[n][m]._Bond__order = o
                if not acceptors:
                    break

        if double_bonded:  # delete quinones
            for n in double_bonded:
                for m in rings.pop(n):
                    rings[m].discard(n)

            for n in [n for n, ms in rings.items() if not ms]:  # imide leads to isolated atoms
                del rings[n]
            if not rings:
                return False
            while True:
                try:
                    n = next(n for n, ms in rings.items() if len(ms) == 1)
                except StopIteration:
                    break
                m = rings.pop(n).pop()
                if n in pyroles:
                    rings[m].discard(n)
                else:
                    pm = rings.pop(m)
                    pm.discard(n)
                    for x in pm:
                        rings[x].discard(m)
        if not rings:
            return False

        n_sssr = sum(len(x) for x in rings.values()) // 2 - len(rings) + len(self._connected_components(rings))
        if not n_sssr:
            return False
        rings = self._sssr(rings, n_sssr)  # search rings again

        seen = set()
        for ring in rings:
            seen.update(ring)
        charges.update(fixed_charges)

        # reset bonds to single
        for ring in tetracycles:
            if seen.issuperset(ring):
                n, *_, m = ring
                bonds[n][m]._Bond__order = 1
                for n, m in zip(ring, ring[1:]):
                    bonds[n][m]._Bond__order = 1

        for ring in rings:
            n, *_, m = ring
            bonds[n][m]._Bond__order = 4
            for n, m in zip(ring, ring[1:]):
                bonds[n][m]._Bond__order = 4

        self.flush_cache()
        for ring in freaks:  # aromatize rule based
            rs = set(ring)
            for q in freak_rules:
                # used low-level API for speedup
                components, closures = q._compiled_query
                if any(q._get_mapping(components[0], closures, atoms, bonds, rs, self.atoms_order)):
                    n, *_, m = ring
                    bonds[n][m]._Bond__order = 4
                    for n, m in zip(ring, ring[1:]):
                        bonds[n][m]._Bond__order = 4
                    break
        if freaks:
            self.flush_cache()  # flush again
        self._fix_stereo()  # check if any stereo centers vanished.
        return True


__all__ = ['Thiele']
