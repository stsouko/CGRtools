# -*- coding: utf-8 -*-
#
#  Copyright 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from typing import Tuple


class Pharmacophore:
    __slots__ = ()

    @cached_property
    def hydrogen_acceptors(self) -> Tuple[int, ...]:
        """
        Atom that could form hydrogen bonds.

        N: R-N(-[H,R,Ar])-[H,R], Ar:N:Ar, Ar-NH2
        O: R-O-[H,R,Ar], [R,Ar]-C(=O)-[H,R,Ar,OH,OR,NR2], R-C(=O)-[O-], Ar-[O-], Ar-OH
        """
        atoms = self._atoms
        bonds = self._bonds
        charges = self._charges
        hybridizations = self._hybridizations
        neighbors = self._neighbors

        out = []
        for n, a in atoms.items():
            if not neighbors[n]:  # not charged and not water or ammonia
                continue
            if not charges[n]:
                if a.atomic_number == 7:  # N
                    if hybridizations[n] == 1:
                        if neighbors[n] == 1:  # ArNH2 or RNH2
                            if all(hybridizations[m] in (1, 4) for m in bonds[n]):
                                out.append(n)
                        elif all(hybridizations[m] not in (2, 3) for m in bonds[n]):  # secondary and third amines
                            if sum(hybridizations[m] == 4 for m in bonds[n]) <= 1:  # only one Aryl neighbor
                                out.append(n)
                    # pyridine, imidazole
                    elif hybridizations[n] == 4:
                        if neighbors[n] == 2:
                            out.append(n)
                elif a.atomic_number == 8:  # O
                    if hybridizations[n] == 2:  # carbonyl
                        out.append(n)
                    elif all(hybridizations[m] not in (2, 3) for m in bonds[n]):  # [R,Ar]-O-[H,R]
                        if sum(hybridizations[m] == 4 for m in bonds[n]) <= 1:  # only one Aryl neighbor
                            out.append(n)
            elif charges[n] == -1:  # R-C(=O)-[O-] or Ar-[O-]
                if a.atomic_number == 8 and not any(charges[m] > 0 for m in bonds[n]):
                    # exclude zwitterions: nitro etc
                    out.append(n)
        return tuple(out)

    @cached_property
    def halogen_donors(self) -> Tuple[Tuple[int, int], ...]:
        """
        Halogen(I)-Carbon pairs.
        """
        atoms = self._atoms
        bonds = self._bonds
        out = []
        for n, a in atoms.items():
            if a.atomic_number in (9, 17, 35, 53):
                env = bonds[n]
                if len(env) == 1:
                    m = next(iter(env))
                    if atoms[m].atomic_number == 6:
                        out.append((n, m))
        return tuple(out)

    @cached_property
    def halogen_acceptors(self) -> Tuple[Tuple[int, int], ...]:
        """
        (O,N,S)-(C,N,P,S) pairs.
        """
        bonds = self._bonds
        atoms = self._atoms
        charges = self._charges
        out = []
        for n, a in atoms.items():
            if a.atomic_number in (7, 8, 16) and not charges[n]:
                env = [m for m in bonds[n] if atoms[m].atomic_number in (6, 7, 15, 16)]
                if len(env) == 1:
                    out.append((n, env[0]))
        return tuple(out)

    @cached_property
    def positive_charged(self) -> Tuple[int, ...]:
        """
        Atoms with positive formal charge, except zwitterions
        """
        bonds = self._bonds
        charges = self._charges
        return tuple(n for n, c in charges.items() if c > 0 and not any(charges[m] < 0 for m in bonds[n]))

    @cached_property
    def negative_charged(self) -> Tuple[int, ...]:
        """
        Atoms with negative formal charge, except zwitterions
        """
        bonds = self._bonds
        charges = self._charges
        return tuple(n for n, c in charges.items() if c < 0 and not any(charges[m] > 0 for m in bonds[n]))

    @cached_property
    def metal_ligands(self):
        """
        Atoms that could possibly be involved in chelating a metal ion.
        O: R-OH, R-O-R, C(=O)-OH, C(=O)-O-R, R-C(=O)H, R-C(=O)-R C(=O)-NR2
        N: N(Ar, not pyrole), R-NR2
        S: R-SH, R-S-R
        """
        return


__all__ = ['Pharmacophore']
