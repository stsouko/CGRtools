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
from itertools import product
from math import sqrt, acos
from typing import Tuple, Optional, List, NamedTuple, Union


class Hydrophobic(NamedTuple):
    source: int
    target: int
    distance: float


class HydrogenDonor(NamedTuple):
    source: int
    target: int
    distance: float
    acceptor_distance: float  # acceptor-H bond length
    donor_distance: float  # donor-H bond length
    angle: float


class HydrogenAcceptor(NamedTuple):
    source: int
    target: int
    distance: float
    acceptor_distance: float  # acceptor-H bond length
    donor_distance: float  # donor-H bond length
    angle: float


def distance(n, m):
    nx, ny, nz = n
    mx, my, mz = m
    return sqrt((nx - mx) ** 2 + (ny - my) ** 2 + (nz - mz) ** 2)


def angle(n, m, k):
    """
    m----n
          \
           k
    """
    nx, ny, nz = n
    mx, my, mz = m
    kx, ky, kz = k
    nmx, nmy, nmz = mx - nx, my - ny, mz - nz
    nkx, nky, nkz = kx - nx, ky - ny, kz - nz
    return acos((nmx * nkx + nmy * nky + nmz * nkz) / sqrt(nmx ** 2 + nmy ** 2 + nmz ** 2))


class Pharmacophore:
    __slots__ = ()

    def find_contacts(self, other: 'Pharmacophore') -> List[Union[Hydrophobic, HydrogenDonor, HydrogenAcceptor]]:
        """
        Find pharmocophore contacts in 3d. Required explicit hydrogens.
        """
        contacts = []
        contacts.extend(self.__hydrophobic_interactions(other))
        return contacts

    def __hydrophobic_interactions(self, other):
        """
        Detection of hydrophobic contacts.

        Definition: All pairs of qualified carbon atoms within a distance of hydroph_dist_max
        """
        config = self._pharmacophore_config
        min_dist = config['min_dist']
        hydroph_dist_max = config['hydroph_dist_max']

        contacts = []
        s_xyz = self._conformers[0]
        o_xyz = other._conformers[0]
        for n, m in product(self.hydrophobic_centers, other.hydrophobic_centers):
            d = distance(s_xyz[n], o_xyz[m])
            if min_dist < d < hydroph_dist_max:
                contacts.append(Hydrophobic(n, m, d))
        return contacts

    def __hydrogen_bond_interactions(self, other):
        """
        Detection of hydrogen bonds between sets of acceptors and donor pairs.

        Definition: All pairs of hydrogen bond acceptor and donors with
        donor hydrogens and acceptor showing a distance within hbond_dist_max
        and donor angles above hbond_don_angle_min
        """
        config = self._pharmacophore_config
        min_dist = config['min_dist']
        hbond_dist_max = config['hbond_dist_max']
        hbond_don_angle_min = config['hbond_don_angle_min']

        seen = set()
        contacts = []
        s_xyz = self._conformers[0]
        o_xyz = other._conformers[0]
        for (n, h), m in product(self.hydrogen_donors, other.hydrogen_acceptors):
            dist_nm = distance(s_xyz[n], o_xyz[m])
            if min_dist < dist_nm < hbond_dist_max:
                a = angle(s_xyz[h], s_xyz[n], o_xyz[m])
                if a >= hbond_don_angle_min:
                    seen.add((n, m))
                    contacts.append(HydrogenDonor(n, m, dist_nm, distance(s_xyz[h], o_xyz[m]),
                                                  distance(s_xyz[n], s_xyz[h]), a))

        for n, (m, h) in product(self.hydrogen_acceptors, other.hydrogen_donors):
            if (n, m) in seen:
                continue
            dist_nm = distance(s_xyz[n], o_xyz[m])
            if min_dist < dist_nm < hbond_dist_max:
                a = angle(o_xyz[h], o_xyz[m], s_xyz[n])
                if a >= hbond_don_angle_min:
                    contacts.append(HydrogenAcceptor(n, m, dist_nm, distance(o_xyz[h], s_xyz[n]),
                                                     distance(o_xyz[m], o_xyz[h]), a))
        return contacts

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
                    # imidazole, pyrroline, guanidine
                    elif hybridizations[n] == 2:
                        out.append(n)
                    # pyridine
                    elif hybridizations[n] == 4:
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
    def hydrogen_donors(self) -> Tuple[Tuple[int, Optional[int]], ...]:
        """
        NH, OH, SH groups.
        """
        atoms = self._atoms
        bonds = self._bonds
        charges = self._charges
        hybridizations = self._hybridizations
        neighbors = self._neighbors

        out = []
        for n, a in atoms.items():
            if hybridizations[n] == 1:  # amines, alcoholes, phenoles, thioles
                if a.atomic_number in (8, 16) and not charges[n] and neighbors[n] == 1 or \
                        a.atomic_number == 7 and (not charges[n] and neighbors[n] in (1, 2) or
                                                  charges[n] == 1 and neighbors[n] in (1, 2, 3)):
                    out.append((n, next((m for m in bonds[n] if atoms[m].atomic_number == 1), None)))
            # imine, guanidine
            elif hybridizations[n] == 2 and a.atomic_number == 7 and not charges[n] and neighbors[n] == 1:
                out.append((n, next((m for m in bonds[n] if atoms[m].atomic_number == 1), None)))
        return tuple(out)

    @cached_property
    def halogen_donors(self) -> Tuple[Tuple[int, int], ...]:
        """
        Carbon - Halogen(I) pairs.
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
                        out.append((m, n))
        return tuple(out)

    @cached_property
    def halogen_acceptors(self) -> Tuple[Tuple[int, int], ...]:
        """
        (C,N,P,S) - Terminal(O,N,S) pairs.
        """
        bonds = self._bonds
        atoms = self._atoms
        charges = self._charges
        out = []
        for n, a in atoms.items():
            if a.atomic_number in (7, 8, 16) and charges[n] <= 0:
                env = [m for m in bonds[n] if atoms[m].atomic_number in (6, 7, 15, 16)]
                if len(env) == 1:
                    out.append((env[0], n))
        return tuple(out)

    @cached_property
    def positive_charged_centers(self) -> Tuple[int, ...]:
        """
        Atoms with positive formal charge, except zwitterions.
        Guanidines-H+ and same ions which have delocalization of charge will be added fully.
        Supported delocalization between S, O, N atoms
        """
        atoms = self._atoms
        bonds = self._bonds
        charges = self._charges
        neighbors = self._neighbors
        hybridizations = self._hybridizations

        out = set()
        for n, c in charges.items():
            if not neighbors[n]:  # skip [NH4-], etc
                continue
            if c > 0 and not any(charges[m] < 0 for m in bonds[n]):
                out.add(n)
                if atoms[n].atomic_number == 7 and hybridizations[n] == 2:
                    m = next(m for m, b in bonds[n].items() if b.order == 2)
                    for x, b in bonds[m].items():
                        if x != n and b.order == 1 and not charges[x] and atoms[x].atomic_number in (7, 8, 16):
                            out.add(x)
        return tuple(out)

    @cached_property
    def negative_charged_centers(self) -> Tuple[int, ...]:
        """
        Atoms with negative formal charge, except zwitterions.

        Carboxyles and same ions which have delocalization of charge will be added fully.
        Supported delocalization between S, O, N atoms
        """
        atoms = self._atoms
        bonds = self._bonds
        charges = self._charges
        neighbors = self._neighbors

        out = set()
        for n, c in charges.items():
            if not neighbors[n]:  # skip [OH-], etc
                continue
            if c < 0 and not any(charges[m] > 0 for m in bonds[n]):
                out.add(n)
                if atoms[n].atomic_number in (8, 16):
                    m = next(iter(bonds[n]))
                    for x, b in bonds[m].items():
                        if x != n and b.order == 2 and not charges[x] and atoms[x].atomic_number in (7, 8, 16):
                            out.add(x)
        return tuple(out)

    @cached_property
    def hydrophobic_centers(self) -> Tuple[int, ...]:
        """
        Hydrophobic atoms.

        Hydrocarbons atoms.
        """
        atoms = self._atoms
        bonds = self._bonds
        charges = self._charges

        out = []
        for n, a in atoms.items():
            if a.atomic_number == 6 and not charges[n] and all(atoms[m].atomic_number in (1, 6) for m in bonds[n]):
                out.append(n)
        return tuple(out)

    @cached_property
    def aromatic_centers(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Aromatic rings atoms.

        Ignored compounds like N1C=CC=C2C=CC=C12.
        """
        atoms = self._atoms
        bonds = self._bonds
        charges = self._charges
        hybridizations = self._hybridizations
        aroma = self.aromatic_rings

        # search pyroles
        pyroles = []
        for ring in self.sssr:
            if len(ring) == 5 and ring not in aroma and all(not charges[n] for n in ring):
                sp2 = 0
                pa = None
                for n in ring:
                    if hybridizations[n] in (2, 4):
                        sp2 += 1
                    elif atoms[n].atomic_number in (7, 8, 15, 16):
                        pa = n
                if sp2 == 4 and pa:
                    pyroles.append(ring)
        if pyroles:
            aroma = aroma + tuple(pyroles)
        return aroma

    @cached_property
    def metal_ligands(self):
        """
        Atoms that could possibly be involved in chelating a metal ion.
        O: R-OH, R-O-R, C(=O)-OH, C(=O)-O-R, R-C(=O)H, R-C(=O)-R C(=O)-NR2
        N: N(Ar, not pyrole), R-NR2
        S: R-SH, R-S-R
        """
        raise NotImplemented

    _pharmacophore_config = {'min_dist': .5, 'hydroph_dist_max': 4., 'hbond_dist_max': 4.1,
                             'hbond_don_angle_min': 1.75}


__all__ = ['Pharmacophore']
