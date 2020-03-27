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
from math import sqrt, acos, pi
from typing import Tuple, Optional, List, NamedTuple, Union

pi2 = pi / 2


class Hydrophobic(NamedTuple):
    source: int
    target: int
    distance: float


class Salts(NamedTuple):
    source: int
    target: int
    distance: float


class HydrogenDonor(NamedTuple):
    source: int
    target: int
    distance: float  # between atoms connected by hydrogen
    acceptor_distance: float  # acceptor-H bond length
    donor_distance: float  # donor-H bond length
    angle: float


class HydrogenAcceptor(NamedTuple):
    source: int
    target: int
    distance: float  # between atoms connected by hydrogen
    acceptor_distance: float  # acceptor-H bond length
    donor_distance: float  # donor-H bond length
    angle: float


class PiStack(NamedTuple):
    source: Tuple[int, ...]
    target: Tuple[int, ...]
    distance: float  # distance between centers
    angle: float  # angle between plains
    t_shaped: bool
    offset: float  # distance of projected centers


class PiCation(NamedTuple):
    source: Tuple[int, ...]
    target: int
    distance: float  # distance between centers
    offset: float  # distance of projected centers


class CationPi(NamedTuple):
    source: int
    target: Tuple[int, ...]
    distance: float  # distance between centers
    offset: float  # distance of projected centers


class HalogenAcceptor(NamedTuple):
    source: int
    target: int
    distance: float
    acceptor_angle: float  # acceptor group to halogen atom angle
    donor_angle: float  # halogen donor group to halogen acceptor atom angle


class HalogenDonor(NamedTuple):
    source: int
    target: int
    distance: float
    acceptor_angle: float  # acceptor group to halogen atom angle
    donor_angle: float  # halogen donor group to halogen acceptor atom angle


def distance(n, m):
    nx, ny, nz = n
    mx, my, mz = m
    return sqrt((nx - mx) ** 2 + (ny - my) ** 2 + (nz - mz) ** 2)


def points_angle(n, m, k):
    # m<---n
    #       \
    #        v
    #        k

    nx, ny, nz = n
    mx, my, mz = m
    kx, ky, kz = k
    return vectors_angle((mx - nx, my - ny, mz - nz), (kx - nx, ky - ny, kz - nz))


def vectors_angle(n, m):
    # n<---0
    #       \
    #        v
    #        m
    nmx, nmy, nmz = n
    nkx, nky, nkz = m

    nmd = sqrt(nmx ** 2 + nmy ** 2 + nmz ** 2)
    nkd = sqrt(nkx ** 2 + nky ** 2 + nkz ** 2)

    # normalization
    nmx /= nmd
    nmy /= nmd
    nmz /= nmd
    nkx /= nkd
    nky /= nkd
    nkz /= nkd
    return acos(nmx * nkx + nmy * nky + nmz * nkz)


def ring_math(ring, xyz):
    # ring center
    cx, cy, cz = xyz[ring[0]]
    for n in ring[1:]:
        x, y, z = xyz[n]
        cx += x
        cy += y
        cz += z
    lr = len(ring)
    cx /= lr
    cy /= lr
    cz /= lr

    # ring vector
    r1x, r1y, r1z = x - cx, y - cy, z - cz
    x, y, z = xyz[ring[0]]
    r2x, r2y, r2z = x - cx, y - cy, z - cz
    nx = r1y * r2z - r1z * r2y
    ny = r1z * r2x - r1x * r2z
    nz = r1x * r2y - r1y * r2x
    return (cx, cy, cz), (nx, ny, nz)


def projected_distance(normal, centroid, point):
    nx, ny, nz = normal

    cx, cy, cz = centroid
    px, py, pz = point

    cpx = px - cx
    cpy = py - cy
    cpz = pz - cz

    sn = cpx * nx + cpy * ny + cpz * nz
    sd = nx ** 2 + ny ** 2 + nz ** 2
    sb = sn / sd
    x, y, z = cpx - nx * sb, cpy - ny * sb, cpz - nz * sb
    return sqrt(x ** 2 + y ** 2 + z ** 2)


class Pharmacophore:
    __slots__ = ()

    def find_contacts(self, other: 'Pharmacophore') -> List[Union[Hydrophobic, Salts, HydrogenDonor, HydrogenAcceptor,
                                                                  PiStack, PiCation, CationPi, HalogenAcceptor,
                                                                  HalogenDonor]]:
        """
        Find pharmocophore contacts in 3d. Required explicit hydrogens.
        """
        contacts = []
        contacts.extend(self.__electrostatic_interactions(other))
        contacts.extend(self.__hydrogen_bond_interactions(other))
        contacts.extend(self.__pi_interactions(other))
        contacts.extend(self.__halogen_bond_interactions(other))
        return contacts

    def __electrostatic_interactions(self, other):
        """
        Detection of electrostatic contacts.
        """
        config = self._pharmacophore_config
        min_dist = config['min_dist']
        hydroph_dist_max = config['hydroph_dist_max']
        salt_bridge_dist_max = config['salt_bridge_dist_max']

        contacts = []
        s_xyz = self._conformers[0]
        o_xyz = other._conformers[0]
        for n, m in product(self.hydrophobic_centers, other.hydrophobic_centers):
            d = distance(s_xyz[n], o_xyz[m])
            if min_dist < d < hydroph_dist_max:
                contacts.append(Hydrophobic(n, m, d))
        for n, m in product(self.positive_charged_centers, other.negative_charged_centers):
            d = distance(s_xyz[n], o_xyz[m])
            if min_dist < d < salt_bridge_dist_max:
                contacts.append(Salts(n, m, d))
        for n, m in product(self.negative_charged_centers, other.positive_charged_centers):
            d = distance(s_xyz[n], o_xyz[m])
            if min_dist < d < salt_bridge_dist_max:
                contacts.append(Salts(n, m, d))
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
        for (n, h), m in product(self.hydrogen_donor_centers, other.hydrogen_acceptor_centers):
            dist_nm = distance(s_xyz[n], o_xyz[m])
            if min_dist < dist_nm < hbond_dist_max:
                a = points_angle(s_xyz[h], s_xyz[n], o_xyz[m])
                if a >= hbond_don_angle_min:
                    seen.add((n, m))
                    contacts.append(HydrogenDonor(n, m, dist_nm, distance(s_xyz[h], o_xyz[m]),
                                                  distance(s_xyz[n], s_xyz[h]), a))

        for n, (m, h) in product(self.hydrogen_acceptor_centers, other.hydrogen_donor_centers):
            if (n, m) in seen:
                continue
            dist_nm = distance(s_xyz[n], o_xyz[m])
            if min_dist < dist_nm < hbond_dist_max:
                a = points_angle(o_xyz[h], o_xyz[m], s_xyz[n])
                if a > hbond_don_angle_min:
                    contacts.append(HydrogenAcceptor(n, m, dist_nm, distance(o_xyz[h], s_xyz[n]),
                                                     distance(o_xyz[m], o_xyz[h]), a))
        return contacts

    def __pi_interactions(self, other):
        """
        Detection of pi-stackings between aromatic ring systems and pi-Cation interaction between aromatic rings and
        positively charged groups.

        For tertiary and quaternary amines, check also the angle between the ring and the nitrogen.
        """
        config = self._pharmacophore_config
        min_dist = config['min_dist']
        pi_stack_dist_max = config['pi_stack_dist_max']
        pi_stack_ang_dev = config['pi_stack_ang_dev']
        pi_stack_offset_max = config['pi_stack_offset_max']
        pi_cation_dist_max = config['pi_cation_dist_max']
        pi_stack_ang_dev_perp = pi2 - pi_stack_ang_dev

        contacts = []
        s_xyz = self._conformers[0]
        o_xyz = other._conformers[0]

        s_centers = []
        s_normals = []
        for r in self.aromatic_centers:
            c, n = ring_math(r, s_xyz)
            s_centers.append(c)
            s_normals.append(n)
        o_centers = []
        o_normals = []
        for r in other.aromatic_centers:
            c, n = ring_math(r, o_xyz)
            o_centers.append(c)
            o_normals.append(n)

        for (nr, nc, nn), (mr, mc, mn) in product(zip(self.aromatic_centers, s_centers, s_normals),
                                                  zip(other.aromatic_centers, o_centers, o_normals)):
            d = distance(nc, mc)
            if min_dist < d < pi_stack_dist_max:
                a = vectors_angle(nn, mn)
                if a > pi2:
                    a = pi - a
                if a < pi_stack_ang_dev:
                    t_shaped = False
                elif a > pi_stack_ang_dev_perp:
                    t_shaped = True
                else:
                    continue

                offset = min(projected_distance(nn, nc, mc), projected_distance(mn, mc, nc))
                if offset < pi_stack_offset_max:
                    contacts.append(PiStack(nr, mr, d, a, t_shaped, offset))

        for (nr, nc, nn), m in product(zip(self.aromatic_centers, s_centers, s_normals),
                                       other.positive_charged_centers):
            d = distance(nc, o_xyz[m])
            if min_dist < d < pi_cation_dist_max:
                offset = projected_distance(nn, nc, o_xyz[m])
                if offset < pi_stack_offset_max:
                    contacts.append(PiCation(nr, m, d, offset))
        for (mr, mc, mn), n in product(zip(other.aromatic_centers, o_centers, o_normals),
                                       self.positive_charged_centers):
            d = distance(mc, s_xyz[n])
            if min_dist < d < pi_cation_dist_max:
                offset = projected_distance(mn, mc, s_xyz[n])
                if offset < pi_stack_offset_max:
                    contacts.append(CationPi(n, mr, d, offset))
        return contacts

    def __halogen_bond_interactions(self, other):
        """
        Detect all halogen bonds of the type Y-O...X-C
        """
        config = self._pharmacophore_config
        min_dist = config['min_dist']
        halogen_dist_max = config['halogen_dist_max']
        halogen_angle_dev = config['halogen_angle_dev']
        halogen_acc_angle = config['halogen_acc_angle']
        halogen_don_angle = config['halogen_don_angle']
        acc_ap = halogen_acc_angle + halogen_angle_dev
        acc_am = halogen_acc_angle - halogen_angle_dev
        don_ap = halogen_don_angle + halogen_angle_dev
        don_am = halogen_don_angle - halogen_angle_dev

        contacts = []
        s_xyz = self._conformers[0]
        o_xyz = other._conformers[0]
        for (an, a), (hn, h) in product(self.halogen_acceptor_centers, other.halogen_donor_centers):
            d = distance(s_xyz[a], o_xyz[h])
            if min_dist < d < halogen_dist_max:
                acc_angle = points_angle(s_xyz[a], s_xyz[an], o_xyz[h])
                if acc_am < acc_angle < acc_ap:
                    don_angle = points_angle(o_xyz[h], s_xyz[a], o_xyz[hn])
                    if don_am < don_angle < don_ap:
                        contacts.append(HalogenAcceptor(a, h, d, acc_angle, don_angle))
        for (an, a), (hn, h) in product(other.halogen_acceptor_centers, self.halogen_donor_centers):
            d = distance(s_xyz[a], o_xyz[h])
            if min_dist < d < halogen_dist_max:
                acc_angle = points_angle(s_xyz[a], s_xyz[an], o_xyz[h])
                if acc_am < acc_angle < acc_ap:
                    don_angle = points_angle(o_xyz[h], s_xyz[a], o_xyz[hn])
                    if don_am < don_angle < don_ap:
                        contacts.append(HalogenDonor(h, a, d, acc_angle, don_angle))
        return contacts

    def __metal_complex_interactions(self, other):
        """Find all metal complexes between metals and appropriate groups in both protein and ligand, as well as water"""
        data = namedtuple('metal_complex', 'metal metal_orig_idx metal_type target target_orig_idx target_type '
                                           'coordination_num distance resnr restype '
                                           'reschain  restype_l reschain_l resnr_l location rms, geometry num_partners complexnum')
        pairings_dict = {}
        pairings = []
        # #@todo Refactor
        metal_to_id = {}
        metal_to_orig_atom = {}
        for metal, target in itertools.product(metals, metal_binding_lig + metal_binding_bs):
            distance = euclidean3d(metal.m.coords, target.atom.coords)
            if not distance < config.METAL_DIST_MAX:
                continue
            if metal.m not in pairings_dict:
                pairings_dict[metal.m] = [(target, distance), ]
                metal_to_id[metal.m] = metal.m_orig_idx
                metal_to_orig_atom[metal.m] = metal.orig_m
            else:
                pairings_dict[metal.m].append((target, distance))
        for cnum, metal in enumerate(pairings_dict):
            rms = 0.0
            excluded = []
            # cnum +1 being the complex number
            contact_pairs = pairings_dict[metal]
            num_targets = len(contact_pairs)
            vectors_dict = defaultdict(list)
            for contact_pair in contact_pairs:
                target, distance = contact_pair
                vectors_dict[target.atom.idx].append(vector(metal.coords, target.atom.coords))

            # Listing of coordination numbers and their geometries
            configs = {2: ['linear', ],
                       3: ['trigonal.planar', 'trigonal.pyramidal'],
                       4: ['tetrahedral', 'square.planar'],
                       5: ['trigonal.bipyramidal', 'square.pyramidal'],
                       6: ['octahedral', ]}

            # Angle signatures for each geometry (as seen from each target atom)
            ideal_angles = {'linear': [[180.0]] * 2,
                            'trigonal.planar': [[120.0, 120.0]] * 3,
                            'trigonal.pyramidal': [[109.5, 109.5]] * 3,
                            'tetrahedral': [[109.5, 109.5, 109.5, 109.5]] * 4,
                            'square.planar': [[90.0, 90.0, 90.0, 90.0]] * 4,
                            'trigonal.bipyramidal': [[120.0, 120.0, 90.0, 90.0]] * 3 + [[90.0, 90.0, 90.0, 180.0]] * 2,
                            'square.pyramidal': [[90.0, 90.0, 90.0, 180.0]] * 4 + [[90.0, 90.0, 90.0, 90.0]],
                            'octahedral': [[90.0, 90.0, 90.0, 90.0, 180.0]] * 6}
            angles_dict = {}

            for target in vectors_dict:
                cur_vector = vectors_dict[target]
                other_vectors = []
                for t in vectors_dict:
                    if not t == target:
                        [other_vectors.append(x) for x in vectors_dict[t]]
                angles = [vecangle(pair[0], pair[1]) for pair in itertools.product(cur_vector, other_vectors)]
                angles_dict[target] = angles

            all_total = []  # Record fit information for each geometry tested
            gdata = namedtuple('gdata', 'geometry rms coordination excluded diff_targets')  # Geometry Data
            # Can't specify geometry with only one target
            if num_targets == 1:
                final_geom = 'NA'
                final_coo = 1
                excluded = []
                rms = 0.0
            else:
                for coo in sorted(configs, reverse=True):  # Start with highest coordination number
                    geometries = configs[coo]
                    for geometry in geometries:
                        signature = ideal_angles[geometry]  # Set of ideal angles for geometry, from each perspective
                        geometry_total = 0
                        geometry_scores = []  # All scores for one geometry (from all subsignatures)
                        used_up_targets = []  # Use each target just once for a subsignature
                        not_used = []
                        coo_diff = num_targets - coo  # How many more observed targets are there?

                        # Find best match for each subsignature
                        for subsignature in signature:  # Ideal angles from one perspective
                            best_target = None  # There's one best-matching target for each subsignature
                            best_target_score = 999

                            for k, target in enumerate(angles_dict):
                                if target not in used_up_targets:
                                    observed_angles = angles_dict[
                                        target]  # Observed angles from perspective of one target
                                    single_target_scores = []
                                    used_up_observed_angles = []
                                    for i, ideal_angle in enumerate(subsignature):
                                        # For each angle in the signature, find the best-matching observed angle
                                        best_match = None
                                        best_match_diff = 999
                                        for j, observed_angle in enumerate(observed_angles):
                                            if j not in used_up_observed_angles:
                                                diff = abs(ideal_angle - observed_angle)
                                                if diff < best_match_diff:
                                                    best_match_diff = diff
                                                    best_match = j
                                        if best_match is not None:
                                            used_up_observed_angles.append(best_match)
                                            single_target_scores.append(best_match_diff)
                                    # Calculate RMS for target angles
                                    target_total = sum(
                                        [x ** 2 for x in single_target_scores]) ** 0.5  # Tot. score targ/sig
                                    if target_total < best_target_score:
                                        best_target_score = target_total
                                        best_target = target

                            used_up_targets.append(best_target)
                            geometry_scores.append(best_target_score)
                            # Total score is mean of RMS values
                            geometry_total = np.mean(geometry_scores)
                        # Record the targets not used for excluding them when deciding for a final geometry
                        [not_used.append(target) for target in angles_dict if target not in used_up_targets]
                        all_total.append(gdata(geometry=geometry, rms=geometry_total, coordination=coo,
                                               excluded=not_used, diff_targets=coo_diff))

            # Make a decision here. Starting with the geometry with lowest difference in ideal and observed partners ...
            # Check if the difference between the RMS to the next best solution is not larger than 0.5
            if not num_targets == 1:  # Can't decide for any geoemtry in that case
                all_total = sorted(all_total, key=lambda x: abs(x.diff_targets))
                for i, total in enumerate(all_total):
                    next_total = all_total[i + 1]
                    this_rms, next_rms = total.rms, next_total.rms
                    diff_to_next = next_rms - this_rms
                    if diff_to_next > 0.5:
                        final_geom, final_coo, rms, excluded = total.geometry, total.coordination, total.rms, total.excluded
                        break
                    elif next_total.rms < 3.5:
                        final_geom, final_coo, = next_total.geometry, next_total.coordination
                        rms, excluded = next_total.rms, next_total.excluded
                        break
                    elif i == len(all_total) - 2:
                        final_geom, final_coo, rms, excluded = "NA", "NA", float('nan'), []
                        break

            # Record all contact pairing, excluding those with targets superfluous for chosen geometry
            only_water = set([x[0].location for x in contact_pairs]) == {'water'}
            if not only_water:  # No complex if just with water as targets
                write_message("Metal ion %s complexed with %s geometry (coo. number %r/ %i observed).\n"
                              % (metal.type, final_geom, final_coo, num_targets), indent=True)
                for contact_pair in contact_pairs:
                    target, distance = contact_pair
                    if target.atom.idx not in excluded:
                        metal_orig_atom = metal_to_orig_atom[metal]
                        restype_l, reschain_l, resnr_l = whichrestype(metal_orig_atom), whichchain(
                            metal_orig_atom), whichresnumber(metal_orig_atom)
                        contact = data(metal=metal, metal_orig_idx=metal_to_id[metal], metal_type=metal.type,
                                       target=target, target_orig_idx=target.atom_orig_idx, target_type=target.type,
                                       coordination_num=final_coo, distance=distance, resnr=target.resnr,
                                       restype=target.restype, reschain=target.reschain, location=target.location,
                                       rms=rms, geometry=final_geom, num_partners=num_targets, complexnum=cnum + 1,
                                       resnr_l=resnr_l, restype_l=restype_l, reschain_l=reschain_l)
                        pairings.append(contact)
        return filter_contacts(pairings)

    @cached_property
    def hydrogen_acceptor_centers(self) -> Tuple[int, ...]:
        """
        Atom that could form hydrogen bonds.

        N: R-N(-[H,R,Ar])-[H,R], Ar:N:Ar, Ar-NH2, [R,H]-N=R
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
    def hydrogen_donor_centers(self) -> Tuple[Tuple[int, Optional[int]], ...]:
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
    def halogen_donor_centers(self) -> Tuple[Tuple[int, int], ...]:
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
    def halogen_acceptor_centers(self) -> Tuple[Tuple[int, int], ...]:
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
            if not neighbors[n]:  # skip [NH4+], etc
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
    def metal_ligands_centers(self):
        """
        Atoms that could possibly be involved in chelating a metal ion.
        O: R-OH, R-O-R, C(=O)-OH, C(=O)-O-R, R-C(=O)H, R-C(=O)-R C(=O)-NR2
        N: N(Ar, not pyrole), R-NR2
        S: R-SH, R-S-R
        """
        self.hydrogen_acceptor_centers
        raise NotImplemented

    @cached_property
    def metal_centers(self):
        METAL_IONS = ['Ca', 'Co', 'Mg', 'Mn', 'FE', 'CU', 'ZN', 'FE2', 'FE3', 'FE4', 'LI', 'NA', 'K', 'RB', 'SR', 'CS',
                      'BA',
                      'CR', 'NI', 'FE1', 'NI', 'RU', 'RU1', 'RH', 'RH1', 'PD', 'AG', 'CD', 'LA', 'W', 'W1', 'OS', 'IR',
                      'PT',
                      'PT1', 'AU', 'HG', 'CE', 'PR', 'SM', 'EU', 'GD', 'TB', 'YB', 'LU', 'AL', 'GA', 'IN', 'SB', 'TL',
                      'PB']



    _pharmacophore_config = {'min_dist': .5, 'hydroph_dist_max': 4., 'hbond_dist_max': 4.1,
                             'hbond_don_angle_min': 1.75, 'pi_stack_dist_max': 5.5, 'pi_stack_ang_dev': 0.52,
                             'pi_stack_offset_max': 2.0, 'salt_bridge_dist_max': 5.5, 'halogen_dist_max': 4.0,
                             'halogen_angle_dev': 0.52, 'halogen_acc_angle': 2.1, 'halogen_don_angle': 2.88,
                             'pi_cation_dist_max': 6.}


__all__ = ['Pharmacophore']
