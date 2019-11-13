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
from itertools import product, permutations
from math import hypot, pi, acos, cos, sin, ceil, atan2


def rotate_vector(x1, y1, x2, y2):
    """
    rotate x,y vector over x2-x1, y2-y1 angle
    """
    angle = atan2(y2, x2)
    cos_rad = cos(angle)
    sin_rad = sin(angle)
    return cos_rad * x1 - sin_rad * y1, sin_rad * x1 + cos_rad * y1


def crosses(a, b, c, d):
    ax, ay = a
    bx, by = b
    cx, cy = c
    dx, dy = d

    x, y = None, None
    a1 = ay - by
    b1 = bx - ax
    c1 = ax * by - bx * ay
    a2 = cy - dy
    b2 = dx - cx
    c2 = cx * dy - dx * cy
    um1 = a1 * b1
    if a1 - .001 <= 0 <= a1 + .001:
        a1 = 0
    if b1 - .001 <= 0 <= b1 + .001:
        b1 = 0
    if a2 - .001 <= 0 <= a2 + .001:
        a2 = 0
    if b2 - .001 <= 0 <= b2 + .001:
        b2 = 0
    if not a1 and a2 or not b2 and b1:
        x = -c2 / a2
        y = -c1 / b1
    elif not a2 and a1 or not b1 and b2:
        x = -c1 / a1
        y = -c2 / b2
    elif not a1 and not a2 or not b1 and b2 or um1 - .001 <= a2 * b2 <= um1 + .001:
        return False

    #   calculate intersection point
    if x is None and y is None:
        x = (b2 * c1 - c2 * b1) / (b1 * a2 - b2 * a1)
        if not b1:
            b1 = .001
        y = (-c1 - a1 * x) / b1
    for px, py in (a, b, c, d):
        if x - .001 <= px <= x + .001 and y - .001 <= py <= y + .001:
            return False

    #   check found point in ab and in cd
    vr1 = hypot(ax - x, ay - y)
    vr2 = hypot(bx - x, by - y)
    ab = hypot(bx - ax, by - ay)
    vr3 = hypot(cx - x, cy - y)
    vr4 = hypot(dx - x, dy - y)
    cd = hypot(dx - cx, dy - cy)
    return ab == vr1 + vr2 and cd == vr3 + vr4


def superposition(def_dict):
    out = False
    for k, v in def_dict:
        if out:
            break
        kx, ky = v
        for key, value in def_dict:
            if out:
                break
            if k == key:
                continue
            key_x, key_y = value
            d_x, d_y = abs(key_x - kx), abs(key_y - ky)
            if d_x < .1 and d_y < .1:
                out = True
    return out


class Calculate2D:
    __slots__ = ()

    @staticmethod
    def towards(bond1, bond2):
        return bond1 == 1 and bond2 == 3 or bond1 == 3 and bond2 == 1 or bond1 == 2 and bond2 == 2

    def clean2d(self):
        plane = self._plane
        atoms = self._atoms
        bonds = self._bonds
        atoms_order = self.atoms_order

        for_cycles = defaultdict(set)
        non_standard_cycle = defaultdict(set)
        for cycle in self.sssr:
            lc = len(cycle)
            if lc != 6:
                for_cycles[lc].union(set(cycle))
            for elem in cycle:
                non_standard_cycle[elem].add(lc)

        size = len(atoms)
        dx, dy = n_dist / 2, n_dist * sin(pi / 3)
        components = {k: set(v) for k, v in bonds.items()}
        start, neigh = min((x for x in components.items()), key=lambda x: len(x[1]))
        stack = [[(next_atom, start, bonds[start][next_atom].order, None, 0, 0)] for next_atom in neigh]

        path = []
        hashed_path = set()
        while stack:
            atom, prev_atom, bond1, current_coordinates, v, _ = stack[-1].pop()
            if atom not in hashed_path:
                bond2 = bonds[prev_atom][atom].order
                if current_coordinates is None:
                    current_coordinates = (dx, dy)
                    plane[atom], plane[prev_atom] = current_coordinates, (0, 0)
                    vx, vy = mx, my = dx, dy
                    path.append((start, (0, 0)))
                    path.append((atom, start, current_coordinates))
                    hashed_path.add(start)
                    hashed_path.add(atom)
                else:
                    x, y = v
                    nx, ny = plane[prev_atom]
                    if self.towards(bond1, bond2):
                        current_coordinates = (x + nx, y + ny)
                    mx, my = plane[atom] = current_coordinates
                    vx, vy = mx - nx, my - ny
                    path.append((atom, prev_atom, current_coordinates))
                    hashed_path.add(atom)

            if len(path) == size:
                bad_points = superposition(plane.items())
                cross = False
                long = False

                # for vector ab
                if not bad_points:
                    for a, v in bonds.items():
                        ax, ay = plane[a]
                        for b in v:
                            bx, by = plane[b]
                            if not .800 < hypot(bx - ax, by - ay) < .850:
                                long = True

                if not bad_points and not long:
                    for i, tpl in enumerate(path[1:], 1):
                        if cross:
                            break
                        n, m, crdnts = tpl
                        nx, ny = plane[n]
                        mx, my = plane[m]
                        for _n, _m, _crdnts in path[i+1:]:
                            if cross:
                                break
                            if n != _n and m != _m:
                                _nx, _ny = plane[_n]
                                _mx, _my = plane[_m]
                                cross = crosses((nx, ny), (mx, my), (_nx, _ny), (_mx, _my))

                if not bad_points and not cross and not long:
                    yield path
                del stack[-1]
                if stack:
                    path = path[:stack[-1][-1][-1]]
                    hashed_path = {x for x, *_ in path}

            elif atom != start:  # we finished. next step is final closure
                for_stack = []
                closures = []
                for next_atom in components[atom]:
                    if next_atom == prev_atom:  # only forward. behind us is the homeland
                        continue
                    elif next_atom in hashed_path:  # closure found
                        closures.append(next_atom)
                    else:
                        for_stack.append(next_atom)

                if len(for_stack) == 1:  # easy path grow. next bond double or include single for pyroles
                    lp = len(path)
                    next_atom = for_stack[0]
                    opposite = stack[-1].copy()
                    vx1, vy1 = rotate_vector(dx, dy, vx, vy)
                    vx2, vy2 = rotate_vector(dx, -dy, vx, vy)
                    new_coords1 = (vx1 + mx, vy1 + my)
                    new_coords2 = (vx2 + mx, vy2 + my)
                    new_coords = False
                    for hash_atom in hashed_path:
                        v_x, v_y = plane[hash_atom]
                        n1_x, n1_y = new_coords1
                        n2_x, n2_y = new_coords2
                        d1_x, d2_x = abs(n1_x - v_x), abs(n2_x - v_x)
                        d1_y, d2_y = abs(n1_y - v_y), abs(n2_y - v_y)
                        if d1_x < .3 and d1_y < .3:
                            new_coords = new_coords2
                            if d2_x < .3 and d2_y < .3:
                                new_coords = None
                        elif d2_x < .3 and d2_y < .3:
                            new_coords = new_coords1
                            if d1_x < .3 and d1_y < .3:
                                new_coords = None
                    if new_coords is None:
                        del stack[-1]
                        if stack:
                            path = path[:stack[-1][-1][-1]]
                            hashed_path = {x for x, *_ in path}
                    elif new_coords:
                        stack[-1].append((next_atom, atom, bond2, new_coords, (vx, vy), lp))
                    else:
                        stack[-1].append((next_atom, atom, bond2, new_coords1, (vx, vy), lp))
                        if not self.towards(bond2, bonds[atom][next_atom].order):
                            opposite.append((next_atom, atom, bond2, new_coords2, (vx, vy), lp))
                            stack.append(opposite)

                elif len(for_stack) == 2:  # fork
                    lp = len(path)
                    next_atom1, next_atom2 = for_stack
                    vx1, vy1 = rotate_vector(dx, dy, vx, vy)
                    vx2, vy2 = rotate_vector(dx, -dy, vx, vy)
                    new_coords1 = (vx1 + mx, vy1 + my)
                    new_coords2 = (vx2 + mx, vy2 + my)

                    new_coords = False
                    for hash_atom in hashed_path:
                        v_x, v_y = plane[hash_atom]
                        n1_x, n1_y = new_coords1
                        n2_x, n2_y = new_coords2
                        d1_x, d2_x = abs(n1_x - v_x), abs(n2_x - v_x)
                        d1_y, d2_y = abs(n1_y - v_y), abs(n2_y - v_y)
                        if d1_x < .1 and d1_y < .1 or d2_x < .1 and d2_y < .1:
                            new_coords = None

                    if new_coords is None:
                        del stack[-1]
                        if stack:
                            path = path[:stack[-1][-1][-1]]
                            hashed_path = {x for x, *_ in path}
                    else:
                        opposite = stack[-1].copy()
                        stack[-1].append((next_atom1, atom, bond2, new_coords1, (vx1, vy1), None))
                        stack[-1].append((next_atom2, atom, bond2, new_coords2, (vx2, vy2), lp))
                        if not self.towards(bond2, bonds[atom][next_atom1].order) and \
                                not self.towards(bond2, bonds[atom][next_atom2].order) and \
                                atoms_order[next_atom1] != atoms_order[next_atom2]:
                            opposite.append((next_atom1, atom, bond2, new_coords2, (vx2, vy2), None))
                            opposite.append((next_atom2, atom, bond2, new_coords1, (vx1, vy1), lp))
                            stack.append(opposite)

                elif len(for_stack) == 3:
                    lp = len(path)
                    dx1, dx2, dx3 = 0, .825, 0
                    dy1, dy2, dy3 = .825, 0, -.825
                    combinations = permutations(for_stack)
                    vx1, vy1 = rotate_vector(dx1, dy1, vx, vy)
                    vx2, vy2 = rotate_vector(dx2, dy2, vx, vy)
                    vx3, vy3 = rotate_vector(dx3, dy3, vx, vy)
                    new_coords1 = (vx1 + mx, vy1 + my)
                    new_coords2 = (vx2 + mx, vy2 + my)
                    new_coords3 = (vx3 + mx, vy3 + my)

                    new_coords = False
                    for hash_atom in hashed_path:
                        v_x, v_y = plane[hash_atom]
                        n1_x, n1_y = new_coords1
                        n2_x, n2_y = new_coords2
                        n3_x, n3_y = new_coords3
                        d1_x, d2_x, d3_x = abs(n1_x - v_x), abs(n2_x - v_x), abs(n3_x - v_x)
                        d1_y, d2_y, d3_y = abs(n1_y - v_y), abs(n2_y - v_y), abs(n3_y - v_y)

                        if d1_x < .1 and d1_y < .1 or d2_x < .1 and d2_y < .1 or d3_x < .1 and d3_y < .1:
                            new_coords = None

                    if new_coords is None:
                        del stack[-1]
                        if stack:
                            path = path[:stack[-1][-1][-1]]
                            hashed_path = {x for x, *_ in path}
                    else:
                        opposite = stack[-1].copy()
                        opposite1 = stack[-1].copy()
                        next_atom1, next_atom2, next_atom3 = next(combinations)
                        stack[-1].append((next_atom1, atom, bond2, (vx1 + mx, vy1 + my), (vx1, vy1), None))
                        stack[-1].append((next_atom2, atom, bond2, (vx2 + mx, vy2 + my), (vx2, vy2), None))
                        stack[-1].append((next_atom3, atom, bond2, (vx3 + mx, vy3 + my), (vx3, vy3), lp))

                        for nxt_atm1, nxt_atm2, nxt_atm3 in combinations:
                            if atoms_order[next_atom1] == atoms_order[nxt_atm1] \
                                    or atoms_order[next_atom1] == atoms_order[nxt_atm2] \
                                    or atoms_order[next_atom1] == atoms_order[nxt_atm3]:
                                continue
                            opposite.append((nxt_atm1, atom, bond2, (vx1 + mx, vy1 + my), (vx1, vy1), None))
                            opposite.append((nxt_atm2, atom, bond2, (vx2 + mx, vy2 + my), (vx2, vy2), None))
                            opposite.append((nxt_atm3, atom, bond2, (vx3 + mx, vy3 + my), (vx3, vy3), lp))
                            stack.append(opposite)
                            opposite = opposite1.copy()

                elif len(for_stack) == 4:
                    lp = len(path)
                    angle = pi / 2.5
                    dx1, dx2 = -.825 * cos(angle), .825 * cos(angle/2)
                    dy1, dy2 = .825 * sin(angle), .825 * sin(angle/2)
                    combinations = permutations(for_stack)
                    vx1, vy1 = rotate_vector(dx1, dy1, vx, vy)
                    vx2, vy2 = rotate_vector(dx2, dy2, vx, vy)
                    vx3, vy3 = rotate_vector(dx2, -dy2, vx, vy)
                    vx4, vy4 = rotate_vector(dx1, -dy1, vx, vy)
                    new_coords1 = (vx1 + mx, vy1 + my)
                    new_coords2 = (vx2 + mx, vy2 + my)
                    new_coords3 = (vx3 + mx, vy3 + my)
                    new_coords4 = (vx4 + mx, vy4 + my)

                    new_coords = False
                    for hash_atom in hashed_path:
                        v_x, v_y = plane[hash_atom]
                        n1_x, n1_y = new_coords1
                        n2_x, n2_y = new_coords2
                        n3_x, n3_y = new_coords3
                        n4_x, n4_y = new_coords4
                        d1_x, d2_x, d3_x, d4_x = abs(n1_x - v_x), abs(n2_x - v_x), abs(n3_x - v_x), abs(n4_x - v_x)
                        d1_y, d2_y, d3_y, d4_y = abs(n1_y - v_y), abs(n2_y - v_y), abs(n3_y - v_y), abs(n4_y - v_y)

                        if d1_x < .1 and d1_y < .1 or d2_x < .1 and d2_y < .1 or d3_x < .1 and d3_y < .1 \
                                or d4_x < .1 and d4_y < .1:
                            new_coords = None

                    if new_coords is None:
                        del stack[-1]
                        if stack:
                            path = path[:stack[-1][-1][-1]]
                            hashed_path = {x for x, *_ in path}
                    else:
                        opposite = stack[-1].copy()
                        opposite1 = stack[-1].copy()
                        next_atom1, next_atom2, next_atom3, next_atom4 = next(combinations)
                        stack[-1].append((next_atom1, atom, bond2, (vx1 + mx, vy1 + my), (vx1, vy1), None))
                        stack[-1].append((next_atom2, atom, bond2, (vx2 + mx, vy2 + my), (vx2, vy2), None))
                        stack[-1].append((next_atom3, atom, bond2, (vx3 + mx, vy3 + my), (vx3, vy3), None))
                        stack[-1].append((next_atom4, atom, bond2, (vx4 + mx, vy4 + my), (vx4, vy4), lp))

                        for nxt_atm1, nxt_atm2, nxt_atm3, nxt_atm4 in combinations:
                            if atoms_order[next_atom1] == atoms_order[nxt_atm1] \
                                    or atoms_order[next_atom1] == atoms_order[nxt_atm2] \
                                    or atoms_order[next_atom1] == atoms_order[nxt_atm3] \
                                    or atoms_order[next_atom1] == atoms_order[nxt_atm4]:
                                continue
                            opposite.append((nxt_atm1, atom, bond2, (vx1 + mx, vy1 + my), (vx1, vy1), None))
                            opposite.append((nxt_atm2, atom, bond2, (vx2 + mx, vy2 + my), (vx2, vy2), None))
                            opposite.append((nxt_atm3, atom, bond2, (vx3 + mx, vy3 + my), (vx3, vy3), None))
                            opposite.append((nxt_atm4, atom, bond2, (vx4 + mx, vy4 + my), (vx4, vy4), lp))
                            stack.append(opposite)
                            opposite = opposite1.copy()

                elif len(for_stack) == 5:
                    lp = len(path)
                    combinations = permutations(for_stack)
                    vx1, vy1 = rotate_vector(-dx, dy, vx, vy)
                    vx2, vy2 = rotate_vector(dx, dy, vx, vy)
                    vx3, vy3 = rotate_vector(.825, 0, vx, vy)
                    vx4, vy4 = rotate_vector(dx, -dy, vx, vy)
                    vx5, vy5 = rotate_vector(-dx, -dy, vx, vy)
                    new_coords1 = (vx1 + mx, vy1 + my)
                    new_coords2 = (vx2 + mx, vy2 + my)
                    new_coords3 = (vx3 + mx, vy3 + my)
                    new_coords4 = (vx4 + mx, vy4 + my)
                    new_coords5 = (vx5 + mx, vy5 + my)

                    new_coords = False
                    for hash_atom in hashed_path:
                        v_x, v_y = plane[hash_atom]
                        n1_x, n1_y = new_coords1
                        n2_x, n2_y = new_coords2
                        n3_x, n3_y = new_coords3
                        n4_x, n4_y = new_coords4
                        n5_x, n5_y = new_coords5
                        d1_x, d2_x, d3_x, d4_x, d5_x = abs(n1_x - v_x), abs(n2_x - v_x), abs(n3_x - v_x), abs(n4_x - v_x), abs(n5_x - v_x)
                        d1_y, d2_y, d3_y, d4_y, d5_y = abs(n1_y - v_y), abs(n2_y - v_y), abs(n3_y - v_y), abs(n4_y - v_y), abs(n5_y - v_y)

                        if d1_x < .1 and d1_y < .1 or d2_x < .1 and d2_y < .1 or d3_x < .1 and d3_y < .1 \
                                or d4_x < .1 and d4_y < .1 or d5_x < .1 and d5_y < .1:
                            new_coords = None

                    if new_coords is None:
                        del stack[-1]
                        if stack:
                            path = path[:stack[-1][-1][-1]]
                            hashed_path = {x for x, *_ in path}
                    else:
                        opposite = stack[-1].copy()
                        opposite1 = stack[-1].copy()
                        next_atom1, next_atom2, next_atom3, next_atom4, next_atom5 = next(combinations)
                        stack[-1].append((next_atom1, atom, bond2, (vx1 + mx, vy1 + my), (vx1, vy1), None))
                        stack[-1].append((next_atom2, atom, bond2, (vx2 + mx, vy2 + my), (vx2, vy2), None))
                        stack[-1].append((next_atom3, atom, bond2, (vx3 + mx, vy3 + my), (vx3, vy3), None))
                        stack[-1].append((next_atom4, atom, bond2, (vx4 + mx, vy4 + my), (vx4, vy4), None))
                        stack[-1].append((next_atom5, atom, bond2, (vx5 + mx, vy5 + my), (vx5, vy5), lp))

                        for nxt_atm1, nxt_atm2, nxt_atm3, nxt_atm4, nxt_atm5 in combinations:
                            if atoms_order[next_atom1] == atoms_order[nxt_atm1] \
                                    or atoms_order[next_atom1] == atoms_order[nxt_atm2] \
                                    or atoms_order[next_atom1] == atoms_order[nxt_atm3] \
                                    or atoms_order[next_atom1] == atoms_order[nxt_atm4]\
                                    or atoms_order[next_atom1] == atoms_order[nxt_atm5]:
                                continue
                            opposite.append((nxt_atm1, atom, bond2, (vx1 + mx, vy1 + my), (vx1, vy1), None))
                            opposite.append((nxt_atm2, atom, bond2, (vx2 + mx, vy2 + my), (vx2, vy2), None))
                            opposite.append((nxt_atm3, atom, bond2, (vx3 + mx, vy3 + my), (vx3, vy3), None))
                            opposite.append((nxt_atm4, atom, bond2, (vx4 + mx, vy4 + my), (vx4, vy4), None))
                            opposite.append((nxt_atm5, atom, bond2, (vx5 + mx, vy5 + my), (vx5, vy5), lp))
                            stack.append(opposite)
                            opposite = opposite1.copy()

                elif closures:
                    atom_x, atom_y = plane[atom]
                    for closure in closures:
                        clos_x, clos_y = plane[closure]
                        vect_x, vect_y = clos_x - atom_x, clos_y - atom_y
                        if not .800 < hypot(vect_x, vect_y) < .850:
                            del stack[-1]
                            if stack:
                                path = path[:stack[-1][-1][-1]]
                                hashed_path = {x for x, *_ in path}
                            break


n_dist = .825

__all__ = ['Calculate2D']
