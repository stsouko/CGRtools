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
# from ..algorithms.depict import rotate_vector


def rotate_vector(x1, y1, x2, y2):
    """
    rotate x,y vector over x2-x1, y2-y1 angle
    """
    angle = atan2(y2, x2)
    cos_rad = cos(angle)
    sin_rad = sin(angle)
    return cos_rad * x1 - sin_rad * y1, sin_rad * x1 + cos_rad * y1


def rotate_vector2(x1, y1, angle):
    """
    rotate x,y vector over angle
    """
    cos_rad = cos(angle)
    sin_rad = sin(angle)
    return cos_rad * x1 + sin_rad * y1, -sin_rad * x1 + cos_rad * y1


def dfs_paths(graph, start, goal):
    stack = [(start, [start])]
    while stack:
        (vertex, path) = stack.pop()
        for next in set(graph[vertex]) - set(path):
            if next == goal:
                yield path + [next]
            else:
                stack.append((next, path + [next]))


class Calculate2D:
    __slots__ = ()

    def calculate_cycles(self):
        cycles = self.sssr
        plane = self._plane
        cycle = cycles[0]
        lc = len(cycle)
        positive, negative = 1, -1

        if lc == 3:
            angle = pi / 3
            a, b, c = cycle
            x, y = .825 / 2, sin(angle) * .825
            positive, negative = (x, y), (x, -y)
            direction = negative
            plane[a], plane[b], plane[c] = (0, 0), (.825, 0), direction

        elif lc in (4, 5, 6, 7, 8):
            angle = 2 * pi / lc
            a, b = cycle[0], cycle[1]
            plane[a], plane[b] = (0, 0), (.825, 0)
            seen = {a, b}
            direction = positive * angle

            stack = [(2, (.825, 0), lc - 2)]
            while stack:
                index, coords, count = stack.pop(0)
                atom = cycle[index]
                count -= 1
                x1, y1 = coords
                if atom not in seen:
                    x2, y2 = rotate_vector2(x1, y1, direction)
                    plane[atom] = (x2 + .825, y2)
                    if count:
                        stack.append((index + 1, (x2 + .825, y2), count))

        else:
            a = tuple(range(9, 300, 4))
            b = tuple(range(10, 300, 4))
            c = tuple(range(11, 300, 4))
            d = tuple(range(12, 301, 4))
            # angle = 2 * pi / 3
            dy, dx = .825 / 2, sin(pi / 3) * .825
            n, m, o = (dx, -dy / 2), (0, dy / 2), (-dx, -dy / 2)
            if lc in b or lc in c:
                n, m, o = (dx, dy / 2), (0, -dy / 2), (-dx, dy / 2)
            plane[cycle[0]], plane[cycle[1]], plane[cycle[2]] = n, m, o
            mid = ceil(lc / 2)
            stack = [(3, o, lc - 3)]
            while stack:
                index, coords, count = stack.pop(0)
                x1, y1 = coords
                count -= 1
                atom = cycle[index]

                if lc in a:
                    if index % 2:
                        x2 = x1 - dx
                        y2 = y1 + dy
                    else:
                        x2 = x1 - dx
                        y2 = y1 - dy
                    if count < mid:
                        if index % 2:
                            x2 = x1 + dx
                            y2 = y1 + dy
                        else:
                            x2 = x1 + dx
                            y2 = y1 - dy

                        if count == mid - 1:
                            x2 = x1
                            y2 = y1 + .825

                elif lc in b:
                    if index % 2:
                        x2 = x1 - dx
                        y2 = y1 - dy
                    else:
                        x2 = x1 - dx
                        y2 = y1 + dy

                    if count < mid:
                        if index % 2:
                            x2 = x1 + dx
                            y2 = y1 - dy
                        else:
                            x2 = x1 + dx
                            y2 = y1 + dy

                        if count == mid - 1:
                            x2 = x1
                            y2 = y1 + .825

                elif lc in c:
                    if index % 2:
                        x2 = x1 - dx
                        y2 = y1 - dy
                    else:
                        x2 = x1 - dx
                        y2 = y1 + dy

                    if count < mid:
                        if index % 2:
                            x2 = x1 + dx
                            y2 = y1 - dy
                        else:
                            x2 = x1 + dx
                            y2 = y1 + dy

                        if count == mid - 1:
                            x2 = x1
                            y2 = y1 + .825

                        if not count:
                            x2 = x1 - dx
                            y2 = y1 - dy

                elif lc in d:
                    if index % 2:
                        x2 = x1 - dx
                        y2 = y1 + dy
                    else:
                        x2 = x1 - dx
                        y2 = y1 - dy

                    if count < mid:
                        if index % 2:
                            x2 = x1 + dx
                            y2 = y1 + dy
                        else:
                            x2 = x1 + dx
                            y2 = y1 - dy

                        if count == mid - 1:
                            x2 = x1
                            y2 = y1 + .825

                plane[atom] = (x2, y2)
                if count:
                    stack.append((index + 1, (x2, y2), count))

    @cached_property
    def __compiled(self):
        atoms = self._atoms
        bonds = self._bonds
        atoms_order = self.atoms_order

        closures = defaultdict(list)
        components = []
        seen = set()
        while len(seen) < len(atoms):
            start = max(atoms.keys() - seen, key=lambda x: atoms_order[x])
            seen.add(start)
            stack = [(n, start, bond) for n, bond in sorted(bonds[start].items(),
                                                            key=lambda x: atoms_order[x[0]])]
            order = [start]
            components.append(order)

            while stack:
                front, back, *_ = atom = stack.pop()
                if front not in seen:
                    order.append(atom)
                    for n, bond in sorted(bonds[front].items(), key=lambda x: atoms_order[x[0]]):
                        if n != back:
                            if n not in seen:
                                stack.append((n, front, bond))
                            else:
                                closures[front].append((n, bond))
                    seen.add(front)
        return components, closures

    @staticmethod
    def vector(n, m, plane):
        nx, ny = plane[n]
        mx, my = plane[m]
        return mx - nx, my - ny

    @staticmethod
    def towards(bond1, bond2):
        return bond1 == 1 and bond2 == 3 or bond1 == 3 and bond2 == 1 or bond1 == 2 and bond2 == 2

    def clean2d(self):
        plane = self._plane
        atoms = self._atoms
        bonds = self._bonds
        atoms_order = self.atoms_order

        dx, dy = n_dist / 2, n_dist * sin(pi / 3)
        components = {k: set(v) for k, v in bonds.items()}
        start = next(n for n, ms in components.items() if len(ms) == 1)
        stack = [[(next_atom, start, bonds[start][next_atom].order, None, None, 0)] for next_atom in components[start]]

        size = len(atoms) - 1
        path = []
        hashed_path = set()
        while stack:
            atom, prev_atom, bond1, current_coordinates, v, _ = stack[-1].pop()
            bond2 = bonds[prev_atom][atom].order
            if current_coordinates is None:
                plane[atom], plane[prev_atom] = (dx, dy), (0, 0)
                vx, vy = mx, my = dx, dy
            else:
                x, y = v
                nx, ny = plane[prev_atom]
                if self.towards(bond1, bond2):
                    current_coordinates = (x + nx, y + ny)
                mx, my = plane[atom] = current_coordinates
                vx, vy = mx - nx, my - ny
            path.append((atom, prev_atom, bond1, current_coordinates))
            hashed_path.add(atom)

            if len(path) == size:
                bad_points = False
                cross = False
                for k, v in plane.items():
                    kx, ky = v
                    for key, value in plane.items():
                        if k == key:
                            continue
                        key_x, key_y = value
                        d_x, d_y = abs(kx - key_x), abs(ky - key_y)
                        if d_x < .1 and d_y < .1:
                            bad_points = True

                # _bonds = list(self.bonds())
                # for i, triple in enumerate(_bonds):
                #     n, m, bond = triple
                #     nx, ny = plane[n]
                #     mx, my = plane[m]
                #     for _triple in _bonds[i:]:
                #         _n, _m, _bond = _triple
                #         _nx, _ny = plane[_n]
                #         _mx, _my = plane[_m]

                if not bad_points and not cross:
                    yield path
                del stack[-1]
                if stack:
                    nnn = stack[-1][-1][-1]
                    path = path[:nnn]
                    hashed_path = {x for x, *_ in path}

            elif atom != start:  # we finished. next step is final closure
                for_stack = []
                closures = []
                loop = 0
                for next_atom in components[atom]:
                    if next_atom == prev_atom:  # only forward. behind us is the homeland
                        continue
                    elif next_atom == start:
                        loop = next_atom
                    elif next_atom in hashed_path:  # closure found
                        closures.append(next_atom)
                    else:
                        for_stack.append(next_atom)

                if loop:
                    for_stack.remove(loop)
                    sx, sy = plane[start]
                    lp = len(path)
                    vx1, vy1 = rotate_vector(dx, dy, vx, vy)
                    vx2, vy2 = rotate_vector(dx, -dy, vx, vy)
                    nx1, ny1 = vx1 + mx, vy1 + my
                    mx2, my2 = vx2 + mx, vy2 + my
                    d1_x, d2_x = abs(nx1 - sx), abs(ny1 - sy)
                    d1_y, d2_y = abs(mx2 - sx), abs(my2 - sy)
                    if d1_x < .001 and d1_y < .001:
                        stack[-1].append((loop, atom, bond2, (nx1, ny1), (vx, vy), lp))
                    elif d2_x < .001 and d2_y < .001:
                        stack[-1].append((loop, atom, bond2, (mx2, my2), (vx, vy), lp))
                    else:
                        del stack[-1]
                        if stack:
                            path = path[:stack[-1][-1][-1]]
                            hashed_path = {x for x, *_ in path}

                elif len(for_stack) == 1:  # easy path grow. next bond double or include single for pyroles
                    lp = len(path)
                    next_atom = for_stack[0]
                    opposite = stack[-1].copy()
                    vx1, vy1 = rotate_vector(dx, dy, vx, vy)
                    vx2, vy2 = rotate_vector(dx, -dy, vx, vy)
                    new_coords1 = (vx1 + mx, vy1 + my)
                    new_coords2 = (vx2 + mx, vy2 + my)
                    new_coords = False
                    # for hash_atom in hashed_path:
                    #     v_x, v_y = plane[hash_atom]
                    for n, v in plane.items():
                        v_x, v_y = v
                        n1_x, n1_y = new_coords1
                        n2_x, n2_y = new_coords2
                        d1_x, d2_x = abs(n1_x - v_x), abs(n2_x - v_x)
                        d1_y, d2_y = abs(n1_y - v_y), abs(n2_y - v_y)
                        if d1_x < .1 and d1_y < .1:
                            new_coords = new_coords2
                            if d2_x < .1 and d2_y < .1:
                                new_coords = None
                        elif d2_x < .1 and d2_y < .1:
                            new_coords = new_coords1
                            if d1_x < .1 and d1_y < .1:
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
                    opposite = stack[-1].copy()
                    vx1, vy1 = rotate_vector(dx, dy, vx, vy)
                    vx2, vy2 = rotate_vector(dx, -dy, vx, vy)
                    new_coords1 = (vx1 + mx, vy1 + my)
                    new_coords2 = (vx2 + mx, vy2 + my)

                    new_coords = False
                    for n, v in plane.items():
                        v_x, v_y = v
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
                        stack[-1].append((next_atom1, atom, bond2, new_coords1, (vx1, vy1), None))
                        stack[-1].append((next_atom2, atom, bond2, new_coords2, (vx2, vy2), lp))
                        if not self.towards(bond2, bonds[atom][next_atom1].order) and \
                                not self.towards(bond2, bonds[atom][next_atom2].order) and \
                                atoms_order[next_atom1] != atoms_order[next_atom2]:
                            opposite.append((next_atom1, atom, bond2, (vx2 + mx, vy2 + my), (vx2, vy2), None))
                            opposite.append((next_atom2, atom, bond2, (vx1 + mx, vy1 + my), (vx1, vy1), lp))
                            stack.append(opposite)

                elif len(for_stack) == 3:
                    dx1, dx2, dx3 = 0, .825, 0
                    dy1, dy2, dy3 = .825, 0, -.825
                    combinations = permutations(for_stack)
                    vx1, vy1 = rotate_vector(dx1, dy1, vx, vy)
                    vx2, vy2 = rotate_vector(dx2, dy2, vx, vy)
                    vx3, vy3 = rotate_vector(dx3, dy3, vx, vy)

                    opposite = stack[-1].copy()
                    opposite1 = stack[-1].copy()
                    lp = len(path)
                    next_atom1, next_atom2, next_atom3 = next(combinations)
                    stack[-1].append((next_atom1, atom, bond2, (vx1 + mx, vy1 + my), (vx1, vy1), None))
                    stack[-1].append((next_atom2, atom, bond2, (vx2 + mx, vy2 + my), (vx2, vy2), None))
                    stack[-1].append((next_atom3, atom, bond2, (vx3 + mx, vy3 + my), (vx3, vy3), lp))
                    for combination in combinations:
                        nxt_atm1, nxt_atm2, nxt_atm3 = combination
                        if atoms_order[next_atom1] == atoms_order[nxt_atm1] \
                                or atoms_order[next_atom1] == atoms_order[nxt_atm2] \
                                or atoms_order[next_atom1] == atoms_order[nxt_atm3]:
                            continue
                        opposite.append((nxt_atm1, atom, bond2, (vx1 + mx, vy1 + my), (vx1, vy1), None))
                        opposite.append((nxt_atm2, atom, bond2, (vx2 + mx, vy2 + my), (vx2, vy2), None))
                        opposite.append((nxt_atm3, atom, bond2, (vx3 + mx, vy3 + my), (vx3, vy3), lp))
                        stack.append(opposite)
                        opposite = opposite1.copy()

                elif closures:
                    atom_x, atom_y = plane[atom]
                    clos_x, clos_y = plane[closures.pop()]
                    vect_x, vect_y = clos_x - atom_x, clos_y - atom_y
                    if not .800 < hypot(vect_x, vect_y) < .850:
                        del stack[-1]
                        if stack:
                            path = path[:stack[-1][-1][-1]]
                            hashed_path = {x for x, *_ in path}


n_dist = .825
dist_3_sm = n_dist / 2
dist_3_la = n_dist*cos(pi/6)
dist_4 = n_dist*cos(pi/4)
templates_1 = [(-dist_3_la, dist_3_sm), (dist_3_la, dist_3_sm), (dist_3_la, -dist_3_sm), (-dist_3_la, -dist_3_sm)]
templates_2 = [((-dist_3_la, dist_3_sm), (dist_3_la, dist_3_sm)),
               ((-dist_3_la, -dist_3_sm), (dist_3_la, -dist_3_sm)),
               ((-dist_3_sm, dist_3_la), (-dist_3_sm, -dist_3_la)),
               ((dist_3_sm, dist_3_la), (dist_3_sm, -dist_3_la))]
templates_3 = [((0, n_dist), (dist_3_la, -dist_3_sm), (-dist_3_la, -dist_3_sm)),
               ((0, -n_dist), (-dist_3_la, dist_3_sm), (dist_3_la, dist_3_sm)),
               ((-n_dist, 0), (dist_3_sm, dist_3_la), (dist_3_sm, -dist_3_la)),
               ((n_dist, 0), (-dist_3_sm, -dist_3_la), (-dist_3_sm, dist_3_la))]
templates_4 = [((-n_dist, 0), (0, n_dist), (n_dist, 0), (0, -n_dist)),
               ((-dist_4, dist_4), (dist_4, dist_4), (dist_4, -dist_4), (-dist_4, -dist_4))]
normal_angles = {1: {1: False, 2: False, 3: True, 4: True},
                 2: {1: False, 2: True, 3: True, 4: True},
                 3: {1: True, 2: True, 3: True, 4: True},
                 4: {1: True, 2: True, 3: True, 4: True}}

__all__ = ['Calculate2D']
