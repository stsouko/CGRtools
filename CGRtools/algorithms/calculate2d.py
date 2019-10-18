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
from itertools import product
from math import hypot, pi, acos, cos, sin, ceil, atan2
from ..algorithms.depict import rotate_vector


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

    def calculate_chains(self):
        plane = self._plane
        atoms = self._atoms
        bonds = self._bonds

        dx, dy = n_dist * sin(pi / 3), n_dist / 2
        components = {k: set(v) for k, v in bonds.items()}
        start = next(n for n, ms in components.items() if len(ms) == 1)
        stack = [[(next_atom, start, bonds[start][next_atom].order, 'up1', None, 0)] for next_atom in components[start]]

        size = len(atoms) - 1
        good_path = []
        path = []
        hashed_path = set()
        while stack:
            atom, prev_atom, bond, direct, v, _ = stack[-1].pop()
            if prev_atom == start:
                plane[atom], plane[prev_atom] = (dx, dy), (0, 0)
            else:
                # hypot can`t for vx == 0 and vy in (-.825, 825)
                old_x, old_y = plane[prev_atom]
                if self.towards(bond, bonds[prev_atom][atom].order):
                    vx, vy = v
                    ang = hypot(vx, vy)
                    if ang > 0:
                        if ang == vy:
                            plane[atom] = (old_x, old_y + .825)
                            direct = 'down1'
                        else:
                            plane[atom] = (old_x + dx, old_y + dy)
                    else:
                        if ang == -vy:
                            plane[atom] = (old_x, old_y - .825)
                            direct = 'up1'
                        else:
                            plane[atom] = (old_x + dx, old_y - dy)
                elif direct == 'up1':
                    plane[atom] = (old_x + dx, old_y + dy)
                elif direct == 'up2':
                    plane[atom] = (old_x, old_y + .825)
                elif direct == 'down1':
                    plane[atom] = (old_x + dx, old_y - dy)
                elif direct == 'down2':
                    plane[atom] = (old_x, old_y - .825)

            path.append((atom, prev_atom, bond))
            hashed_path.add(atom)

            if len(path) == size:
                good_path.append(path)
                del stack[-1]
                if stack:
                    path = path[:stack[-1][-1][-1]]
                    hashed_path = {x for x, *_ in path}

            elif atom != start:
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

                if loop:  # we found starting point.
                    if direct == 'up1':  # finish should be single bonded
                        del stack[-1]
                        if stack:
                            path = path[:stack[-1][-1][-1]]
                            hashed_path = {x for x, *_ in path}
                        continue
                    else:  # finish should be double bonded
                        vx, vy = self.vector(prev_atom, atom, plane)
                        stack[-1].insert(0, (loop, atom, bonds[prev_atom][atom].order, 'down1', (vx, vy), None))
                        direct = 'up1'
                # if direct == 'up1':  # double in - single out. quinone has two single bonds
                #     # for next_atom in closures:
                #     #     path.append((next_atom, atom, bonds[prev_atom][atom].order))  # closures always single-bonded
                #     #     stack[-1].remove((atom, next_atom, 1, None))  # remove fork from stack
                #     for next_atom in for_stack:
                #         vx, vy = self.vector(prev_atom, atom, plane)
                #         stack[-1].append((next_atom, atom, bonds[prev_atom][atom].order, 'down1', (vx, vy), None))
                # elif direct == 'up2':
                #     for next_atom in for_stack:
                #         vx, vy = self.vector(prev_atom, atom, plane)
                #         stack[-1].append((next_atom, atom, bonds[prev_atom][atom].order, 'up1', (vx, vy), None))
                if len(for_stack) == 1:  # easy path grow. next bond double or include single for pyroles
                    next_atom = for_stack[0]
                    vx, vy = self.vector(prev_atom, atom, plane)
                    if direct == 'down1':
                        d = 'up1'
                    else:
                        d = 'down1'
                    stack[-1].append((next_atom, atom, bonds[prev_atom][atom].order, d, (vx, vy), None))
                    # if closures:
                    #     next_atom = closures[0]
                    #     path.append((next_atom, atom, 1))  # closures always single-bonded
                    #     stack[-1].remove((atom, next_atom, 1, None))  # remove fork from stack
                elif for_stack:  # fork
                    next_atom1, next_atom2 = for_stack
                    # if next_atom1 in double_bonded:  # quinone next from fork
                    #     if next_atom2 in double_bonded:  # bad path
                    #         del stack[-1]
                    #         if stack:
                    #             path = path[:stack[-1][-1][-1]]
                    #             hashed_path = {x for x, *_ in path}
                    #     else:
                    vx, vy = self.vector(prev_atom, atom, plane)
                    if direct == 'up1':
                        stack[-1].append((next_atom1, atom, bonds[prev_atom][atom].order, 'up2', (vx, vy), None))
                        stack[-1].append((next_atom2, atom, bonds[prev_atom][atom].order, 'down1', (vx, vy), None))
                    elif direct == 'down1':
                        stack[-1].append((next_atom1, atom, bonds[prev_atom][atom].order, 'up1', (vx, vy), None))
                        stack[-1].append((next_atom2, atom, bonds[prev_atom][atom].order, 'down2', (vx, vy), None))
                    # elif next_atom2 in double_bonded:  # quinone next from fork
                    #     stack[-1].append((next_atom2, atom, 1, None))
                    #     stack[-1].append((next_atom1, atom, 2, None))
                    # else:  # new path
                    #     opposite = stack[-1].copy()
                    #     stack[-1].append((next_atom1, atom, 1, None))
                    #     stack[-1].append((next_atom2, atom, 2, len(path)))
                    #     opposite.append((next_atom2, atom, 1, None))
                    #     opposite.append((next_atom1, atom, 2, None))
                    #     stack.append(opposite)

    def clean2d(self):
        cycles = self.sssr
        plane = self._plane

        if cycles:
            self.calculate_cycles()
        else:
            self.calculate_chains()


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
