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
from math import hypot, pi, acos, cos, sin, ceil
from ..algorithms.depict import rotate_vector


def rotate_vector2(x1, y1, angle):
    """
    rotate x,y vector over angle
    """
    cos_rad = cos(angle)
    sin_rad = sin(angle)
    return cos_rad * x1 + sin_rad * y1, -sin_rad * x1 + cos_rad * y1


# def optimal(temps, a, b, c):
#     ax, ay = a
#     bx, by = b
#     cx, cy = c
#     min_distance = 10000
#     fa = fc = 0
#     for temp in temps:
#         t1, t2 = temp
#         tx1, ty1 = t1
#         tx2, ty2 = t2
#
#
#         tx1, ty1, tx2, ty2 = tx1 + bx, ty1 + by, tx2 + bx, ty2 + by
#         fa_1, fc_1 = (tx1 - ax, ty1 - ay), (tx2 - cx, ty2 - cy)
#         fa_2, fc_2 = (tx1 - cx, ty1 - cy), (tx2 - ax, ty2 - ay)
#         sum_1 = hypot(*fa_1) + hypot(*fc_1)
#         sum_2 = hypot(*fa_2) + hypot(*fc_2)
#         summ = min(sum_1, sum_2)
#         if summ < min_distance:
#             min_distance = summ
#             if summ == sum_1:
#                 fa, fc = fa_1, fc_1
#             else:
#                 fa, fc = fa_2, fc_2
#     return fa, fc


class Calculate2D:
    __slots__ = ()

    @cached_property
    def __angles(self):
        bonds = self._bonds
        sssr = self.sssr
        angles = []
        for n, m_bond in bonds.items():
            neighbors = tuple(m_bond)
            ln = len(m_bond)
            cycle_angle = None
            for cycle in sssr:
                lc = len(cycle)
                if n in cycle:
                    if lc == 3:
                        cycle_angle = pi / 3
                    else:
                        cycle_angle = pi - 2 * pi / lc
            if ln == 2:
                a, c = neighbors
                ange = normal_angles[bonds[a][n].order][bonds[n][c].order]
                if cycle_angle:
                    ange = cycle_angle
                angles.append((a, n, c, ange))
            # elif ln == 4:
            #     print(m_bond)
            #     print(neighbors)
            #     angles.append((neighbors[0], n, neighbors[3], pi))
            #     angles.append((neighbors[1], n, neighbors[2], pi))
            #     angles.append((neighbors[0], n, neighbors[2], pi / 2))
            #     angles.append((neighbors[2], n, neighbors[3], pi / 2))
            #     angles.append((neighbors[3], n, neighbors[1], pi / 2))
            #     angles.append((neighbors[1], n, neighbors[0], pi / 2))
            #     print(angles)
            #     # for v, w in zip(neighbors, neighbors[1:]):
            #     #     angles.append((v, n, w, pi / 2))
            # elif ln != 1:
            #     ange = 2 * pi / ln
            #     if cycle_angle:
            #         ange = cycle_angle
            #     angles.append((neighbors[-1], n, neighbors[0], ange))
            #     for v, w in zip(neighbors, neighbors[1:]):
            #         angles.append((v, n, w, ange))
        return angles

    def __get_dist_forces(self):
        plane = self._plane

        # for hypot(x, y)
        for n, (x, y) in plane.items():
            if not x and not y:
                plane[n] = (.0001, .0001)
        force_fields = {k: (.00001, .00001) for k in self._atoms}

        # distance forces
        max_diff = 0.1
        for n, m, bond in self.bonds():
            nx, ny = plane[n]
            mx, my = plane[m]
            x, y = nx - mx, ny - my
            c_dist = hypot(x, y)

            diff = abs(n_dist - c_dist)
            if diff > max_diff:
                max_diff = diff

            elong = n_dist / c_dist / 2 - .5
            dx, dy = x * elong, y * elong
            x, y = force_fields[n]
            force_fields[n] = (x + dx, y + dy)
            x, y = force_fields[m]
            force_fields[m] = (x - dx, y - dy)

        return force_fields, max_diff

    def __get_angle_forces(self):
        # angle forces
        #
        # a
        # |
        # b---c
        #
        plane = self._plane
        force_fields = {k: (.00001, .00001) for k in self._atoms}

        max_diff = 0.1
        print(self.__angles)
        for a, b, c, opt in self.__angles:
            ax, ay = plane[a]
            bx, by = plane[b]
            cx, cy = plane[c]
            bax, bay = ax - bx, ay - by
            bcx, bcy = cx - bx, cy - by

            l_ba, l_bc = hypot(bay, bax), hypot(bcy, bcx)
            bax /= l_ba
            bay /= l_ba
            bcx /= l_bc
            bcy /= l_bc

            # bis_x = bax + bcx
            # bis_y = bay + bcy
            # bis_l = hypot(bis_y, bis_x)
            angle = acos(bax * bcx + bay * bcy)
            if angle < .01:
                dx, dy = rotate_vector(0, .2, bcx, bcy)
                fax, fay = force_fields[a]
                force_fields[a] = fax + dx, fay + dy
                fcx, fcy = force_fields[c]
                force_fields[c] = fcx - dx, fcy - dy

            diff = angle - opt
            adiff = abs(diff)
            if adiff > max_diff:
                max_diff = adiff

            rx, ry = rotate_vector2(bcx, bcy, diff)
            rx, ry = rx + bx, ry + by
            fx, fy = rx - cx, ry - cy
            # force = 2 * (opt - angle) * adiff
            # ratio = force / bis_l
            # bis_x *= ratio
            # bis_y *= ratio

            kx, ky = force_fields[b]
            force_fields[b] = (kx + fx, ky + fy)
            # for n, m_bond in bonds.items():
            #     ln = len(m_bond)
            #     nx, ny = plane[n]
            #     if ln == 1:
            #         m, = m_bond
            #         mx, my = plane[m]
            #         min_distance = 10000
            #         fx = fy = 0
            #         for temp in templates_1:
            #             tx, ty = temp
            #             tx, ty = tx + nx, ty + ny
            #             fm = (tx - mx, ty - my)
            #             summ = hypot(*fm)
            #             if summ < min_distance:
            #                 min_distance = summ
            #                 fx, fy = fm
            #         fmx, fmy = force_fields[m]
            #         force_fields[m] = fmx + fx, fmy + fy
        return force_fields, max_diff

    def _changes(self, forces):
        """
        changes in coordinates
        :param forces: dict of atoms forces
        :return:
        """
        plane = self._plane

        stack = 1
        steps = 1
        while stack:
            # for x in range(1):
            stack = 0
            force = max(hypot(x, y) for x, y in forces.values())
            ratio = .2 / force
            if ratio > 1:
                ratio = 1
            for atom, (x2, y2) in forces.items():
                x1, y1 = plane[atom]
                x2, y2 = x2 * ratio, y2 * ratio
                plane[atom] = (x1 + x2, y1 + y2)
            if force > .05:
                stack = 1
            if steps >= 200:
                break
            steps += 1

    def clean2d(self):
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
            print(cycle)
            a = tuple(range(9, 300, 4))
            b = tuple(range(10, 300, 4))
            d = tuple(range(11, 300, 4))
            e = tuple(range(12, 301, 4))
            if lc in a:
                pass
            elif lc in b:
                pass
            elif lc in d:
                pass
            elif lc in e:
                pass
            else:
                raise Exception
            angle = 2 * pi / 3
            dx, dy = .825 / 2, sin(pi / 3) * .825
            c = cycle[2]
            plane[cycle[0]], plane[cycle[1]], plane[c] = (.825, 0), (0, 0), (-dx, -dy)
            mid = ceil(lc / 2)
            stack = [(3, plane[c], lc - 3)]
            while stack:
                index, coords, count = stack.pop(0)
                x1, y1 = coords
                count -= 1
                atom = cycle[index]

                if lc in a:
                    if index % 2:
                        x2 = x1 - .825
                        y2 = y1
                    else:
                        x2 = x1 - dx
                        y2 = y1 - dy
                    if count < mid:
                        if index % 2:
                            x2 = x1 + dx
                            y2 = y1 + dy
                        else:
                            x2 = x1 + .825
                            y2 = y1

                        if count == mid - 1:
                            x2 = x1 - dx
                            y2 = y1 + dy

                elif lc in b:
                    if index % 2:
                        x2 = x1 - .825
                        y2 = y1
                    else:
                        x2 = x1 - dx
                        y2 = y1 + dy

                    if count < mid:
                        if index % 2:
                            x2 = x1 + dx
                            y2 = y1 - dy
                        else:
                            x2 = x1 + .825
                            y2 = y1

                        if count == mid - 1:
                            x2 = x1 + dx
                            y2 = y1 + dy

                # if count == 6:
                #     break
                # x2, y2 = rotate_vector2(x1, y1, direction)
                # x2 += .825
                plane[atom] = (x2, y2)
                print(atom)
                print(x2, y2)
                print(count)
                if count:
                    stack.append((index + 1, (x2, y2), count))
            # angle = 2 * pi / 3
            # alternate = False
            # if angle < pi / 6:
            #     alternate = True
            # chain = [(0, 0), (.825, 0)]
            # lc -= 1
            # a, b = (.825, 0), (2 * .825, 0)
            # ax, ay = a
            # bx, by = b
            # abx, aby = bx - ax, by - ay
            # abx, aby = rotate_vector2(abx, aby, angle)
            # stack = [(abx, aby)]
            # abx, aby = rotate_vector2(abx, aby, -angle)
            # stack.append((abx, aby))

            # while stack:
            #     x1, y1 = stack[1]
            #     x2, y2 = stack[2]
            #     stack = []
            #     lc -= 1
            #     chain.append((x1, y1))
            #     if lc:
            #         x3, y3 = x1 * 2, y1 * 2
            #         stack.append(3)


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
normal_angles = {1: {1: 2/3*pi, 2: 2/3*pi, 3: pi, 4: pi},
                 2: {1: 2/3*pi, 2: pi, 3: pi, 4: pi},
                 3: {1: pi, 2: pi, 3: pi, 4: pi},
                 4: {1: pi, 2: pi, 3: pi, 4: pi}}

__all__ = ['Calculate2D']
