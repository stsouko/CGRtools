# -*- coding: utf-8 -*-
#
#  Copyright 2019 Ramil Nugmanov <stsouko@live.ru>
#  Copyright 2019 Dinar Batyrshin <batyrshin-dinar@mail.ru>
#  This file is part of CGRtools.
#
from numpy import dot, linalg, zeros
from math import hypot, pi, atan2
from ..algorithms.depict import rotate_vector

distance_0 = 0.825
normal_angles = {1: {1: 2/3*pi, 2: 2/3*pi, 3: pi},
                 2: {1: 2/3*pi, 2: pi, 3: pi},
                 3: {1: pi, 2: pi, 3: pi}}


def angle(v1, v2):
    '''
    :params v1, v2: input vectors
    :return: angle between vectors v1, v2 in radians
    '''
    return atan2(linalg.det([v1, v2]), dot(v1, v2))


class Clean2D:
    def clean(self):
        self.matrix_of_force_fields()
    
    def matrix_of_force_fields(self):
        plane = self._plane
        atoms = self._atoms
        bonds = self._bonds
        size = len(atoms)
        matrix = zeros((2, size))
        bond = (0, 0, 0)
        for n in atoms.keys():
            nx, ny = plane[n]
            bond_1 = bonds[n]
            for m in atoms.keys():
                if not n == m:
                    mx, my = plane[m]
                    x, y = mx-nx, my-ny
                    distance = hypot(x, y)
                    matrix[0][n-1], matrix[1][n-1] = rotate_vector(0, distance_0 - distance, x, y)
                    if m in bond_1:
                        if bond[1] == n:
                            o = bond[0]
                            ox, oy = plane[o]
                            ang = angle((nx-ox, ny-oy), (x, y))
                            sign = normal_angles[bond[2].order][bond_1[m].order] - ang
                        bond = (n, m, bond_1[m])
