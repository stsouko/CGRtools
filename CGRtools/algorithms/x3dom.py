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
from math import acos, sqrt, cos, sin, ceil


def plane_normal(nmx, nmy, nmz, nox, noy, noz):
    # return normal to plane of two vectors nm and no
    # m <--- n
    #         \
    #          v
    #          o
    return nmy * noz - nmz * noy, nox * nmz - nmx * noz, nmx * noy - nmy * nox


def unit_vector(nmx, nmy, nmz):
    nmd = sqrt(nmx ** 2 + nmy ** 2 + nmz ** 2)
    return nmx / nmd, nmy / nmd, nmz / nmd


# @njit(f8(f8, f8, f8, f8, f8, f8), cache=True)
def get_angle(nx, ny, nz, mx, my, mz):
    ch = nx * mx + ny * my + nz * mz
    zn = (nx ** 2 + ny ** 2 + nz ** 2) ** .5 * (mx ** 2 + my ** 2 + mz ** 2) ** .5
    if zn < .0001:
        angle = .0
    else:
        angle = acos(ch/zn)
    return angle


def vector_normal(nmx, nmy, nmz):
    # return normal to vector nm
    if not -.0001 < nmx < .0001:
        return (- nmy - nmz) / nmx, 1, 1
    elif not -.0001 < nmy < .0001:
        return 1, (- nmx - nmz) / nmy, 1
    else:
        return 1, 1, (- nmx - nmy) / nmz


class JupyterWidget:
    def __init__(self, xml, width, height):
        self.xml = xml
        self.width = width
        self.height = height

    def _repr_html_(self):
        return ("<script type='text/javascript' src='https://www.x3dom.org/download/x3dom.js'></script>"
                "<link rel='stylesheet' type='text/css' href='https://www.x3dom.org/download/x3dom.css'>"
                f'<div style="width: {self.width}; height: {self.height}">{self.xml}</div>')

    def __html__(self):
        return self._repr_html_()


class X3dom:
    __slots__ = ()

    def depict3d(self, index: int = 0) -> str:
        """Get X3DOM XML string.

        :param index: index of conformer
        """
        xyz = self._conformers[index]
        mx = sum(x for x, _, _ in xyz.values()) / len(xyz)
        my = sum(y for _, y, _ in xyz.values()) / len(xyz)
        mz = sum(z for _, _, z in xyz.values()) / len(xyz)
        xyz = {n: (x - mx, y - my, z - mz) for n, (x, y, z) in xyz.items()}
        atoms = self.__render_atoms(xyz)
        bonds = self._render_3d_bonds(xyz)
        return f'<x3d width=100% height=100%>\n  <scene>\n{atoms}{bonds}  </scene>\n</x3d>'

    def view3d(self, index: int = 0, width='600px', height='400px'):
        """
        Jupyter widget for 3D visualization.

        :param index: index of conformer
        :param width: widget width
        :param height: widget height
        """
        return JupyterWidget(self.depict3d(index), width, height)

    def __render_atoms(self, xyz):
        config = self._render_config
        colors = config['atoms_colors']
        mapping_color = config['mapping_color']
        carbon = config['carbon']
        font = config['font_size']
        radius = config['atom_radius'] = -.1
        if radius < 0:
            multiplier = -radius
            radius = 0
        elif not radius:
            multiplier = .2

        atoms = []
        if carbon:
            for n, a in self._atoms.items():
                r = radius or a.atomic_radius * multiplier
                fr = r * 0.71
                atoms.append(f"    <transform translation='{' '.join(format(x, '.2f') for x in xyz[n])}'>\n"
                             "      <billboard axisOfRotation='0 0 0'>\n        <group>\n"
                             f"          <transform translation='{fr:.2f} {fr:.2f} 0'>\n"
                             '            <shape>\n              <appearance>\n'
                             f"                <material diffuseColor='{mapping_color}'/>\n"
                             f"              </appearance>\n              <text string='{a.atomic_symbol}'>\n"
                             f"                <fontstyle family='sans' size='{font:.2f}' justify='first'/>\n"
                             "              </text>\n            </shape>\n          </transform>\n"
                             "          <shape>\n            <appearance>\n"
                             f"              <material diffuseColor='{colors[a.atomic_number - 1]}'/>\n"
                             f"            </appearance>\n            <sphere radius='{r:.2f}'/>\n"
                             "          </shape>\n        </group>\n      </billboard>\n    </transform>\n")
        else:
            for n, a in self._atoms.items():
                r = radius or a.atomic_radius * multiplier
                atoms.append(f"    <transform translation='{' '.join(format(x, '.2f') for x in xyz[n])}'>\n"
                             "      <shape>\n        <appearance>\n"
                             f"          <material diffuseColor='{colors[a.atomic_number - 1]}'/>\n"
                             f"        </appearance>\n        <sphere radius='{r:.2f}'/>\n"
                             "      </shape>\n    </transform>\n")
        return ''.join(atoms)


class X3domMolecule(X3dom):
    __slots__ = ()

    def _render_3d_bonds(self, xyz):
        bonds = self._bonds
        lengths = self.lengths
        doubles = self.doubles
        config = self._render_config

        double_space = config['double_space']
        triple_space = .25
        dash1, dash2 = config['dashes']
        bond_color = '#808080'
        radius = .02
        xml = []
        for n, m, bond in self.bonds():
            order = bond.order
            # order = 3
            nx, ny, nz = xyz[n]
            mx, my, mz = xyz[m]

            nmx, nmy, nmz = mx - nx, my - ny, mz - nz
            length = sqrt(nmx ** 2 + nmy ** 2 + nmz ** 2)
            if length < .001:
                continue

            rotation_angle = acos(nmy / length)
            lengths[(n, m)] = lengths[(m, n)] = (length, rotation_angle)
            x, y, z = nx + nmx / 2, ny + nmy / 2, nz + nmz / 2
            if order in (1, 4):
                xml.append(f"    <transform translation='{x:.2f} {y:.2f} {z:.2f}' rotation='{nmz:.2f} 0 "
                           f"{-nmx:.2f} {rotation_angle:.2f}'>\n      <shape>\n        <appearance>\n"
                           f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                           f"       </appearance>\n        <cylinder radius='{radius}' height='{length:.2f}'>\n"
                           "        </cylinder>\n      </shape>\n    </transform>\n")
            elif order == 2:
                if n in doubles:
                    # normal for plane n m o
                    norm_x, norm_y, norm_z = plane_normal(nmx, nmy, nmz, *doubles[n])
                elif m in doubles:
                    # normal for plane n m o
                    norm_x, norm_y, norm_z = plane_normal(nmx, nmy, nmz, *doubles[m])
                else:
                    third = sorted([x for x in bonds[n] if x != m], key=lambda a: len(bonds[a]))
                    if third:
                        ox, oy, oz = xyz[third[0]]
                        nox, noy, noz = ox - nx, oy - ny, oz - nz
                    else:
                        nox, noy, noz = vector_normal(nmx, nmy, nmz)

                    # normal for plane n m o
                    normx, normy, normz = unit_vector(*plane_normal(nmx, nmy, nmz, nox, noy, noz))

                    # normal for plane n m normal
                    norm_x, norm_y, norm_z = plane_normal(nmx, nmy, nmz, normx, normy, normz)

                doubles[n] = doubles[m] = (norm_x, norm_y, norm_z)
                norm_dist = sqrt(norm_x ** 2 + norm_y ** 2 + norm_z ** 2)

                if norm_dist < .0001:
                    coef = 1
                else:
                    coef = double_space / norm_dist

                dx, dy, dz = norm_x * coef, norm_y * coef, norm_z * coef
                xml.append(
                    f"    <transform translation='{x + dx:.2f} {y + dy:.2f} {z + dz:.2f}' rotation='{nmz:.2f} 0 "
                    f"{-nmx:.2f} {rotation_angle:.2f}'>\n      <shape>\n        <appearance>\n"
                    f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                    f"       </appearance>\n        <cylinder radius='{radius}' height='{length:.2f}'>\n"
                    "        </cylinder>\n      </shape>\n    </transform>\n")
                xml.append(
                    f"    <transform translation='{x - dx:.2f} {y - dy:.2f} {z - dz:.2f}' rotation='{nmz:.2f} 0 "
                    f"{-nmx:.2f} {rotation_angle:.2f}'>\n      <shape>\n        <appearance>\n"
                    f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                    f"       </appearance>\n        <cylinder radius='{radius}' height='{length:.2f}'>\n"
                    "        </cylinder>\n      </shape>\n    </transform>\n")
            elif order == 3:
                d = triple_space / 2
                r1 = triple_space * sqrt(3) / 3
                r2 = triple_space * sqrt(3) / 6
                nox, noy, noz = vector_normal(nmx, nmy, nmz)

                # normal for plane n m o
                normx, normy, normz = unit_vector(*plane_normal(nmx, nmy, nmz, nox, noy, noz))
                vecRx, vecRy, vecRz = normx * r1, normy * r1, normz * r1

                # normal for plane n m normal
                norm_x, norm_y, norm_z = unit_vector(*plane_normal(nmx, nmy, nmz, normx, normy, normz))
                vecx, vecy, vecz = norm_x * d, norm_y * d, norm_z * d

                xml.append(f"    <transform translation='{x + vecRx:.2f} {y + vecRy:.2f} {z + vecRz:.2f}'"
                           f" rotation='{nmz:.2f} 0 {-nmx:.2f} {rotation_angle:.2f}'>\n      <shape>\n"
                           f"        <appearance>\n          <material diffusecolor='{bond_color}'>\n"
                           f"          </material>\n       </appearance>\n        <cylinder radius='{radius}'"
                           f" height='{length:.2f}'>\n        </cylinder>\n      </shape>\n    </transform>\n")

                xx, yy, zz = x - normx * r2, y - normy * r2, z - normz * r2
                xml.append(f"    <transform translation='{xx - vecx:.2f} {yy - vecy:.2f} {zz - vecz:.2f}'"
                           f" rotation='{nmz:.2f} 0 {-nmx:.2f} {rotation_angle:.2f}'>\n      <shape>\n"
                           f"        <appearance>\n          <material diffusecolor='{bond_color}'>\n"
                           f"          </material>\n       </appearance>\n        <cylinder radius='{radius}'"
                           f" height='{length:.2f}'>\n        </cylinder>\n      </shape>\n    </transform>\n")
                xml.append(f"    <transform translation='{xx + vecx:.2f} {yy + vecy:.2f} {zz + vecz:.2f}'"
                           f" rotation='{nmz:.2f} 0 {-nmx:.2f} {rotation_angle:.2f}'>\n      <shape>\n"
                           f"        <appearance>\n          <material diffusecolor='{bond_color}'>\n"
                           f"          </material>\n       </appearance>\n        <cylinder radius='{radius}'"
                           f" height='{length:.2f}'>\n        </cylinder>\n      </shape>\n    </transform>\n")
            else:
                dash = dash1 + dash2
                d = dash / length
                dx, dy, dz = nmx * d, nmy * d, nmz * d

                b = int((length - dash1) // dash)
                t = (length - (b * dash)) / length
                nx, ny, nz = nx + nmx * t / 2, ny + nmy * t / 2, nz + nmz * t / 2
                for _ in range(b):
                    xml.append(f"    <transform translation='{nx:.2f} {ny:.2f} {nz:.2f}' rotation='{nmz:.2f} 0 "
                               f"{-nmx:.2f} {rotation_angle:.2f}'>\n      <shape>\n        <appearance>\n"
                               f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                               f"       </appearance>\n        <cylinder radius='{radius}' height='{dash1:.2f}'>\n"
                               "        </cylinder>\n      </shape>\n    </transform>\n")
                    nx += dx
                    ny += dy
                    nz += dz
                xml.append(f"    <transform translation='{nx:.2f} {ny:.2f} {nz:.2f}' rotation='{nmz:.2f} 0 "
                           f"{-nmx:.2f} {rotation_angle:.2f}'>\n      <shape>\n        <appearance>\n"
                           f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                           f"       </appearance>\n        <cylinder radius='{radius}' height='{dash1:.2f}'>\n"
                           "        </cylinder>\n      </shape>\n    </transform>\n")

        # lengths = self.lengths
        # for ring in self.connected_rings:
        #     cx = sum(xyz[n][0] for n in ring) / len(ring)
        #     cy = sum(xyz[n][1] for n in ring) / len(ring)
        #     cz = sum(xyz[n][2] for n in ring) / len(ring)
        #
        #     for n, m in zip(ring, ring[1:]):
        #         nx, ny, nz = xyz[n]
        #         mx, my, mz = xyz[m]
        #         nmx, nmy, nmz = mx - nx, my - ny, mz - nz
        #         if (n, m) in lengths:
        #             length, rotation_angle = lengths[(n, m)]
        #         else:
        #             length = sqrt(nmx ** 2 + nmy ** 2 + nmz ** 2)
        #             rotation_angle = acos(nmy / length)
        #         if length < .001:
        #             continue
        #
        #         aromatic = self.__render_aromatic_bond(nx, ny, nz, mx, my, mz, cx, cy, cz, length)
        #         if aromatic:
        #             ax, ay, az, bx, by, bz = aromatic
        #             abx, aby, abz = bx - ax, by - ay, bz - az
        #             for _ in range(1, 5):
        #                 ax += abx * .15
        #                 ay += aby * .15
        #                 nz += nmz * .15
        #                 xml.append(
        #                     f"    <transform translation='{ax:.2f} {ay:.2f} {nz:.2f}' rotation='{nmz:.2f} 0 "
        #                     f"{-nmx:.2f} {rotation_angle:.2f}'>\n      <shape>\n        <appearance>\n"
        #                     f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
        #                     f"       </appearance>\n        <cylinder radius='{radius}' height='{dash2:.2f}'>\n"
        #                     "        </cylinder>\n      </shape>\n    </transform>\n")
        #
        #     i, j = ring[-1], ring[0]
        #     nx, ny, nz = xyz[i]
        #     mx, my, mz = xyz[j]
        #     nmx, nmy, nmz = mx - nx, my - ny, mz - nz
        #     if (i, j) in lengths:
        #         length, rotation_angle = lengths[(i, j)]
        #     else:
        #         length = sqrt(nmx ** 2 + nmy ** 2 + nmz ** 2)
        #         rotation_angle = acos(nmy / length)
        #     if length < .001:
        #         continue
        #
        #     aromatic = self.__render_aromatic_bond(nx, ny, nz, mx, my, mz, cx, cy, cz, length)
        #     if aromatic:
        #         ax, ay, az, bx, by, bz = aromatic
        #         abx, aby, abz = bx - ax, by - ay, bz - az
        #         for _ in range(1, 5):
        #             ax += abx * .15
        #             ay += aby * .15
        #             nz += nmz * .15
        #             xml.append(f"    <transform translation='{ax:.2f} {ay:.2f} {nz:.2f}' rotation='{nmz:.2f} {0} "
        #                        f"{-nmx:.2f} {rotation_angle:.2f}'>\n      <shape>\n        <appearance>\n"
        #                        f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
        #                        f"       </appearance>\n        <cylinder radius='{radius}' height='{dash2:.2f}'>\n"
        #                        "        </cylinder>\n      </shape>\n    </transform>\n")
        return ''.join(xml)

    def __render_aromatic_bond(self, n_x, n_y, n_z, m_x, m_y, m_z, c_x, c_y, c_z, nm_ln):
        config = self._render_config

        aromatic_space = config['aromatic_space'] = 2
        dash3, dash4 = config['aromatic_dashes']
        # n aligned xyz
        nm_x, nm_y, nm_z, nc_x, nc_y, nc_z = m_x - n_x, m_y - n_y, m_z - n_z, c_x - n_x, c_y - n_y, c_z - n_y
        nc_ln = sqrt(nc_x ** 2 + nc_y ** 2 + nc_z ** 2)
        angle = get_angle(nm_x, nm_y, nm_z, nc_x, nc_y, nc_z)

        # nm reoriented xy
        nmr_x, nmr_y = nm_ln, 0
        cnr_x, cnr_y = cos(angle) * nc_ln, sin(angle) * nc_ln

        if -.0001 < cnr_y < .0001:
            coef = 1
        else:
            coef = aromatic_space / cnr_y

        a_x, a_y, a_z = nm_x * coef, nm_y * coef, nm_z * coef
        b_x, b_y, b_z = (c_x - m_x + nm_ln) * coef, (c_y - m_y) * coef, (c_z - m_z) * coef

        # if cnr_y and aromatic_space / cnr_y < .65:
        #     if cnr_y > 0:
        #         r_y = aromatic_space
        #     else:
        #         r_y = -aromatic_space
        #         cr_y = -cnr_y
        #
        #     ar_x = aromatic_space * cr_x / cr_y
        #     br_x = mr_x - aromatic_space * (mr_x - cr_x) / cr_y
        #
        #     # backward reorienting
        #     an_x, an_y = rotate_vector(ar_x, r_y, mn_x, mn_y)
        #     bn_x, bn_y = rotate_vector(br_x, r_y, mn_x, mn_y)
        #     a_x, a_y = n_x + an_x, n_y + an_y
        #     b_x, b_y = n_x + bn_x, n_y + bn_y

        return a_x, a_y, a_z, b_x, b_y, b_z

    doubles = {}
    lengths = {}


class X3domCGR(X3dom):
    __slots__ = ()

    def _render_3d_bonds(self, xyz):
        config = self._render_config
        return ''


__all__ = ['X3domMolecule', 'X3domCGR']
