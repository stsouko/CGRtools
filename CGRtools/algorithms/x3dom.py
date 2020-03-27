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
from math import acos, sqrt, hypot
from .depict import rotate_vector


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
        font = config['font']
        radius = config['atom_radius']
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
        config = self._render_config

        double_space = config['double_space']
        triple_space = .07
        dash1, dash2 = config['dashes']
        bond_color = '#808080'
        radius = .02
        xml = []
        for n, m, bond in self.bonds():
            order = bond.order
            nx, ny, nz = xyz[n]
            mx, my, mz = xyz[m]

            nmx, nmy, nmz = mx - nx, my - ny, mz - nz
            # norm_x, norm_y, norm_z = 0, 1, 0
            dev = sqrt(nmx ** 2 + nmy ** 2 + nmz ** 2)
            if dev < .001:
                continue

            angle = acos(nmy / dev)
            x, y, z = nx + nmx / 2, ny + nmy / 2, nz + nmz / 2
            length = sqrt(nmx ** 2 + nmy ** 2 + nmz ** 2)

            if order in (1, 4):
                xml.append(f"    <transform translation='{x:.2f} {y:.2f} {z:.2f}' rotation='{nmz:.2f} {0} {-nmx:.2f} "
                           f"{angle:.2f}'>\n      <shape>\n        <appearance>\n"
                           f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                           f"       </appearance>\n        <cylinder radius='{radius}' height='{length:.2f}'>\n"
                           "        </cylinder>\n      </shape>\n    </transform>\n")
            elif order == 2:
                dx, dy = rotate_vector(0, double_space, nmx, nmy)
                xml.append(
                    f"    <transform translation='{x + dx:.2f} {y - dy:.2f} {z:.2f}' rotation='{nmz:.2f} {0} {-nmx:.2f} "
                    f"{angle:.2f}'>\n      <shape>\n        <appearance>\n"
                    f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                    f"       </appearance>\n        <cylinder radius='{radius}' height='{length:.2f}'>\n"
                    "        </cylinder>\n      </shape>\n    </transform>\n")
                xml.append(
                    f"    <transform translation='{x - dx:.2f} {y + dy:.2f} {z:.2f}' rotation='{nmz:.2f} {0} {-nmx:.2f} "
                    f"{angle:.2f}'>\n      <shape>\n        <appearance>\n"
                    f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                    f"       </appearance>\n        <cylinder radius='{radius}' height='{length:.2f}'>\n"
                    "        </cylinder>\n      </shape>\n    </transform>\n")
            elif order == 3:
                dx, dy = rotate_vector(0, triple_space, nmx, nmy)
                dz = hypot(dx, dy)
                xml.append(
                    f"    <transform translation='{x + dx:.2f} {y - dy:.2f} {z + dz:.2f}' rotation='{nmz:.2f} {0} {-nmx:.2f} "
                    f"{angle:.2f}'>\n      <shape>\n        <appearance>\n"
                    f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                    f"       </appearance>\n        <cylinder radius='{radius}' height='{length:.2f}'>\n"
                    "        </cylinder>\n      </shape>\n    </transform>\n")
                xml.append(
                    f"    <transform translation='{x:.2f} {y:.2f} {z - dz:.2f}' rotation='{nmz:.2f} {0} {-nmx:.2f} "
                    f"{angle:.2f}'>\n      <shape>\n        <appearance>\n"
                    f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                    f"       </appearance>\n        <cylinder radius='{radius}' height='{length:.2f}'>\n"
                    "        </cylinder>\n      </shape>\n    </transform>\n")
                xml.append(
                    f"    <transform translation='{x - dx:.2f} {y + dy:.2f} {z + dz:.2f}' rotation='{nmz:.2f} {0} {-nmx:.2f} "
                    f"{angle:.2f}'>\n      <shape>\n        <appearance>\n"
                    f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                    f"       </appearance>\n        <cylinder radius='{radius}' height='{length:.2f}'>\n"
                    "        </cylinder>\n      </shape>\n    </transform>\n")
            else:
                for _ in range(1, 5):
                    nx += nmx * .2
                    ny += nmy * .2
                    nz += nmz * .2
                    xml.append(f"    <transform translation='{nx:.2f} {ny:.2f} {nz:.2f}' rotation='{nmz:.2f} {0} "
                               f"{-nmx:.2f} {angle:.2f}'>\n      <shape>\n        <appearance>\n"
                               f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                               f"       </appearance>\n        <cylinder radius='{radius}' height='{dash1:.2f}'>\n"
                               "        </cylinder>\n      </shape>\n    </transform>\n")

        for ring in self.connected_rings:
            cx = sum(xyz[n][0] for n in ring) / len(ring)
            cy = sum(xyz[n][1] for n in ring) / len(ring)

            for n, m in zip(ring, ring[1:]):
                nx, ny, nz = xyz[n]
                mx, my, mz = xyz[m]
                nmx, nmy, nmz = mx - nx, my - ny, mz - nz
                dev = sqrt(nmx ** 2 + nmy ** 2 + nmz ** 2)
                if dev < .001:
                    continue

                angle = acos(nmy / dev)
                aromatic = self.__render_aromatic_bond(nx, ny, nz, mx, my, mz, cx, cy)
                if aromatic:
                    ax, ay, bx, by = aromatic
                    abx, aby = bx - ax, by - ay
                    for _ in range(1, 5):
                        ax += abx * .2
                        ay += aby * .2
                        nz += nmz * .2
                        xml.append(
                            f"    <transform translation='{ax:.2f} {ay:.2f} {nz:.2f}' rotation='{nmz:.2f} {0} "
                            f"{-nmx:.2f} {angle:.2f}'>\n      <shape>\n        <appearance>\n"
                            f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                            f"       </appearance>\n        <cylinder radius='{radius}' height='{dash2:.2f}'>\n"
                            "        </cylinder>\n      </shape>\n    </transform>\n")

            nx, ny, nz = xyz[ring[-1]]
            mx, my, mz = xyz[ring[0]]
            nmx, nmy, nmz = mx - nx, my - ny, mz - nz
            dev = sqrt(nmx ** 2 + nmy ** 2 + nmz ** 2)
            if dev < .001:
                continue

            angle = acos(nmy / dev)
            aromatic = self.__render_aromatic_bond(nx, ny, nz, mx, my, mz, cx, cy)
            if aromatic:
                ax, ay, bx, by = aromatic
                abx, aby = bx - ax, by - ay
                for _ in range(1, 5):
                    ax += abx * .2
                    ay += aby * .2
                    nz += nmz * .2
                    xml.append(f"    <transform translation='{ax:.2f} {ay:.2f} {nz:.2f}' rotation='{nmz:.2f} {0} "
                               f"{-nmx:.2f} {angle:.2f}'>\n      <shape>\n        <appearance>\n"
                               f"          <material diffusecolor='{bond_color}'>\n          </material>\n"
                               f"       </appearance>\n        <cylinder radius='{radius}' height='{dash2:.2f}'>\n"
                               "        </cylinder>\n      </shape>\n    </transform>\n")
        return ''.join(xml)

    def __render_aromatic_bond(self, n_x, n_y, n_z, m_x, m_y, m_z, c_x, c_y):
        config = self._render_config

        aromatic_space = config['aromatic_space']
        dash3, dash4 = config['aromatic_dashes']
        # n aligned xy
        mn_x, mn_y, mn_z, cn_x, cn_y = m_x - n_x, m_y - n_y, m_z - n_z, c_x - n_x, c_y - n_y

        # nm reoriented xy
        mr_x, mr_y = hypot(mn_x, mn_y), 0
        cr_x, cr_y = rotate_vector(cn_x, cn_y, mn_x, -mn_y)

        if cr_y and aromatic_space / cr_y < .65:
            if cr_y > 0:
                r_y = aromatic_space
            else:
                r_y = -aromatic_space
                cr_y = -cr_y

            ar_x = aromatic_space * cr_x / cr_y
            br_x = mr_x - aromatic_space * (mr_x - cr_x) / cr_y

            # backward reorienting
            an_x, an_y = rotate_vector(ar_x, r_y, mn_x, mn_y)
            bn_x, bn_y = rotate_vector(br_x, r_y, mn_x, mn_y)
            a_x, a_y = n_x + an_x, n_y + an_y
            b_x, b_y = n_x + bn_x, n_y + bn_y

            return a_x, a_y, b_x, b_y


class X3domCGR(X3dom):
    __slots__ = ()

    def _render_3d_bonds(self, xyz):
        config = self._render_config
        return ''


__all__ = ['X3domMolecule', 'X3domCGR']
