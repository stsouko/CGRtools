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
from math import acos, sqrt


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
        bonds = []
        for n, m, bond in self.bonds():
            nx, ny, nz = xyz[n]
            mx, my, mz = xyz[m]
            _nx, _ny, _nz = 0, 0, 0

            nmx, nmy, nmz = mx - nx, my - ny, mz - nz
            ncx, ncy, ncz = 0, 0, 1
            angle = acos((nmx * ncx + nmy * ncy + nmz * ncz) / sqrt(nmx ** 2 + nmy ** 2 + nmz ** 2))

            bonds.append(f'    <transform translation="{nx:.2f} {ny:.2f} {nz:.2f}" rotation="{nmy:.2f} {- nmx:.2f} {0} '
                         f'{angle:.2f}">\n      <shape\n>        <appearance\n>          <material diffusecolor="0 1 1"'
                         f' ambientintensity="0.2" emissivecolor="0,0,0" shininess="0.2" specularcolor="0,0,0"\n>'
                         f'          </material\n>        </appearance\n>        <cylinder radius=".04" height=".825"\n>'
                         f'        </cylinder\n>      </shape\n>    </transform>\n')
        return ''.join(bonds)


class X3domCGR(X3dom):
    __slots__ = ()

    def _render_3d_bonds(self, xyz):
        config = self._render_config
        return ''


__all__ = ['X3domMolecule', 'X3domCGR']
