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


class X3dom:
    __slots__ = ()

    def depict3d(self, index: int = 0) -> str:
        """Get X3DOM XML string.
        :param index: index of conformer
        """
        return f'<x3d width=100% height=100%>\n  <scene>\n{self.__render_atoms(index)}{self._render_3d_bonds(index)}' \
               '  </scene>\n</x3d>'

    def __render_atoms(self, index):
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

        xyz = self._conformers[index]
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

    def _render_3d_bonds(self, index):
        config = self._render_config
        xyz = self._conformers[index]
        return ''


class X3domCGR(X3dom):
    __slots__ = ()

    def _render_3d_bonds(self, index):
        config = self._render_config
        xyz = self._conformers[index]
        return ''


__all__ = ['X3domMolecule', 'X3domCGR']
