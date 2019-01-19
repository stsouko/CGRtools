# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from collections import defaultdict
from math import atan2, sin, cos
from ..cache import cached_method
from ..periodictable import cpk


class Depict:
    def depict(self, carbon=False, colors=None, font=.4, double_space=.04, triple_space=.07, aromatic_space=.08,
               dashes=(.2, .1), embedding=False):
        if colors is None:
            colors = cpk
        min_x = min(x.x for x in self._node.values())
        max_x = max(x.x for x in self._node.values())
        min_y = min(x.y for x in self._node.values())
        max_y = max(x.y for x in self._node.values())

        svg = []
        radius = {}
        sup_font = .75 * font
        up_font = -.5 * font
        for n, atom in self.atoms():
            tmp, radius[n] = self._render_atom(atom, colors[atom.element], font, sup_font, up_font,
                                               carbon or not bool(self._adj[n]))
            svg.extend(tmp)
        if svg:
            svg.insert(0, '  <g font-family="sans-serif">')
            svg.append('  </g>')

        svg.append('  <g fill="none" stroke="black" stroke-width=".03">')
        for n, m, bond in self.bonds():
            na, ma = self._node[n], self._node[m]
            nx, ny, mx, my = na.x, na.y, ma.x, ma.y

            # indent bond from atom
            if radius[n]:
                dx, dy = self._rotate_vector(radius[n], 0, mx, my, nx, ny)
                nx += dx
                ny -= dy
                rn = True
            else:
                rn = False
            if radius[m]:
                dx, dy = self._rotate_vector(radius[m], 0, mx, my, nx, ny)
                mx -= dx
                my += dy
                rm = True
            else:
                rm = False
            svg.extend(self._render_bond(bond, nx, ny, mx, my, rn, rm,
                                         double_space, triple_space, aromatic_space, dashes))
        svg.append('  </g>')
        if not embedding:
            width = max_x - min_x + 2.5 * font
            height = max_y - min_y + 2.5 * font
            svg.insert(0, f'<svg width="{width:.2f}cm" height="{height:.2f}cm" '
                       f'viewBox="{min_x - 1.25 * font:.2f} {-max_y - 1.25 * font:.2f} {width:.2f} {height:.2f}" '
                       f'xmlns="http://www.w3.org/2000/svg" version="1.1">')
            svg.append('</svg>')
        if embedding:
            return '\n'.join(svg), max_x, max_y
        return '\n'.join(svg)

    @cached_method
    def _repr_svg_(self):
        return self.depict()

    @staticmethod
    def _rotate_vector(x, y, x2, y2, x1, y1):
        """
        rotate x,y vector over x2-x1, y2-y1 angle
        """
        angle = atan2(y2 - y1, x2 - x1)
        cos_rad = cos(angle)
        sin_rad = sin(angle)
        return cos_rad * x + sin_rad * y, -sin_rad * x + cos_rad * y


class DepictMolecule(Depict):
    @staticmethod
    def _render_atom(atom, color, font, sup_font, up_font, carbon):
        svg = []
        if atom.element != 'C' or carbon or atom.charge or atom.multiplicity or atom.isotope != atom.common_isotope:
            x_shift = -shifts[atom.element] * font
            y_shift = .35 * font
            radius = -1.5 * x_shift
            svg.append(f'    <g fill="{color}">')
            svg.append(f'      <text x="{atom.x + x_shift:.2f}" y="{y_shift - atom.y:.2f}" font-size="{font:.2f}">'
                       f'{atom.element}</text>')
            if atom.charge:
                svg.append(f'      <text x="{atom.x - x_shift:.2f}" y="{-y_shift - atom.y:.2f}" '
                           f'font-size="{sup_font:.2f}">{charge_str[atom.charge]}</text>')
            if atom.multiplicity:
                svg.append(f'      <text x="{atom.x + x_shift:.2f}" y="{up_font - atom.y:.2f}" '
                           f'font-size="{sup_font:.2f}">{multiplicity_str[atom.multiplicity]}</text>')
            if atom.isotope != atom.common_isotope:
                svg.append(f'      <text x="{atom.x - font:.2f}" y="{-y_shift - atom.y:.2f}" '
                           f'font-size="{sup_font:.2f}">{atom.isotope}</text>')
            svg.append('    </g>')
        else:
            radius = 0
        return svg, radius

    def _render_bond(self, bond, nx, ny, mx, my, rn, rm, double_space, triple_space, aromatic_space, dashes):
        if bond.order == 1:
            return [f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" />']
        elif bond.order == 2:
            dx, dy = self._rotate_vector(0, double_space, mx, my, nx, ny)
            return [f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" x2="{mx + dx:.2f}" y2="{-my + dy:.2f}" />',
                    f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" />']
        elif bond.order == 3:
            dx, dy = self._rotate_vector(0, triple_space, mx, my, nx, ny)
            return [f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" x2="{mx + dx:.2f}" y2="{-my + dy:.2f}" />',
                    f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" />',
                    f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" />']
        elif bond.order == 4:
            dx, dy = self._rotate_vector(0, aromatic_space, mx, my, nx, ny)
            return [f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" />',
                    f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" '
                    f'stroke-dasharray="{dashes[0]:.2f} {dashes[1]:.2f}" />']
        return [f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" '
                f'stroke-dasharray="{dashes[0]:.2f} {dashes[1]:.2f}" />']

    def __rings(self):
        rings = defaultdict(list)
        for n, x in enumerate(self.sssr, start=1):
            for y in x:
                rings[y].append(n)
        centers = {}
        for n, x in enumerate(self.sssr, start=1):
            atoms = [self._node[x] for x in x]
            centers[n] = (sum(x.x for x in atoms) / len(x), sum(x.y for x in atoms) / len(x))


class DepictReaction:
    def depict(self, carbon=False, colors=None, font=.4, double_space=.04, triple_space=.07, aromatic_space=.08,
               dashes=(.2, .1)):
        if not self._arrow:
            self.fix_positions()

        svg = ['  <defs><marker id="arrow" markerWidth="10" markerHeight="10" refX="0" refY="3" orient="auto">'
               '<path d="M0,0 L0,6 L9,3 z" /></marker></defs>',
               f'  <line x1="{self._arrow[0]:.2f}" y1="-1" x2="{self._arrow[1]:.2f}" y2="-1" fill="none" '
               'stroke="black" stroke-width=".04" marker-end="url(#arrow)" />']

        r_max_x = r_max_y = 0
        for ml in (self.reagents, self.reactants, self.products):
            for m in ml:
                tmp, max_x, max_y = m.depict(carbon, colors, font, double_space, triple_space,
                                             aromatic_space, dashes, True)
                svg.append(tmp)
                if max_x > r_max_x:
                    r_max_x = max_x
                if max_y > r_max_y:
                    r_max_y = max_y
        width = r_max_x + 2.5 * font
        height = r_max_y + 2.5 * font

        svg.insert(0, f'<svg width="{width:.2f}cm" height="{height:.2f}cm" '
                   f'viewBox="{-1.25 * font:.2f} {-r_max_y - 1.25 * font:.2f} {width:.2f} {height:.2f}" '
                   'xmlns="http://www.w3.org/2000/svg" version="1.1">')
        svg.append('</svg>')
        return '\n'.join(svg)

    @cached_method
    def _repr_svg_(self):
        return self.depict()


shifts = {'H': .35, 'He': .35, 'Li': .35, 'Be': .35, 'B': .35, 'C': .35, 'N': .35, 'O': .35,
          'F': .35, 'Ne': .35, 'Na': .35, 'Mg': .35, 'Al': .35, 'Si': .35, 'P': .35, 'S': .35,
          'Cl': .35, 'Ar': .35, 'K': .35, 'Ca': .35, 'Sc': .35, 'Ti': .35, 'V': .35, 'Cr': .35,
          'Mn': .35, 'Fe': .35, 'Co': .35, 'Ni': .35, 'Cu': .35, 'Zn': .35, 'Ga': .35, 'Ge': .35,
          'As': .35, 'Se': .35, 'Br': .35, 'Kr': .35, 'Rb': .35, 'Sr': .35, 'Y': .35, 'Zr': .35,
          'Nb': .35, 'Mo': .35, 'Tc': .35, 'Ru': .35, 'Rh': .35, 'Pd': .35, 'Ag': .35, 'Cd': .35,
          'In': .35, 'Sn': .35, 'Sb': .35, 'Te': .35, 'I': .35, 'Xe': .35, 'Cs': .35, 'Ba': .35,
          'La': .35, 'Ce': .35, 'Pr': .35, 'Nd': .35, 'Pm': .35, 'Sm': .35, 'Eu': .35, 'Gd': .35,
          'Tb': .35, 'Dy': .35, 'Ho': .35, 'Er': .35, 'Tm': .35, 'Yb': .35, 'Lu': .35, 'Hf': .35,
          'Ta': .35, 'W': .35, 'Re': .35, 'Os': .35, 'Ir': .35, 'Pt': .35, 'Au': .35, 'Hg': .35,
          'Tl': .35, 'Pb': .35, 'Bi': .35, 'Po': .35, 'At': .35, 'Rn': .35, 'Fr': .35, 'Ra': .35,
          'Ac': .35, 'Th': .35, 'Pa': .35, 'U': .35, 'Np': .35, 'Pu': .35, 'Am': .35, 'Cm': .35,
          'Bk': .35, 'Cf': .35, 'Es': .35, 'Fm': .35, 'Md': .35, 'No': .35, 'Lr': .35, 'Rf': .35,
          'Db': .35, 'Sg': .35, 'Bh': .35, 'Hs': .35, 'Mt': .35, 'Ds': .35, 'Rg': .35, 'Cn': .35,
          'Nh': .35, 'Fl': .35, 'Mc': .35, 'Lv': .35, 'Ts': .35, 'Og': .35}

charge_str = {-3: '3⁃', -2: '2⁃', -1: '⁃', 1: '+', 2: '2+', 3: '3+'}
multiplicity_str = {1: '↑↓', 2: '↑', 3: '↑↑'}

__all__ = ['DepictMolecule', 'DepictReaction']
