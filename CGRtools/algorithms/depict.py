# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
#  Copyright 2019 Dinar Batyrshin <batyrshin-dinar@mail.ru>
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
from math import atan2, sin, cos, hypot
from uuid import uuid4
from ..cache import cached_method
from ..periodictable import cpk


def rotate_vector(x1, y1, x2, y2):
    """
    rotate x,y vector over x2-x1, y2-y1 angle
    """
    angle = atan2(y2, x2)
    cos_rad = cos(angle)
    sin_rad = sin(angle)
    return cos_rad * x1 + sin_rad * y1, -sin_rad * x1 + cos_rad * y1


class Depict:
    def depict(self, *, embedding=False):
        min_x = min(x.x for x in self._atoms.values())
        max_x = max(x.x for x in self._atoms.values())
        min_y = min(x.y for x in self._atoms.values())
        max_y = max(x.y for x in self._atoms.values())

        bonds = self._render_bonds()
        atoms, masks = self._render_atoms()

        if embedding:
            return atoms, bonds, masks, max_x, max_y

        width = max_x - min_x + 2.5 * self._render_font
        height = max_y - min_y + 2.5 * self._render_font
        viewbox_x = min_x - 1.25 * self._render_font
        viewbox_y = -max_y - 1.25 * self._render_font

        svg = [f'<svg width="{width:.2f}cm" height="{height:.2f}cm" '
               f'viewBox="{viewbox_x:.2f} {viewbox_y:.2f} {width:.2f} '
               f'{height:.2f}" xmlns="http://www.w3.org/2000/svg" version="1.1">']

        if bonds:
            if masks:
                uid = str(uuid4())
                svg.append(f'  <defs>\n    <mask id="mask-{uid}">\n'
                           f'      <rect x="{min_x - 1.25 * self._render_font:.2f}" '
                           f'y="{-max_y - 1.25 * self._render_font:.2f}" '
                           f'width="{width:.2f}" height="{height:.2f}" fill="white"/>')
                svg.extend(masks)
                svg.append(f'    </mask>\n  </defs>\n  <g fill="none" stroke="{self._render_bonds_color}" '
                           f'stroke-width=".03"  mask="url(#mask-{uid})">')
            else:
                svg.append(f'  <g fill="none" stroke="{self._render_bonds_color}" stroke-width=".03">')
            svg.append(f'    <line x1="{viewbox_x:.2f}" y1="{viewbox_y:.2f}" '
                       f'x2="{viewbox_x + width:.2f}" y2="{viewbox_y:.2f}" stroke="none"/>')
            svg.extend(bonds)
            svg.append('  </g>')

        if atoms:
            svg.append('  <g font-family="sans-serif">')
            svg.extend(atoms)
            svg.append('  </g>')

        svg.append('</svg>')
        return '\n'.join(svg)

    @cached_method
    def _repr_svg_(self):
        return self.depict()

    _render_carbon = False
    _render_atoms_colors = cpk
    _render_bonds_color = 'black'
    _render_font = .4
    _render_sup_font = .75 * _render_font
    _render_up_font = -.5 * _render_font
    _render_double_space = .04
    _render_triple_space = .07
    _render_aromatic_space = .08
    _render_dashes = (.2, .1)


class DepictMolecule(Depict):
    def _render_bonds(self):
        svg = []
        nodes = self._atoms
        for n, m, bond in self.bonds():
            n, m = nodes[n], nodes[m]
            nx, ny, mx, my = n.x, n.y, m.x, m.y
            if bond.order in (1, 4):
                svg.append(f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}"/>')
            elif bond.order == 2:
                dx, dy = rotate_vector(0, self._render_double_space, mx - nx, my - ny)
                svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" x2="{mx + dx:.2f}" y2="{-my + dy:.2f}"/>')
                svg.append(f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" x2="{mx - dx:.2f}" y2="{-my - dy:.2f}"/>')
            elif bond.order == 3:
                dx, dy = rotate_vector(0, self._render_triple_space, mx - nx, my - ny)
                svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" x2="{mx + dx:.2f}" y2="{-my + dy:.2f}"/>')
                svg.append(f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}"/>')
                svg.append(f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" x2="{mx - dx:.2f}" y2="{-my - dy:.2f}"/>')
            else:
                svg.append(f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" '
                           f'stroke-dasharray="{self._render_dashes[0]:.2f} {self._render_dashes[1]:.2f}"/>')

        for ring in self.aromatic_rings:
            c_x = sum(nodes[x].x for x in ring) / len(ring)
            c_y = sum(nodes[y].y for y in ring) / len(ring)

            for n, m in zip(ring, ring[1:]):
                n, m = nodes[n], nodes[m]
                aromatic = self.__render_aromatic_bond(n.x, n.y, m.x, m.y, c_x, c_y)
                if aromatic:
                    svg.append(aromatic)

            n, m = nodes[ring[-1]], nodes[ring[0]]
            aromatic = self.__render_aromatic_bond(n.x, n.y, m.x, m.y, c_x, c_y)
            if aromatic:
                svg.append(aromatic)
        return svg

    def __render_aromatic_bond(self, n_x, n_y, m_x, m_y, c_x, c_y):
        # n aligned xy
        mn_x, mn_y, cn_x, cn_y = m_x - n_x, m_y - n_y, c_x - n_x, c_y - n_y

        # nm reoriented xy
        mr_x, mr_y = hypot(mn_x, mn_y), 0
        cr_x, cr_y = rotate_vector(cn_x, cn_y, mn_x, mn_y)

        if self._render_aromatic_space / cr_y < .65:
            if cr_y > 0:
                ar_y = br_y = self._render_aromatic_space
            else:
                ar_y = br_y = -self._render_aromatic_space

            ar_x = self._render_aromatic_space * cr_x / abs(cr_y)
            br_x = ((abs(cr_y) - self._render_aromatic_space) * mr_x + self._render_aromatic_space * cr_x) / abs(cr_y)

            # backward reorienting
            an_x, an_y = rotate_vector(ar_x, ar_y, mn_x, -mn_y)
            bn_x, bn_y = rotate_vector(br_x, br_y, mn_x, -mn_y)
            a_x, a_y = an_x + n_x, an_y + n_y
            b_x, b_y = bn_x + n_x, bn_y + n_y

            return f'    <line x1="{a_x:.2f}" y1="{-a_y:.2f}" x2="{b_x:.2f}" y2="{-b_y:.2f}" ' \
                f'stroke-dasharray="{self._render_dashes[0]:.2f} {self._render_dashes[1]:.2f}"/>'

    def _render_atoms(self):
        svg = []
        mask = []
        for n, atom in self.atoms():
            single = not self._bonds[n]
            if single or atom.element != 'C' or self._render_carbon or atom.charge or atom.multiplicity or \
                    atom.isotope != atom.common_isotope:
                x_shift = -shifts[atom.element] * self._render_font
                y_shift = .35 * self._render_font
                radius = -1.5 * x_shift
                svg.append(f'    <g fill="{self._render_atoms_colors[atom.element]}">')
                svg.append(f'      <text x="{atom.x + x_shift:.2f}" y="{y_shift - atom.y:.2f}" '
                           f'font-size="{self._render_font:.2f}">{atom.element}</text>')
                if atom.charge:
                    svg.append(f'      <text x="{atom.x - x_shift:.2f}" y="{-y_shift - atom.y:.2f}" '
                               f'font-size="{self._render_sup_font:.2f}">{_render_charge_str[atom.charge]}</text>')
                if atom.multiplicity:
                    svg.append(f'      <text x="{atom.x + x_shift:.2f}" y="{self._render_up_font - atom.y:.2f}" '
                               f'font-size="{self._render_sup_font:.2f}">'
                               f'{_render_multiplicity_str[atom.multiplicity]}</text>')
                if atom.isotope != atom.common_isotope:
                    svg.append(f'      <text x="{atom.x - self._render_font:.2f}" y="{-y_shift - atom.y:.2f}" '
                               f'font-size="{self._render_sup_font:.2f}">{atom.isotope}</text>')
                svg.append('    </g>')
                mask.append(f'      <circle cx="{atom.x:.2f}" cy="{-atom.y:.2f}" r="{radius:.2f}" fill="black"/>')
        return svg, mask


class DepictReaction:
    def depict(self):
        if not self._arrow:
            self.fix_positions()

        r_atoms = []
        r_bonds = []
        r_masks = []

        r_max_x = r_max_y = 0
        for ml in (self.reactants, self.reagents, self.products):
            for m in ml:
                atoms, bonds, masks, max_x, max_y = m.depict(embedding=True)
                r_atoms.extend(atoms)
                r_bonds.extend(bonds)
                r_masks.extend(masks)
                if max_x > r_max_x:
                    r_max_x = max_x
                if max_y > r_max_y:
                    r_max_y = max_y

        width = r_max_x + 2.5 * self._render_font
        height = r_max_y + 2.5 * self._render_font
        viewbox_x = -1.25 * self._render_font
        viewbox_y = -r_max_y - 1.25 * self._render_font

        svg = [f'<svg width="{width:.2f}cm" height="{height:.2f}cm" '
               f'viewBox="{viewbox_x:.2f} {viewbox_y:.2f} {width:.2f} '
               f'{height:.2f}" xmlns="http://www.w3.org/2000/svg" version="1.1">']

        if r_bonds:
            svg.append('  <defs>\n    <marker id="arrow" markerWidth="10" markerHeight="10" '
                       'refX="0" refY="3" orient="auto">\n      <path d="M0,0 L0,6 L9,3"/>\n    </marker>')
            if r_masks:
                uid = str(uuid4())
                svg.append(f'    <mask id="mask-{uid}">\n'
                           f'      <rect x="{-1.25 * self._render_font:.2f}" '
                           f'y="{-r_max_y - 1.25 * self._render_font:.2f}" '
                           f'width="{width:.2f}" height="{height:.2f}" fill="white"/>')
                svg.extend(r_masks)
                svg.append('    </mask>\n  </defs>\n'
                           f'  <line x1="{self._arrow[0]:.2f}" y1="-1" x2="{self._arrow[1]:.2f}" y2="-1" fill="none" '
                           'stroke="black" stroke-width=".04" marker-end="url(#arrow)"/>\n'
                           f'  <g fill="none" stroke="{self._render_bonds_color}" '
                           f'stroke-width=".03" mask="url(#mask-{uid})">')
            else:
                svg.append('  </defs>')
                svg.append(f'  <line x1="{self._arrow[0]:.2f}" y1="-1" x2="{self._arrow[1]:.2f}" y2="-1" fill="none" '
                           'stroke="black" stroke-width=".04" marker-end="url(#arrow)"/>')
                svg.append(f'  <g fill="none" stroke="{self._render_bonds_color}" stroke-width=".03">')
            svg.append(f'    <line x1="{viewbox_x:.2f}" y1="{viewbox_y:.2f}" '
                       f'x2="{viewbox_x + width:.2f}" y2="{viewbox_y:.2f}" stroke="none"/>')
            svg.extend(r_bonds)
            svg.append('  </g>')

        if r_atoms:
            svg.append('  <g font-family="sans-serif">')
            svg.extend(r_atoms)
            svg.append('  </g>')

        svg.append('</svg>')
        return '\n'.join(svg)

    @cached_method
    def _repr_svg_(self):
        return self.depict()

    _render_font = .4
    _render_bonds_color = 'black'


class DepictCGR(Depict):
    def _render_bonds(self):
        svg = []
        nodes = self._atoms
        _ar_bond_colors = defaultdict(dict)
        for n, m, bond in self.bonds():
            order, p_order = bond.order, bond.p_order
            na, ma = nodes[n], nodes[m]
            nx, ny, mx, my = na.x, na.y, ma.x, ma.y
            rv = lambda s: rotate_vector(0, s, mx - nx, my - ny)
            if order == 1:
                if p_order == 1:
                    svg.append(f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}"/>')
                elif p_order == 4:
                    _ar_bond_colors[n][m] = 'green'
                    _ar_bond_colors[m][n] = 'green'
                    svg.append(f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}"/>')
                elif p_order == 2:
                    dx, dy = rv(self._render_double_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" stroke="green"/>')
                elif p_order == 3:
                    dx, dy = rv(self._render_triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}"/>\n'
                               f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" '
                               f'x2="{mx:.2f}" y2="{-my:.2f}" stroke="green"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" stroke="green"/>')
                elif p_order is None:
                    svg.append(f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" '
                               'stroke="red"/>')
                else:
                    dx, dy = rv(self._render_double_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my - dy:.2f}" stroke-dasharray="{self._render_dashes[0]:.2f} '
                               f'{self._render_dashes[1]:.2f}" stroke="green"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" stroke="red"/>')
            elif order == 4:
                if p_order == 4:
                    _ar_bond_colors[n][m] = 'black'
                    _ar_bond_colors[m][n] = 'black'
                    svg.append(f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}"/>')
                elif p_order == 1:
                    _ar_bond_colors[n][m] = 'red'
                    _ar_bond_colors[m][n] = 'red'
                    svg.append(f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}"/>')
                elif p_order == 2:
                    _ar_bond_colors[n][m] = 'red'
                    _ar_bond_colors[m][n] = 'red'
                    dx, dy = rv(self._render_double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{-ny + dx:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}"/>\n'
                               f'      <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" stroke="green"/>')
                elif p_order == 3:
                    _ar_bond_colors[n][m] = 'red'
                    _ar_bond_colors[m][n] = 'red'
                    dx, dy = rv(self._render_triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}"  stroke="green"/>\n'
                               f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" stroke="green"/>')
                elif p_order is None:
                    _ar_bond_colors[n][m] = 'red'
                    _ar_bond_colors[m][n] = 'red'
                    svg.append(f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" '
                               'stroke="red"/>')
                else:
                    _ar_bond_colors[n][m] = 'black'
                    _ar_bond_colors[m][n] = 'black'
                    svg.append(f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" stroke="red"/>')
            elif order == 2:
                if p_order == 2:
                    dx, dy = rv(self._render_double_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(self._render_double_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" stroke="red"/>')
                elif p_order == 4:
                    _ar_bond_colors[n][m] = 'green'
                    _ar_bond_colors[m][n] = 'green'
                    dx, dy = rv(self._render_double_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" stroke="red"/>')
                elif p_order == 3:
                    dx, dy = rv(self._render_triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}"/>\n'
                               f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" stroke="green"/>')
                elif p_order is None:
                    dx, dy = rv(self._render_double_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}" stroke="red"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" stroke="red"/>')
                elif p_order == 5:
                    dx, dy = rv(self._render_triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}" stroke-dasharray="{self._render_dashes[0]:.2f} '
                               f'{self._render_dashes[1]:.2f}" stroke="green"/>\n'
                               f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" '
                               f'x2="{mx:.2f}" y2="{-my:.2f}" stroke="red"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" stroke="red"/>')
            elif order == 3:
                if p_order == 3:
                    dx, dy = rv(self._render_triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}"/>\n'
                               f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(self._render_triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}"/>\n'
                               f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" stroke="red"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" 'f'stroke="red"/>')
                elif p_order == 4:
                    _ar_bond_colors[n][m] = 'green'
                    _ar_bond_colors[m][n] = 'green'
                    dx, dy = rv(self._render_triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" x2="{mx + dx:.2f}" '
                               f'y2="{-my + dy:.2f}" stroke="red"/>\n'
                               f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" x2="{mx - dx:.2f}" '
                               f'y2="{-my - dy:.2f}" stroke="red"/>')
                elif p_order == 2:
                    dx, dy = rv(self._render_triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}"/>\n'
                               f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" stroke="red"/>')
                elif p_order is None:
                    dx, dy = rv(self._render_triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}" stroke="red"/>\n'
                               f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" '
                               f'x2="{mx:.2f}" y2="{-my:.2f}" stroke="red"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" stroke="red"/>')
                else:
                    dx, dy = rv(self._render_double_space)
                    svg.append(f'    <line x1="{nx + 3 * dx:.2f}" y1="{-ny + 3 * dy:.2f}" x2="{mx + 3 * dx:.2f}" '
                               f'y2="{-my + 3 * dy:.2f}" stroke-dasharray="{self._render_dashes[0]:.2f} '
                               f'{self._render_dashes[1]:.2f}" stroke="green"/>\n'
                               f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}" stroke="red"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" stroke="red"/>\n'
                               f'    <line x1="{nx - 3 * dx:.2f}" y1="{-ny - 3 * dy:.2f}" x2="{mx - 3 * dx:.2f}" '
                               f'y2="{-my - 3 * dy:.2f}" stroke="red"/>')
            elif order is None:
                if p_order == 1:
                    svg.append(f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" '
                               'stroke="green"/>')
                elif p_order == 4:
                    _ar_bond_colors[n][m] = 'green'
                    _ar_bond_colors[m][n] = 'green'
                    svg.append(f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" '
                               'stroke="green"/>')
                elif p_order == 2:
                    dx, dy = rv(self._render_double_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" x2="{mx + dx:.2f}" '
                               f'y2="{-my + dy:.2f}" stroke="green"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" x2="{mx - dx:.2f}" '
                               f'y2="{-my - dy:.2f}" stroke="green"/>')
                elif p_order == 3:
                    dx, dy = rv(self._render_triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}" stroke="green"/>\n'
                               f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" '
                               f'x2="{mx:.2f}" y2="{-my:.2f}" stroke="green"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" stroke="green"/>')
                else:
                    svg.append(f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" '
                               f'stroke-dasharray="{self._render_dashes[0]:.2f} '
                               f'{self._render_dashes[1]:.2f}" stroke="green"/>')
            else:
                if p_order == 5:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" '
                               f'stroke-dasharray="{self._render_dashes[0]:.2f} {self._render_dashes[1]:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(self._render_double_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}" stroke-dasharray="{self._render_dashes[0]:.2f} '
                               f'{self._render_dashes[1]:.2f}"stroke="red"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" stroke="green"/>')
                elif p_order == 4:
                    _ar_bond_colors[n][m] = 'black'
                    svg.append(f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" '
                               'stroke="green"/>')
                elif p_order == 2:
                    dx, dy = rv(self._render_triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}"/>'
                               f'stroke-dasharray="{self._render_dashes[0]:.2f} '
                               f'{self._render_dashes[1]:.2f}" stroke="red"/>\n'
                               f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" '
                               f'x2="{mx:.2f}" y2="{-my:.2f}" stroke="green"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" stroke="green"/>')
                elif p_order == 3:
                    dx, dy = rv(self._render_double_space)
                    svg.append(f'    <line x1="{nx + 3 * dx:.2f}" y1="{-ny + 3 * dy:.2f}" '
                               f'x2="{mx + 3 * dx:.2f}" y2="{-my + 3 * dy:.2f}" '
                               f'stroke-dasharray="{self._render_dashes[0]:.2f} '
                               f'{self._render_dashes[1]:.2f}" stroke="red"/>\n'
                               f'    <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{-my + dy:.2f}" stroke="green"/>\n'
                               f'    <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" stroke="green"/>\n'
                               f'    <line x1="{nx - 3 * dx:.2f}" y1="{-ny - 3 * dy:.2f}" '
                               f'x2="{mx - 3 * dx:.2f}" y2="{-my - 3 * dy:.2f}" stroke="green"/>')
                else:
                    svg.append(f'    <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" '
                               f'stroke-dasharray="{self._render_dashes[0]:.2f} '
                               f'{self._render_dashes[1]:.2f}" stroke="red"/>')

        for ring in self.aromatic_rings:
            c_x = sum(nodes[x].x for x in ring) / len(ring)
            c_y = sum(nodes[y].y for y in ring) / len(ring)

            for n, m in zip(ring, ring[1:]):
                na, ma = nodes[n], nodes[m]
                aromatic = self.__render_aromatic_bond(na.x, na.y, ma.x, ma.y, c_x, c_y)
                if aromatic:
                    a_x, a_y, b_x, b_y = aromatic

                    svg.append(f'    <line x1="{a_x:.2f}" y1="{-a_y:.2f}" x2="{b_x:.2f}" y2="{-b_y:.2f}" '
                               f'stroke-dasharray="{self._render_dashes[0]:.2f} {self._render_dashes[1]:.2f}" '
                               f'stroke="{_ar_bond_colors[n][m]}"/>')

            n, m = ring[-1], ring[0]
            na, ma = nodes[n], nodes[m]
            aromatic = self.__render_aromatic_bond(na.x, na.y, ma.x, ma.y, c_x, c_y)
            if aromatic:
                a_x, a_y, b_x, b_y = aromatic

                svg.append(f'    <line x1="{a_x:.2f}" y1="{-a_y:.2f}" x2="{b_x:.2f}" y2="{-b_y:.2f}" '
                           f'stroke-dasharray="{self._render_dashes[0]:.2f} {self._render_dashes[1]:.2f}" '
                           f'stroke="{_ar_bond_colors[n][m]}"/>')

        return svg

    def __render_aromatic_bond(self, n_x, n_y, m_x, m_y, c_x, c_y):
        # n aligned xy
        mn_x, mn_y, cn_x, cn_y = m_x - n_x, m_y - n_y, c_x - n_x, c_y - n_y

        # nm reoriented xy
        mr_x, mr_y = hypot(mn_x, mn_y), 0
        cr_x, cr_y = rotate_vector(cn_x, cn_y, mn_x, mn_y)

        if self._render_aromatic_space / cr_y < .65:
            if cr_y > 0:
                ar_y = br_y = self._render_aromatic_space
            else:
                ar_y = br_y = -self._render_aromatic_space

            ar_x = self._render_aromatic_space * cr_x / abs(cr_y)
            br_x = ((abs(cr_y) - self._render_aromatic_space) * mr_x + self._render_aromatic_space * cr_x) / abs(cr_y)

            # backward reorienting
            an_x, an_y = rotate_vector(ar_x, ar_y, mn_x, -mn_y)
            bn_x, bn_y = rotate_vector(br_x, br_y, mn_x, -mn_y)

            return an_x + n_x, an_y + n_y, bn_x + n_x, bn_y + n_y

    def _render_atoms(self):
        svg = []
        mask = []
        for n, atom in self.atoms():
            single = not self._bonds[n]
            if single or atom.element != 'C' or self._render_carbon or atom.charge or atom.multiplicity or \
                    atom.isotope != atom.common_isotope:
                x_shift = -shifts[atom.element] * self._render_font
                y_shift = .35 * self._render_font
                radius = -1.5 * x_shift
                svg.append(f'    <g fill="{self._render_atoms_colors[atom.element]}">')
                svg.append(f'      <text x="{atom.x + x_shift:.2f}" y="{y_shift - atom.y:.2f}" '
                           f'font-size="{self._render_font:.2f}" '
                           f'>{atom.element}</text>')
                if atom.charge != atom.p_charge:
                    svg.append(f'      <text x="{atom.x - x_shift:.2f}" y="{-y_shift - atom.y:.2f}" '
                               f'font-size="{self._render_sup_font:.2f}">'
                               f'{_render_p_charge_str[atom.charge][atom.p_charge]}</text>')
                elif atom.charge:
                    svg.append(f'      <text x="{atom.x - x_shift:.2f}" y="{-y_shift - atom.y:.2f}" '
                               f'font-size="{self._render_sup_font:.2f}">{_render_charge_str[atom.charge]}</text>')
                if atom.multiplicity:
                    svg.append(f'      <text x="{atom.x + x_shift:.2f}" y="{self._render_up_font - atom.y:.2f}" '
                               f'font-size="{self._render_sup_font:.2f}">'
                               f'{_render_multiplicity_str[atom.multiplicity]}</text>')
                if atom.isotope != atom.common_isotope:
                    svg.append(f'      <text x="{atom.x - self._render_font:.2f}" y="{-y_shift - atom.y:.2f}" '
                               f'font-size="{self._render_sup_font:.2f}">{atom.isotope}</text>')
                svg.append('    </g>')
                mask.append(f'      <circle cx="{atom.x}" cy="{-atom.y}" r="{radius}" fill="black"/>')
        return svg, mask

    _render_aromatic_space = .14


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

_render_charge_str = {-3: '3⁃', -2: '2⁃', -1: '⁃', 1: '+', 2: '2+', 3: '3+'}
_render_p_charge_str = {-3: {-2: '-3»-2', -1: '-3»-', 0: '-3»0', 1: '-3»+', 2: '-3»2', 3: '-3»3'},
                        -2: {-3: '-2»-3', -1: '-2»-', 0: '-2»0', 1: '-2»+', 2: '-2»2', 3: '-2»3'},
                        -1: {-3: '-»-3', -2: '-»-2', 0: '-»0', 1: '-»+', 2: '-»2', 3: '-»3'},
                        0: {-3: '0»-3', -2: '0»-2', -1: '0»-', 1: '0»+', 2: '0»2', 3: '0»3'},
                        1: {-3: '+»-3', -2: '+»-2', -1: '+»-', 0: '+»0', 2: '+»2', 3: '+»3'},
                        2: {-3: '2»-3', -2: '2»-2', -1: '2»-', 0: '2»0', 1: '2»+', 3: '2»3'},
                        3: {-3: '3»-3', -2: '3»-2', -1: '3»-', 0: '3»0', 1: '3»+', 2: '3»2'}}
_render_multiplicity_str = {1: '↑↓', 2: '↑', 3: '↑↑'}

__all__ = ['DepictMolecule', 'DepictReaction', 'DepictCGR']
