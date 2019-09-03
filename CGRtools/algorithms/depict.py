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
from CachedMethods import cached_method
from collections import defaultdict
from functools import partial
from math import atan2, sin, cos, hypot
from uuid import uuid4
from ..periodictable.cpk import cpk


def rotate_vector(x1, y1, x2, y2):
    """
    rotate x,y vector over x2-x1, y2-y1 angle
    """
    angle = atan2(y2, x2)
    cos_rad = cos(angle)
    sin_rad = sin(angle)
    return cos_rad * x1 + sin_rad * y1, -sin_rad * x1 + cos_rad * y1


class Depict:
    __slots__ = ()

    def depict(self, *, embedding=False):
        min_x = min(x for x, _ in self._plane.values())
        max_x = max(x for x, _ in self._plane.values())
        min_y = min(y for _, y in self._plane.values())
        max_y = max(y for _, y in self._plane.values())

        bonds = self._render_bonds()
        atoms, masks = self._render_atoms()

        if embedding:
            return atoms, bonds, masks, max_x, max_y

        font = self._render_config['font']
        font125 = 1.25 * font
        width = max_x - min_x + 3.0 * font
        height = max_y - min_y + 2.5 * font
        viewbox_x = min_x - font125
        viewbox_y = -max_y - font125

        svg = [f'<svg width="{width:.2f}cm" height="{height:.2f}cm" '
               f'viewBox="{viewbox_x:.2f} {viewbox_y:.2f} {width:.2f} '
               f'{height:.2f}" xmlns="http://www.w3.org/2000/svg" version="1.1">']

        if bonds:
            if masks:
                uid = str(uuid4())
                svg.append(f'  <defs>\n    <mask id="mask-{uid}">\n'
                           f'      <rect x="{viewbox_x:.2f}" y="{viewbox_y:.2f}" '
                           f'width="{width:.2f}" height="{height:.2f}" fill="white"/>\n      <g fill="black">')
                svg.extend(masks)
                svg.append(f'      </g>\n    </mask>\n  </defs>\n  <g fill="none" '
                           f'stroke="{self._render_config["bond_color"]}" '
                           f'stroke-width="{self._render_config["bond_width"]:.2f}"  mask="url(#mask-{uid})">')
                if len(bonds) == 1:  # SVG BUG adhoc
                    svg.append(f'    <line x1="{viewbox_x:.2f}" y1="{viewbox_y:.2f}" '
                               f'x2="{viewbox_x + width:.2f}" y2="{viewbox_y:.2f}" stroke="none"/>')
            else:
                svg.append(f'  <g fill="none" stroke="{self._render_config["bond_color"]}" '
                           f'stroke-width="{self._render_config["bond_width"]:.2f}">')
            svg.extend(bonds)
            svg.append('  </g>')

        if atoms:
            svg.append('  <g font-family="monospace">')
            svg.extend(atoms)
            svg.append('  </g>')

        svg.append('</svg>')
        return '\n'.join(svg)

    @cached_method
    def _repr_svg_(self):
        return self.depict()

    _render_config = {'carbon': False, 'atoms_colors': cpk, 'bond_color': 'black', 'font': .25, 'dashes': (.2, .1),
                      'aromatic_space': .08, 'triple_space': .07, 'double_space': .04, 'mapping': True,
                      'mapping_color': '#788CFF', 'bond_width': .03, 'query_color': '#5D8AA8'}


class DepictMolecule(Depict):
    __slots__ = ()

    def _render_bonds(self):
        svg = []
        plane = self._plane
        double_space = self._render_config['double_space']
        triple_space = self._render_config['triple_space']
        dash1, dash2 = self._render_config['dashes']
        for n, m, bond in self.bonds():
            order = bond.order
            nx, ny = plane[n]
            mx, my = plane[m]
            ny, my = -ny, -my
            if order in (1, 4):
                svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
            elif order == 2:
                dx, dy = rotate_vector(0, double_space, mx - nx, ny - my)
                svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" x2="{mx + dx:.2f}" y2="{my + dy:.2f}"/>')
                svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" x2="{mx - dx:.2f}" y2="{my - dy:.2f}"/>')
            elif order == 3:
                dx, dy = rotate_vector(0, triple_space, mx - nx, ny - my)
                svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" x2="{mx + dx:.2f}" y2="{my + dy:.2f}"/>')
                svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" x2="{mx - dx:.2f}" y2="{my - dy:.2f}"/>')
            else:
                svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                           f'stroke-dasharray="{dash1:.2f} {dash2:.2f}"/>')

        for ring in self.aromatic_rings:
            cx = sum(plane[n][0] for n in ring) / len(ring)
            cy = sum(plane[n][1] for n in ring) / len(ring)

            for n, m in zip(ring, ring[1:]):
                nx, ny = plane[n]
                mx, my = plane[m]
                aromatic = self.__render_aromatic_bond(nx, ny, mx, my, cx, cy)
                if aromatic:
                    svg.append(aromatic)

            nx, ny = plane[ring[-1]]
            mx, my = plane[ring[0]]
            aromatic = self.__render_aromatic_bond(nx, ny, mx, my, cx, cy)
            if aromatic:
                svg.append(aromatic)
        return svg

    def __render_aromatic_bond(self, n_x, n_y, m_x, m_y, c_x, c_y):
        aromatic_space = self._render_config['aromatic_space']
        dash1, dash2 = self._render_config['dashes']
        # n aligned xy
        mn_x, mn_y, cn_x, cn_y = m_x - n_x, m_y - n_y, c_x - n_x, c_y - n_y

        # nm reoriented xy
        mr_x, mr_y = hypot(mn_x, mn_y), 0
        cr_x, cr_y = rotate_vector(cn_x, cn_y, mn_x, mn_y)

        if aromatic_space / cr_y < .65:
            if cr_y > 0:
                ar_y = br_y = aromatic_space
            else:
                ar_y = br_y = -aromatic_space

            ar_x = aromatic_space * cr_x / abs(cr_y)
            br_x = ((abs(cr_y) - aromatic_space) * mr_x + aromatic_space * cr_x) / abs(cr_y)

            # backward reorienting
            an_x, an_y = rotate_vector(ar_x, ar_y, mn_x, -mn_y)
            bn_x, bn_y = rotate_vector(br_x, br_y, mn_x, -mn_y)
            a_x, a_y = an_x + n_x, an_y + n_y
            b_x, b_y = bn_x + n_x, bn_y + n_y

            return f'    <line x1="{a_x:.2f}" y1="{-a_y:.2f}" x2="{b_x:.2f}" y2="{-b_y:.2f}" ' \
                f'stroke-dasharray="{dash1:.2f} {dash2:.2f}"/>'

    def _render_atoms(self):
        bonds = self._bonds
        plane = self._plane
        hydrogens = self._hydrogens
        charges = self._charges
        radicals = self._radicals

        mapping = self._render_config['mapping']
        carbon = self._render_config['carbon']
        atoms_colors = self._render_config['atoms_colors']
        font = self._render_config['font']
        font2 = .2 * font
        font3 = .3 * font
        font5 = .5 * font
        font6 = .6 * font
        font8 = .8 * font

        svg = []
        maps = []
        mask = []
        for n, atom in self._atoms.items():
            x, y = plane[n]
            y = -y
            symbol = atom.atomic_symbol
            if not bonds[n] or symbol != 'C' or carbon or atom.charge or atom.is_radical or atom.isotope:
                h = hydrogens[n]
                if h == 1:
                    h = 'H'
                elif h:
                    h = f'H<tspan  dy="{font5:.2f}" font-size="{font8:.2f}">{h}</tspan>'
                else:
                    h = ''
                if mapping:
                    maps.append(f'      <text x="{x:.2f}" y="{y:.2f}" dx="-{font3:.2f}" dy="{font8:.2f}" '
                                f'text-anchor="end">{n}</text>')

                svg.append(f'    <g fill="{atoms_colors[atom.atomic_number - 1]}">')
                svg.append(f'      <text x="{x:.2f}" y="{y:.2f}" dx="-{font3:.2f}" dy="{font3:.2f}" '
                           f'font-size="{font:.2f}">{symbol}{h}</text>')
                if charges[n]:
                    svg.append(f'      <text x="{x:.2f}" y="{y:.2f}" dx="{font2:.2f}" dy="-{font5:.2f}" '
                               f'font-size="{font6:.2f}">{_render_charge[charges[n]]}{"↑" if radicals[n] else ""}'
                               f'</text>')
                elif radicals[n]:
                    svg.append(f'      <text x="{x:.2f}" y="{y:.2f}" dx="{font2:.2f}" dy="-{font5:.2f}" '
                               f'font-size="{font6:.2f}">↑</text>')
                if atom.isotope:
                    svg.append(f'      <text x="{x:.2f}" y="{y:.2f}" dx="-{font3:.2f}" dy="-{font2:.2f}" '
                               f'font-size="{font6:.2f}" text-anchor="end">{atom.isotope}</text>')
                svg.append('    </g>')
                mask.append(f'        <circle cx="{x:.2f}" cy="{y:.2f}" r="{font5:.2f}"/>')
            elif mapping:
                maps.append(f'      <text x="{x:.2f}" y="{y:.2f}" text-anchor="middle">{n}</text>')
        if mapping:
            svg.append(f'    <g fill="{self._render_config["mapping_color"]}" font-size="{font8:.2f}">')
            svg.extend(maps)
            svg.append('    </g>')
        return svg, mask


class DepictReaction:
    __slots__ = ()

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

        config = Depict._render_config
        font = config['font']
        font125 = 1.25 * font
        width = r_max_x + 3.0 * font
        height = r_max_y + 2.5 * font
        viewbox_x = -font125
        viewbox_y = -r_max_y - font125

        svg = [f'<svg width="{width:.2f}cm" height="{height:.2f}cm" '
               f'viewBox="{viewbox_x:.2f} {viewbox_y:.2f} {width:.2f} '
               f'{height:.2f}" xmlns="http://www.w3.org/2000/svg" version="1.1">']

        if r_bonds:
            svg.append('  <defs>\n    <marker id="arrow" markerWidth="10" markerHeight="10" '
                       'refX="0" refY="3" orient="auto">\n      <path d="M0,0 L0,6 L9,3"/>\n    </marker>')
            if r_masks:
                uid = str(uuid4())
                svg.append(f'    <mask id="mask-{uid}">\n'
                           f'      <rect x="{viewbox_x:.2f}" y="{viewbox_y:.2f}" '
                           f'width="{width:.2f}" height="{height:.2f}" fill="white"/>\n      <g fill="black">')
                svg.extend(r_masks)
                svg.append('      </g>\n    </mask>\n  </defs>\n'
                           f'  <line x1="{self._arrow[0]:.2f}" y1="-1" x2="{self._arrow[1]:.2f}" y2="-1" fill="none" '
                           'stroke="black" stroke-width=".04" marker-end="url(#arrow)"/>\n'
                           f'  <g fill="none" stroke="{config["bond_color"]}" '
                           f'stroke-width="{config["bond_width"]:.2f}" mask="url(#mask-{uid})">')
                if len(r_bonds) == 1:  # SVG BUG adhoc
                    svg.append(f'    <line x1="{viewbox_x:.2f}" y1="{viewbox_y:.2f}" '
                               f'x2="{viewbox_x + width:.2f}" y2="{viewbox_y:.2f}" stroke="none"/>')
            else:
                svg.append('  </defs>')
                svg.append(f'  <line x1="{self._arrow[0]:.2f}" y1="-1" x2="{self._arrow[1]:.2f}" y2="-1" fill="none" '
                           'stroke="black" stroke-width=".04" marker-end="url(#arrow)"/>')
                svg.append(f'  <g fill="none" stroke="{config["bonds_color"]}" '
                           f'stroke-width="{config["bond_width"]:.2f}">')
            svg.extend(r_bonds)
            svg.append('  </g>')

        if r_atoms:
            svg.append('  <g font-family="monospace">')
            svg.extend(r_atoms)
            svg.append('  </g>')

        svg.append('</svg>')
        return '\n'.join(svg)

    @cached_method
    def _repr_svg_(self):
        return self.depict()


class DepictCGR(Depict):
    def _render_bonds(self):
        svg = []
        plane = self._plane
        config = self._render_config
        double_space = config['double_space']
        triple_space = config['triple_space']
        dash1, dash2 = config['dashes']
        broken = config['broken_color']
        formed = config['formed_color']

        ar_bond_colors = defaultdict(dict)
        for n, m, bond in self.bonds():
            order, p_order = bond.order, bond.p_order
            nx, ny = plane[n]
            mx, my = plane[m]
            ny, my = -ny, -my
            rv = partial(rotate_vector, 0, x2=mx - nx, y2=ny - my)
            if order == 1:
                if p_order == 1:
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my + dy:.2f}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" stroke="{formed}"/>')
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my + dy:.2f}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"/>')
                elif p_order is None:
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" stroke="{broken}"/>')
                else:
                    dx, dy = rv(double_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} '
                               f'{dash2:.2f}" stroke="{formed}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
            elif order == 4:
                if p_order == 4:
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 1:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 2:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    dx, dy = rv(double_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my + dy:.2f}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"/>')
                elif p_order == 3:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    dx, dy = rv(triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my + dy:.2f}"  stroke="{formed}"/>')
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"/>')
                elif p_order is None:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" stroke="{broken}"/>')
                else:
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" stroke="{broken}"/>')
            elif order == 2:
                if p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my + dy:.2f}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(double_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my + dy:.2f}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    dx, dy = rv(double_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my + dy:.2f}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my + dy:.2f}"/>')
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"/>')
                elif p_order is None:
                    dx, dy = rv(double_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                else:
                    dx, dy = rv(triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my + dy:.2f}" stroke-dasharray="{dash1:.2f} '
                               f'{dash2:.2f}" stroke="{formed}"/>')
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" stroke="{broken}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
            elif order == 3:
                if p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my + dy:.2f}"/>')
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my + dy:.2f}"/>')
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" stroke="{broken}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" 'f'stroke="{broken}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    dx, dy = rv(triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" x2="{mx + dx:.2f}" '
                               f'y2="{my + dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" x2="{mx - dx:.2f}" '
                               f'y2="{my - dy:.2f}" stroke="{broken}"/>')
                elif p_order == 2:
                    dx, dy = rv(triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my + dy:.2f}"/>')
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                elif p_order is None:
                    dx, dy = rv(triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" '
                               f'x2="{mx:.2f}" y2="{my:.2f}" stroke="{broken}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                else:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'    <line x1="{nx + dx3:.2f}" y1="{ny + dy3:.2f}" x2="{mx + dx3:.2f}" '
                               f'y2="{my + dy3:.2f}" stroke-dasharray="{dash1:.2f} '
                               f'{dash2:.2f}" stroke="{formed}"/>')
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'    <line x1="{nx - dx3:.2f}" y1="{ny - dy3:.2f}" x2="{mx - dx3:.2f}" '
                               f'y2="{my - dy3:.2f}" stroke="{broken}"/>')
            elif order is None:
                if p_order == 1:
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" stroke="{formed}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" stroke="{formed}"/>')
                elif p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" x2="{mx + dx:.2f}" '
                               f'y2="{my + dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" x2="{mx - dx:.2f}" '
                               f'y2="{my - dy:.2f}" stroke="{formed}"/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" '
                               f'x2="{mx:.2f}" y2="{my:.2f}" stroke="{formed}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"/>')
                else:
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"/>')
            else:
                if p_order == 5:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(double_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" x2="{mx + dx:.2f}" y2="{my + dy:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"/>')
                elif p_order == 4:
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" stroke="{formed}"/>')
                elif p_order == 2:
                    dx, dy = rv(triple_space)
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" x2="{mx + dx:.2f}" y2="{my + dy:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" stroke="{formed}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"/>')
                elif p_order == 3:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'    <line x1="{nx + dx3:.2f}" y1="{ny + dy3:.2f}" '
                               f'x2="{mx + dx3:.2f}" y2="{my + dy3:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')
                    svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'    <line x1="{nx - dx3:.2f}" y1="{ny - dy3:.2f}" '
                               f'x2="{mx - dx3:.2f}" y2="{my - dy3:.2f}" stroke="{formed}"/>')
                else:
                    svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')

        for ring in self.aromatic_rings:
            cx = sum(plane[x][0] for x in ring) / len(ring)
            cy = sum(plane[x][1] for x in ring) / len(ring)

            for n, m in zip(ring, ring[1:]):
                nx, ny = plane[n]
                mx, my = plane[m]
                aromatic = self.__render_aromatic_bond(nx, ny, mx, my, cx, cy, ar_bond_colors[n].get(m))
                if aromatic:
                    svg.append(aromatic)

            n, m = ring[-1], ring[0]
            nx, ny = plane[n]
            mx, my = plane[m]
            aromatic = self.__render_aromatic_bond(nx, ny, mx, my, cx, cy, ar_bond_colors[n].get(m))
            if aromatic:
                svg.append(aromatic)
        return svg

    def __render_aromatic_bond(self, n_x, n_y, m_x, m_y, c_x, c_y, color):
        aromatic_space = self._render_config['aromatic_space']
        dash1, dash2 = self._render_config['dashes']
        # n aligned xy
        mn_x, mn_y, cn_x, cn_y = m_x - n_x, m_y - n_y, c_x - n_x, c_y - n_y

        # nm reoriented xy
        mr_x, mr_y = hypot(mn_x, mn_y), 0
        cr_x, cr_y = rotate_vector(cn_x, cn_y, mn_x, mn_y)

        if aromatic_space / cr_y < .65:
            if cr_y > 0:
                ar_y = br_y = aromatic_space
            else:
                ar_y = br_y = -aromatic_space

            ar_x = aromatic_space * cr_x / abs(cr_y)
            br_x = ((abs(cr_y) - aromatic_space) * mr_x + aromatic_space * cr_x) / abs(cr_y)

            # backward reorienting
            an_x, an_y = rotate_vector(ar_x, ar_y, mn_x, -mn_y)
            bn_x, bn_y = rotate_vector(br_x, br_y, mn_x, -mn_y)
            if color:
                return f'    <line x1="{an_x + n_x:.2f}" y1="{-an_y - n_y:.2f}"' \
                    f' x2="{bn_x + n_x:.2f}" y2="{-bn_y - n_y:.2f}" ' \
                    f'stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{color}"/>'
            return f'    <line x1="{an_x + n_x:.2f}" y1="{-an_y - n_y:.2f}"' \
                f' x2="{bn_x + n_x:.2f}" y2="{-bn_y - n_y:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}"/>'

    def _render_atoms(self):
        bonds = self._bonds
        plane = self._plane
        charges = self._charges
        p_charges = self._p_charges
        radicals = self._radicals
        p_radicals = self._p_radicals

        carbon = self._render_config['carbon']
        atoms_colors = self._render_config['atoms_colors']
        font = self._render_config['font']
        font2 = .2 * font
        font3 = .3 * font
        font5 = .5 * font
        font6 = .6 * font

        svg = []
        mask = []
        for n, atom in self._atoms.items():
            x, y = plane[n]
            y = -y
            symbol = atom.atomic_symbol
            if not bonds[n] or symbol != 'C' or carbon or atom.charge or atom.is_radical or atom.isotope:
                svg.append(f'    <g fill="{atoms_colors[atom.atomic_number - 1]}">')
                svg.append(f'      <text x="{x:.2f}" y="{y:.2f}" dx="-{font3:.2f}" dy="{font3:.2f}" '
                           f'font-size="{font:.2f}">{symbol}</text>')

                if radicals[n]:
                    r = '↑' if p_radicals[n] else '↑↓'
                elif p_radicals[n]:
                    r = '↓↑'
                else:
                    r = ''

                if charges[n] != p_charges[n]:
                    svg.append(f'      <text x="{x:.2f}" y="{y:.2f}" dx="{font2:.2f}" dy="-{font5:.2f}" '
                               f'font-size="{font6:.2f}">{_render_p_charge[charges[n]][p_charges[n]]}{r}</text>')
                if charges[n]:
                    svg.append(f'      <text x="{x:.2f}" y="{y:.2f}" dx="{font2:.2f}" dy="-{font5:.2f}" '
                               f'font-size="{font6:.2f}">{_render_charge[charges[n]]}{r}</text>')
                elif r:
                    svg.append(f'      <text x="{x:.2f}" y="{y:.2f}" dx="{font2:.2f}" dy="-{font5:.2f}" '
                               f'font-size="{font6:.2f}">{r}</text>')

                if atom.isotope:
                    svg.append(f'      <text x="{x:.2f}" y="{y:.2f}" dx="-{font3:.2f}" dy="-{font2:.2f}" '
                               f'font-size="{font6:.2f}" text-anchor="end">{atom.isotope}</text>')
                svg.append('    </g>')
                mask.append(f'        <circle cx="{x:.2f}" cy="{y:.2f}" r="{font5:.2f}"/>')
        return svg, mask

    _render_config = {'broken_color': 'red', 'formed_color': 'green', **Depict._render_config, 'aromatic_space': .14}


class DepictQuery(Depict):

    def _render_bonds(self):
        svg = []
        plane = self._plane
        double_space = self._render_config['double_space']
        triple_space = self._render_config['triple_space']
        dash1, dash2 = self._render_config['dashes']
        for n, m, bond in self.bonds():
            order = bond.order
            nx, ny = plane[n]
            mx, my = plane[m]
            ny, my = -ny, -my
            if order == 1:
                svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
            elif order == 4:
                dx, dy = rotate_vector(0, double_space, mx - nx, ny - my)
                svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" x2="{mx + dx:.2f}" y2="{my + dy:.2f}"/>')
                svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" x2="{mx - dx:.2f}" y2="{my - dy:.2f}" '
                           f'stroke-dasharray="{dash1:.2f} {dash2:.2f}"/>')
            elif order == 2:
                dx, dy = rotate_vector(0, double_space, mx - nx, ny - my)
                svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" x2="{mx + dx:.2f}" y2="{my + dy:.2f}"/>')
                svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" x2="{mx - dx:.2f}" y2="{my - dy:.2f}"/>')
            elif order == 3:
                dx, dy = rotate_vector(0, triple_space, mx - nx, ny - my)
                svg.append(f'    <line x1="{nx + dx:.2f}" y1="{ny + dy:.2f}" x2="{mx + dx:.2f}" y2="{my + dy:.2f}"/>')
                svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                svg.append(f'    <line x1="{nx - dx:.2f}" y1="{ny - dy:.2f}" x2="{mx - dx:.2f}" y2="{my - dy:.2f}"/>')
            else:
                svg.append(f'    <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                           f'stroke-dasharray="{dash1:.2f} {dash2:.2f}"/>')
        return svg

    def _render_atoms(self):
        plane = self._plane
        mapping = self._render_config['mapping']
        carbon = self._render_config['carbon']
        atoms_colors = self._render_config['atoms_colors']
        font = self._render_config['font']
        font2 = .2 * font
        font3 = .3 * font
        font4 = .4 * font
        font6 = .6 * font
        font7 = .7 * font
        font8 = .8 * font

        svg = []
        mask = []
        maps = []
        nghbrs = []
        hbrdztns = []
        for n, atom in self._atoms.items():
            x, y = plane[n]
            y = -y
            single = not self._bonds[n]
            symbol = atom.atomic_symbol
            if single or symbol != 'C' or carbon or atom.charge or atom.is_radical:
                svg.append(f'    <g fill="{atoms_colors[atom.atomic_number - 1]}">')
                svg.append(f'      <text x="{x - font4:.2f}" y="{font3 + y:.2f}" '
                           f'font-size="{font:.2f}">{symbol}</text>')
                if mapping:
                    maps.append(f'      <text x="{x:.2f}" y="{y:.2f}" dx="0" dy="{font8:.2f}" '
                                f'text-anchor="end">{n}</text>')

                if atom.charge:
                    svg.append(f'      <text x="{x:.2f}" y="{y:.2f}" dx="{font2:.2f}" dy="-{font4:.2f}" '
                               f'font-size="{font6:.2f}">{_render_charge[atom.charge]}'
                               f'{"↑" if atom.is_radical else ""}</text>')
                elif atom.is_radical:
                    svg.append(f'      <text x="{x:.2f}" y="{y:.2f}" dx="{font2:.2f}" dy="-{font4:.2f}" '
                               f'font-size="{font6:.2f}">↑</text>')
                svg.append('    </g>')
                mask.append(f'      <circle cx="{x:.2f}" cy="{y:.2f}" r="{font7:.2f}"/>')
            elif mapping:
                maps.append(f'      <text x="{x:.2f}" y="{y:.2f}" dx="0" dy="{font8:.2f}" '
                            f'text-anchor="end">{n}</text>')

            level = 1
            if atom.neighbors:
                level = 1.6
                nn = [str(x) for x in atom.neighbors]
                nghbrs.append(f'      <text x="{x:.2f}" y="{y:.2f}" dx="{0:.2f}" dy="{font8:.2f}" '
                              f'text-anchor="start">{"".join(nn)}</text>')

            if atom.hybridization:
                hh = [_render_hybridization[x] for x in atom.hybridization]
                hbrdztns.append(f'      <text x="{x:.2f}" y="{y:.2f}" dx="{0:.2f}" dy="{level * font8:.2f}" '
                                f'text-anchor="start">{"".join(hh)}</text>')
        hybs = True
        if nghbrs:
            svg.append(f'    <g fill="{self._render_config["query_color"]}" font-size="{font7:.2f}">')
            svg.extend(nghbrs)
            if hbrdztns:
                svg.extend(hbrdztns)
                hybs = False
            svg.append('    </g>')

        if hbrdztns and hybs:
            svg.append(f'    <g fill="{self._render_config["query_color"]}" font-size="{font7:.2f}">')
            svg.extend(hbrdztns)
            svg.append('    </g>')

        if mapping:
            svg.append(f'    <g fill="{self._render_config["mapping_color"]}" font-size="{font7:.2f}">')
            svg.extend(maps)
            svg.append('    </g>')

        return svg, mask


_render_hybridization = {1: 's', 2: 'd', 3: 't', 4: 'a'}
_render_charge = {-3: '3⁃', -2: '2⁃', -1: '⁃', 1: '+', 2: '2+', 3: '3+'}
_render_p_charge = {-3: {-2: '-3»-2', -1: '-3»-', 0: '-3»0', 1: '-3»+', 2: '-3»2', 3: '-3»3'},
                    -2: {-3: '-2»-3', -1: '-2»-', 0: '-2»0', 1: '-2»+', 2: '-2»2', 3: '-2»3'},
                    -1: {-3: '-»-3', -2: '-»-2', 0: '-»0', 1: '-»+', 2: '-»2', 3: '-»3'},
                    0: {-3: '0»-3', -2: '0»-2', -1: '0»-', 1: '0»+', 2: '0»2', 3: '0»3'},
                    1: {-3: '+»-3', -2: '+»-2', -1: '+»-', 0: '+»0', 2: '+»2', 3: '+»3'},
                    2: {-3: '2»-3', -2: '2»-2', -1: '2»-', 0: '2»0', 1: '2»+', 3: '2»3'},
                    3: {-3: '3»-3', -2: '3»-2', -1: '3»-', 0: '3»0', 1: '3»+', 2: '3»2'}}


__all__ = ['DepictMolecule', 'DepictReaction', 'DepictCGR', 'DepictQuery']
