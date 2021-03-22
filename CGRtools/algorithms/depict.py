# -*- coding: utf-8 -*-
#
#  Copyright 2018-2021 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019-2020 Dinar Batyrshin <batyrshin-dinar@mail.ru>
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
from CachedMethods import cached_method
from collections import defaultdict
from functools import partial
from math import atan2, sin, cos, hypot
from uuid import uuid4
from typing import Tuple


cpk = tuple('''
 #909090                                                                                         #D9FFFF
 #CC80FF #C2FF00                                         #FFB5B5 #101010 #3050F8 #FF0D0D #90E050 #B3E3F5
 #AB5CF2 #8AFF00                                         #BFA6A6 #F0C8A0 #FF8000 #C6C600 #1FF01F #80D1E3
 #8F40D4 #3DFF00 #E6E6E6 #BFC2C7 #A6A6AB #8A99C7 #9C7AC7 
                 #E06633 #F090A0 #50D050 #C88033 #7D80B0 #C28F8F #668F8F #BD80E3 #FFA100 #A62929 #5CB8D1
 #702EB0 #00FF00 #94FFFF #94E0E0 #73C2C9 #54B5B5 #3B9E9E 
                 #248F8F #0A7D8C #006985 #C0C0C0 #FFD98F #A67573 #668080 #9E63B5 #D47A00 #940094 #429EB0
 #57178F #00C900 #70D4FF 
                 #FFFFC7 #D9FFC7 #C7FFC7 #A3FFC7 #8FFFC7 #61FFC7 #45FFC7
                 #30FFC7 #1FFFC7 #00FF9C #00E675 #00D452 #00BF38 #00AB24 
                         #4DC2FF #4DA6FF #2194D6 #267DAB
                 #266696 #175487 #D0D0E0 #FFD123 #B8B8D0 #A6544D #575961 #9E4FB5 #AB5C00 #754F45 #428296
 #420066 #007D00 #70ABFA
                 #00BAFF #00A1FF #008FFF #0080FF #006BFF #545CF2 #785CE3
                 #8A4FE3 #A136D4 #B31FD4 #B31FBA #B30DA6 #BD0D87 #C70066
                         #CC0059 #D1004F #D90045 #E00038 
                 #E6002E #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026 #EB0026
'''.split())
_font1 = .1 * .5
_font3 = .3 * .5
_font4 = .4 * .5


def rotate_vector(x1, y1, x2, y2):
    """
    rotate x,y vector over x2-x1, y2-y1 angle
    """
    angle = atan2(y2, x2)
    cos_rad = cos(angle)
    sin_rad = sin(angle)
    return cos_rad * x1 - sin_rad * y1, sin_rad * x1 + cos_rad * y1


class Depict:
    __slots__ = ()

    def depict(self, *, embedding=False):
        values = self._plane.values()
        min_x = min(x for x, _ in values)
        max_x = max(x for x, _ in values)
        min_y = min(y for _, y in values)
        max_y = max(y for _, y in values)

        config = self._render_config
        bonds = self._render_bonds()
        atoms, masks = self._render_atoms()
        if embedding:
            return atoms, bonds, masks, min_x, min_y, max_x, max_y

        font_size = config['font_size']
        font125 = 1.25 * font_size
        width = max_x - min_x + 3.0 * font_size
        height = max_y - min_y + 2.5 * font_size
        viewbox_x = min_x - font125
        viewbox_y = -max_y - font125

        svg = [f'<svg width="{width:.2f}cm" height="{height:.2f}cm" '
               f'viewBox="{viewbox_x:.2f} {viewbox_y:.2f} {width:.2f} '
               f'{height:.2f}" xmlns="http://www.w3.org/2000/svg" version="1.1">']
        svg.extend(self._graph_svg(atoms, bonds, masks, viewbox_x, viewbox_y, width, height))
        svg.append('</svg>')
        return '\n'.join(svg)

    @classmethod
    def _graph_svg(cls, atoms, bonds, masks, viewbox_x, viewbox_y, width, height):
        config = cls._render_config
        svg = ['  <g>']
        if bonds:
            if masks:
                uid = str(uuid4())
                svg.append(f'    <defs>\n      <mask id="mask-{uid}">\n'
                           f'        <rect x="{viewbox_x:.2f}" y="{viewbox_y:.2f}" '
                           f'width="{width:.2f}" height="{height:.2f}" fill="white"/>')
                svg.extend(cls._masks_svg(masks))
                svg.append('      </mask>\n    </defs>\n'
                           f'    <g fill="none" stroke="{config["bond_color"]}" '
                           f'stroke-width="{config["bond_width"]:.2f}"  mask="url(#mask-{uid})">')
                if len(bonds) == 1:  # SVG BUG adhoc
                    svg.append(f'      <line x1="{viewbox_x:.2f}" y1="{viewbox_y:.2f}" '
                               f'x2="{viewbox_x + width:.2f}" y2="{viewbox_y:.2f}" stroke="none"/>')
            else:
                svg.append(f'    <g fill="none" stroke="{config["bond_color"]}" '
                           f'stroke-width="{config["bond_width"]:.2f}">')
            svg.extend(bonds)
            svg.append('    </g>')

        if atoms:
            svg.append('    <g font-family="monospace">')
            svg.extend(atoms)
            svg.append('    </g>')
        svg.append('  </g>')
        return svg

    @classmethod
    def _masks_svg(cls, masks):
        config = cls._render_config

        font_size = config['font_size']
        other_size = config['other_size']
        mapping_size = config['mapping_size']
        other_font_style = config['other_font_style']
        symbols_font_style = config['symbols_font_style']

        svg = []
        stroke_width_s = font_size * .1
        stroke_width_o = other_size * .1
        stroke_width_m = mapping_size * .1

        if 'center' in masks:
            svg.append('        <g fill="black">')
            svg.extend(masks['center'])
            svg.append('        </g>')

        svg.append(f'        <g font-family="monospace" stroke="black">')
        if 'symbols' in masks:
            svg.append(f'          <g font-family="{symbols_font_style}" font-size="{font_size:.2f}" '
                       f'stroke-width="{stroke_width_s:.2f}">')
            svg.extend(masks['symbols'])
            svg.append('          </g>')
        if 'aam' in masks:
            svg.append(f'          <g font-size="{mapping_size:.2f}" stroke-width="{stroke_width_m:.2f}">')
            svg.extend(masks['aam'])
            svg.append('          </g>')
        if 'other' in masks:
            svg.append(f'          <g font-family="{other_font_style}" font-size="{other_size}" '
                       f'stroke-width="{stroke_width_o:.2f}">')
            svg.extend(masks['other'])
            svg.append('          </g>')
        if 'span' in masks:
            svg.append(f'          <g font-family="{symbols_font_style}" font-size="{font_size:.2f}" '
                       f'stroke-width="{stroke_width_o:.2f}">')
            svg.extend(masks['span'])
            svg.append('          </g>')
        svg.append('        </g>')
        return svg

    @classmethod
    def depict_settings(cls, *, carbon: bool = False, bond_color: str = 'black', font_size: float = .5,
                        aam: bool = True, aam_color: str = '#0305A7', bond_width: float = .04,
                        dashes: Tuple[float, float] = (.2, .1), query_color: str = '#5D8AA8', atoms_colors: tuple = cpk,
                        dx_ci: float = _font1, dy_ci: float = _font4, triple_space: float = .13,
                        aromatic_dashes: Tuple[float, float] = (.15, .05), dy_nh: float = _font3,
                        formed_color: str = 'green', monochrome: bool = False,  atom_radius: float = .2,
                        dy_m: float = _font4, symbols_font_style: str = 'sans-serif', other_size: float = .6 * .5,
                        double_space: float = .06, dx_m: float = _font1, span_dy: float = _font3,
                        span_size: float = .7 * .5, dx_nh: float = .15 * .5, other_font_style: str = 'monospace',
                        cgr_aromatic_space: float = .18, aam_size: float = .5 * .5, other_color: str = 'black',
                        broken_color: str = 'red', aromatic_space: float = .14, bond_radius=.02):
        """
        Settings for depict of chemical structures

        :param carbon: if True, depict atom C
        :param font_size: font size
        :param aam_size: atom-to-atom mapping font size
        :param span_size: font size for hydrogen count
        :param other_size: isotope, radical, charges, neighbors and hybridization symbols size
        :param bond_width: bond width
        :param bond_color: color of bonds
        :param aam_color: atom-to-atom mapping color
        :param query_color: hybridization and neighbors color
        :param atoms_colors: atom colors where key is atomic number - 1, value is atom color (str)
        :param broken_color: only CGRContainer: color of broken bond
        :param formed_color: only CGRContainer: color of formed bond
        :param other_color: color for charges, radicals, isotopes
        :param symbols_font_style: font style for atom symbols
        :param other_font_style: font style for mapping, charges, radicals, isotopes, hybridization and neighbors
        :param aam: if True, depict mapping
        :param monochrome: if True, colors of items in molecule not used
        :param dashes: first value is long of visible line, second is long of invisible line
        :param aromatic_space: space between simple and aromatic bonds
        :param triple_space: space between simple and triple bonds
        :param double_space: space between simple and double bonds
        :param cgr_aromatic_space: only CGRContainer: space between simple and aromatic bonds
        :param aromatic_dashes: first value is long of visible line, second is long of invisible line
        :param atom_radius: radius of atoms spheres in depict3d. if negative is multiplier to covalent radii
        :param bond_radius: radius of bonds spheres in depict3d
        :param dx_ci: x-axis offset relative to the center of the atom symbol for radical, charges, isotope
        :param dy_ci: y-axis offset relative to the center of the atom symbol for radical, charges, isotope
        :param dx_m: x-axis offset relative to the center of the atom symbol for atom-to-atom mapping
        :param dy_m: y-axis offset relative to the center of the atom symbol for atom-to-atom mapping
        :param dx_nh: x-axis offset relative to the center of the atom symbol for neighbors and hybridization
        :param dy_nh: y-axis offset relative to the center of the atom symbol for neighbors and hybridization
        :param span_dy: y-axis offset relative to the center of the atom symbol for hydrogen count
        """

        config = cls._render_config
        config['carbon'] = carbon
        config['dashes'] = dashes
        config['span_dy'] = span_dy
        config['mapping'] = aam
        config['font_size'] = font_size
        config['span_size'] = span_size
        config['other_size'] = other_size
        config['monochrome'] = monochrome
        config['bond_color'] = bond_color
        config['bond_width'] = bond_width
        config['query_color'] = query_color
        config['other_color'] = other_color
        config['bond_radius'] = bond_radius
        config['atom_radius'] = -atom_radius
        config['mapping_size'] = aam_size
        config['atoms_colors'] = atoms_colors
        config['triple_space'] = triple_space
        config['double_space'] = double_space
        config['broken_color'] = broken_color
        config['formed_color'] = formed_color
        config['mapping_color'] = aam_color
        config['aromatic_space'] = aromatic_space
        config['aromatic_dashes'] = aromatic_dashes
        config['dx_m'], config['dy_m'] = dx_m, dy_m
        config['other_font_style'] = other_font_style
        config['dx_ci'], config['dy_ci'] = dx_ci, dy_ci
        config['dx_nh'], config['dy_nh'] = dx_nh, dy_nh
        config['cgr_aromatic_space'] = cgr_aromatic_space
        config['symbols_font_style'] = symbols_font_style

    @cached_method
    def _repr_svg_(self):
        return self.depict()

    _render_config = {'carbon': False, 'atoms_colors': cpk, 'bond_color': 'black', 'font_size': .5, 'dashes': (.2, .1),
                      'aromatic_space': .14, 'triple_space': .13, 'double_space': .06, 'mapping': True, 'dx_m': _font1,
                      'mapping_color': '#0305A7', 'bond_width': .04, 'query_color': '#5D8AA8', 'broken_color': 'red',
                      'formed_color': 'green', 'aromatic_dashes': (.15, .05), 'atom_radius': -.2, 'monochrome': False,
                      'mapping_size': .5 * .5, 'symbols_font_style': 'sans-serif', 'other_size': .6 * .5,
                      'dy_m': _font4, 'span_dy': _font3, 'span_size': .7 * .5, 'dx_ci': _font1, 'dy_ci': _font4,
                      'cgr_aromatic_space': .18, 'dx_nh': .15 * .5, 'dy_nh': _font3, 'other_font_style': 'monospace',
                      'other_color': 'black', 'bond_radius': .02}


class DepictMolecule(Depict):
    __slots__ = ()

    def _render_bonds(self):
        svg = []
        plane = self._plane
        config = self._render_config

        double_space = config['double_space']
        triple_space = config['triple_space']
        dash1, dash2 = config['dashes']
        for n, m, bond in self.bonds():
            order = bond.order
            nx, ny = plane[n]
            mx, my = plane[m]
            ny, my = -ny, -my
            if order in (1, 4):
                svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
            elif order == 2:
                dx, dy = rotate_vector(0, double_space, mx - nx, ny - my)
                svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
            elif order == 3:
                dx, dy = rotate_vector(0, triple_space, mx - nx, ny - my)
                svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
            else:
                svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
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
        config = self._render_config

        aromatic_space = config['aromatic_space']
        dash3, dash4 = config['aromatic_dashes']
        # n aligned xy
        mn_x, mn_y, cn_x, cn_y = m_x - n_x, m_y - n_y, c_x - n_x, c_y - n_y

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

            return f'      <line x1="{a_x:.2f}" y1="{-a_y:.2f}" x2="{b_x:.2f}" y2="{-b_y:.2f}" ' \
                   f'stroke-dasharray="{dash3:.2f} {dash4:.2f}"/>'

    def _render_atoms(self):
        bonds = self._bonds
        plane = self._plane
        charges = self._charges
        radicals = self._radicals
        hydrogens = self._hydrogens
        config = self._render_config

        carbon = config['carbon']
        mapping = config['mapping']
        span_size = config['span_size']
        font_size = config['font_size']
        monochrome = config['monochrome']
        other_size = config['other_size']
        atoms_colors = config['atoms_colors']
        mapping_font = config['mapping_size']
        dx_m, dy_m = config['dx_m'], config['dy_m']
        dx_ci, dy_ci = config['dx_ci'], config['dy_ci']
        symbols_font_style = config['symbols_font_style']

        # for cumulenes
        cumulenes = {y for x in self._cumulenes(heteroatoms=True) if len(x) > 2 for y in x[1:-1]}

        if monochrome:
            map_fill = other_fill = 'black'
        else:
            map_fill = config['mapping_color']
            other_fill = config['other_color']

        svg = []
        maps = []
        others = []
        font2 = .2 * font_size
        font3 = .3 * font_size
        font4 = .4 * font_size
        font5 = .5 * font_size
        font6 = .6 * font_size
        font7 = .7 * font_size
        font15 = .15 * font_size
        font25 = .25 * font_size
        mask = defaultdict(list)
        for n, atom in self._atoms.items():
            x, y = plane[n]
            y = -y
            symbol = atom.atomic_symbol
            if not bonds[n] or symbol != 'C' or carbon or atom.charge or atom.is_radical or atom.isotope \
                    or n in cumulenes:
                h = hydrogens[n]
                if h == 1:
                    h = 'H'
                    span = ''
                elif h:
                    span = f'<tspan  dy="{config["span_dy"]:.2f}" font-size="{span_size:.2f}">{h}</tspan>'
                    h = 'H'
                else:
                    h = ''
                    span = ''

                if charges[n]:
                    t = f'{_render_charge[charges[n]]}{"↑" if radicals[n] else ""}'
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                  f'{t}</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'{t}</text>')
                elif radicals[n]:
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">↑</text>')
                    mask['other'].append(f'            <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}"'
                                         f' dy="-{dy_ci:.2f}">↑</text>')
                if atom.isotope:
                    t = atom.isotope
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'text-anchor="end">{t}</text>')
                    mask['other'].append(f'            <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_ci:.2f}"'
                                         f' dy="-{dy_ci:.2f}" text-anchor="end">{t}</text>')
                svg.append(f'      <g fill="{"black" if monochrome else atoms_colors[atom.atomic_number - 1]}" '
                           f'font-family="{symbols_font_style }">')
                if len(symbol) > 1:
                    dx = font7
                    dx_mm = dx_m + font5
                    if symbol[-1] in ('l', 'i', 'r', 't'):
                        rx = font6
                        ax = font25
                    else:
                        rx = font7
                        ax = font15
                    mask['center'].append(f'          <ellipse cx="{x - ax:.2f}" cy="{y:.2f}" rx="{rx}" ry="{font4}"/>')
                else:
                    if symbol == 'I':
                        dx = font15
                        dx_mm = dx_m
                    else:
                        dx = font4
                        dx_mm = dx_m + font2
                    mask['center'].append(f'          <circle cx="{x:.2f}" cy="{y:.2f}" r="{font4:.2f}"/>')
                svg.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" dy="{font4:.2f}" '
                           f'font-size="{font_size:.2f}">{symbol}{h}{span}</text>')
                mask['symbols'].append(f'            <text x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" '
                                       f'dy="{font4:.2f}">{symbol}{h}</text>')
                if span:
                    mask['span'].append(f'            <text x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" dy="{font4:.2f}">'
                                        f'{symbol}{h}{span}</text>')
                svg.append('      </g>')

                if mapping:
                    maps.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_mm:.2f}" dy="{dy_m + font3:.2f}" '
                                f'text-anchor="end">{n}</text>')
                    mask['aam'].append(f'            <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_mm:.2f}" '
                                       f'dy="{dy_m + font3:.2f}" text-anchor="end">{n}</text>')

            elif mapping:
                maps.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_m:.2f}" dy="{dy_m:.2f}" '
                            f'text-anchor="end">{n}</text>')
                mask['aam'].append(f'            <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_m:.2f}" dy="{dy_m:.2f}" '
                                   f'text-anchor="end">{n}</text>')
        if others:
            svg.append(f'      <g font-family="{config["other_font_style"]}" fill="{other_fill}" '
                       f'font-size="{other_size:.2f}">')
            svg.extend(others)
            svg.append('      </g>')
        if mapping:
            svg.append(f'      <g fill="{map_fill}" font-size="{mapping_font:.2f}">')
            svg.extend(maps)
            svg.append('      </g>')
        return svg, mask


class DepictReaction:
    __slots__ = ()

    def depict(self):
        if not self._arrow:
            self.fix_positions()

        r_atoms = []
        r_bonds = []
        r_masks = []
        r_max_x = r_max_y = r_min_y = 0
        for m in self.molecules():
            atoms, bonds, masks, min_x, min_y, max_x, max_y = m.depict(embedding=True)
            r_atoms.append(atoms)
            r_bonds.append(bonds)
            r_masks.append(masks)
            if max_x > r_max_x:
                r_max_x = max_x
            if max_y > r_max_y:
                r_max_y = max_y
            if min_y < r_min_y:
                r_min_y = min_y

        config = Depict._render_config
        font_size = config['font_size']
        font125 = 1.25 * font_size
        width = r_max_x + 3.0 * font_size
        height = r_max_y - r_min_y + 2.5 * font_size
        viewbox_x = -font125
        viewbox_y = -r_max_y - font125

        svg = [f'<svg width="{width:.2f}cm" height="{height:.2f}cm" '
               f'viewBox="{viewbox_x:.2f} {viewbox_y:.2f} {width:.2f} '
               f'{height:.2f}" xmlns="http://www.w3.org/2000/svg" version="1.1">\n'
               '  <defs>\n    <marker id="arrow" markerWidth="10" markerHeight="10" '
               'refX="0" refY="3" orient="auto">\n      <path d="M0,0 L0,6 L9,3"/>\n    </marker>\n  </defs>\n'
               f'  <line x1="{self._arrow[0]:.2f}" y1="0" x2="{self._arrow[1]:.2f}" y2="0" '
               'fill="none" stroke="black" stroke-width=".04" marker-end="url(#arrow)"/>']

        sings_plus = self._signs
        if sings_plus:
            svg.append(f'  <g fill="none" stroke="black" stroke-width=".04">')
            for x in sings_plus:
                svg.append(f'    <line x1="{x + .35:.2f}" y1="0" x2="{x + .65:.2f}" y2="0"/>')
                svg.append(f'    <line x1="{x + .5:.2f}" y1="0.15" x2="{x + .5:.2f}" y2="-0.15"/>')
            svg.append('  </g>')

        for atoms, bonds, masks in zip(r_atoms, r_bonds, r_masks):
            svg.extend(Depict._graph_svg(atoms, bonds, masks, viewbox_x, viewbox_y, width, height))
        svg.append('</svg>')
        return '\n'.join(svg)

    @staticmethod
    def depict_settings(**kwargs):
        """Settings for depict of chemical structures"""
        Depict.depict_settings(**kwargs)

    @cached_method
    def _repr_svg_(self):
        return self.depict()


class DepictCGR(Depict):
    __slots__ = ()

    def _render_bonds(self):
        plane = self._plane
        config = self._render_config

        broken = config['broken_color']
        formed = config['formed_color']
        dash1, dash2 = config['dashes']
        double_space = config['double_space']
        triple_space = config['triple_space']

        svg = []
        ar_bond_colors = defaultdict(dict)
        for n, m, bond in self.bonds():
            order, p_order = bond.order, bond.p_order
            nx, ny = plane[n]
            mx, my = plane[m]
            ny, my = -ny, -my
            rv = partial(rotate_vector, 0, x2=mx - nx, y2=ny - my)
            if order == 1:
                if p_order == 1:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order is None:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
                else:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
            elif order == 4:
                if p_order == 4:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 1:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 2:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order == 3:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"  stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order is None:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = broken
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
                else:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = None
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
            elif order == 2:
                if p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order is None:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                else:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
            elif order == 3:
                if p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" 'f'stroke="{broken}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" '
                               f'y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" '
                               f'y2="{my + dy:.2f}" stroke="{broken}"/>')
                elif p_order == 2:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                elif p_order is None:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" '
                               f'x2="{mx:.2f}" y2="{my:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                else:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'      <line x1="{nx + dx3:.2f}" y1="{ny - dy3:.2f}" x2="{mx + dx3:.2f}" '
                               f'y2="{my - dy3:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx3:.2f}" y1="{ny + dy3:.2f}" x2="{mx - dx3:.2f}" '
                               f'y2="{my + dy3:.2f}" stroke="{broken}"/>')
            elif order is None:
                if p_order == 1:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = formed
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                elif p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" '
                               f'y2="{my - dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" '
                               f'y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                else:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"/>')
            else:
                if p_order == 8:
                    svg.append(f'        <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order == 4:
                    ar_bond_colors[n][m] = ar_bond_colors[m][n] = None
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                elif p_order == 2:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order == 3:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'      <line x1="{nx + dx3:.2f}" y1="{ny - dy3:.2f}" x2="{mx + dx3:.2f}" '
                               f'y2="{my - dy3:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx3:.2f}" y1="{ny + dy3:.2f}" '
                               f'x2="{mx - dx3:.2f}" y2="{my + dy3:.2f}" stroke="{formed}"/>')
                else:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
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
        config = self._render_config

        dash1, dash2 = config['dashes']
        dash3, dash4 = config['aromatic_dashes']
        aromatic_space = config['cgr_aromatic_space']
        # n aligned xy
        mn_x, mn_y, cn_x, cn_y = m_x - n_x, m_y - n_y, c_x - n_x, c_y - n_y

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
            if color:
                return f'      <line x1="{an_x + n_x:.2f}" y1="{-an_y - n_y:.2f}" x2="{bn_x + n_x:.2f}" ' \
                       f'y2="{-bn_y - n_y:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{color}"/>'
            elif color is None:
                dash3, dash4 = dash1, dash2
            return f'      <line x1="{an_x + n_x:.2f}" y1="{-an_y - n_y:.2f}"' \
                   f' x2="{bn_x + n_x:.2f}" y2="{-bn_y - n_y:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}"/>'

    def _render_atoms(self):
        bonds = self._bonds
        plane = self._plane
        charges = self._charges
        radicals = self._radicals
        p_charges = self._p_charges
        config = self._render_config
        p_radicals = self._p_radicals

        carbon = config['carbon']
        font_size = config['font_size']
        other_size = config['other_size']
        monochrome = config['monochrome']
        atoms_colors = config['atoms_colors']
        dx_ci, dy_ci = config['dx_ci'], config['dy_ci']
        symbols_font_style = config['symbols_font_style']
        font4 = .4 * font_size
        font6 = .6 * font_size
        font7 = .7 * font_size
        font15 = .15 * font_size
        font25 = .25 * font_size
        din_charges = {m[0]: m[1] != n[1] for m, n in zip(charges.items(), p_charges.items())}
        din_radicals = {m[0]: m[1] != n[1] for m, n in zip(radicals.items(), p_radicals.items())}

        if monochrome:
            other_fill = 'black'
        else:
            other_fill = config['other_color']

        svg = []
        others = []
        mask = defaultdict(list)
        for n, atom in self._atoms.items():
            x, y = plane[n]
            y = -y
            symbol = atom.atomic_symbol
            if not bonds[n] or symbol != 'C' or carbon or atom.charge or atom.is_radical or atom.isotope \
                    or din_charges[n] or din_radicals[n]:

                if radicals[n]:
                    r = '↑' if p_radicals[n] else '↑↓'
                elif p_radicals[n]:
                    r = '↓↑'
                else:
                    r = ''

                if charges[n] != p_charges[n]:
                    t = _render_p_charge[charges[n]][p_charges[n]]
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}">{t}{r}</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'{t}{r}</text>')
                if charges[n]:
                    t = _render_charge[charges[n]]
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}">{t}{r}</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'{t}{r}</text>')
                elif r:
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}">{r}</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'{r}</text>')

                if atom.isotope:
                    t = atom.isotope
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}" text-anchor="end">{t}</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_ci:.2f}"'
                                         f' dy="-{dy_ci:.2f}">{t}</text>')

                svg.append(f'      <g fill="{"black" if monochrome else atoms_colors[atom.atomic_number - 1]}" '
                           f'font-family="{symbols_font_style}">')
                if len(symbol) > 1:
                    dx = font7
                    if symbol[-1] in ('l', 'i', 'r', 't'):
                        rx = font6
                        ax = font25
                    else:
                        rx = font7
                        ax = font15
                    mask['center'].append(f'          <ellipse cx="{x - ax:.2f}" cy="{y:.2f}" rx="{rx}" ry="{font4}"/>')
                else:
                    if symbol == 'I':
                        dx = font15
                    else:
                        dx = font4
                    mask['center'].append(f'          <circle cx="{x:.2f}" cy="{y:.2f}" r="{font4:.2f}"/>')
                svg.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" dy="{font4:.2f}" '
                           f'font-size="{font_size:.2f}">{symbol}</text>')
                mask['symbols'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" '
                                       f'dy="{font4:.2f}">{symbol}</text>')
                svg.append('      </g>')

        if others:
            svg.append(f'      <g font-family="{config["other_font_style"]}" fill="{other_fill}" '
                       f'font-size="{other_size:.2f}">')
            svg.extend(others)
            svg.append('      </g>')
        return svg, mask


class DepictQuery(Depict):
    __slots__ = ()

    def _render_bonds(self):
        svg = []
        plane = self._plane
        config = self._render_config

        dash1, dash2 = config['dashes']
        double_space = config['double_space']
        triple_space = config['triple_space']
        dash3, dash4 = config['aromatic_dashes']
        for n, m, bond in self.bonds():
            nx, ny = plane[n]
            mx, my = plane[m]
            ny, my = -ny, -my
            if bond == 1:
                svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
            elif bond == 4:
                dx, dy = rotate_vector(0, double_space, mx - nx, ny - my)
                svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" y2="{my + dy:.2f}" '
                           f'stroke-dasharray="{dash3:.2f} {dash4:.2f}"/>')
            elif bond == 2:
                dx, dy = rotate_vector(0, double_space, mx - nx, ny - my)
                svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
            elif bond == 3:
                dx, dy = rotate_vector(0, triple_space, mx - nx, ny - my)
                svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
            else:  # other query bonds
                svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                           f'stroke-dasharray="{dash1:.2f} {dash2:.2f}"/>')
        return svg

    def _render_atoms(self):
        plane = self._plane
        config = self._render_config

        carbon = config['carbon']
        mapping = config['mapping']
        font_size = config['font_size']
        other_size = config['other_size']
        monochrome = config['monochrome']
        atoms_colors = config['atoms_colors']
        mapping_size = config['mapping_size']
        dx_m, dy_m = config['dx_m'], config['dy_m']
        other_font_style = config['other_font_style']
        dx_ci, dy_ci = config['dx_ci'], config['dy_ci']
        dx_nh, dy_nh = config['dx_nh'], config['dy_nh']
        symbols_font_style = config['symbols_font_style']

        font2 = .2 * font_size
        font3 = .3 * font_size
        font4 = .4 * font_size
        font5 = .5 * font_size
        font6 = .6 * font_size
        font7 = .7 * font_size
        font15 = .15 * font_size
        font25 = .25 * font_size

        # for cumulenes
        cumulenes = {y for x in self._cumulenes(heteroatoms=True) if len(x) > 2 for y in x[1:-1]}

        if monochrome:
            map_fill = query_fill = other_fill = 'black'
        else:
            other_fill = config['other_color']
            map_fill = config['mapping_color']
            query_fill = config['query_color']

        svg = []
        maps = []
        others = []
        nghbrs = []
        hbrdztns = []
        mask = defaultdict(list)
        for n, atom in self._atoms.items():
            x, y = plane[n]
            y = -y
            single = not self._bonds[n]
            symbol = atom.atomic_symbol
            if single or symbol != 'C' or carbon or atom.charge or atom.is_radical or n in cumulenes:
                if atom.charge:
                    t = _render_charge[atom.charge]
                    tt = "↑" if atom.is_radical else ""
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}">{t}{tt}</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'{t}{tt}</text>')
                elif atom.is_radical:
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}">↑</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'↑</text>')

                svg.append(f'      <g fill="{"black" if monochrome else atoms_colors[atom.atomic_number - 1]}" '
                           f'font-family="{symbols_font_style}">')
                if len(symbol) > 1:
                    dx = font7
                    dx_mm = dx_m + font5
                    dx_nhh = dx_nh + font5
                    if symbol[-1] in ('l', 'i', 'r', 't'):
                        rx = font6
                        ax = font25
                    else:
                        rx = font7
                        ax = font15
                    mask['center'].append(f'          <ellipse cx="{x - ax:.2f}" cy="{y:.2f}" rx="{rx}" ry="{font4}"/>')
                else:
                    if symbol == 'I':
                        dx = font15
                        dx_mm = dx_m
                        dx_nhh = dx_nh
                    else:
                        dx = font4
                        dx_mm = dx_m + font2
                        dx_nhh = dx_nh + font2
                    mask['center'].append(f'          <circle cx="{x:.2f}" cy="{y:.2f}" r="{font4:.2f}"/>')
                svg.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" dy="{font4:.2f}" '
                           f'font-size="{font_size:.2f}">{symbol}</text>')
                mask['symbols'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" '
                                       f'dy="{font4:.2f}">{symbol}</text>')
                svg.append('      </g>')
                dx_nnh = dx_nhh

                if mapping:
                    maps.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_mm:.2f}" dy="{dy_m + font3:.2f}" '
                                f'text-anchor="end">{n}</text>')
                    mask['aam'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_mm:.2f}" '
                                       f'dy="{dy_m + font3:.2f}" text-anchor="end">{n}</text>')
            elif mapping:
                dx_nnh = dx_nh
                maps.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_m}" dy="{dy_m:.2f}" '
                            f'text-anchor="end">{n}</text>')
                mask['aam'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="-{dx_m}" dy="{dy_m:.2f}" '
                                   f'text-anchor="end">{n}</text>')
            else:
                dx_nnh = dx_nh

            if atom.neighbors:
                level = .6 * other_size
                nn = "".join(str(x) for x in atom.neighbors)
                nghbrs.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_nnh}" dy="{dy_nh:.2f}" '
                              f'text-anchor="start">{nn}</text>')
                mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_nnh:.2f}" dy="{dy_nh:.2f}">'
                                     f'{nn}</text>')
            else:
                level = 0

            if atom.hybridization:
                hh = "".join(_render_hybridization[x] for x in atom.hybridization)
                hbrdztns.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_nnh}" dy="{level + dy_nh:.2f}" '
                                f'text-anchor="start">{hh}</text>')
                mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_nnh:.2f}" '
                                     f'dy="{level + dy_nh:.2f}">{hh}</text>')

        if nghbrs:
            svg.append(f'      <g fill="{query_fill}" font-family="{other_font_style}" font-size="{other_size:.2f}">')
            svg.extend(nghbrs)
            if hbrdztns:
                svg.extend(hbrdztns)
            svg.append('      </g>')
        elif hbrdztns:
            svg.append(f'      <g fill="{query_fill}" font-family="{other_font_style}" font-size="{other_size:.2f}">')
            svg.extend(hbrdztns)
            svg.append('    </g>')

        if others:
            svg.append(f'      <g fill="{other_fill}" font-family="{other_font_style}" font-size="{other_size:.2f}">')
            svg.extend(others)
            svg.append('      </g>')

        if mapping:
            svg.append(f'      <g fill="{map_fill}" font-size="{mapping_size:.2f}">')
            svg.extend(maps)
            svg.append('      </g>')

        return svg, mask


class DepictQueryCGR(Depict):
    __slots__ = ()

    def _render_bonds(self):
        svg = []
        plane = self._plane
        config = self._render_config

        dash1, dash2 = config['dashes']
        broken = config['broken_color']
        formed = config['formed_color']
        double_space = config['double_space']
        triple_space = config['triple_space']
        dash3, dash4 = config['aromatic_dashes']

        for n, m, bond in self.bonds():
            order, p_order = bond.order, bond.p_order
            nx, ny = plane[n]
            mx, my = plane[m]
            ny, my = -ny, -my
            rv = partial(rotate_vector, 0, x2=mx - nx, y2=ny - my)
            if order == 1:
                if p_order == 1:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                elif p_order == 4:
                    dx, dy = rv(double_space)
                    svg.append(
                        f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{formed}"/>')
                elif p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order is None:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
                else:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} '
                               f'{dash2:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
            elif order == 4:
                if p_order == 4:
                    dx, dy = rv(double_space)
                    svg.append(
                        f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(double_space)
                    svg.append(
                        f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{broken}"/>')
                elif p_order == 2:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"  stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{broken}"/>')
                elif p_order == 3:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx3:.2f}" y1="{ny + dy3:.2f}" x2="{mx - dx3:.2f}" '
                               f'y2="{my + dy3:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx + dx3:.2f}" y1="{ny - dy3:.2f}" x2="{mx + dx3:.2f}" '
                               f'y2="{my - dy3:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{broken}"/>')
                elif p_order is None:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" '
                               f'y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{broken}"/>')
                else:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{broken}"/>')
            elif order == 2:
                if p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                elif p_order == 4:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{formed}"/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order is None:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                else:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} '
                               f'{dash2:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
            elif order == 3:
                if p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"  stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" 'f'stroke="{broken}"/>')
                elif p_order == 4:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx3:.2f}" y1="{ny + dy3:.2f}" '
                               f'x2="{mx - dx3:.2f}" y2="{my + dy3:.2f}"/>')
                    svg.append(f'      <line x1="{nx + dx3:.2f}" y1="{ny - dy3:.2f}" x2="{mx + dx3:.2f}" '
                               f'y2="{my - dy3:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{formed}"/>')
                elif p_order == 2:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                elif p_order is None:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" '
                               f'x2="{mx:.2f}" y2="{my:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                else:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'      <line x1="{nx + dx3:.2f}" y1="{ny - dy3:.2f}" x2="{mx + dx3:.2f}" '
                               f'y2="{my - dy3:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx3:.2f}" y1="{ny + dy3:.2f}" x2="{mx - dx3:.2f}" '
                               f'y2="{my + dy3:.2f}" stroke="{broken}"/>')
            elif order is None:
                if p_order == 1:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                elif p_order == 4:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{formed}"/>')
                elif p_order == 2:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}" '
                               f'y2="{my - dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}" '
                               f'y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order == 3:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                else:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{formed}"/>')
            else:
                if p_order == 8:
                    svg.append(f'        <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}"/>')
                elif p_order == 1:
                    dx, dy = rv(double_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order == 4:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" x2="{mx - dx:.2f}"'
                               f' y2="{my + dy:.2f}" stroke-dasharray="{dash3:.2f} {dash4:.2f}" stroke="{formed}"/>')
                elif p_order == 2:
                    dx, dy = rv(triple_space)
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" x2="{mx + dx:.2f}"'
                               f' y2="{my - dy:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}"'
                               f' stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                elif p_order == 3:
                    dx, dy = rv(double_space)
                    dx3 = 3 * dx
                    dy3 = 3 * dy
                    svg.append(f'      <line x1="{nx + dx3:.2f}" y1="{ny - dy3:.2f}" x2="{mx + dx3:.2f}" '
                               f'y2="{my - dy3:.2f}" stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')
                    svg.append(f'      <line x1="{nx + dx:.2f}" y1="{ny - dy:.2f}" '
                               f'x2="{mx + dx:.2f}" y2="{my - dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx:.2f}" y1="{ny + dy:.2f}" '
                               f'x2="{mx - dx:.2f}" y2="{my + dy:.2f}" stroke="{formed}"/>')
                    svg.append(f'      <line x1="{nx - dx3:.2f}" y1="{ny + dy3:.2f}" '
                               f'x2="{mx - dx3:.2f}" y2="{my + dy3:.2f}" stroke="{formed}"/>')
                else:
                    svg.append(f'      <line x1="{nx:.2f}" y1="{ny:.2f}" x2="{mx:.2f}" y2="{my:.2f}" '
                               f'stroke-dasharray="{dash1:.2f} {dash2:.2f}" stroke="{broken}"/>')
        return svg

    def _render_atoms(self):
        plane = self._plane
        config = self._render_config

        carbon = config['carbon']
        font_size = config['font_size']
        other_size = config['other_size']
        monochrome = config['monochrome']
        atoms_colors = config['atoms_colors']
        other_font_style = config['other_font_style']
        dx_ci, dy_ci = config['dx_ci'], config['dy_ci']
        dx_nh, dy_nh = config['dx_nh'], config['dy_nh']
        symbols_font_style = config['symbols_font_style']

        font2 = .2 * font_size
        font4 = .4 * font_size
        font5 = .5 * font_size
        font6 = .6 * font_size
        font7 = .7 * font_size
        font15 = .15 * font_size
        font25 = .25 * font_size

        if monochrome:
            query_fill = other_fill = 'black'
        else:
            other_fill = config['other_color']
            query_fill = config['query_color']

        svg = []
        others = []
        nghbrs = []
        hbrdztns = []
        mask = defaultdict(list)
        for n, atom in self._atoms.items():
            x, y = plane[n]
            y = -y
            single = not self._bonds[n]
            symbol = atom.atomic_symbol
            if single or symbol != 'C' or carbon or atom.charge or atom.is_radical:
                if atom.charge:
                    t = _render_charge[atom.charge]
                    tt = "↑" if atom.is_radical else ""
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}">{t}{tt}</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'{t}{tt}</text>')
                elif atom.is_radical:
                    others.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}" '
                                  f'font-size="{other_size:.2f}">↑</text>')
                    mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_ci:.2f}" dy="-{dy_ci:.2f}">'
                                         f'↑</text>')

                svg.append(f'      <g fill="{"black" if monochrome else atoms_colors[atom.atomic_number - 1]}" '
                           f'font-family="{symbols_font_style}">')
                if len(symbol) > 1:
                    dx = font7
                    dx_nhh = dx_nh + font5
                    if symbol[-1] in ('l', 'i', 'r', 't'):
                        rx = font6
                        ax = font25
                    else:
                        rx = font7
                        ax = font15
                    mask['center'].append(f'          <ellipse cx="{x - ax:.2f}" cy="{y:.2f}" rx="{rx}" ry="{font4}"/>')
                else:
                    if symbol == 'I':
                        dx = font15
                        dx_nhh = dx_nh
                    else:
                        dx = font4
                        dx_nhh = dx_nh + font2
                    mask['center'].append(f'          <circle cx="{x:.2f}" cy="{y:.2f}" r="{font4:.2f}"/>')
                svg.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" dy="{font4:.2f}" '
                           f'font-size="{font_size:.2f}">{symbol}</text>')
                mask['symbols'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="-{dx:.2f}" '
                                       f'dy="{font4:.2f}">{symbol}</text>')
                svg.append('      </g>')
                dx_nnh = dx_nhh
            else:
                dx_nnh = dx_nh

            if atom.neighbors:
                level = .6 * other_size
                nn = ''.join(str(x) for x in atom.neighbors)
                if atom.p_neighbors:
                    pn = ''.join(str(x) for x in atom.p_neighbors)
                else:
                    pn = '0'
                nghbrs.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_nnh:.2f}" dy="{dy_nh:.2f}" '
                              f'text-anchor="start">{nn}»{pn}</text>')
                mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_nnh:.2f}" dy="{dy_nh:.2f}">'
                                     f'{nn}»{pn}</text>')
            else:
                level = 0

            if atom.hybridization:
                hh = ''.join(_render_hybridization[x] for x in atom.hybridization)
                if atom.p_hybridization:
                    ph = ''.join(_render_hybridization[x] for x in atom.p_hybridization)
                else:
                    ph = '0'
                hbrdztns.append(f'        <text x="{x:.2f}" y="{y:.2f}" dx="{dx_nnh:.2f}" dy="{level + dy_nh:.2f}" '
                                f'text-anchor="start">{hh}»{ph}</text>')
                mask['other'].append(f'           <text x="{x:.2f}" y="{y:.2f}" dx="{dx_nnh:.2f}" '
                                     f'dy="{level + dy_nh:.2f}">{hh}»{ph}</text>')

        if nghbrs:
            svg.append(f'      <g fill="{query_fill}" font-family="{other_font_style}" font-size="{other_size:.2f}">')
            svg.extend(nghbrs)
            if hbrdztns:
                svg.extend(hbrdztns)
            svg.append('      </g>')
        elif hbrdztns:
            svg.append(f'      <g fill="{query_fill}" font-family="{other_font_style}" font-size="{other_size:.2f}">')
            svg.extend(hbrdztns)
            svg.append('      </g>')

        if others:
            svg.append(f'      <g fill="{other_fill}" font-family="{other_font_style}" font-size="{other_size:.2f}">')
            svg.extend(others)
            svg.append('      </g>')

        return svg, mask


_render_hybridization = {1: 's', 2: 'd', 3: 't', 4: 'a'}
_render_charge = {-3: '3-', -2: '2-', -1: '-', 1: '+', 2: '2+', 3: '3+'}
_render_p_charge = {-3: {-2: '-3»-2', -1: '-3»-', 0: '-3»0', 1: '-3»+', 2: '-3»2', 3: '-3»3'},
                    -2: {-3: '-2»-3', -1: '-2»-', 0: '-2»0', 1: '-2»+', 2: '-2»2', 3: '-2»3'},
                    -1: {-3: '-»-3', -2: '-»-2', 0: '-»0', 1: '-»+', 2: '-»2', 3: '-»3'},
                    0: {-3: '0»-3', -2: '0»-2', -1: '0»-', 1: '0»+', 2: '0»2', 3: '0»3'},
                    1: {-3: '+»-3', -2: '+»-2', -1: '+»-', 0: '+»0', 2: '+»2', 3: '+»3'},
                    2: {-3: '2»-3', -2: '2»-2', -1: '2»-', 0: '2»0', 1: '2»+', 3: '2»3'},
                    3: {-3: '3»-3', -2: '3»-2', -1: '3»-', 0: '3»0', 1: '3»+', 2: '3»2'}}


__all__ = ['DepictMolecule', 'DepictReaction', 'DepictCGR', 'DepictQuery', 'DepictQueryCGR']
