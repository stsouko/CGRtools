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
    def __bond(self):
        svg = []
        bonds = defaultdict(set)
        nodes = self._node
        for n, m, bond in self.bonds():
            na, ma = nodes[n], nodes[m]
            bond, aromatic = self._render_bond(bond, na.x, na.y, ma.x, ma.y)
            if aromatic:
                bonds[n].add(m)
                bonds[m].add(n)
            svg.append(bond)

        for ring in self.aromatic_rings:
            c_x = sum(nodes[x].x for x in ring) / len(ring)
            c_y = sum(nodes[y].y for y in ring) / len(ring)

            for n, m in zip(ring, ring[1:]):
                if m in bonds[n]:
                    n, m = nodes[n], nodes[m]
                    svg.append(self.__render_aromatic_bond(n.x, n.y, m.x, m.y, c_x, c_y))

            n, m = ring[-1], ring[0]
            if m in bonds[n]:
                n, m = nodes[n], nodes[m]
                svg.append(self.__render_aromatic_bond(n.x, n.y, m.x, m.y, c_x, c_y))
        return svg

    def __render_aromatic_bond(self, n_x, n_y, m_x, m_y, c_x, c_y):
        # n aligned xy
        mn_x, mn_y, cn_x, cn_y = m_x - n_x, m_y - n_y, c_x - n_x, c_y - n_y

        # nm reoriented xy
        mr_x, mr_y = hypot(mn_x, mn_y), 0
        cr_x, cr_y = rotate_vector(cn_x, cn_y, mn_x, mn_y)

        if self.aromatic_space / cr_y < .65:
            if cr_y > 0:
                ar_y = br_y = self.aromatic_space
            else:
                ar_y = br_y = -self.aromatic_space

            ar_x = self.aromatic_space * cr_x / abs(cr_y)
            br_x = ((abs(cr_y) - self.aromatic_space) * mr_x + self.aromatic_space * cr_x) / abs(cr_y)

            # backward reorienting
            an_x, an_y = rotate_vector(ar_x, ar_y, mn_x, -mn_y)
            bn_x, bn_y = rotate_vector(br_x, br_y, mn_x, -mn_y)
            a_x, a_y = an_x + n_x, an_y + n_y
            b_x, b_y = bn_x + n_x, bn_y + n_y

            return f'      <line x1="{a_x:.2f}" y1="{-a_y:.2f}" x2="{b_x:.2f}" y2="{-b_y:.2f}" ' \
                f'stroke-dasharray="{self.dashes[0]:.2f} {self.dashes[1]:.2f}" />'

    def depict(self, embedding=False):
        min_x = min(x.x for x in self._node.values())
        max_x = max(x.x for x in self._node.values())
        min_y = min(x.y for x in self._node.values())
        max_y = max(x.y for x in self._node.values())

        bonds = ['    <g fill="none" stroke="black" stroke-width=".03"  mask="url(#mask)">']
        bonds.extend(self.__bond())
        bonds.append('    </g>')

        width = max_x - min_x + 2.5 * self.font
        height = max_y - min_y + 2.5 * self.font
        atoms = []
        svg = ['    <defs>', '      <mask id="mask" maskContentUnits="userSpaceOnUse">'
               f'<rect x="{min_x - 1.25 * self.font:.2f}" y="{-max_y - 1.25 * self.font:.2f}" '
               f'width="{width:.2f}" height="{height:.2f}" fill="white" />']
        for n, atom in self.atoms():
            tmp, mask = self._render_atom(atom, not bool(self._adj[n]))
            atoms.extend(tmp)
            svg.extend(mask)
        svg.append('       </mask>')
        svg.append('    </defs>')

        if atoms:
            atoms.insert(0, '     <g font-family="sans-serif">')
            atoms.append('    </g>')

        if not embedding:
            svg.insert(0, f'<svg width="{width:.2f}cm" height="{height:.2f}cm" '
                          f'viewBox="{min_x - 1.25 * self.font:.2f} {-max_y - 1.25 * self.font:.2f} {width:.2f} '
                          f'{height:.2f}" '
                          f'xmlns="http://www.w3.org/2000/svg" version="1.1">')
            atoms.append('</svg>')
        svg.extend(bonds)
        svg.extend(atoms)
        if embedding:
            return '\n'.join(svg), max_x, max_y
        return '\n'.join(svg)

    @cached_method
    def _repr_svg_(self):
        return self.depict()

    carbon = False
    colors = cpk
    font = .4
    sup_font = .75 * font
    up_font = -.5 * font
    double_space = .04
    triple_space = .07
    aromatic_space = .08
    dashes = (.2, .1)


class DepictMolecule(Depict):
    def _render_atom(self, atom, single):
        svg = []
        mask = []
        if single or atom.element != 'C' or self.carbon or atom.charge or atom.multiplicity or \
                atom.isotope != atom.common_isotope:
            x_shift = -shifts[atom.element] * self.font
            y_shift = .35 * self.font
            radius = -1.5 * x_shift
            svg.append(f'      <g fill="{self.colors[atom.element]}">')
            svg.append(f'          <text x="{atom.x + x_shift:.2f}" y="{y_shift - atom.y:.2f}" '
                       f'font-size="{self.font:.2f}" '
                       f'>{atom.element}</text>')
            if atom.charge:
                svg.append(f'          <text x="{atom.x - x_shift:.2f}" y="{-y_shift - atom.y:.2f}" '
                           f'font-size="{self.sup_font:.2f}">{charge_str[atom.charge]}</text>')
            if atom.multiplicity:
                svg.append(f'          <text x="{atom.x + x_shift:.2f}" y="{self.up_font - atom.y:.2f}" '
                           f'font-size="{self.sup_font:.2f}">{multiplicity_str[atom.multiplicity]}</text>')
            if atom.isotope != atom.common_isotope:
                svg.append(f'          <text x="{atom.x - self.font:.2f}" y="{-y_shift - atom.y:.2f}" '
                           f'font-size="{self.sup_font:.2f}">{atom.isotope}</text>')
            svg.append('      </g>')
            mask = [f'          <circle cx="{atom.x}" cy="{-atom.y}" r="{radius}" fill="black"/>']
        return svg, mask

    def _render_bond(self, bond, nx, ny, mx, my):
        if bond.order == 1:
            return f'      <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" />', False
        elif bond.order == 2:
            dx, dy = rotate_vector(0, self.double_space, mx - nx, my - ny)
            return f'      <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" x2="{mx + dx:.2f}" y2="{-my + dy:.2f}" />\n' \
                   f'      <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" />', \
                   False
        elif bond.order == 3:
            dx, dy = rotate_vector(0, self.triple_space, mx - nx, my - ny)
            return f'      <line x1="{nx + dx:.2f}" y1="{-ny + dy:.2f}" x2="{mx + dx:.2f}" y2="{-my + dy:.2f}" />\n' \
                   f'      <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" />\n' \
                   f'      <line x1="{nx - dx:.2f}" y1="{-ny - dy:.2f}" x2="{mx - dx:.2f}" y2="{-my - dy:.2f}" />', \
                   False
        elif bond.order == 4:
            return f'      <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" />', True
        return f'      <line x1="{nx:.2f}" y1="{-ny:.2f}" x2="{mx:.2f}" y2="{-my:.2f}" \n' \
               f'stroke-dasharray="{self.dashes[0]:.2f} {self.dashes[1]:.2f}" />', False


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
        for ml in (self.reactants, self.reagents, self.products):
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
