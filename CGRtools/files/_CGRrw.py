# -*- coding: utf-8 -*-
#
#  Copyright 2014-2018 Ramil Nugmanov <stsouko@live.ru>
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
from abc import abstractmethod
from collections import defaultdict
from itertools import count
from io import StringIO, BytesIO, TextIOWrapper
from logging import warning
from pathlib import Path
from ..containers import ReactionContainer, MoleculeContainer, CGRContainer, QueryContainer, QueryCGRContainer
from ..exceptions import MappingError
from ..periodictable import elements_set


elements_set = elements_set.copy()
elements_set.discard('A')
cgr_keys = dict(extrabond=2, dynbond=2, dynatom=1, isotope=1, atomhyb=1, atomneighbors=1, dynatomhyb=1,
                dynatomneighbors=1)


def _pyramid_volume(n, u, v, w):
    nx, ny, nz = n
    ux, uy, uz = u
    vx, vy, vz = v
    wx, wy, wz = w
    ux -= nx
    uy -= ny
    uz -= nz
    vx -= nx
    vy -= ny
    vz -= nz
    wx -= nx
    wy -= ny
    wz -= nz
    return ux * (vy * wz - vz * wy) + uy * (vz * wx - vx * wz) + uz * (vx * wy - vy * wx)


class WithMixin:
    def __init__(self, file, mode='r'):
        if mode not in ('r', 'w', 'rb'):
            raise ValueError('invalid mode')
        if not file:
            raise ValueError('invalid file')

        if isinstance(file, str):
            self._file = open(file, mode)
            self.__is_buffer = False
        elif isinstance(file, Path):
            self._file = file.open(mode)
            self.__is_buffer = False
        elif isinstance(file, StringIO) and mode in 'rw':
            self._file = file
        elif isinstance(file, TextIOWrapper) and mode == 'r':
            self._file = file
        elif isinstance(file, BytesIO) and mode == 'rb':
            self._file = file
        elif hasattr(file, 'read') and file.mode == mode:  # check if file is open(filename, mode)
            self._file = file
        else:
            raise TypeError('invalid file')
        self.__write = mode == 'w'

    def close(self, force=False):
        if self.__write:
            self.write = self.__write_adhoc
            self.__write = False

        if not self.__is_buffer or force:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    @staticmethod
    def __write_adhoc(_):
        raise ValueError('I/O operation on closed writer')

    __is_buffer = True


class CGRread:
    def __init__(self, remap=True, ignore=False, colors=None):
        self.__remap = remap
        self._ignore = ignore
        self.__colors = colors

    def _convert_reaction(self, reaction):
        if not (reaction['reagents'] or reaction['products'] or reaction['reactants']):
            raise ValueError('empty reaction')
        maps = {'reagents': [], 'products': [], 'reactants': []}
        for i, tmp in maps.items():
            for molecule in reaction[i]:
                used = set()
                for atom in molecule['atoms']:
                    m = atom['mapping']
                    if m:
                        if m in used:
                            if not self._ignore:
                                raise MappingError('mapping in molecules should be unique')
                            warning(f'non-unique mapping in molecule: {m}')
                        else:
                            used.add(m)
                    tmp.append(m)

        length = count(max(max(maps['products'], default=0), max(maps['reagents'], default=0),
                           max(maps['reactants'], default=0)) + 1)

        ''' map unmapped atoms.
        '''
        for i, tmp in maps.items():
            used = set()
            maps[i] = remap = []
            for m in tmp:
                if not m:
                    remap.append(next(length))
                elif m in used:
                    if not self._ignore:
                        raise MappingError('mapping in reagents or products or reactants should be unique')
                    # force remap non unique atoms in molecules.
                    remap.append(next(length))
                    warning(f'mapping changed: {m} to {remap[-1]}')
                else:
                    remap.append(m)
                    used.add(m)

        if maps['reactants']:
            tmp = (set(maps['reagents']) | set(maps['products'])) & set(maps['reactants'])
            if tmp:
                if not self._ignore:
                    raise MappingError(f'reactants has map intersection with reagents or products: {tmp}')
                warning(f'reactants has map intersection with reagents or products: {tmp}')
                maps['reactants'] = [x if x not in tmp else next(length) for x in maps['reactants']]

        ''' find breaks in map. e.g. 1,2,5,6. 3,4 - skipped
        '''
        if self.__remap:
            lose = sorted(set(range(1, next(length))) - set(maps['reagents']) - set(maps['products']) -
                          set(maps['reactants']), reverse=True)
            if lose:
                for i, tmp in maps.items():
                    if not tmp:
                        continue
                    for j in lose:
                        maps[i] = [x if x < j else x - 1 for x in tmp]
        ''' end
        '''
        rc = ReactionContainer(meta=reaction['meta'])
        for i, tmp in maps.items():
            shift = 0
            for j in reaction[i]:
                atom_len = len(j['atoms'])
                remapped = {x: y for x, y in enumerate(tmp[shift: atom_len + shift])}
                shift += atom_len
                g = self.__convert_structure(j, remapped, {})
                rc[i].append(g)
        return rc

    def _convert_structure(self, molecule):
        if self.__remap:
            remapped = {n: k for n, k in enumerate(range(1, len(molecule['atoms']) + 1))}
        else:
            length = count(max(x['mapping'] for x in molecule['atoms']) + 1)
            remapped, used = {}, set()
            for n, atom in enumerate(molecule['atoms']):
                m = atom['mapping']
                if not m:
                    remapped[n] = next(length)
                elif m in used:
                    if not self._ignore:
                        raise MappingError('mapping in molecules should be unique')
                    remapped[n] = next(length)
                    warning(f'mapping in molecule changed: {m} to {remapped[n]}')
                else:
                    remapped[n] = m
                    used.add(m)

        if self.__colors and self.__colors in molecule['meta']:
            colors = self.__parse_colors(molecule['meta'][self.__colors])
        else:
            colors = {}
        g = self.__convert_structure(molecule, remapped, colors)
        g.meta.update(molecule['meta'])
        return g

    @staticmethod
    def __parse_dyn(value):
        value = [x.split('>') for x in value.split(',')]
        if len(value) > 1:
            r_value, p_value = [], []
            for r, p in value:
                r_value.append(int(r))
                p_value.append(int(p))
        else:
            r_value = int(value[0][0])
            p_value = int(value[0][1])
        return r_value, p_value

    @classmethod
    def __parse_cgr_atom(cls, value, base, name):
        p_name = f'p_{name}'
        diff = [base + int(x) for x in value[1:].split(',')]
        if len(diff) > 1:
            if diff[0] == base:  # for list of charges c0,1,-2,+3...
                return {name: diff}, True, False
            elif len(diff) % 2 == 1:  # for dyn charges c1,-1,0... -1,0 is dyn group relatively to base
                return {name: [base, *diff[1::2]], p_name: diff[::2]}, True, True
        else:
            return {name: base, p_name: diff[0]}, False, True

    @staticmethod
    def __parse_colors(colors):
        adhoc, before, res = [], [], defaultdict(dict)
        for x in colors:
            if (len(x) == 81 or len(x) == 75 and not before) and x[-1] == '+':
                before.append(x[:-1])
            else:
                before.append(x)
                adhoc.append(''.join(before).split())
                before = []

        for population, *keys in adhoc:
            p = int(population)
            for x in keys:
                atom, val = x.split(':')
                a = int(atom)
                res[a][p] = val
        return res

    def __convert_structure(self, molecule, mapping, colors):
        atom_data, bond_data, stereo = defaultdict(dict), defaultdict(dict), []
        atoms = molecule['atoms']
        bonds = molecule['bonds']
        is_cgr = False
        is_query = any(x['element'] in ('A', '*') for x in atoms)

        for e_atom, e_type, e_value in molecule['extra']:
            if e_type == 'atomlist':
                if set(e_value) - elements_set:
                    raise ValueError('atom list contain invalid atoms')
                if not is_query:
                    is_query = True
                atom_data[e_atom]['element'] = e_value
            else:
                if set(e_value) - elements_set:
                    raise ValueError('atom not_list contain invalid atoms')
                if not is_query:
                    is_query = True
                atom_data[e_atom]['element'] = elements_set - set(e_value)

        for c_atoms, c_type, c_value in molecule['cgr']:
            n = c_atoms[0]
            if c_type == 'dynatom':
                key = c_value[0]
                if key == 'x':
                    x, y, z = (float(x) for x in c_value[1:].split(','))
                    atom_data[n].update(p_x=x, p_y=y, p_z=z)
                    if not is_cgr:
                        is_cgr = True
                elif key == 'c':
                    value, iq, ic = self.__parse_cgr_atom(c_value, atoms[n]['charge'], 'charge')
                    atom_data[n].update(value)
                    if iq and not is_query:
                        is_query = True
                    if ic and not is_cgr:
                        is_cgr = True
                elif key == 'r':
                    value, iq, ic = self.__parse_cgr_atom(c_value, atoms[n]['multiplicity'] or 0, 'multiplicity')
                    atom_data[n].update(value)
                    if iq and not is_query:
                        is_query = True
                    if ic and not is_cgr:
                        is_cgr = True
                else:
                    raise ValueError('unknown dynatom')
            elif c_type == 'extrabond':
                value = [self.__bondlabels[x] for x in c_value.split(',')]
                if len(value) > 1:
                    if not is_query:
                        is_query = True
                    value = {'order': value}
                else:
                    value = {'order': value[0]}
                m = c_atoms[1]
                bond_data[n][m] = bond_data[m][n] = value
            elif c_type == 'dynbond':
                value = [x.split('>') for x in c_value.split(',')]
                if not is_cgr:
                    is_cgr = True
                if len(value) > 1:
                    if not is_query:
                        is_query = True
                    r_value, p_value = [], []
                    for r, p in value:
                        r_value.append(self.__bondlabels[r])
                        p_value.append(self.__bondlabels[p])
                    value = {'order': r_value, 'p_order': p_value}
                else:
                    value = {'order': self.__bondlabels[value[0][0]], 'p_order': self.__bondlabels[value[0][1]]}
                m = c_atoms[1]
                bond_data[n][m] = bond_data[m][n] = value
            elif c_type == 'isotope':
                atom_data[n]['isotope'] = [int(x) for x in c_value.split(',')]
                if not is_query:
                    is_query = True
            elif c_type == 'atomhyb':
                atom_data[n]['hybridization'] = [int(x) for x in c_value.split(',')]
                if not is_query:
                    is_query = True
            elif c_type == 'atomneighbors':
                atom_data[n]['neighbors'] = [int(x) for x in c_value.split(',')]
                if not is_query:
                    is_query = True
            elif c_type == 'dynatomhyb':
                atom_data[n]['hybridization'], atom_data[c_atoms[1]]['p_hybridization'] = \
                    self.__parse_dyn(c_value)
                if not is_cgr:
                    is_cgr = True
                if not is_query:
                    is_query = True
            elif c_type == 'dynatomneighbors':
                atom_data[n]['neighbors'], atom_data[c_atoms[1]]['p_neighbors'] = self.__parse_dyn(c_value)
                if not is_cgr:
                    is_cgr = True
                if not is_query:
                    is_query = True

        if is_cgr:
            for k, v in atom_data.items():
                atoms[k].update(v)
            for atom in atoms:
                if 'p_charge' not in atom:
                    atom['p_charge'] = atom['charge']
                if 'p_multiplicity' not in atom:
                    atom['p_multiplicity'] = atom['multiplicity']
                if 'p_x' not in atom:
                    atom.update(p_x=atom['x'], p_y=atom['y'], p_z=atom['z'])

            for n, m, bond, _ in bonds:
                if n in bond_data and m in bond_data[n]:
                    if bond['order'] != 8:
                        raise ValueError('invalid CGR spec')
                    bond.update(bond_data[n][m])
                if 'p_order' not in bond:
                    bond['p_order'] = bond['order']

            if is_query:
                for atom in atoms:
                    del atom['mark']
                    del atom['mapping']
                for k, v in atom_data.items():
                    if 'hybridization' in v and 'p_hybridization' not in v:
                        atoms[k]['p_hybridization'] = v['hybridization']
                    if 'neighbors' in v and 'p_neighbors' not in v:
                        atoms[k]['p_neighbors'] = v['neighbors']
                g = QueryCGRContainer()
            else:
                g = CGRContainer()
        elif is_query:
            for k, v in atom_data.items():
                atoms[k].update(v)
            for atom in atoms:
                del atom['mark']
                del atom['mapping']
            if bond_data:
                for n, m, bond, _ in bonds:
                    if n in bond_data and m in bond_data[n]:
                        if bond['order'] != 8:
                            raise ValueError('invalid CGR spec')
                        bond.update(bond_data[n][m])
            g = QueryContainer()
        else:
            if not any(x['z'] for x in atoms) and any(x['x'] for x in atoms) and any(x['y'] for x in atoms):
                # 0d, 1d and 3d stereo unsupported
                for n, m, bond, s in bonds:
                    if s:
                        if bond['order'] not in (1, 4):
                            raise ValueError('invalid wedge stereo spec')
                        stereo.append((mapping[n], mapping[m], s))
            g = MoleculeContainer()

        for n, atom in enumerate(atoms):
            element = atom['element']
            if element == 'D':
                atom['element'] = 'H'
                atom['isotope'] = 2
            elif element == 'T':
                atom['element'] = 'H'
                atom['isotope'] = 3
            elif element == '*':
                atom['element'] = 'A'
            if n in colors:
                atom['color'] = colors[n]
            g.add_atom(atom, mapping[n])

        for n, m, b, _ in bonds:
            n_map, m_map = mapping[n], mapping[m]
            g.add_bond(n_map, m_map, b)
        for n, m, s in stereo:
            self.__add_stereo(g, n, m, s)
        return g

    @classmethod
    def __add_stereo(cls, g, a1, a2, mark):
        adj = g._adj[a1]
        atom = g._node[a1]
        implicit = g.atom_implicit_h(a1)
        total = implicit + len(adj)

        if total == 4:  # tetrahedron
            order = ((g._node[x], mark if x == a2 else 0) for x in adj)
            if implicit:
                vol = _pyramid_volume((atom.x, atom.y, 0), *((atom.x, atom.y, z) for atom, z in order))
            else:
                vol = _pyramid_volume(*((atom.x, atom.y, z) for atom, z in order))
            atom.stereo = 1 if vol > 0 else -1
        else:
            warning('unsupported stereo or stereo impossible. tetrahedron only supported')

    __marks = {mark: ('s_%s' % mark, 'p_%s' % mark, 'sp_%s' % mark) for mark in ('neighbors', 'hyb', 'bond')}
    __bondlabels = {'0': None, '1': 1, '2': 2, '3': 3, '4': 4, '9': 9, 'n': None, 's': 9}
    __stereolabels = {'0': None, '1': 1, '-1': -1, '+1': 1, 'n': None, 1: 1, 6: -1, -1: -1}


class CGRwrite:
    def __init__(self, extralabels=False, mark_to_map=False, xyz=False, fix_position=True, colors='ISIDA_COLORS'):
        self.__xyz = xyz
        self.__mark_to_map = mark_to_map
        self.__extralabels = extralabels
        self._fix_position = fix_position
        self.__colors = colors

    def _convert_structure(self, g, shift=0):
        data = {'colors': {}}
        if isinstance(g, CGRContainer):
            format_atom = self.__format_dyn_atom
            format_bond = self.__format_dyn_bond
            stereo_map = {}
        elif isinstance(g, QueryCGRContainer):
            format_atom = self.__format_dyn_query_atom
            format_bond = self.__format_dyn_query_bond
            stereo_map = {}
        elif isinstance(g, QueryContainer):
            format_atom = self.__format_query_atom
            format_bond = self.__format_query_bond
            stereo_map = {}
        else:
            format_atom = self.__format_atom
            format_bond = self.__format_bond
            stereo_map = self.__wedge_map(g)

            #  colors supported only for molecules due to ambiguity for other types
            colors = defaultdict(list)
            for n, atom in enumerate(g._node.values(), start=1):
                if atom.color:
                    for part, val in atom.color.items():
                        colors[part].append(f'{n}:{val}')
            if colors:
                data['colors'][self.__colors] = '\n'.join('%s %s' % (x, ' '.join(y)) for x, y in colors.items())

        list_x, list_y = zip(*((x.x, x.y) for x in g._node.values()))
        min_x, max_x = min(list_x), max(list_x)
        data['y_shift'] = y_shift = -(max(list_y) + min(list_y)) / 2
        x_shift = shift - min_x
        if self._fix_position:
            data.update(max_x=max_x + x_shift, min_x=shift)
        else:
            data.update(max_x=max_x, min_x=min_x)

        cgr_data, atoms, bonds = [], [], []
        renum = {}
        for n, (i, atom) in enumerate(g._node.items(), start=1):
            renum[i] = n
            cgr, symbol, isotope, charge, multiplicity, mark, elements = format_atom(atom, n)
            cgr_data.extend(cgr)

            x, y = atom.x, atom.y
            if self._fix_position:
                x += x_shift
                y += y_shift
            atoms.append({'mapping': (atom.mark or 0) if self.__mark_to_map else i, 'symbol': symbol,
                          'isotope': isotope, 'charge': charge, 'multiplicity': multiplicity, 'mark': mark,
                          'x': x, 'y': y, 'z': atom.z, 'elements': elements, 'id': n})

        seen = set()
        for n, m_bond in g._adj.items():
            seen.add(n)
            for m, bond in m_bond.items():
                if m in seen:
                    continue
                if stereo_map:
                    stereo = stereo_map.get((n, m))
                    if stereo:
                        n_map, m_map = renum[n], renum[m]
                    else:
                        stereo = stereo_map.get((m, n))
                        if stereo:
                            n_map, m_map = renum[m], renum[n]
                        else:
                            n_map, m_map = renum[n], renum[m]
                else:
                    n_map, m_map = renum[n], renum[m]
                    stereo = None

                cgr, order = format_bond(bond, n_map, m_map)
                cgr_data.extend(cgr)
                bonds.append((n_map, m_map, order, self._stereo_map(stereo)))

        data['structure'] = (atoms, bonds, cgr_data)
        return data

    def __format_atom(self, atom, n):
        isotope = self._isotope_map(atom.isotope if atom.isotope != atom.common_isotope else None, n)
        return (), atom.element, isotope, self._charge_map(atom.charge), \
            self._multiplicity_map(atom.multiplicity, n), self._mark_map(atom.mark), self._atom_list_map(None, n)

    def __format_dyn_atom(self, atom, n):
        cgr = []
        nt = (n,)
        if atom.charge != atom.p_charge:
            cgr.append((nt, 'dynatom', f'c{atom.p_charge - atom.charge:+d}'))
        if atom.multiplicity != atom.p_multiplicity:
            cgr.append((nt, 'dynatom', f'r{(atom.p_multiplicity or 0) - (atom.multiplicity or 0):+d}'))

        if self.__xyz:
            dx, dy, dz = atom.p_x - atom.x, atom.p_y - atom.y, atom.p_z - atom.z
            if abs(dx) > .0001 or abs(dy) > .0001 or abs(dz) > .0001:
                cgr.append((nt, 'dynatom', f'x{dx:.4f},{dy:.4f},{dz:.4f}'))
        isotope = self._isotope_map(atom.isotope if atom.isotope != atom.common_isotope else None, n)
        return cgr, atom.element, isotope, self._charge_map(atom.charge), \
            self._multiplicity_map(atom.multiplicity, n), self._mark_map(atom.mark), self._atom_list_map(None, n)

    def __format_query_atom(self, atom, n):
        charge = atom.charge[0]
        cgr = []
        nt = (n,)
        if len(atom.element) > self._half_table:
            symbol = 'L'
            elements = self._atom_not_list_map(elements_set.difference(atom.element), n)
        elif len(atom.element) > 1:
            symbol = 'L'
            elements = self._atom_list_map(atom.element, n)
        else:
            symbol = atom.element[0]
            elements = self._atom_list_map(None, n)

        if len(atom.charge) > 1:
            cgr.append((nt, 'dynatom', 'c0,' + ','.join(f'{x - charge:+d}' for x in atom.charge[1:])))

        if len(atom.multiplicity) > 1:
            multiplicity = atom.multiplicity[0]
            cgr.append((nt, 'dynatom', 'r0,' + ','.join(f'{x - multiplicity:+d}' for x in atom.multiplicity[1:])))
        elif atom.multiplicity:
            multiplicity = atom.multiplicity[0]
        else:
            multiplicity = None

        if len(atom.isotope) > 1:
            isotope = None
            cgr.append((nt, 'isotope', ','.join(str(x) for x in atom.isotope)))
        elif atom.isotope:
            isotope = atom.isotope[0]
        else:
            isotope = None

        if atom.hybridization:
            cgr.append((nt, 'atomhyb', ','.join(str(x) for x in atom.hybridization)))
        if atom.neighbors:
            cgr.append((nt, 'atomneighbors', ','.join(str(x) for x in atom.neighbors)))
        return cgr, symbol, self._isotope_map(isotope, n), self._charge_map(charge), \
            self._multiplicity_map(multiplicity, n), self._mark_map(None), elements

    def __format_dyn_query_atom(self, atom, n):
        charge = atom.charge[0]
        cgr = []
        nt = (n,)
        if len(atom.element) > self._half_table:
            symbol = 'L'
            elements = self._atom_not_list_map(elements_set.difference(atom.element), n)
        elif len(atom.element) > 1:
            symbol = 'L'
            elements = self._atom_list_map(atom.element, n)
        else:
            symbol = atom.element[0]
            elements = self._atom_list_map(None, n)

        if atom.charge != atom.p_charge:
            if len(atom.charge) > 1:
                cgr.append((nt, 'dynatom', f'c{atom.p_charge[0] - charge:+d},' +
                                           ','.join(f'{x - charge:+d},{y - charge:+d}'
                                                    for x, y in zip(atom.charge[1:], atom.p_charge[1:]))))
            else:
                cgr.append((nt, 'dynatom', f'c{atom.p_charge[0] - charge:+d},'))
        elif len(atom.charge) > 1:
            cgr.append((nt, 'dynatom', 'c0,' + ','.join(f'{x - charge:+d}' for x in atom.charge[1:])))

        if atom.multiplicity != atom.p_multiplicity:
            multiplicity = atom.multiplicity[0]
            if len(atom.multiplicity) > 1:
                cgr.append((nt, 'dynatom', f'r{atom.p_multiplicity[0] - multiplicity:+d},' +
                                           ','.join(f'{x - multiplicity:+d},{y - multiplicity:+d}'
                                                    for x, y in zip(atom.multiplicity[1:], atom.p_multiplicity[1:]))))
            else:
                cgr.append((nt, 'dynatom', f'r{atom.p_multiplicity[0] - multiplicity:+d}'))
        elif len(atom.multiplicity) > 1:
            multiplicity = atom.multiplicity[0]
            cgr.append((nt, 'dynatom', 'r0,' + ','.join(f'{x - multiplicity:+d}' for x in atom.multiplicity[1:])))
        elif atom.multiplicity:
            multiplicity = atom.multiplicity[0]
        else:
            multiplicity = None

        if len(atom.isotope) > 1:
            isotope = None
            cgr.append((nt, 'isotope', ','.join(str(x) for x in atom.isotope)))
        elif atom.isotope:
            isotope = atom.isotope[0]
        else:
            isotope = None

        if atom.hybridization != atom.p_hybridization:
            cgr.append((nt, 'dynatomhyb',
                        ','.join(f'{x}>{y}' for x, y in zip(atom.hybridization, atom.p_hybridization))))
        elif atom.hybridization:
            cgr.append((nt, 'atomhyb', ','.join(str(x) for x in atom.hybridization)))

        if atom.neighbors != atom.p_neighbors:
            cgr.append((nt, 'dynatomneighbors', ','.join(f'{x}>{y}' for x, y in zip(atom.neighbors, atom.p_neighbors))))
        elif atom.neighbors:
            cgr.append((nt, 'atomneighbors', ','.join(str(x) for x in atom.neighbors)))

        if self.__xyz:
            dx, dy, dz = atom.p_x - atom.x, atom.p_y - atom.y, atom.p_z - atom.z
            if abs(dx) > .0001 or abs(dy) > .0001 or abs(dz) > .0001:
                cgr.append((nt, 'dynatom', f'x{dx:.4f},{dy:.4f},{dz:.4f}'))
        return cgr, symbol, self._isotope_map(isotope, n), self._charge_map(charge), \
            self._multiplicity_map(multiplicity, n), self._mark_map(None), elements

    def __format_bond(self, bond, n, m):
        if bond.order == 9:
            return (((n, m), 'extrabond', self._bond_map(9)),), self._bond_map(8)
        return (), self._bond_map(bond.order)

    def __format_dyn_bond(self, bond, n, m):
        if bond.order != bond.p_order:
            return (((n, m), 'dynbond',
                     f'{self._bond_map(bond.order)}>{self._bond_map(bond.p_order)}'),), self._bond_map(8)
        elif bond.order == 9:
            return (((n, m), 'extrabond', self._bond_map(9)),), self._bond_map(8)
        return (), self._bond_map(bond.order)

    def __format_query_bond(self, bond, n, m):
        if len(bond.order) > 1:
            return (((n, m), 'extrabond', ','.join(self._bond_map(x) for x in bond.order)),), self._bond_map(8)
        elif bond.order == (9,):
            return (((n, m), 'extrabond', self._bond_map(9)),), self._bond_map(8)
        return (), self._bond_map(bond.order[0])

    def __format_dyn_query_bond(self, bond, n, m):
        if bond.order != bond.p_order:
            return (((n, m), 'dynbond', ','.join(f'{self._bond_map(x)}>{self._bond_map(y)}'
                                                 for x, y in zip(bond.order, bond.p_order))),), self._bond_map(8)
        elif len(bond.order) > 1:
            return (((n, m), 'extrabond', ','.join(self._bond_map(x) for x in bond.order)),), self._bond_map(8)
        elif bond.order == (9,):
            return (((n, m), 'extrabond', self._bond_map(9)),), self._bond_map(8)
        return (), self._bond_map(bond.order[0])

    @staticmethod
    def __wedge_map(g):
        stereo_map = {}
        nodes = list(g._node.items())
        while True:
            failed = []
            for i, tmp in enumerate(nodes, start=1):
                n, atom = tmp
                s = atom.stereo
                if not s:
                    continue
                neighbors = list(g._adj[n])
                len_n = len(neighbors)
                if len_n in (3, 4):  # tetrahedron
                    for _ in range(len_n):
                        if (neighbors[0], n) in stereo_map:
                            neighbors.append(neighbors.pop(0))
                        else:
                            failed.append(tmp)
                            break
                    else:
                        failed.insert(0, tmp)
                        failed.extend(nodes[i:])
                        stereo_map = None
                        nodes = failed
                        break

                    if len_n == 4:
                        zero = g._node[neighbors[0]]
                        zero = (zero.x, zero.y, 1)
                        first = g._node[neighbors[1]]
                        first = (first.x, first.y, 0)
                    else:
                        zero = (atom.x, atom.y, 0)
                        first = g._node[neighbors[0]]
                        first = (first.x, first.y, 1)

                    second = g._node[neighbors[-2]]
                    third = g._node[neighbors[-1]]
                    vol = _pyramid_volume(zero, first, (second.x, second.y, 0), (third.x, third.y, 0))

                    stereo_map[(n, neighbors[0])] = 1 if vol > 0 and s == 1 or vol < 0 and s == -1 else -1
            else:
                return stereo_map

    @staticmethod
    @abstractmethod
    def _stereo_map(x):
        pass

    @staticmethod
    @abstractmethod
    def _charge_map(x):
        pass

    @staticmethod
    @abstractmethod
    def _isotope_map(x, n):
        pass

    @staticmethod
    @abstractmethod
    def _mark_map(x):
        pass

    @staticmethod
    @abstractmethod
    def _multiplicity_map(x, n):
        pass

    @staticmethod
    @abstractmethod
    def _bond_map(x):
        pass

    @staticmethod
    @abstractmethod
    def _atom_list_map(x, n):
        pass

    @staticmethod
    @abstractmethod
    def _atom_not_list_map(x, n):
        pass

    _half_table = len(elements_set) // 2
