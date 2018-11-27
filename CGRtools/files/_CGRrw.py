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
from itertools import count, chain
from io import StringIO, BytesIO, TextIOWrapper
from pathlib import Path
from warnings import warn
from ..containers import ReactionContainer, MoleculeContainer, CGRContainer, QueryContainer, QueryCGRContainer
from ..exceptions import InvalidStereo, InvalidAtom, InvalidConfig, MapError, InvalidData
from ..periodictable import elements_set
from ..periodictable import common_isotopes


fromMDL = (0, 3, 2, 1, 0, -1, -2, -3)


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
            raise InvalidConfig('invalid mode')
        if not file:
            raise InvalidConfig('invalid file')
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
            raise InvalidConfig('invalid file')
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
    def __init__(self, remap=True, ignore=False, is_template=False, colors=None):
        self.__remap = remap
        self.__ignore = ignore
        self.__is_template = is_template
        self.__colors = colors

    def _get_reaction(self, reaction):
        spmaps = None
        maps = {}
        for i in ('reagents', 'products', 'reactants'):
            maps[i] = tmp = []
            for m in reaction[i]:
                mm = []
                for a in m['atoms']:
                    am = a['map']
                    if not self.__ignore and am and am in mm:
                        raise MapError('mapping in molecules should be unique')
                    else:
                        mm.append(am)
                    tmp.append(am)

        length = count(max(max(maps['products'], default=0), max(maps['reagents'], default=0),
                           max(maps['reactants'], default=0)) + 1)

        ''' map unmapped atoms.
        '''
        for i in ('reagents', 'products', 'reactants'):
            mi, used = [], []
            for m in maps[i]:
                k = m or next(length)
                if k in used:
                    if not self.__ignore:
                        raise MapError('mapping in reagents or products or reactants should be unique')
                    # force remap non unique atoms in molecules.
                    mi.append((next(length), k))
                else:
                    mi.append((k, k))
                    used.append(k)

            maps[i] = mi

        if maps['reactants']:
            spmaps = set(x for x, _ in chain(maps['reagents'], maps['products']))
            d = set(x for x, _ in maps['reactants']).intersection(spmaps)
            if d:
                if not self.__ignore:
                    raise MapError('reactants has map intersection with reagents or products')
                maps['reactants'] = [(x if x not in d else next(length), y) for x, y in maps['reactants']]

        ''' find breaks in map. e.g. 1,2,5,6. 3,4 - skipped
        '''
        if self.__remap:
            lose = sorted(set(range(1, next(length)))
                          .difference(spmaps or (x for x, _ in chain(maps['reagents'], maps['products'])))
                          .difference(x for x, _ in maps['reactants']), reverse=True)
            if lose:
                for i in ('reagents', 'products', 'reactants'):
                    if not maps[i]:
                        continue
                    for j in lose:
                        maps[i] = [(x if x < j else x - 1, y) for x, y in maps[i]]
        ''' end
        '''
        rc = ReactionContainer(meta={x: '\n'.join(y) for x, y in reaction['meta'].items()})

        colors = {}
        if self.__colors:
            for k, v in reaction['meta'].items():
                if k.startswith(self.__colors):
                    try:
                        mol_num = k.split('.', 1)[1]
                        colors[int(mol_num)] = v
                    except (ValueError, IndexError):
                        warn(f'coloring data invalid: {k}, {v}', ResourceWarning)

        counter = count(1)
        for i in ('reagents', 'products', 'reactants'):
            shift = 0
            for j in reaction[i]:
                atom_len = len(j['atoms'])
                remapped = {x: y for x, y in enumerate(maps[i][shift: atom_len + shift], start=1)}
                shift += atom_len
                g = self.__parse_molecule(j, remapped, colors=colors.get(next(counter)))
                rc[i].append(g)
        return rc

    def _get_molecule(self, molecule):
        if self.__remap:
            remapped = {k: (k, k) for k in range(1, len(molecule['atoms']) + 1)}
        else:
            length = count(max(x['map'] for x in molecule['atoms']) + 1)
            remapped, used = {}, []
            for k, l in enumerate(molecule['atoms'], start=1):
                m = l['map'] or next(length)
                if m in used:
                    if not self.__ignore:
                        raise MapError('mapping in molecules should be unique')

                    remapped[k] = (next(length), m)
                else:
                    remapped[k] = (m, m)
                    used.append(m)

        colors = molecule['meta'].get(self.__colors) if self.__colors else None
        return self.__parse_molecule(molecule, remapped, meta=molecule['meta'], colors=colors)

    @classmethod
    def __parsedyn(cls, name, value):
        tmp0, tmp1 = [], []
        for x in value.split(','):
            s, *_, p = x.split('>')
            if s != p:
                if name == 'order':
                    tmp0.append(cls.__bondlabels[s])
                    tmp1.append(cls.__bondlabels[p])
                else:
                    tmp0.append(s != 'n' and int(s) or None)
                    tmp1.append(p != 'n' and int(p) or None)

        p_name = f'p_{name}'
        if len(tmp0) == 1:
            out = {}
            if tmp0[0] is not None:
                out[name] = tmp0[0]
            if tmp1[0] is not None:
                out[p_name] = tmp1[0]
            return out, False
        elif len(tmp0) > 1:
            tmp2 = list(zip(tmp0, tmp1))
            if len(set(tmp2)) == len(tmp2):
                return {name: tmp0, p_name: tmp1}, True
        warn(f'CGR data invalid: {name}, {value}', ResourceWarning)

    @classmethod
    def __parselist(cls, name, value):
        if name == 'order':
            tmp = [cls.__bondlabels[x] for x in value.split(',')]
        else:
            tmp = [(x != 'n' and int(x) or None) for x in value.split(',')]

        p_name = f'p_{name}'
        if len(tmp) == 1:
            if tmp[0] is not None:
                return {name: tmp[0], p_name: tmp[0]}, False
        elif len(set(tmp)) == len(tmp) and None not in tmp:
            return {name: tmp, p_name: tmp}, True
        warn(f'CGR data invalid: {name}, {value}', ResourceWarning)

    @classmethod
    def __parse_cgr_atom(cls, value, base, name):
        p_name = f'p_{name}'
        diff = [int(x) for x in value]
        if len(diff) > 1:
            if diff[0] == 0:  # for list of charges c0,1,-2,+3...
                tmp = [base + x for x in diff]
                if len(set(tmp)) == len(tmp):
                    return {name: tmp, p_name: tmp}, True
            elif len(diff) % 2 == 1:  # for dyn charges c1,-1,0... -1,0 is dyn group relatively to base
                s = [base] + [base + x for x in diff[1::2]]
                p = [base + x for x in diff[::2]]
                tmp = list(zip(s, p))
                if len(set(tmp)) == len(tmp):
                    return {name: s, p_name: p}, True
        else:
            tmp = diff[0]
            if tmp:
                p = base + tmp
                return {name: base, p_name: p}, False

    @classmethod
    def __cgr_dat(cls, _type, value):
        if _type == 'atomlist':
            if set(value).difference(elements_set):
                raise InvalidAtom('atom list contain invalid atoms')
            return dict(element=value), True

        elif _type == 'atomnotlist':
            if set(value).difference(elements_set):
                raise InvalidAtom('atom not list contain invalid atoms')
            return dict(element=list(elements_set.difference(value))), True

        elif _type == 'atomhyb':
            return cls.__parselist('hyb', value)

        elif _type == 'atomneighbors':
            return cls.__parselist('neighbors', value)

        elif _type == 'dynatomhyb':
            return cls.__parsedyn('hyb', value)

        elif _type == 'dynatomneighbors':
            return cls.__parsedyn('neighbors', value)

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

    def __parse_molecule(self, molecule, mapping, meta=None, colors=None):
        meta = meta and {x: '\n'.join(y) for x, y in meta.items()} or None
        cgr_dat_atom, cgr_dat_bond, cgr_dat_stereo, isotope_dat = defaultdict(dict), defaultdict(dict), {}, {}
        any_bonds, normal_stereo = [], []
        charge_dat = {k['atoms'][0]: int(k['value']) for k in molecule['CGR_DAT'] if k['type'] == 'charge'}
        radical_dat = {k['atoms'][0]: int(k['value']) for k in molecule['CGR_DAT'] if k['type'] == 'radical'}
        colors = self.__parse_colors(colors) if colors else {}
        is_query = self.__is_template

        for k in molecule['CGR_DAT']:
            k_atoms = k['atoms']
            k_type = k['type']
            if len(k_atoms) != self._cgr_keys.get(k_type):
                raise InvalidData('CGR data invalid. invalid number of atoms')
            if not k_atoms[0] or len(k_atoms) == 2 and not k_atoms[1]:
                raise InvalidData('CGR data invalid. contains zero atoms')

            if k_type in ('charge', 'radical'):
                continue
            elif k_type == 'dynatom':
                key = k['value'][0]
                value = k['value'][1:].split(',')
                a1 = k_atoms[0]
                if key == 'x':
                    x, y, z = (float(x) for x in value)
                    cgr_dat_atom[a1].update(p_x=x, p_y=y, p_z=z)
                else:
                    if key == 'c':
                        val = self.__parse_cgr_atom(value, charge_dat.get(a1, molecule['atoms'][a1 - 1]['charge']),
                                                    'charge')
                        if val:
                            cgr_dat_atom[a1].update(val[0])
                            if val[1] and not is_query:
                                is_query = True
                        else:
                            warn(f'dynatom charge invalid data. skipped: {value}', ResourceWarning)
                    elif key == 'r':
                        val = self.__parse_cgr_atom(value, radical_dat.get(a1, 0), 'multiplicity')
                        if val:
                            cgr_dat_atom[a1].update(val[0])
                            if val[1] and not is_query:
                                is_query = True
                        else:
                            warn(f'dynatom radical invalid data. skipped: {value}', ResourceWarning)
            elif k_type == 'dynstereo':
                s_stereo, p_stereo = k['value'].split('>')
                cgr_dat_stereo[k_atoms] = (self.__stereolabels[s_stereo], self.__stereolabels[p_stereo])
            elif k_type == 'extrabond':
                val = self.__parselist('order', k['value'])
                if val:
                    a1, a2 = k_atoms
                    cgr_dat_bond[a1][a2] = cgr_dat_bond[a2][a1] = val[0]
                    if val[1] and not is_query:
                        is_query = True
            elif k_type == 'dynbond':
                val = self.__parsedyn('order', k['value'])
                if val:
                    a1, a2 = k_atoms
                    cgr_dat_bond[a1][a2] = cgr_dat_bond[a2][a1] = val[0]
                    if val[1] and not is_query:
                        is_query = True
            elif k_type == 'isotope':
                tmp = k['value'].split(',')
                if len(tmp) > 1:
                    cgr_dat_atom[k_atoms[0]]['isotope'] = [int(x) for x in tmp]
                    if not is_query:
                        is_query = True
                else:
                    isotope_dat[k_atoms[0]] = int(tmp[0])
            else:
                val = self.__cgr_dat(k_type, k['value'])
                if val:
                    cgr_dat_atom[k_atoms[0]].update(val[0])
                    if val[1] and not is_query:
                        is_query = True

        if not is_query:
            is_query = any(x['element'] in ('A', '*') for x in molecule['atoms'])
        is_cgr = is_query or bool(cgr_dat_atom or cgr_dat_bond or cgr_dat_stereo)

        if is_query:
            g = QueryContainer()
        elif is_cgr:
            g = CGRContainer()
        else:
            g = MoleculeContainer()
        if meta:
            g.meta.update(meta)

        for k, l in enumerate(molecule['atoms'], start=1):
            atom_map, parsed_map = mapping[k]
            element = l['element']
            if element == 'D':
                element = 'H'
                isotope_dat[k] = 2
            elif element == 'T':
                element = 'H'
                isotope_dat[k] = 3

            if is_cgr:
                if k in cgr_dat_atom:
                    atom_dat = cgr_dat_atom[k]
                    if 'element' not in atom_dat:
                        atom_dat['element'] = element
                    if 'charge' not in atom_dat:
                        lc = charge_dat.get(k, l['charge'])
                        atom_dat.update(charge=lc, p_charge=lc)
                    if 'multiplicity' not in atom_dat:
                        lr = radical_dat.get(k)
                        atom_dat.update(multiplicity=lr, p_multiplicity=lr)
                    if 'p_x' in atom_dat:
                        atom_dat['p_x'] += l['x']
                        atom_dat['p_y'] += l['y']
                        atom_dat['p_z'] += l['z']
                    else:
                        atom_dat.update(p_x=l['x'], p_y=l['y'], p_z=l['z'])
                    if k in colors:
                        atom_dat['color'] = colors[k]
                    if 'isotope' not in atom_dat:
                        if k in isotope_dat:
                            atom_dat['isotope'] = isotope_dat[k]
                        elif l['isotope']:
                            atom_dat['isotope'] = common_isotopes[element] + l['isotope']

                    g.add_node(atom_map, mark=l['mark'], mapping=parsed_map, x=l['x'], y=l['y'], z=l['z'], **atom_dat)
                else:
                    lc = charge_dat.get(k, l['charge'])
                    lr = radical_dat.get(k)
                    g.add_node(atom_map, element=element, mark=l['mark'], mapping=parsed_map,
                               isotope=isotope_dat.get(k) or common_isotopes[element] + l['isotope'],
                               charge=lc, p_charge=lc, multiplicity=lr, p_multiplicity=lr,
                               x=l['x'], y=l['y'], z=l['z'], p_x=l['x'], p_y=l['y'], p_z=l['z'], color=colors.get(k))
            else:
                g.add_node(atom_map, element=element, mark=l['mark'], mapping=parsed_map,
                           isotope=isotope_dat.get(k) or common_isotopes[element] + l['isotope'],
                           color=colors.get(k), charge=charge_dat.get(k, l['charge']),
                           multiplicity=radical_dat.get(k), x=l['x'], y=l['y'], z=l['z'])

        for k, l, m, s in molecule['bonds']:
            k_map, l_map = mapping[k][0], mapping[l][0]
            if is_cgr:
                if m == 8 and k in cgr_dat_bond and l in cgr_dat_bond[k]:  # dynbond only for any bond accept.
                    g.add_edge(k_map, l_map, **cgr_dat_bond[k][l])
                else:
                    g.add_edge(k_map, l_map, order=m, p_order=m)
                if m == 8:
                    any_bonds.append((k, l))
                    any_bonds.append((l, k))
            else:
                g.add_edge(k_map, l_map, order=m)

            if s in (1, 6) and m in (1, 4):
                s_mark = self.__stereolabels[s]
                normal_stereo.append((k_map, l_map, s_mark))

        while True:
            failed_cgr_dat_stereo = {}
            for (k, l), (s_mark, p_mark) in cgr_dat_stereo.items():
                if (k, l) in any_bonds:  # remove dynamic stereo if bond not any [8]. INVALID SPEC.
                    k_map, l_map = mapping[k][0], mapping[l][0]
                    try:
                        g.add_stereo(k_map, l_map, s_mark, p_mark)
                    except InvalidStereo as e:
                        warn(str(e), ResourceWarning)
                        failed_cgr_dat_stereo[(k, l)] = (s_mark, p_mark)
                    except InvalidAtom as e:
                        if not self.__ignore:
                            raise
                        warn(str(e), ResourceWarning)

            failed_normal_stereo = []
            for k_map, l_map, s_mark in normal_stereo:
                try:
                    if is_cgr:
                        g.add_stereo(k_map, l_map, s_mark, s_mark)
                    else:
                        g.add_stereo(k_map, l_map, s_mark)
                except InvalidStereo as e:
                    warn(str(e), ResourceWarning)
                    failed_normal_stereo.append((k_map, l_map, s_mark))
                except InvalidAtom as e:
                    if not self.__ignore:
                        raise
                    warn(str(e), ResourceWarning)

            if failed_cgr_dat_stereo and len(cgr_dat_stereo) > len(failed_cgr_dat_stereo) or \
                    failed_normal_stereo and len(normal_stereo) > len(failed_normal_stereo):
                cgr_dat_stereo = failed_cgr_dat_stereo
                normal_stereo = failed_normal_stereo
                continue
            break

        return g

    __marks = {mark: ('s_%s' % mark, 'p_%s' % mark, 'sp_%s' % mark) for mark in ('neighbors', 'hyb', 'bond')}
    __bondlabels = {'0': None, '1': 1, '2': 2, '3': 3, '4': 4, '9': 9, 'n': None, 's': 9}
    __stereolabels = {'0': None, '1': 1, '-1': -1, '+1': 1, 'n': None, 1: 1, 6: -1, -1: -1}
    _cgr_keys = dict(extrabond=2, dynbond=2, dynatom=1, atomnotlist=1, atomlist=1, isotope=1, charge=1, radical=1,
                     atomhyb=1, atomneighbors=1, dynatomhyb=1, dynatomneighbors=1, dynstereo=2)


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

        cgr_data, extra_data, atoms, bonds = [], [], [], []
        renum = {}
        for n, (i, atom) in enumerate(g._node.items(), start=1):
            renum[i] = n
            cgr, extra, charge, element = format_atom(atom, n)
            cgr_data.extend(cgr)
            extra_data.extend(extra)

            x, y = atom.x, atom.y
            if self._fix_position:
                x += x_shift
                y += y_shift
            atoms.append({'map': (atom.mark or 0) if self.__mark_to_map else i, 'charge': charge, 'element': element,
                          'mark': atom.mark or 0, 'x': x, 'y': y, 'z': atom.z})

        seen = set()
        for i, j_bond in g._adj.items():
            seen.add(i)
            for j, bond in j_bond.items():
                if j in seen:
                    continue

                stereo = stereo_map.get((i, j))
                if not stereo:
                    stereo = stereo_map.get((j, i))
                    if stereo:
                        i, j = j, i

                cgr, order = format_bond(bond, renum[i], renum[j])
                cgr_data.extend(cgr)
                bonds.append((renum[i], renum[j], order, self._stereo_map[stereo]))

        data['structure'] = (atoms, bonds, extra_data, cgr_data)
        return data

    def __format_atom(self, atom, n):
        extra = []
        if atom.multiplicity:
            extra.append((n, 'radical', self._radical_map[atom.multiplicity]))
        if atom.isotope != atom.common_isotope:
            extra.append((n, 'isotope', atom.isotope))
        return (), extra, self._charge_map[atom.charge], atom.element

    def __format_dyn_atom(self, atom, n):
        cgr = []
        extra = []
        nt = (n,)
        if atom.charge != atom.p_charge:
            cgr.append((nt, 'dynatom', f'c{atom.p_charge - atom.charge:+d}'))
        if atom.multiplicity != atom.p_multiplicity:
            cgr.append((nt, 'dynatom', f'r{atom.p_radical - atom.radical:+d}'))
        if atom.multiplicity:
            extra.append((n, 'radical', self._radical_map[atom.multiplicity]))
        if atom.isotope != atom.common_isotope:
            extra.append((n, 'isotope', atom.isotope))

        if self.__xyz:
            dx, dy, dz = atom.p_x - atom.x, atom.p_y - atom.y, atom.p_z - atom.z
            if abs(dx) > .0001 or abs(dy) > .0001 or abs(dz) > .0001:
                cgr.append((nt, 'dynatom', f'x{dx:.4f},{dy:.4f},{dz:.4f}'))
        return cgr, extra, self._charge_map[atom.charge], atom.element

    def __format_query_atom(self, atom, n):
        charge = atom.charge[0]
        cgr = []
        extra = []
        nt = (n,)
        if len(atom.element) > 1:
            extra.append((n, 'atomlist', atom.element))
            element = 'L'
        else:
            element = atom.element[0]

        if len(atom.charge) > 1:
            cgr.append((nt, 'dynatom', 'c0,' + ','.join(f'{x - charge:+d}' for x in atom.charge[1:])))
        if len(atom.multiplicity) > 1:
            radical = atom.radical[0]
            cgr.append((nt, 'dynatom', 'r0,' + ','.join(f'{x - radical:+d}' for x in atom.radical[1:])))
        if atom.multiplicity:
            extra.append((n, 'radical', self._radical_map[atom.multiplicity[0]]))

        if atom.isotope:
            extra.append((n, 'isotope', ','.join(str(x) for x in atom.isotope)))

        if atom.hybridization:
            cgr.append((nt, 'atomhyb', ','.join(str(x) for x in atom.hybridization)))
        if atom.neighbors:
            cgr.append((nt, 'atomneighbors', ','.join(str(x) for x in atom.neighbors)))
        return cgr, extra, self._charge_map[charge], element

    def __format_dyn_query_atom(self, atom, n):
        charge = atom.charge[0]
        cgr = []
        extra = []
        nt = (n,)
        if len(atom.element) > 1:
            extra.append(dict(atom=n, value=atom.element, type='atomlist'))
            element = 'L'
        else:
            element = atom.element[0]

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
            radical = atom.radical[0]
            if len(atom.multiplicity) > 1:
                cgr.append((nt, 'dynatom', f'r{atom.p_radical[0] - radical:+d},' +
                                           ','.join(f'{x - radical:+d},{y - radical:+d}'
                                                    for x, y in zip(atom.radical[1:], atom.p_radical[1:]))))
            else:
                cgr.append((nt, 'dynatom', f'r{atom.p_radical[0] - radical:+d}'))
        elif len(atom.multiplicity) > 1:
            radical = atom.radical[0]
            cgr.append((nt, 'dynatom', 'r0,' + ','.join(f'{x - radical:+d}' for x in atom.radical[1:])))
        if atom.multiplicity:
            extra.append((n, 'radical', self._radical_map[atom.multiplicity[0]]))

        if atom.isotope:
            extra.append((n, 'isotope', ','.join(str(x) for x in atom.isotope)))

        if atom.hybridization != atom.p_hybridization:
            cgr.append((nt, 'dynatomhyb',
                        ','.join(f'{x}>{y}' for x, y in zip(atom.hybridization, atom.p_hybridization))))
        elif atom.hybridization:
            cgr.append((nt, 'atomhyb', ','.join(str(x) for x in atom.hybridization)))

        if atom.neighbors != atom.p_neighbors:
            cgr.append((nt, 'dynatomneighbors', ','.join(f'{x}>{y}' for x, y in zip(atom.neighbors, atom.p_neighbors))))
        elif atom.neighbors:
            cgr.append((nt, 'atomneighbors', ','.join(str(x) for x in atom.neighbors)))
        return cgr, extra, self._charge_map[charge], element

    def __format_bond(self, bond, n, m):
        if bond.order == 9:
            return (((n, m), 'extrabond', self._bond_map[9]),), self._bond_map[8]
        return (), self._bond_map[bond.order]

    def __format_dyn_bond(self, bond, n, m):
        if bond.order != bond.p_order:
            return (((n, m), 'dynbond',
                     f'{self._bond_map[bond.order]}>{self._bond_map[bond.p_order]}'),), self._bond_map[8]
        elif bond.order == 9:
            return (((n, m), 'extrabond', self._bond_map[9]),), self._bond_map[8]
        return (), self._bond_map[bond.order]

    def __format_query_bond(self, bond, n, m):
        if len(bond.order) > 1:
            return (((n, m), 'extrabond', ','.join(self._bond_map[x] for x in bond.order)),), self._bond_map[8]
        elif bond.order == (9,):
            return (((n, m), 'extrabond', self._bond_map[9]),), self._bond_map[8]
        return (), self._bond_map[bond.order[0]]

    def __format_dyn_query_bond(self, bond, n, m):
        if bond.order != bond.p_order:
            return (((n, m), 'dynbond', ','.join(f'{self._bond_map[x]}>{self._bond_map[y]}'
                                                 for x, y in zip(bond.order, bond.p_order))),), self._bond_map[8]
        elif len(bond.order) > 1:
            return (((n, m), 'extrabond', ','.join(self._bond_map[x] for x in bond.order)),), self._bond_map[8]
        elif bond.order == (9,):
            return (((n, m), 'extrabond', self._bond_map[9]),), self._bond_map[8]
        return (), self._bond_map[bond.order[0]]

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

    @property
    @abstractmethod
    def _stereo_map(self):
        pass

    @property
    @abstractmethod
    def _charge_map(self):
        pass

    @property
    @abstractmethod
    def _radical_map(self):
        pass

    @property
    @abstractmethod
    def _bond_map(self):
        pass

    _half_table = len(elements_set) // 2
