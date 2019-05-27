# -*- coding: utf-8 -*-
#
#  Copyright 2014-2019 Ramil Nugmanov <stsouko@live.ru>
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
from io import StringIO, BytesIO, TextIOWrapper, BufferedIOBase, BufferedReader
from logging import warning
from pathlib import Path
from ..containers import ReactionContainer, MoleculeContainer, CGRContainer, QueryContainer, QueryCGRContainer
from ..exceptions import MappingError
from ..periodictable import elements_set


elements_set = elements_set.copy()
elements_set.discard('A')
cgr_keys = dict(extrabond=2, dynbond=2, dynatom=1, atomhyb=1, atomneighbors=1, dynatomhyb=1, dynatomneighbors=1)


class WithMixin:
    def __init__(self, file, mode='r'):
        if mode not in ('r', 'w', 'rb'):
            raise ValueError('invalid mode')
        if not file:
            raise ValueError('invalid file')

        if isinstance(file, str):
            self._file = open(file, mode)
            self._is_buffer = False
        elif isinstance(file, Path):
            self._file = file.open(mode)
            self._is_buffer = False
        elif isinstance(file, (TextIOWrapper, StringIO)) and mode in ('r', 'w'):
            self._file = file
        elif isinstance(file, (BytesIO, BufferedReader, BufferedIOBase)) and mode == 'rb':
            self._file = file
        else:
            raise TypeError('invalid file. '
                            'TextIOWrapper, StringIO, BytesIO, BufferedReader and BufferedIOBase subclasses possible')
        self.__write = mode == 'w'

    def close(self, force=False):
        """
        close opened file

        :param force: force closing of externally opened file or buffer
        """
        if self.__write:
            self.write = self.__write_adhoc
            self.__write = False

        if not self._is_buffer or force:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    @staticmethod
    def __write_adhoc(_):
        raise ValueError('I/O operation on closed writer')

    _is_buffer = True


class CGRread:
    def __init__(self, remap=True, ignore=False):
        self.__remap = remap
        self._ignore = ignore

    def _convert_reaction(self, reaction):
        if not (reaction['reactants'] or reaction['products'] or reaction['reagents']):
            raise ValueError('empty reaction')
        maps = {'reactants': [], 'products': [], 'reagents': []}
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

        length = count(max(max(maps['products'], default=0), max(maps['reactants'], default=0),
                           max(maps['reagents'], default=0)) + 1)

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

        if maps['reagents']:
            tmp = (set(maps['reactants']) | set(maps['products'])) & set(maps['reagents'])
            if tmp:
                e = f'reagents has map intersection with reactants or products: {tmp}'
                if not self._ignore:
                    raise MappingError(e)
                warning(e)
                maps['reagents'] = [x if x not in tmp else next(length) for x in maps['reagents']]

        ''' find breaks in map. e.g. 1,2,5,6. 3,4 - skipped
        '''
        if self.__remap:
            lose = sorted(set(range(1, next(length))) - set(maps['reactants']) - set(maps['products']) -
                          set(maps['reagents']), reverse=True)
            if lose:
                for i, tmp in maps.items():
                    if not tmp:
                        continue
                    for j in lose:
                        maps[i] = tmp = [x if x < j else x - 1 for x in tmp]
        ''' end
        '''
        rc = ReactionContainer(meta=reaction['meta'])
        for i, tmp in maps.items():
            shift = 0
            for j in reaction[i]:
                atom_len = len(j['atoms'])
                remapped = {x: y for x, y in enumerate(tmp[shift: atom_len + shift])}
                shift += atom_len
                g = self.__convert_structure(j, remapped)
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

        g = self.__convert_structure(molecule, remapped)
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

    def __convert_structure(self, molecule, mapping):
        atom_data, bond_data, stereo = defaultdict(dict), defaultdict(dict), []
        atoms = molecule['atoms']
        bonds = molecule['bonds']
        is_cgr = False
        is_query = any(x['element'] in ('A', '*') for x in atoms)

        for e_atom, e_type, e_value in molecule['extra']:
            if e_type == 'atomlist':
                if set(e_value) - elements_set:
                    raise ValueError('atom list contain invalid atoms')
                atom_data[e_atom]['element'] = e_value
                if not is_query:
                    is_query = True
            else:
                if set(e_value) - elements_set:
                    raise ValueError('atom not_list contain invalid atoms')
                atom_data[e_atom]['element'] = elements_set - set(e_value)
                if not is_query:
                    is_query = True

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
                    atom_data[n]['p_charge'] = atoms[n]['charge'] + int(c_value[1:])
                    if not is_cgr:
                        is_cgr = True
                elif key == 'r':
                    atom_data[n]['p_multiplicity'] = (atoms[n]['multiplicity'] or 0) + int(c_value[1:])
                    if not is_cgr:
                        is_cgr = True
                else:
                    raise ValueError('unknown dynatom')
            elif c_type == 'extrabond':
                value = {'order': self.__bondlabels[c_value]}
                m = c_atoms[1]
                bond_data[n][m] = bond_data[m][n] = value
            elif c_type == 'dynbond':
                bond, p_bond = c_value.split('>')
                value = {'order': self.__bondlabels[bond], 'p_order': self.__bondlabels[p_bond]}
                m = c_atoms[1]
                bond_data[n][m] = bond_data[m][n] = value
                if not is_cgr:
                    is_cgr = True
            elif c_type == 'atomhyb':
                atom_data[n]['hybridization'] = [int(x) for x in c_value.split(',')]
                if not is_query:
                    is_query = True
            elif c_type == 'atomneighbors':
                atom_data[n]['neighbors'] = [int(x) for x in c_value.split(',')]
                if not is_query:
                    is_query = True
            elif c_type == 'dynatomhyb':
                atom_data[n]['hybridization'], atom_data[n]['p_hybridization'] = self.__parse_dyn(c_value)
                if not is_cgr:
                    is_cgr = True
                if not is_query:
                    is_query = True
            elif c_type == 'dynatomneighbors':
                atom_data[n]['neighbors'], atom_data[n]['p_neighbors'] = self.__parse_dyn(c_value)
                if not is_cgr:
                    is_cgr = True
                if not is_query:
                    is_query = True

        prepared_bonds = []
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

            for n, m, bond in bonds:
                if n in bond_data and m in bond_data[n]:
                    if bond != 8:
                        raise ValueError('invalid CGR spec')
                    bond = bond_data[n][m]
                    if 'p_order' not in bond:
                        bond['p_order'] = bond['order']
                    prepared_bonds.append((n, m, bond))
                else:
                    prepared_bonds.append((n, m, {'order': bond, 'p_order': bond}))

            if is_query:
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
            for n, m, bond in bonds:
                if n in bond_data and m in bond_data[n]:
                    if bond != 8:
                        raise ValueError('invalid CGR spec')
                    prepared_bonds.append((n, m, bond_data[n][m]))
                else:
                    prepared_bonds.append((n, m, {'order': bond}))
            g = QueryContainer()
        else:
            for n, m, bond in bonds:
                prepared_bonds.append((n, m, {'order': bond}))
            g = MoleculeContainer()

        parsed_mapping = [x.pop('mapping') for x in atoms]
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
            g.add_atom(atom, mapping[n])

        if not is_query and not is_cgr:
            for n, m in enumerate(parsed_mapping):
                g.atom(mapping[n])._parsed_mapping = m

        for n, m, b in prepared_bonds:
            n_map, m_map = mapping[n], mapping[m]
            g.add_bond(n_map, m_map, b)
        return g

    __bondlabels = {'0': None, '1': 1, '2': 2, '3': 3, '4': 4, '5': 5, '9': 5, 'n': None, 's': 5}


class CGRwrite:
    def __init__(self, xyz=False):
        self.__xyz = xyz

    def _convert_structure(self, g):
        if isinstance(g, CGRContainer):
            format_atom = self.__format_dyn_atom
            format_bond = self.__format_dyn_bond
            stereo_map = None
        elif isinstance(g, MoleculeContainer):
            format_atom = self.__format_atom
            format_bond = self.__format_bond
            stereo_map = None
        else:
            raise TypeError('queries is read-only')

        cgr_data, atoms, bonds = [], [], []
        renum = {}
        for n, (i, atom) in enumerate(g._node.items(), start=1):
            renum[i] = n
            cgr, symbol, isotope, charge, multiplicity = format_atom(atom, n)
            cgr_data.extend(cgr)
            atoms.append({'mapping': i, 'symbol': symbol, 'isotope': isotope, 'charge': charge,
                          'multiplicity': multiplicity, 'x': atom.x, 'y': atom.y, 'z': atom.z, 'id': n})

        if stereo_map:  # wedge stereo
            for n, m, bond in g.bonds():
                stereo = stereo_map.get((n, m))
                if stereo:
                    n_map, m_map = renum[n], renum[m]
                else:
                    stereo = stereo_map.get((m, n))
                    if stereo:
                        n_map, m_map = renum[m], renum[n]
                    else:
                        n_map, m_map = renum[n], renum[m]

                cgr, order = format_bond(bond, n_map, m_map)
                cgr_data.extend(cgr)
                bonds.append((n_map, m_map, order, self._stereo_map(stereo)))
        else:
            for n, m, bond in g.bonds():
                n_map, m_map = renum[n], renum[m]
                cgr, order = format_bond(bond, n_map, m_map)
                cgr_data.extend(cgr)
                bonds.append((n_map, m_map, order, self._stereo_map(None)))

        return atoms, bonds, cgr_data

    def __format_atom(self, atom, n):
        isotope = self._isotope_map(atom.isotope if atom.isotope != atom.common_isotope else None, n)
        return (), atom.element, isotope, self._charge_map(atom.charge), self._multiplicity_map(atom.multiplicity, n)

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
        return cgr, atom.element, isotope, self._charge_map(atom.charge), self._multiplicity_map(atom.multiplicity, n)

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
    def _multiplicity_map(x, n):
        pass

    @staticmethod
    @abstractmethod
    def _bond_map(x):
        pass
