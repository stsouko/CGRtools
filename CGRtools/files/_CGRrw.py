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
from collections import defaultdict
from itertools import count
from logging import warning
from ..containers import ReactionContainer, MoleculeContainer, CGRContainer, QueryContainer
from ..containers.cgr import DynamicBond
from ..exceptions import MappingError, NotChiral, IsChiral
from ..periodictable import Element, DynamicElement, QueryElement


cgr_keys = {'dynbond': 2, 'dynatom': 1}
query_keys = {'atomhyb': 'hybridization', 'hybridization': 'hybridization', 'hyb': 'hybridization',
              'atomneighbors': 'neighbors', 'neighbors': 'neighbors'}
common_isotopes = {'H': 1, 'He': 4, 'Li': 7, 'Be': 9, 'B': 11, 'C': 12, 'N': 14, 'O': 16, 'F': 19, 'Ne': 20, 'Na': 23,
                   'Mg': 24, 'Al': 27, 'Si': 28, 'P': 31, 'S': 32, 'Cl': 35, 'Ar': 40, 'K': 39, 'Ca': 40, 'Sc': 45,
                   'Ti': 48, 'V': 51, 'Cr': 52, 'Mn': 55, 'Fe': 56, 'Co': 59, 'Ni': 59, 'Cu': 64, 'Zn': 65, 'Ga': 70,
                   'Ge': 73, 'As': 75, 'Se': 79, 'Br': 80, 'Kr': 84, 'Rb': 85, 'Sr': 88, 'Y': 89, 'Zr': 91, 'Nb': 93,
                   'Mo': 96, 'Tc': 98, 'Ru': 101, 'Rh': 103, 'Pd': 106, 'Ag': 108, 'Cd': 112, 'In': 115, 'Sn': 119,
                   'Sb': 122, 'Te': 128, 'I': 127, 'Xe': 131, 'Cs': 133, 'Ba': 137, 'La': 139, 'Ce': 140, 'Pr': 141,
                   'Nd': 144, 'Pm': 145, 'Sm': 150, 'Eu': 152, 'Gd': 157, 'Tb': 159, 'Dy': 163, 'Ho': 165, 'Er': 167,
                   'Tm': 169, 'Yb': 173, 'Lu': 175, 'Hf': 178, 'Ta': 181, 'W': 184, 'Re': 186, 'Os': 190, 'Ir': 192,
                   'Pt': 195, 'Au': 197, 'Hg': 201, 'Tl': 204, 'Pb': 207, 'Bi': 209, 'Po': 209, 'At': 210, 'Rn': 222,
                   'Fr': 223, 'Ra': 226, 'Ac': 227, 'Th': 232, 'Pa': 231, 'U': 238, 'Np': 237, 'Pu': 244, 'Am': 243,
                   'Cm': 247, 'Bk': 247, 'Cf': 251, 'Es': 252, 'Fm': 257, 'Md': 258, 'No': 259, 'Lr': 260, 'Rf': 261,
                   'Db': 270, 'Sg': 269, 'Bh': 270, 'Hs': 270, 'Mt': 278, 'Ds': 281, 'Rg': 281, 'Cn': 285, 'Nh': 278,
                   'Fl': 289, 'Mc': 289, 'Lv': 293, 'Ts': 297, 'Og': 294}
elements_set = set(common_isotopes)
elements_list = list(common_isotopes)


class CGRRead:
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

        # map unmapped atoms.
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

        # find breaks in map. e.g. 1,2,5,6. 3,4 - skipped
        if self.__remap:
            lose = sorted(set(range(1, next(length))) - set(maps['reactants']) - set(maps['products']) -
                          set(maps['reagents']), reverse=True)
            if lose:
                for i, tmp in maps.items():
                    if not tmp:
                        continue
                    for j in lose:
                        maps[i] = tmp = [x if x < j else x - 1 for x in tmp]

        rc = {'reactants': [], 'products': [], 'reagents': []}
        for i, tmp in maps.items():
            shift = 0
            for j in reaction[i]:
                atom_len = len(j['atoms'])
                remapped = {x: y for x, y in enumerate(tmp[shift: atom_len + shift])}
                shift += atom_len
                g = self.__prepare_structure(j, remapped)
                rc[i].append(g)
        return ReactionContainer(meta=reaction['meta'], **rc)

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

        g = self.__prepare_structure(molecule, remapped)
        g.meta.update(molecule['meta'])
        return g

    @staticmethod
    def _convert_molecule(molecule, mapping):
        g = MoleculeContainer()
        pm = g._parsed_mapping
        for n, atom in enumerate(molecule['atoms']):
            n = g.add_atom(Element.from_symbol(atom['element'])(atom['isotope']), mapping[n],
                           charge=atom['charge'], is_radical=atom['is_radical'], xy=(atom['x'], atom['y']))
            pm[n] = atom['mapping']
        for n, m, b in molecule['bonds']:
            g.add_bond(mapping[n], mapping[m], b)
        if any(a['z'] for a in molecule['atoms']):
            g._conformers.append({mapping[n]: (a['x'], a['y'], a['z']) for n, a in enumerate(molecule['atoms'])})

        stereo = [(mapping[n], mapping[m], s) for n, m, s in molecule['stereo']]
        old_stereo = 0
        while len(stereo) != old_stereo:
            fail_stereo = []
            old_stereo = len(stereo)
            for n, m, s in stereo:
                try:
                    g.add_wedge(n, m, s)
                except NotChiral:
                    fail_stereo.append((n, m, s))
                except IsChiral:
                    warning(f'wedge {{{n}, {m}}} on already chiral atom')
            stereo = fail_stereo

        return g

    @staticmethod
    def _convert_cgr(molecule, mapping):
        atoms = molecule['atoms']
        bonds = defaultdict(dict)

        for nm, _type, value in molecule['cgr']:
            n = nm[0]
            if _type == 'dynatom':
                if value == 'r1':  # only dynatom = r1 acceptable. this means change of radical state
                    atoms[n]['p_is_radical'] = not atoms[n]['is_radical']
                elif value[0] == 'c':
                    atoms[n]['p_charge'] = atoms[n]['charge'] + int(value[1:])
                else:
                    raise ValueError('unknown dynatom')
            else:
                bond, p_bond = value.split('>')
                m = nm[1]
                bonds[n][m] = bonds[m][n] = DynamicBond(int(bond) or None, int(p_bond) or None)

        g = CGRContainer()
        pm = g._parsed_mapping
        for n, atom in enumerate(molecule['atoms']):
            n = g.add_atom(DynamicElement.from_symbol(atom['element'])(atom['isotope']), mapping[n],
                           charge=atom['charge'], is_radical=atom['is_radical'],
                           p_charge=atom.get('p_charge', atom['charge']),
                           p_is_radical=atom.get('p_is_radical', atom['is_radical']),
                           xy=(atom['x'], atom['y']))
            pm[n] = atom['mapping']
        for n, m, b in molecule['bonds']:
            if m in bonds[n] and b != 8:
                raise ValueError('CGR spec invalid')
            g.add_bond(mapping[n], mapping[m], bonds[n].get(m, b))
        return g

    @staticmethod
    def _convert_query(molecule, mapping):
        atoms = molecule['atoms']
        if molecule['atoms_lists']:
            warning('atoms lists not supported yet')

        for n, _type, value in molecule['query']:
            atoms[n][_type] = [int(x) for x in value.split(',')]

        g = QueryContainer()
        pm = g._parsed_mapping
        for n, atom in enumerate(molecule['atoms']):
            n = g.add_atom(QueryElement.from_symbol(atom['element'])(atom['isotope']), mapping[n],
                           charge=atom['charge'], is_radical=atom['is_radical'],
                           neighbors=atom.get('neighbors'), hybridization=atom.get('hybridization'),
                           xy=(atom['x'], atom['y']))
            pm[n] = atom['mapping']
        for n, m, b in molecule['bonds']:
            g.add_bond(mapping[n], mapping[m], b)
        return g

    def __prepare_structure(self, molecule, mapping):
        if molecule['atoms_lists'] or molecule['query']:
            if molecule['cgr']:
                raise ValueError('QueryCGR parsing not supported. try to compose Queries')
            return self._convert_query(molecule, mapping)
        elif molecule['cgr']:
            return self._convert_cgr(molecule, mapping)
        else:
            return self._convert_molecule(molecule, mapping)
