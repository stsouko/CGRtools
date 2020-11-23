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
from ...exceptions import EmptyMolecule


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

query_keys = {'atomhyb': 'hybridization', 'hybridization': 'hybridization', 'hyb': 'hybridization',
              'atomneighbors': 'neighbors', 'neighbors': 'neighbors'}


class MOLRead:
    def __init__(self, line, log_buffer=None):
        self.__atoms_count = int(line[0:3])
        self.__bonds_count = int(line[3:6])
        self.__cgr = {}
        self.__query = []
        self.__atoms = []
        self.__bonds = []
        self.__stereo = []
        if log_buffer is None:
            log_buffer = []
        self.__log_buffer = log_buffer

    def __new__(cls, line, log_buffer=None):
        if line.startswith('  0'):
            raise EmptyMolecule
        return super().__new__(cls)

    def getvalue(self):
        if self.__mend:
            mol = {'atoms': self.__atoms, 'bonds': self.__bonds, 'stereo': self.__stereo}
            if self.__cgr:
                mol['cgr'] = self.__cgr
            if self.__query:
                mol['query'] = self.__query
            return mol
        raise ValueError('molecule not complete')

    def __call__(self, line):
        if self.__mend:
            raise ValueError('parser closed')
        elif len(self.__atoms) < self.__atoms_count:
            try:
                charge = self.__charge_map[line[36:39]]
            except KeyError:
                raise ValueError('invalid charge')
            element = line[31:34].strip()
            isotope = line[34:36]

            if element == 'A':
                self.__query.append((len(self.__atoms), 'element', 'A'))
                if isotope != ' 0':
                    raise ValueError('isotope on query atom')
                isotope = None
            elif element == 'L':
                raise ValueError('list of atoms not supported')
            elif element == 'D':
                element = 'H'
                if isotope != ' 0':
                    raise ValueError('isotope on deuterium atom')
                isotope = 2
            elif isotope != ' 0':
                try:
                    isotope = common_isotopes[element] + int(isotope)
                except KeyError:
                    raise ValueError('invalid element symbol')
            else:
                isotope = None

            mapping = line[60:63]
            self.__atoms.append({'element': element, 'charge': charge, 'isotope': isotope, 'is_radical': False,
                                 'mapping': int(mapping) if mapping else 0,
                                 'x': float(line[0:10]), 'y': float(line[10:20]), 'z': float(line[20:30])})

        elif len(self.__bonds) < self.__bonds_count:
            a1, a2 = int(line[0:3]) - 1, int(line[3:6]) - 1
            s = line[9:12]
            if s == '  1':
                self.__stereo.append((a1, a2, 1))
            elif s == '  6':
                self.__stereo.append((a1, a2, -1))
            elif s != '  0':
                self.__log_buffer.append('unsupported or invalid stereo')
            b = int(line[6:9])  # added ad-hoc for bond type 9
            self.__bonds.append((a1, a2, b if b != 9 else 8))
        elif line.startswith('M  END'):
            cgr = []
            for a in self.__atoms:
                if a['is_radical']:
                    a['is_radical'] = True
            for x in self.__cgr.values():
                try:
                    _type = x['type']
                    if _type in query_keys:
                        atoms = x['atoms']
                        value = x['value']
                        if len(atoms) != 1 or atoms[0] == -1 or not value:
                            raise ValueError(f'CGR Query spec invalid {x}')
                        self.__query.append((atoms[0], query_keys[_type], [int(x) for x in value.split(',')]))
                    elif _type == 'dynatom':
                        atoms = x['atoms']
                        value = x['value']
                        if len(atoms) != 1 or atoms[0] == -1 or not value:
                            raise ValueError(f'CGR spec invalid {x}')
                        if value == 'r1':  # only dynatom = r1 acceptable. this means change of radical state
                            cgr.append((atoms[0], 'radical', None))
                        elif value[0] == 'c':
                            cgr.append((atoms[0], 'charge', int(value[1:])))
                        else:
                            raise ValueError('unknown dynatom')
                    elif _type == 'dynbond':
                        atoms = x['atoms']
                        value = x['value']
                        if len(atoms) != 2 or -1 in atoms or not value:
                            raise ValueError(f'CGR spec invalid {x}')
                        bond, p_bond = x['value'].split('>')
                        cgr.append((atoms, 'bond', (int(bond) or None, int(p_bond) or None)))
                    else:
                        self.__log_buffer.append(f'ignored data: {x}')
                except KeyError:
                    raise ValueError(f'CGR spec invalid {x}')
            self.__cgr = cgr
            self.__mend = True
            return True
        else:
            self.__collect(line)

    def __collect(self, line):
        if line.startswith('M  ALS'):
            raise ValueError('list of atoms not supported')
        elif line.startswith(('M  ISO', 'M  RAD', 'M  CHG')):
            _type = self.__ctf_data[line[3]]
            for i in range(int(line[6:9])):
                i8 = i * 8
                atom = int(line[10 + i8:13 + i8])
                if not atom or atom > len(self.__atoms):
                    raise ValueError('invalid atoms number')
                atom = self.__atoms[atom - 1]
                atom[_type] = int(line[14 + i8:17 + i8])

        elif line.startswith('M  STY'):
            for i in range(int(line[6:9])):
                i8 = i * 8
                if 'DAT' == line[14 + i8:17 + i8]:
                    self.__cgr[int(line[10 + i8:13 + i8])] = {}
        elif line.startswith('M  SAL'):
            i = int(line[7:10])
            if i in self.__cgr:
                self.__cgr[i]['atoms'] = tuple(int(line[14 + 4 * i:17 + 4 * i]) - 1 for i in range(int(line[10:13])))
        elif line.startswith('M  SDT'):
            i = int(line[7:10])
            if i in self.__cgr:
                self.__cgr[i]['type'] = line.split()[-1].lower()
        elif line.startswith('M  SED'):
            i = int(line[7:10])
            if i in self.__cgr:
                self.__cgr[i]['value'] = line[10:].strip().replace('/', '').lower()

    __ctf_data = {'R': 'is_radical', 'C': 'charge', 'I': 'isotope'}
    __charge_map = {'  0': 0, '  1': 3, '  2': 2, '  3': 1, '  4': 0, '  5': -1, '  6': -2, '  7': -3}
    __mend = False


__all__ = ['MOLRead', 'common_isotopes']
