# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
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
from ctypes import c_char, c_double, c_short, c_long, create_string_buffer, POINTER, Structure, cdll, byref
from logging import warning
from re import split
from sys import platform
from traceback import format_exc
from warnings import warn
from ._CGRrw import CGRread, WithMixin
from . import __path__ as files_path
from ..periodictable import common_isotopes


class INCHIread(CGRread, WithMixin):
    """
    INCHI separated per lines files reader. works similar to opened file object. support `with` context manager.
    on initialization accept opened in text mode file, string path to file,
    pathlib.Path object or another buffered reader object.
    line should be start with INCHI string and
    optionally continues with space/tab separated list of key:value [or key=value] data.
    AuxInfo also can be stored in data.

    example:
    InChI=1S/C2H5/c1-2/h1H2,2H3/q+1 AuxInfo=1/0/N:1,2/CRV:1+1/rA:2C+C/rB:s1;/rC:-8,4583,1,1250,0;-7,1247,1,8950,0;
    id:123 key=value\n
    """
    def __init__(self, file, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(CGRread, self).__init__(file)
        self.__data = self.__reader()

    def read(self):
        """
        parse whole file

        :return: list of parsed molecules
        """
        return list(self.__data)

    def __iter__(self):
        return self.__data

    def __next__(self):
        return next(self.__data)

    def __reader(self):
        for line in self._file:
            inchi, *data = line.split()
            meta = {}
            aux = None
            for x in data:
                if x.startswith('AuxInfo='):
                    aux = x
                else:
                    try:
                        k, v = split('[=:]', x, 1)
                        meta[k.strip()] = v.strip()
                    except ValueError:
                        warning(f'invalid metadata entry: {x}')
            try:
                record = self.__parse_aux(aux) if aux else self.__parse_inchi(inchi)
            except ValueError:
                warning(f'line: {inchi} {aux}\nconsist errors:\n{format_exc()}')
                continue

            record['meta'] = meta
            try:
                yield self._convert_structure(record)
            except ValueError:
                warning(f'record consist errors:\n{format_exc()}')

    @staticmethod
    def __parse_inchi(string):
        structure = INCHIStructure()
        lib.GetStructFromINCHI(byref(InputINCHI(string)), byref(structure))

        atoms, bonds = [], []
        seen = set()
        for n in range(structure.num_atoms):
            seen.add(n)
            atom = structure.atom[n]
            element = atom.elname.decode()

            isotope = atom.isotopic_mass
            if isotope in (0, 10000):
                isotope = None
            elif isotope > 10000:
                isotope = isotope - 10000 + common_isotopes[element]

            atoms.append({'element': element, 'charge': int.from_bytes(atom.charge, byteorder='big', signed=True),
                          'mapping': 0, 'x': atom.x, 'y': atom.y, 'z': atom.z, 'isotope': isotope,
                          'multiplicity': int.from_bytes(atom.radical, byteorder='big') or None})

            for k in range(atom.num_bonds):
                m = atom.neighbor[k]
                if m in seen:
                    continue
                order = atom.bond_type[k]
                if order:
                    bonds.append((n, m, order))
        lib.FreeStructFromINCHI(byref(structure))
        return {'atoms': atoms, 'extra': [], 'cgr': [], 'bonds': bonds}

    @staticmethod
    def __parse_aux(string):
        atoms, bonds = [], []
        for n in range(structure.num_atoms):
            atom = structure.atom[n]
            print(atom.x)
            atoms.append({'element': atom.elname, 'charge': atom.charge, 'mapping': 0, 'x': atom.x, 'y': atom.y,
                          'z': atom.z, 'mark': None, 'isotope': atom.isotopic_mass, 'multiplicity': atom.radical})
            # (b['a0'], b['a1'], {'order': self.__bond_map[b['order']]}, None)
        print(atoms)
        return {'atoms': atoms, 'extra': [], 'cgr': [], 'bonds': bonds}


class InputINCHI(Structure):
    def __init__(self, string, options=None):
        if options is None:
            options = create_string_buffer(1)
        else:
            options = create_string_buffer(' '.join(f'{opt_flag}{x}' for x in options).encode())
        super().__init__(create_string_buffer(string.encode()), options)

    _fields_ = [('szInChI', POINTER(c_char)),  # InChI ASCII string to be converted to a strucure
                ('szOptions', POINTER(c_char))  # InChI options: space-delimited; each is preceded
                                                # by '/' or '-' depending on OS and compiler
                ]


class Atom(Structure):
    _fields_ = [('x', c_double), ('y', c_double), ('z', c_double),  # atom coordinates
                ('neighbor', c_short * 20),  # adjacency list: ordering numbers of the adjacent atoms, >= 0
                ('bond_type', c_char * 20),  # inchi_BondType
                ('bond_stereo', c_char * 20),  # inchi_BondStereo2D; negative if the sharp end points to opposite atom
                ('elname', c_char * 6),  # zero-terminated chemical element name: "H", "Si", etc.
                ('num_bonds', c_short),  # number of neighbors, bond types and bond stereo in the adjacency list
                ('num_iso_H', c_char * 4),  # implicit hydrogen atoms
                                            # [0]: number of implicit non-isotopic H (exception: num_iso_H[0]=-1 means
                                            # INCHI adds implicit H automatically),
                                            # [1]: number of implicit isotopic 1H (protium),
                                            # [2]: number of implicit 2H (deuterium),
                                            # [3]: number of implicit 3H (tritium)
                ('isotopic_mass', c_short),  # 0 => non-isotopic; isotopic mass or 10000 + mass - average atomic mass
                ('radical', c_char),  # inchi_Radical,
                ('charge', c_char)  # positive or negative; 0 => no charge
                ]


class Stereo0D(Structure):
    _fields_ = [('neighbor', c_short * 4),  # 4 atoms always
                ('central_atom', c_short),  # central tetrahedral atom or a central atom of allene; otherwise NO_ATOM
                ('type', c_char),  # inchi_StereoType0D
                ('parity', c_char)  # inchi_StereoParity0D: may be a combination of two parities:
                                    # ParityOfConnected | (ParityOfDisconnected << 3), see Note above
                ]


class INCHIStructure(Structure):
    _fields_ = [('atom', POINTER(Atom)),  # array of num_atoms elements
                ('stereo0D', POINTER(Stereo0D)),  # array of num_stereo0D 0D stereo elements or NULL
                ('num_atoms', c_short),  # number of atoms in the structure
                ('num_stereo0D', c_short),  # number of 0D stereo elements
                ('szMessage', POINTER(c_char)),  # Error/warning ASCII message
                ('szLog', POINTER(c_char)),  # log-file ASCII string, contains a human-readable list
                                             # of recognized options and possibly an Error/warn message
                ('WarningFlags', (c_long * 2) * 2)
                ]


class AUXStructure(Structure):
    pass


__all__ = ['INCHIread']

if platform == 'linux':
    opt_flag = '-'
    lib = cdll.LoadLibrary(f'{files_path[0]}/dll/libinchi.so')

elif platform == 'win32':
    opt_flag = '/'
    lib = cdll.LoadLibrary(f'{files_path[0]}/dll/libinchi.dll')
else:
    warn('unsupported platform', ImportWarning)
    __all__ = []
    del INCHIread
