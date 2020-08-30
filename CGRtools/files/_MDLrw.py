# -*- coding: utf-8 -*-
#
#  Copyright 2017-2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from base64 import urlsafe_b64encode
from collections import defaultdict
from csv import reader
from io import StringIO, TextIOWrapper
from itertools import chain, islice
from os.path import abspath, join
from pathlib import Path
from pickle import dump, load, UnpicklingError
from sys import platform
from tempfile import gettempdir
from ._CGRrw import CGRRead, common_isotopes, parse_error
from ..containers import MoleculeContainer, CGRContainer, QueryContainer, QueryCGRContainer
from ..exceptions import EmptyMolecule, NotChiral, IsChiral, ValenceError


query_keys = {'atomhyb': 'hybridization', 'hybridization': 'hybridization', 'hyb': 'hybridization',
              'atomneighbors': 'neighbors', 'neighbors': 'neighbors'}


class MOLRead:
    def __init__(self, line, log_buffer=None):
        atom_count = int(line[0:3])
        if not atom_count:
            raise EmptyMolecule

        self.__bonds_count = int(line[3:6])
        self.__atoms_count = atom_count
        self.__cgr = {}
        self.__query = []
        self.__atoms = []
        self.__bonds = []
        self.__stereo = []
        if log_buffer is None:
            log_buffer = []
        self.__log_buffer = log_buffer

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
            self.__bonds.append((a1, a2, int(line[6:9])))
        elif line.startswith('M  END'):
            cgr = []
            for a in self.__atoms:
                if a['is_radical']:
                    a['is_radical'] = True
            for x in self.__cgr.values():
                try:
                    atoms = x['atoms']
                    _type = x['type']
                    value = x['value']
                    if _type in query_keys:
                        if len(atoms) != 1 or atoms[0] == -1 or not value:
                            raise ValueError(f'CGR Query spec invalid {x}')
                        self.__query.append((atoms[0], query_keys[_type], [int(x) for x in value.split(',')]))
                    elif _type == 'dynatom':
                        if len(atoms) != 1 or atoms[0] == -1 or not value:
                            raise ValueError(f'CGR spec invalid {x}')
                        if value == 'r1':  # only dynatom = r1 acceptable. this means change of radical state
                            cgr.append((atoms[0], 'radical', None))
                        elif value[0] == 'c':
                            cgr.append((atoms[0], 'charge', int(value[1:])))
                        else:
                            raise ValueError('unknown dynatom')
                    elif _type == 'dynbond':
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


class EMOLRead:
    def __init__(self, log_buffer=None):
        self.__atoms = []
        self.__bonds = []
        self.__atom_map = {}
        self.__stereo = []
        if log_buffer is None:
            log_buffer = []
        self.__log_buffer = log_buffer

    def getvalue(self):
        if self.__in_mol or self.__in_mol is None:
            raise ValueError('molecule not complete')
        return {'atoms': self.__atoms, 'bonds': self.__bonds, 'stereo': self.__stereo}

    def __call__(self, line, lineu=None):
        if lineu is None:
            lineu = line.upper()
        if self.__in_mol:
            if lineu.startswith('M  V30 END CTAB'):
                self.__in_mol = False
                return True
            elif self.__atoms_count:
                if lineu.startswith('M  V30 END'):
                    x = lineu[11:].strip()
                    cp = self.__parser
                    self.__parser = None

                    if x == 'ATOM':
                        if cp == self.__atom_parser and len(self.__atoms) == self.__atoms_count:
                            return
                    elif x == 'BOND':
                        if cp == self.__bond_parser and len(self.__bonds) == self.__bonds_count:
                            return
                    else:
                        return
                    raise ValueError('invalid number of %s records or invalid CTAB' % x)

                elif self.__parser:
                    collected = self.__record_collector(line)
                    if collected:
                        self.__parser(collected)

                elif lineu.startswith('M  V30 BEGIN ATOM'):
                    self.__parser = self.__atom_parser
                elif lineu.startswith('M  V30 BEGIN BOND'):
                    self.__parser = self.__bond_parser
                elif lineu.startswith(('M  V30 BEGIN SGROUP', 'M  V30 BEGIN COLLECTION')):
                    self.__parser = self.__ignored_block_parser
                else:
                    raise ValueError('invalid CTAB')

            else:  # M  V30 COUNTS line expected
                a, b, *_ = line[13:].split()
                atom_count = int(a)
                if not atom_count:
                    raise EmptyMolecule
                self.__bonds_count = int(b)
                self.__atoms_count = atom_count

        elif self.__in_mol is not None:
            raise SyntaxError('invalid usage')
        elif not lineu.startswith('M  V30 BEGIN CTAB'):
            raise ValueError('invalid CTAB')
        else:
            self.__in_mol = True

    def __record_collector(self, line):
        if not line.endswith('-\n'):
            line = line[7:]
            if self.__record:
                line = self.__record + line
                self.__record = None

            return next(reader([line], delimiter=' ', quotechar='"', skipinitialspace=True))

        line = line[7:-2]
        if not self.__record:
            self.__record = line
        else:
            self.__record += line

    def __bond_parser(self, line):
        _, t, a1, a2, *kvs = line
        try:
            self.__bonds.append((self.__atom_map[a1], self.__atom_map[a2], int(t)))
        except KeyError:
            raise ValueError('invalid atoms numbers')
        for kv in kvs:
            k, v = kv.split('=')
            if k in ('CFG', 'cfg', 'Cfg'):
                if v == '1':
                    self.__stereo.append((self.__atom_map[a1], self.__atom_map[a2], 1))
                elif v == '3':
                    self.__stereo.append((self.__atom_map[a1], self.__atom_map[a2], -1))
                else:
                    self.__log_buffer.append('invalid or unsupported stereo')
                break

    def __atom_parser(self, line):
        n, a, x, y, z, m, *kvs = line
        i = None
        c = 0
        r = False
        for kv in kvs:
            k, v = kv.split('=', 1)
            k = k.upper()
            if k == 'CHG':
                c = int(v)
            elif k == 'MASS':
                i = int(v)
            elif k == 'RAD':
                r = True

        self.__atom_map[n] = n = len(self.__atoms)
        if a.startswith(('[', 'NOT', 'not', 'Not')):
            raise ValueError('list of atoms not supported')
        elif a == 'D':
            if i:
                raise ValueError('isotope on deuterium atom')
            a = 'H'
            i = 2

        self.__atoms.append({'element': a, 'isotope': i, 'charge': c, 'is_radical': r,
                             'x': float(x), 'y': float(y), 'z': float(z), 'mapping': int(m)})

    def __ignored_block_parser(self, line):
        return

    __record = __atoms_count = __in_mol = __parser = None


class RXNRead:
    def __init__(self, line, ignore=False, log_buffer=None):
        self.__reactants_count = int(line[:3])
        self.__products_count = int(line[3:6]) + self.__reactants_count
        self.__reagents_count = int(line[6:].rstrip() or 0) + self.__products_count
        self.__molecules = []
        self.__ignore = ignore
        if log_buffer is None:
            log_buffer = []
        self.__log_buffer = log_buffer

    def __call__(self, line):
        if self.__parser:
            if self.__parser(line):
                self.__im = 4
                mol = self.__parser.getvalue()
                if self.__title:
                    mol['title'] = self.__title
                self.__molecules.append(mol)
                self.__parser = None
                if len(self.__molecules) == self.__reagents_count:
                    self.__rend = True
                    return True
        elif self.__empty_skip:
            if not line.startswith('$MOL'):
                return
            self.__empty_skip = False
            self.__im = 3
        elif self.__rend:
            raise ValueError('parser closed')
        elif self.__im == 4:
            if not line.startswith('$MOL'):
                raise ValueError('invalid RXN')
            self.__im = 3
        elif self.__im:
            if self.__im == 3:
                self.__title = line.strip()
            self.__im -= 1
        else:
            try:
                self.__parser = MOLRead(line, self.__log_buffer)
            except EmptyMolecule:
                if not self.__ignore:
                    raise
                self.__log_buffer.append('empty molecule ignored')
                if len(self.__molecules) < self.__reactants_count:
                    self.__reactants_count -= 1
                    self.__products_count -= 1
                    self.__reagents_count -= 1
                elif len(self.__molecules) < self.__products_count:
                    self.__products_count -= 1
                    self.__reagents_count -= 1
                elif len(self.__molecules) < self.__reagents_count:
                    self.__reagents_count -= 1
                    if len(self.__molecules) == self.__reagents_count:  # empty molecule is last in list
                        self.__rend = True
                        return True
                self.__empty_skip = True

    def getvalue(self):
        if self.__rend:
            return {'reactants': self.__molecules[:self.__reactants_count],
                    'products': self.__molecules[self.__reactants_count:self.__products_count],
                    'reagents': self.__molecules[self.__products_count:self.__reagents_count]}
        raise ValueError('reaction not complete')

    __parser = None
    __empty_skip = __rend = False
    __im = 4


class ERXNRead:
    def __init__(self, line, ignore=False, log_buffer=None):
        tmp = line[13:].split()
        self.__reactants_count = int(tmp[0])
        self.__products_count = int(tmp[1])
        self.__reagents_count = int(tmp[2]) if len(tmp) == 3 else 0

        self.__reactants = []
        self.__products = []
        self.__reagents = []
        self.__ignore = ignore
        if log_buffer is None:
            log_buffer = []
        self.__log_buffer = log_buffer

    def __call__(self, line):
        lineu = line.upper()
        if self.__empty_skip:
            if not lineu.startswith('M  V30 END CTAB'):
                return
            self.__empty_skip = False
        elif self.__in_mol:
            try:
                x = self.__parser(line, lineu)
            except EmptyMolecule:
                if not self.__ignore:
                    raise
                self.__empty_skip = True
                self.__in_mol -= 1
                if self.__in_mol:
                    self.__parser = EMOLRead(self.__log_buffer)
                self.__log_buffer.append('empty molecule ignored')
            else:
                if x:
                    x = self.__parser.getvalue()
                    self.__in_mol -= 1
                    if self.__in_mol:
                        self.__parser = EMOLRead(self.__log_buffer)
                    if self.__parser_group == 'REACTANT':
                        self.__reactants.append(x)
                    elif self.__parser_group == 'PRODUCT':
                        self.__products.append(x)
                    elif self.__parser_group == 'AGENT':
                        self.__reagents.append(x)
        elif self.__rend:
            raise ValueError('parser closed')
        elif lineu.startswith('M  V30 END'):
            if self.__parser_group != lineu[11:].strip():
                raise ValueError('invalid CTAB')
        elif lineu.startswith('M  V30 BEGIN'):
            x = lineu[13:].strip()
            if x == 'REACTANT':
                self.__in_mol = self.__reactants_count
            elif x == 'PRODUCT':
                self.__in_mol = self.__products_count
            elif x == 'AGENT':
                self.__in_mol = self.__reagents_count
            else:
                raise ValueError('invalid RXN CTAB')
            self.__parser_group = x
            if self.__in_mol:
                self.__parser = EMOLRead(self.__log_buffer)
        elif lineu.startswith('M  END'):
            self.__rend = True
            return True
        else:
            raise ValueError('invalid CTAB')

    def getvalue(self):
        if self.__rend:
            return {'reactants': self.__reactants, 'products': self.__products, 'reagents': self.__reagents}
        raise ValueError('reaction not complete')

    __parser_group = __parser = None
    __rend = __empty_skip = False
    __in_mol = 0


class MDLStereo(CGRRead):
    def __init__(self, calc_cis_trans=False, ignore_stereo=False, **kwargs):
        super().__init__(**kwargs)
        self.__calc_cis_trans = calc_cis_trans
        self.__ignore_stereo = ignore_stereo

    def _convert_molecule(self, molecule, mapping):
        mol = super()._convert_molecule(molecule, mapping)
        if self.__ignore_stereo:
            return mol

        if self.__calc_cis_trans:
            mol.calculate_cis_trans_from_2d()

        stereo = [(mapping[n], mapping[m], s) for n, m, s in molecule['stereo']]
        while stereo:
            fail_stereo = []
            old_stereo = len(stereo)
            for n, m, s in stereo:
                try:
                    mol.add_wedge(n, m, s, clean_cache=False)
                except NotChiral:
                    fail_stereo.append((n, m, s))
                except IsChiral:
                    pass
                except ValenceError:
                    self._info('structure has errors, stereo data skipped')
                    mol.flush_cache()
                    break
            else:
                stereo = fail_stereo
                if len(stereo) == old_stereo:
                    break
                del mol.__dict__['_MoleculeStereo__chiral_centers']
                if self.__calc_cis_trans:
                    mol.calculate_cis_trans_from_2d(clean_cache=False)
                continue
            break
        return mol


class MDLReadMeta(type):
    def __call__(cls, *args, **kwargs):
        if kwargs.get('indexable'):
            cls = type(cls.__name__, (cls,), {'__len__': lambda x: len(x._shifts) - 1, '__module__': cls.__module__})
        obj = object.__new__(cls)
        obj.__init__(*args, **kwargs)
        return obj


class MDLRead(MDLStereo, metaclass=MDLReadMeta):
    def __init__(self, file, **kwargs):
        if isinstance(file, str):
            self._file = open(file)
            self._is_buffer = False
        elif isinstance(file, Path):
            self._file = file.open()
            self._is_buffer = False
        elif isinstance(file, (TextIOWrapper, StringIO)):
            self._file = file
            self._is_buffer = True
        else:
            raise TypeError('invalid file. TextIOWrapper, StringIO subclasses possible')
        super().__init__(**kwargs)

    def close(self, force=False):
        """
        Close opened file

        :param force: force closing of externally opened file or buffer
        """
        if not self._is_buffer or force:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    def _load_cache(self):
        """
        Load existing cache or create new. Working only for UNIX-like systems and local files (not buffers).
        """
        if platform == 'win32' or self._is_buffer:
            return
        try:
            with open(self.__cache_path, 'rb') as f:
                self._shifts = load(f)
        except FileNotFoundError:  # cache not found
            self.reset_index()
        except IsADirectoryError as e:
            raise IsADirectoryError(f'Please delete {self.__cache_path} directory') from e
        except (UnpicklingError, EOFError) as e:  # invalid file. ask user to check it.
            raise UnpicklingError(f'Invalid cache file {self.__cache_path}. Please delete it') from e

    def reset_index(self):
        """
        Create (rewrite) indexation table. Implemented only for object that
        is a real file (the path to the file is specified) because the external grep utility is used.
        """
        if platform != 'win32' and not self._is_buffer:
            self._shifts = self._get_shifts(self._file.name)
            with open(self.__cache_path, 'wb') as f:
                dump(self._shifts, f)
        else:
            raise self._implement_error

    @property
    def __cache_path(self):
        return abspath(join(gettempdir(), 'cgrtools_' + urlsafe_b64encode(abspath(self._file.name).encode()).decode()))

    def read(self):
        """
        Parse whole file

        :return: list of parsed molecules
        """
        return list(iter(self))

    def __iter__(self):
        return (x for x in self._data if not isinstance(x, parse_error))

    def __next__(self):
        return next(iter(self))

    def __getitem__(self, item):
        """
        Getting the item by index from the original file,
        For slices records with errors skipped.
        For indexed access records with errors returned as error container.
        :return: [Molecule, Reaction]Container or list of [Molecule, Reaction]Containers
        """
        if self._shifts:
            _len = len(self._shifts) - 1
            if isinstance(item, int):
                if item >= _len or item < -_len:
                    raise IndexError('List index out of range')
                if item < 0:
                    item += _len
                self.seek(item)
                return next(self._data)
            elif isinstance(item, slice):
                start, stop, step = item.indices(_len)
                if start == stop:
                    return []
                if step == 1:
                    self.seek(start)
                    records = [x for x in islice(self._data, stop - start) if not isinstance(x, parse_error)]
                else:
                    records = []
                    for index in range(start, stop, step):
                        self.seek(index)
                        record = next(self._data)
                        if not isinstance(record, parse_error):
                            records.append(record)
                return records
            else:
                raise TypeError('Indices must be integers or slices')
        raise self._implement_error

    def _prepare_meta(self, meta):
        new_meta = {}
        for k, v in meta.items():
            if v:
                new_meta[k] = '\n'.join(v)
            else:
                self._info(f'invalid metadata entry: {k}: {v}')
        return new_meta

    _shifts = None
    _implement_error = NotImplementedError('Indexable supported in unix-like o.s. and for files stored on disk')


class MDLWrite:
    def __init__(self, file, *, append: bool = False, write3d: int = 0):
        """
        :param write3d: write for Molecules 3D coordinates instead 2D if exists.
            if 0 - 2D only, 1 - first 3D, 2 - all 3D in sequence.
        """
        if isinstance(file, str):
            self._file = open(file, 'a' if append else 'w')
            self._is_buffer = False
        elif isinstance(file, Path):
            self._file = file.open('a' if append else 'w')
            self._is_buffer = False
        elif isinstance(file, (TextIOWrapper, StringIO)):
            self._file = file
            self._is_buffer = True
        else:
            raise TypeError('invalid file. '
                            'TextIOWrapper, StringIO, BytesIO, BufferedReader and BufferedIOBase subclasses possible')

        if not isinstance(write3d, int):
            raise TypeError('int expected')
        elif write3d not in (0, 1, 2):
            raise ValueError('only 0, 1 and 2 expected')
        self.__write = True
        self.__write3d = write3d

    def close(self, force=False):
        """
        close opened file

        :param force: force closing of externally opened file or buffer
        """
        if self.__write:
            self.write = self.__write_closed
            self.__write = False

        if not self._is_buffer or force:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    @staticmethod
    def __write_closed(_):
        raise ValueError('I/O operation on closed writer')

    def _convert_structure(self, g):
        if isinstance(g, MoleculeContainer):
            bonds = self.__convert_molecule(g)
        elif isinstance(g, CGRContainer):
            bonds = self.__convert_cgr(g)
        elif isinstance(g, QueryContainer):
            bonds = self.__convert_query(g)
        elif isinstance(g, QueryCGRContainer):
            raise TypeError('CGR queries not supported')
        else:
            raise TypeError('Graph expected')

        head = f'{g.name}\n\n\n{g.atoms_count:3d}{g.bonds_count:3d}  0  0  0  0            999 V2000\n'

        gc = g._charges
        gr = g._radicals
        props = []
        for n, (m, a) in enumerate(g._atoms.items(), start=1):
            if a.isotope:
                props.append(f'M  ISO  1 {n:3d} {a.isotope:3d}\n')
            if gr[m]:
                props.append(f'M  RAD  1 {n:3d}   2\n')  # invalid for carbenes
            c = gc[m]
            if c in (-4, 4):
                props.append(f'M  CHG  1 {n:3d} {c:3d}\n')

        if self.__write3d and isinstance(g, MoleculeContainer) and g._conformers:
            if self.__write3d == 2:
                out = [self.__merge(head, self.__convert_atoms3d(g, xyz), bonds, props) for xyz in g._conformers]
            else:
                out = self.__merge(head, self.__convert_atoms3d(g, g._conformers[0]), bonds, props)
        else:
            out = self.__merge(head, self.__convert_atoms2d(g), bonds, props)
        return out

    @staticmethod
    def __merge(head, atoms, bonds, props):
        out = [head]
        out.extend(atoms)
        out.extend(bonds)
        out.extend(props)
        out.append('M  END\n')
        return ''.join(out)

    @classmethod
    def __convert_atoms2d(cls, g):
        gc = g._charges
        gp = g._plane

        out = []
        for n, (m, a) in enumerate(g._atoms.items(), start=1):
            x, y = gp[m]
            c = gc[m]
            if c in (-4, 4):
                out.append(f'{x:10.4f}{y:10.4f}    0.0000 {a.atomic_symbol:3s} 0  0  0  0  0  0  0  0  0{m:3d}  0  0\n')
            else:
                out.append(f'{x:10.4f}{y:10.4f}    0.0000 {a.atomic_symbol:3s} 0{cls.__charge_map[c]}  0  0  0  0'
                           f'  0  0  0{m:3d}  0  0\n')
        return out

    @classmethod
    def __convert_atoms3d(cls, g, xyz):
        gc = g._charges

        out = []
        for n, (m, a) in enumerate(g._atoms.items(), start=1):
            x, y, z = xyz[m]
            c = gc[m]
            if c in (-4, 4):
                out.append(f'{x:10.4f}{y:10.4f}{z:10.4f} {a.atomic_symbol:3s} 0  0  0  0  0  0  0  0  0{m:3d}  0  0\n')
            else:
                out.append(f'{x:10.4f}{y:10.4f}{z:10.4f} {a.atomic_symbol:3s} 0{cls.__charge_map[c]}  0  0  0  0'
                           f'  0  0  0{m:3d}  0  0\n')
        return out

    @classmethod
    def __convert_molecule(cls, g):
        bonds = g._bonds
        atoms = {m: n for n, m in enumerate(g._atoms, start=1)}
        wedge = defaultdict(set)
        out = []
        for n, m, s in g._wedge_map:
            out.append(f'{atoms[n]:3d}{atoms[m]:3d}  {bonds[n][m].order}  {s == 1 and "1" or "6"}  0  0  0\n')
            wedge[n].add(m)
            wedge[m].add(n)
        for n, m, b in g.bonds():
            if m not in wedge[n]:
                out.append(f'{atoms[n]:3d}{atoms[m]:3d}  {b.order}  0  0  0  0\n')
        return out

    @staticmethod
    def __convert_cgr(g):
        gpc = g._p_charges
        gpr = g._p_radicals

        atoms = {m: n for n, m in enumerate(g._atoms, start=1)}
        bonds = []
        props = []
        for n, c in g._charges.items():
            pc = gpc[n]
            if c != pc:
                i = len(props) + 1
                props.append(f'M  SAL {i:3d}  1 {atoms[n]:3d}\n'
                             f'M  SDT {i:3d} dynatom\n'
                             f'M  SDD {i:3d}     0.0000{i / 3:10.4f}    DAU   ALL  0       0\n'
                             f'M  SED {i:3d} c{pc - c:+d}\n')
        for n, r in g._radicals.items():
            if r != gpr[n]:
                i = len(props) + 1
                props.append(f'M  SAL {i:3d}  1 {atoms[n]:3d}\n'
                             f'M  SDT {i:3d} dynatom\n'
                             f'M  SDD {i:3d}     0.0000{i / 3:10.4f}    DAU   ALL  0       0\n'
                             f'M  SED {i:3d} r1\n')

        for n, m, b in g.bonds():
            if b.order != b.p_order:
                i = len(props) + 1
                props.append(f'M  SAL {i:3d}  2 {atoms[n]:3d} {atoms[m]:3d}\n'
                             f'M  SDT {i:3d} dynbond\n'
                             f'M  SDD {i:3d}     0.0000{i / 3:10.4f}    DAU   ALL  0       0\n'
                             f'M  SED {i:3d} {b.order or 0}>{b.p_order or 0}\n')
                bonds.append(f'{atoms[n]:3d}{atoms[m]:3d}  8  0  0  0  0\n')
            else:
                bonds.append(f'{atoms[n]:3d}{atoms[m]:3d}  {b.order}  0  0  0  0\n')

        iterator = iter(range(1, len(props) + 1))
        for first in iterator:
            dat = list(chain((first,), islice(iterator, 7)))
            bonds.append(f'M  STY  {len(dat)}')
            bonds.extend(f'{x:4d} DAT' for x in dat)
            bonds.append('\n')

        bonds.extend(props)
        return bonds

    @classmethod
    def __convert_query(cls, g):
        bonds = g._bonds
        atoms = {m: n for n, m in enumerate(g._atoms, start=1)}
        wedge = defaultdict(set)
        out = []
        for n, m, s in g._wedge_map:
            out.append(f'{atoms[n]:3d}{atoms[m]:3d}  {bonds[n][m].order}  {s == 1 and "1" or "6"}  0  0  0\n')
            wedge[n].add(m)
            wedge[m].add(n)
        for n, m, b in g.bonds():
            if m not in wedge[n]:
                out.append(f'{atoms[n]:3d}{atoms[m]:3d}  {b.order}  0  0  0  0\n')
        props = []

        for n, m in g._neighbors.items():
            if m:
                i = len(props) + 1
                props.append(f'M  SAL {i:3d}  1 {atoms[n]:3d}\n'
                             f'M  SDT {i:3d} neighbors\n'
                             f'M  SDD {i:3d}     0.0000{i / 3:10.4f}    DAU   ALL  0       0\n'
                             f'M  SED {i:3d} {",".join(str(x) for x in m)}\n')
        for n, h in g._hybridizations.items():
            if h:
                i = len(props) + 1
                props.append(f'M  SAL {i:3d}  1 {atoms[n]:3d}\n'
                             f'M  SDT {i:3d} hybridization\n'
                             f'M  SDD {i:3d}     0.0000{i / 3:10.4f}    DAU   ALL  0       0\n'
                             f'M  SED {i:3d} {",".join(str(x) for x in h)}\n')

        iterator = iter(range(1, len(props) + 1))
        for first in iterator:
            dat = list(chain((first,), islice(iterator, 7)))
            out.append(f'M  STY  {len(dat)}')
            out.extend(f'{x:4d} DAT' for x in dat)
            out.append('\n')

        out.extend(props)
        return out

    __charge_map = {-3: '  7', -2: '  6', -1: '  5', 0: '  0', 1: '  3', 2: '  2', 3: '  1'}
