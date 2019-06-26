# -*- coding: utf-8 -*-
#
#  Copyright 2017-2019 Ramil Nugmanov <stsouko@live.ru>
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
from csv import reader
from logging import warning
from io import StringIO, TextIOWrapper
from itertools import chain, islice
from os.path import abspath, join
from pathlib import Path
from pickle import dump, load, UnpicklingError
from tempfile import gettempdir
from ._CGRrw import CGRRead, cgr_keys, query_keys, elements_set, elements_list, common_isotopes
from ..containers import MoleculeContainer, CGRContainer, QueryContainer, QueryCGRContainer
from ..exceptions import EmptyMolecule


class MOLRead:
    def __init__(self, line):
        atom_count = int(line[0:3])
        if not atom_count:
            raise EmptyMolecule

        self.__bonds_count = int(line[3:6])
        self.__atoms_count = atom_count
        self.__atoms_lists = {}
        self.__cgr = {}
        self.__query = []
        self.__atoms = []
        self.__bonds = []
        self.__stereo = []

    def getvalue(self):
        if self.__mend:
            return {'atoms': self.__atoms, 'bonds': self.__bonds, 'atoms_lists': self.__atoms_lists, 'cgr': self.__cgr,
                    'stereo': self.__stereo, 'query': self.__query}
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
                self.__atoms_lists[len(self.__atoms)] = elements_list
                element = 'C'  # replace with carbon
                if isotope != ' 0':
                    raise ValueError('isotope on query atom')
                isotope = None
            elif element == 'L':
                if isotope != ' 0':
                    raise ValueError('isotope on query atom')
                isotope = None
            elif element == 'D':
                element = 'H'
                if isotope != ' 0':
                    raise ValueError('isotope on deuterium atom')
                isotope = 2
            elif element not in elements_set:
                raise ValueError('invalid atom symbol')
            elif isotope != ' 0':
                isotope = common_isotopes[element] + int(isotope)
            else:
                isotope = None

            self.__atoms.append({'element': element, 'charge': charge, 'isotope': isotope, 'is_radical': False,
                                 'mapping': int(line[60:63]),
                                 'x': float(line[0:10]), 'y': float(line[10:20]), 'z': float(line[20:30])})

        elif len(self.__bonds) < self.__bonds_count:
            a1, a2 = int(line[0:3]) - 1, int(line[3:6]) - 1
            try:
                stereo = self.__stereo_map[line[9:12]]
            except KeyError:
                warning('unsupported or invalid stereo')
            else:
                if stereo:
                    self.__stereo.append((a1, a2, stereo))
            self.__bonds.append((a1, a2, int(line[6:9])))
        elif line.startswith('M  END'):
            cgr = []
            if any(a['element'] == 'L' for a in self.__atoms):
                raise ValueError('invalid CTAB')
            for a in self.__atoms:
                if a['is_radical']:
                    a['is_radical'] = True
            for x in self.__cgr.values():
                try:
                    xt = x['type']
                    if xt in query_keys:
                        if len(x['atoms']) != 1 or x['atoms'][0] == -1 or not x['value']:
                            raise ValueError(f'CGR Query spec invalid {x}')
                        self.__query.append((x['atoms'][0], query_keys[xt], x['value']))
                    elif xt in cgr_keys:
                        if cgr_keys[xt] != len(x['atoms']) or -1 in x['atoms'] or not x['value']:
                            raise ValueError(f'CGR spec invalid {x}')
                        cgr.append((x['atoms'], xt, x['value']))
                    else:
                        warning(f'ignored data: {x}')
                except KeyError:
                    raise ValueError(f'CGR spec invalid {x}')
            self.__cgr = cgr
            self.__mend = True
            return True
        else:
            self.__collect(line)

    def __collect(self, line):
        if line.startswith('M  ALS'):
            atom = int(line[7:10]) - 1
            if atom < 0 or atom >= len(self.__atoms) or atom in self.__atoms_lists:
                raise ValueError('invalid atom number')
            value = [line[16 + x * 4: 20 + x * 4].strip() for x in range(int(line[10:13]))]
            if not value or not elements_set.issuperset(value):
                raise ValueError('invalid atoms list')

            if line[14] == 'T':
                value = list(elements_set.difference(value))
            elif line[14] != 'F':
                raise ValueError('invalid atoms list')
            self.__atoms_lists[atom] = value
            self.__atoms[atom]['element'] = value[0]
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
    __stereo_map = {'  0': None, '  1': 1, '  6': -1, '  4': None}
    __mend = False


class EMOLRead:
    def __init__(self):
        self.__atoms_lists = {}
        self.__cgr = []
        self.__query = []
        self.__atoms = []
        self.__bonds = []
        self.__sgroup = 0
        self.__atom_map = {}
        self.__stereo = []

    def getvalue(self):
        if self.__in_mol or self.__in_mol is None:
            raise ValueError('molecule not complete')
        return {'atoms': self.__atoms, 'bonds': self.__bonds, 'atoms_lists': self.__atoms_lists, 'cgr': self.__cgr,
                'stereo': self.__stereo, 'query': self.__query}

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
                    elif x == 'SGROUP':
                        if cp == self.__sgroup_parser and self.__sgroup == self.__sgroup_count:
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
                elif lineu.startswith('M  V30 BEGIN SGROUP'):
                    self.__parser = self.__sgroup_parser
                else:
                    raise ValueError('invalid CTAB')

            else:  # M  V30 COUNTS line expected
                a, b, s, *_ = line[13:].split()
                atom_count = int(a)
                if not atom_count:
                    raise EmptyMolecule
                self.__bonds_count = int(b)
                self.__sgroup_count = int(s)
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
                try:
                    stereo = self.__stereo_map[v]
                except KeyError:
                    warning('invalid or unsupported stereo')
                else:
                    if stereo:
                        self.__stereo.append((self.__atom_map[a1], self.__atom_map[a2], stereo))
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
        if a.startswith('['):
            a = a[1:-1].split(',')
            if not elements_set.issuperset(a):
                raise ValueError('invalid atoms list')
            if i:
                raise ValueError('isotope on query atom')
            self.__atoms_lists[n] = a
            a = a[0]
        elif a.startswith(('NOT', 'not', 'Not')):
            a = a[5:-1].split(',')
            if not elements_set.issuperset(a):
                raise ValueError('invalid atoms list')
            if i:
                raise ValueError('isotope on query atom')
            a = list(elements_set.difference(a))
            self.__atoms_lists[n] = a
            a = a[0]
        elif a == 'D':
            if i:
                raise ValueError('isotope on deuterium atom')
            a = 'H'
            i = 2
        elif a == 'A':
            if i:
                raise ValueError('isotope on query atom')
            self.__atoms_lists[n] = elements_list
            a = 'C'

        self.__atoms.append({'element': a, 'isotope': i, 'charge': c, 'is_radical': r,
                             'x': float(x), 'y': float(y), 'z': float(z), 'mapping': int(m)})

    def __sgroup_parser(self, line):
        self.__sgroup += 1
        if line[1] == 'DAT':
            i = int(line[3][7:])
            try:
                if i == 1:
                    atoms = (self.__atom_map[line[4][:-1]],)
                elif i == 2:
                    atoms = (self.__atom_map[line[4]], self.__atom_map[line[5][:-1]])
                else:
                    return
            except (KeyError, IndexError):
                raise ValueError('invalid atom number')

            _type = value = None
            for kv in line[i + 4:]:
                if '=' not in kv:
                    continue
                k, v = kv.split('=', 1)
                if k == 'FIELDNAME':
                    _type = v.lower()
                elif k == 'FIELDDATA':
                    value = v.lower()
            if _type in cgr_keys:
                if value and cgr_keys[_type] == i:
                    self.__cgr.append((atoms, _type, value))
                else:
                    raise ValueError(f'CGR spec invalid {line}')
            elif _type in query_keys:
                if value and i == 1:
                    self.__query.append((atoms[0], query_keys[_type], value))
                else:
                    raise ValueError(f'CGR spec invalid {line}')

    __record = __atoms_count = __in_mol = __parser = None
    __stereo_map = {'1': 1, '3': -1, '2': None}


class RXNRead:
    def __init__(self, line, ignore=False):
        self.__reactants_count = int(line[:3])
        self.__products_count = int(line[3:6]) + self.__reactants_count
        self.__reagents_count = int(line[6:].rstrip() or 0) + self.__products_count
        self.__molecules = []
        self.__ignore = ignore

    def __call__(self, line):
        if self.__parser:
            if self.__parser(line):
                self.__im = 4
                self.__molecules.append(self.__parser.getvalue())
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
            self.__im -= 1
        else:
            try:
                self.__parser = MOLRead(line)
            except EmptyMolecule:
                if not self.__ignore:
                    raise
                self.__empty_skip = True
                if len(self.__molecules) < self.__reactants_count:
                    self.__reactants_count -= 1
                    self.__products_count -= 1
                    self.__reagents_count -= 1
                elif len(self.__molecules) < self.__products_count:
                    self.__products_count -= 1
                    self.__reagents_count -= 1
                elif len(self.__molecules) < self.__reagents_count:
                    self.__reagents_count -= 1
                warning('empty molecule ignored')

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
    def __init__(self, line, ignore=False):
        tmp = line[13:].split()
        self.__reactants_count = int(tmp[0])
        self.__products_count = int(tmp[1])
        self.__reagents_count = int(tmp[2]) if len(tmp) == 3 else 0

        self.__reactants = []
        self.__products = []
        self.__reagents = []
        self.__ignore = ignore

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
                    self.__parser = EMOLRead()
                warning('empty molecule ignored')
            else:
                if x:
                    x = self.__parser.getvalue()
                    self.__in_mol -= 1
                    if self.__in_mol:
                        self.__parser = EMOLRead()
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
                self.__parser = EMOLRead()
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


class MDLReadMeta(type):
    def __call__(cls, *args, **kwargs):
        if kwargs.get('indexable'):
            cls = type(cls.__name__, (cls,), {'__len__': lambda x: len(x._shifts) - 1, '__module__': cls.__module__})
        obj = object.__new__(cls)
        obj.__init__(*args, **kwargs)
        return obj


class MDLRead(CGRRead, metaclass=MDLReadMeta):
    def __init__(self, file, *args, **kwargs):
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
        super().__init__(*args, **kwargs)

    def close(self, force=False):
        """
        close opened file

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
        the method is implemented for the purpose of optimization, byte positions will not be re-read from a file
        that has already been used, if the content of the file has changed, and the name has been left the same,
        the old version of byte offsets will be loaded
        :return: list of byte offsets from existing file
        """
        try:
            with open(self.__cache_path, 'rb') as f:
                return load(f)
        except FileNotFoundError:
            return
        except IsADirectoryError as e:
            raise IsADirectoryError(f'Please delete {self.__cache_path} directory') from e
        except (UnpicklingError, EOFError) as e:
            raise UnpicklingError(f'Invalid cache file {self.__cache_path}. Please delete it') from e

    @property
    def __cache_path(self):
        return abspath(join(gettempdir(), 'cgrtools_' + urlsafe_b64encode(abspath(self._file.name).encode()).decode()))

    def _dump_cache(self, _shifts):
        """
        _shifts dumps in /tmp directory after reboot it will drop
        """
        with open(self.__cache_path, 'wb') as f:
            dump(_shifts, f)

    def read(self):
        """
        parse whole file

        :return: list of parsed molecules
        """
        return list(iter(self))

    def __iter__(self):
        return (x for x in self._data if x is not None)

    def __next__(self):
        return next(iter(self))

    def __getitem__(self, item):
        """
        getting the item by index from the original file,
        if the required record of the file with an error,
        then only the correct record are returned
        :param item: int or slice
        :return: [Molecule, Reaction]Container or list of [Molecule, Reaction]Containers
        """
        if self._shifts:
            _len = len(self._shifts) - 1
            _current_pos = self.tell()

            if isinstance(item, int):
                if item >= _len or item < -_len:
                    raise IndexError('List index out of range')
                if item < 0:
                    item += _len
                self.seek(item)
                records = next(self._data)
            elif isinstance(item, slice):
                start, stop, step = item.indices(_len)
                if start == stop:
                    return []

                if step == 1:
                    self.seek(start)
                    records = [x for x in islice(self._data, 0, stop - start) if x is not None]
                else:
                    records = []
                    for index in range(start, stop, step):
                        self.seek(index)
                        record = next(self._data)
                        if record:
                            records.append(record)
            else:
                raise TypeError('Indices must be integers or slices')

            self.seek(_current_pos)
            if records is None:
                raise IndexError('Data block with requested index contain errors')
            return records
        raise self._implement_error

    @staticmethod
    def _prepare_meta(meta):
        new_meta = {}
        for k, v in meta.items():
            if v:
                new_meta[k] = '\n'.join(v)
            else:
                warning(f'invalid metadata entry: {k}: {v}')
        return new_meta

    _shifts = None
    _implement_error = NotImplementedError('Indexable supported in unix-like o.s. and for files stored on disk')


class MDLWrite:
    def __init__(self, file):
        if isinstance(file, str):
            self._file = open(file, 'w')
            self._is_buffer = False
        elif isinstance(file, Path):
            self._file = file.open('w')
            self._is_buffer = False
        elif isinstance(file, (TextIOWrapper, StringIO)):
            self._file = file
            self._is_buffer = True
        else:
            raise TypeError('invalid file. '
                            'TextIOWrapper, StringIO, BytesIO, BufferedReader and BufferedIOBase subclasses possible')
        self.__write = True

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

        gp = g._plane
        gc = g._charges
        gr = g._radicals
        props = []
        out = [f'\n\n\n{g.atoms_count:3d}{g.bonds_count:3d}  0  0  0  0            999 V2000\n']
        for n, a in g._atoms.items():
            x, y = gp[n]
            out.append(f'{x:10.4f}{y:10.4f}    0.0000 {a.atomic_symbol:3s} 0{self.__charge_map[gc[n]]}  0  0  0  0  0'
                       f'  0  0{n:3d}  0  0\n')
            if a.isotope:
                props.append(f'M  ISO  1 {n:3d} {a.isotope:3d}\n')
            if gr[n]:
                props.append(f'M  RAD  1 {n:3d}   2\n')  # invalid for carbenes
        out.extend(bonds)
        out.extend(props)
        out.append('M  END\n')
        return ''.join(out)

    @staticmethod
    def __convert_molecule(g):
        atoms = {m: n for n, m in enumerate(g._atoms, start=1)}
        return [f'{atoms[n]:3d}{atoms[m]:3d}  {b.order}  0  0  0  0\n' for n, m, b in g.bonds()]

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

    @staticmethod
    def __convert_query(g):
        atoms = {m: n for n, m in enumerate(g._atoms, start=1)}
        bonds = [f'{atoms[n]:3d}{atoms[m]:3d}  {b.order}  0  0  0  0\n' for n, m, b in g.bonds()]
        props = []

        for n, m in g._neighbors.items():
            if m:
                i = len(props) + 1
                props.append(f'M  SAL {i:3d}  1 {atoms[n]:3d}\n'
                             f'M  SDT {i:3d} neighbors\n'
                             f'M  SDD {i:3d}     0.0000{i / 3:10.4f}    DAU   ALL  0       0\n'
                             f'M  SED {i:3d} {",".join(m)}\n')
        for n, h in g._hybridization.items():
            if h:
                i = len(props) + 1
                props.append(f'M  SAL {i:3d}  1 {atoms[n]:3d}\n'
                             f'M  SDT {i:3d} hybridization\n'
                             f'M  SDD {i:3d}     0.0000{i / 3:10.4f}    DAU   ALL  0       0\n'
                             f'M  SED {i:3d} {",".join(h)}\n')

        iterator = iter(range(1, len(props) + 1))
        for first in iterator:
            dat = list(chain((first,), islice(iterator, 7)))
            bonds.append(f'M  STY  {len(dat)}')
            bonds.extend(f'{x:4d} DAT' for x in dat)
            bonds.append('\n')

        bonds.extend(props)
        return bonds

    __stereo_map = {-1: '6', 1: '1', None: '0'}
    __charge_map = {-3: '  7', -2: '  6', -1: '  5', 0: '  0', 1: '  3', 2: '  2', 3: '  1'}
