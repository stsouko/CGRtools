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
from collections import defaultdict
from importlib.util import find_spec
from io import StringIO, BytesIO, TextIOWrapper, BufferedIOBase, BufferedReader
from itertools import chain, count
from logging import warning
from pathlib import Path
from traceback import format_exc
from warnings import warn
from ._CGRrw import CGRRead, cgr_keys, query_keys, elements_set, elements_list
from ..containers import MoleculeContainer, CGRContainer, QueryContainer, QueryCGRContainer
from ..containers.common import Graph
from ..exceptions import EmptyMolecule


def xml_dict(parent_element, stop_list=None):
    stop_list = set() if stop_list is None else set(stop_list)
    out = {}
    for x, y in parent_element.items():
        y = y.strip()
        if y:
            x = '@%s' % x.strip()
            out[x] = y

    text = []
    if len(parent_element):
        elements_grouped = defaultdict(list)
        for element in parent_element:
            name = QName(element).localname
            if name in stop_list:
                text.append(tostring(element, encoding=str, with_tail=False))
            else:
                elements_grouped[name].append(element)

            if element.tail:
                t = element.tail.strip()
                if t:
                    text.append(t)

        for element_tag, element_group in elements_grouped.items():
            if len(element_group) == 1:
                out[element_tag] = xml_dict(element_group[0], stop_list)
            else:
                out[element_tag] = [xml_dict(x, stop_list) for x in element_group]

    if parent_element.text:
        t = parent_element.text.strip()
        if t:
            text.insert(0, t)
    if text:
        out['$'] = ''.join(text)

    return out


class MRVRead(CGRRead):
    """
    ChemAxon MRV files reader. works similar to opened file object. support `with` context manager.
    on initialization accept opened in binary mode file, string path to file,
    pathlib.Path object or another binary buffered reader object
    """
    def __init__(self, file, *args, **kwargs):
        if isinstance(file, str):
            self.__file = open(file, 'rb')
            self.__is_buffer = False
        elif isinstance(file, Path):
            self.__file = file.open('rb')
            self.__is_buffer = False
        elif isinstance(file, (BytesIO, BufferedReader, BufferedIOBase)):
            self.__file = file
            self.__is_buffer = True
        else:
            raise TypeError('invalid file. BytesIO, BufferedReader and BufferedIOBase subclasses possible')
        super().__init__(*args, **kwargs)
        self._data = self.__reader()

    def close(self, force=False):
        """
        close opened file

        :param force: force closing of externally opened file or buffer
        """
        if not self.__is_buffer or force:
            self.__file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    def read(self):
        """
        parse whole file

        :return: list of parsed molecules or reactions
        """
        return list(iter(self))

    def __iter__(self):
        return (x for x in self._data if x is not None)

    def __next__(self):
        return next(iter(self))

    def __reader(self):
        for _, element in iterparse(self.__file, tag='{*}MChemicalStruct'):
            parsed = xml_dict(element)
            element.clear()
            if 'molecule' in parsed and isinstance(parsed['molecule'], dict):
                parsed = parsed['molecule']
                try:
                    record = self.__parse_molecule(parsed)
                except (KeyError, ValueError):
                    warning(f'record consist errors:\n{format_exc()}')
                    yield None
                else:
                    if 'propertyList' in parsed and 'property' in parsed['propertyList']:
                        record['meta'] = self.__parse_property(parsed['propertyList']['property'])
                    else:
                        record['meta'] = {}
                    try:
                        yield self._convert_structure(record)
                    except ValueError:
                        warning(f'record consist errors:\n{format_exc()}')
                        yield None
            elif 'reaction' in parsed and isinstance(parsed['reaction'], dict):
                parsed = parsed['reaction']
                try:
                    record = self.__parse_reaction(parsed)
                except (KeyError, ValueError):
                    warning(f'record consist errors:\n{format_exc()}')
                    yield None
                else:
                    if 'propertyList' in parsed and 'property' in parsed['propertyList']:
                        record['meta'] = self.__parse_property(parsed['propertyList']['property'])
                    else:
                        record['meta'] = {}
                    try:
                        yield self._convert_reaction(record)
                    except ValueError:
                        warning(f'record consist errors:\n{format_exc()}')
                        yield None
            else:
                warning('invalid MDocument')
                yield None

    def __parse_reaction(self, data):
        reaction = {'reactants': [], 'products': [], 'reagents': []}
        title = data.get('@title')
        if title:
            reaction['title'] = title
        for tag, group in (('reactantList', 'reactants'), ('productList', 'products'), ('agentList', 'reagents')):
            if tag in data and 'molecule' in data[tag]:
                molecule = data[tag]['molecule']
                if isinstance(molecule, dict):
                    molecule = (molecule,)
                for m in molecule:
                    try:
                        reaction[group].append(self.__parse_molecule(m))
                    except EmptyMolecule:
                        if not self._ignore:
                            raise
                        warning('empty molecule ignored')
        return reaction

    @staticmethod
    def __parse_property(data):
        meta = {}
        if isinstance(data, dict):
            key = data['@title']
            val = data['scalar']['$'].strip()
            if key and val:
                meta[key] = val
            else:
                warning(f'invalid metadata entry: {data}')
        else:
            for x in data:
                key = x['@title']
                val = x['scalar']['$'].strip()
                if key and val:
                    meta[key] = val
                else:
                    warning(f'invalid metadata entry: {x}')
        return meta

    def __parse_molecule(self, data):
        atoms, bonds, stereo = [], [], []
        atom_map = {}
        if 'atom' in data['atomArray']:
            da = data['atomArray']['atom']
            if isinstance(da, dict):
                da = (da,)
            for n, atom in enumerate(da):
                atom_map[atom['@id']] = n
                atoms.append({'element': atom['@elementType'],
                              'isotope': int(atom['@isotope']) if '@isotope' in atom else None,
                              'charge': int(atom.get('@formalCharge', 0)),
                              'is_radical': '@radical' in atom,
                              'mapping': int(atom.get('@mrvMap', 0))})
                if '@z3' in atom:
                    atoms[-1].update(x=float(atom['@x3']), y=float(atom['@y3']), z=float(atom['@z3']))
                else:
                    atoms[-1].update(x=float(atom['@x2']) / 2, y=float(atom['@y2']) / 2, z=0.)
                if '@mrvQueryProps' in atom:
                    raise ValueError('queries unsupported')
        else:
            atom = data['atomArray']
            for n, (_id, e) in enumerate(zip(atom['@atomID'].split(), atom['@elementType'].split())):
                atom_map[_id] = n
                atoms.append({'element': e, 'charge': 0, 'mapping': 0, 'isotope': None, 'is_radical': False})
            if '@z3' in atom:
                for a, x, y, z in zip(atoms, atom['@x3'].split(), atom['@y3'].split(), atom['@z3'].split()):
                    a['x'] = float(x)
                    a['y'] = float(y)
                    a['z'] = float(z)
            else:
                for a, x, y in zip(atoms, atom['@x2'].split(), atom['@y2'].split()):
                    a['x'] = float(x) / 2
                    a['y'] = float(y) / 2
                    a['z'] = 0.
            if '@isotope' in atom:
                for a, x in zip(atoms, atom['@isotope'].split()):
                    if x != '0':
                        a['isotope'] = int(x)
            if '@formalCharge' in atom:
                for a, x in zip(atoms, atom['@formalCharge'].split()):
                    if x != '0':
                        a['charge'] = int(x)
            if '@mrvMap' in atom:
                for a, x in zip(atoms, atom['@mrvMap'].split()):
                    if x != '0':
                        a['mapping'] = int(x)
            if '@radical' in atom:
                for a, x in zip(atoms, atom['@radical'].split()):
                    if x != '0':
                        a['is_radical'] = True
            if '@mrvQueryProps' in atom:
                raise ValueError('queries unsupported')
        if not atoms:
            raise EmptyMolecule

        if 'bond' in data['bondArray']:
            db = data['bondArray']['bond']
            if isinstance(db, dict):
                db = (db,)
            for bond in db:
                order = self.__bond_map[bond['@queryType' if '@queryType' in bond else '@order']]
                a1, a2 = bond['@atomRefs2'].split()
                if 'bondStereo' in bond:
                    if '$' in bond['bondStereo']:
                        s = bond['bondStereo']['$']
                        if s == 'H':
                            s = -1
                        elif s == 'W':
                            s = 1
                        else:
                            warning('invalid or unsupported stereo')
                        stereo.append((atom_map[a1], atom_map[a2], s))
                    else:
                        warning('incorrect bondStereo tag')
                bonds.append((atom_map[a1], atom_map[a2], order))

        mol = {'atoms': atoms, 'bonds': bonds, 'stereo': stereo}
        if '@title' in data:
            mol['title'] = data['@title']
        return mol

    __bond_map = {'Any': 8, 'any': 8, 'A': 4, 'a': 4, '1': 1, '2': 2, '3': 3}
    __radical_map = {'monovalent': 2, 'divalent': 1, 'divalent1': 1, 'divalent3': 3}


class MRVWrite:
    """
    ChemAxon MRV files writer. works similar to opened for writing file object. support `with` context manager.
    on initialization accept opened for writing in text mode file, string path to file,
    pathlib.Path object or another buffered writer object
    """
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
        self.__writable = True

    def close(self, force=False):
        """
        write close tag of MRV file and close opened file

        :param force: force closing of externally opened file or buffer
        """
        if not self.__finalized:
            self._file.write('</cml>\n')
            self.__finalized = True
        if self.__writable:
            self.write = self.__write_closed
            self.__writable = False

        if not self._is_buffer or force:
            self._file.close()

    def __enter__(self):
        return self

    def __exit__(self, _type, value, traceback):
        self.close()

    @staticmethod
    def __write_closed(_):
        raise ValueError('I/O operation on closed writer')

    def write(self, data):
        """
        write single molecule or reaction into file
        """
        self._file.write('<cml>\n')
        self.__write(data)
        self.write = self.__write

    def __write(self, data):
        self._file.write('<MDocument><MChemicalStruct>')

        if isinstance(data, Graph):
            if data.name:
                self._file.write(f'<molecule title="{data.name}">')
            else:
                self._file.write('<molecule>')

            if data.meta:
                self._file.write('<propertyList>')
                for k, v in data.meta.items():
                    if isinstance(v, str):
                        v = f'<![CDATA[{v}]]>'
                    self._file.write(f'<property title="{k}"><scalar>{v}</scalar></property>')
                self._file.write('</propertyList>')
            self._file.write(self.__convert_structure(data))
            self._file.write('</molecule>')
        else:
            if not data._arrow:
                data.fix_positions()

            if data.name:
                self._file.write(f'<reaction title="{data.name}">')
            else:
                self._file.write('<reaction>')

            if data.meta:
                self._file.write('<propertyList>')
                for k, v in data.meta.items():
                    if isinstance(v, str):
                        v = f'<![CDATA[{v}]]>'
                    self._file.write(f'<property title="{k}"><scalar>{v}</scalar></property>')
                self._file.write('</propertyList>')
            c = count(1)
            for i, j in ((data.reactants, 'reactantList'), (data.products, 'productList'),
                         (data.reagents, 'agentList')):
                if not i:
                    continue
                self._file.write(f'<{j}>')
                for n, m in zip(c, i):
                    if m.name:
                        self._file.write(f'<molecule title="{m.name}" molID="m{n}">')
                    else:
                        self._file.write(f'<molecule molID="m{n}">')
                    self._file.write(self.__convert_structure(m))
                    self._file.write('</molecule>')
                self._file.write(f'</{j}>')

            self._file.write(f'<arrow type="DEFAULT" x1="{data._arrow[0] * 2:.4f}" y1="0" '
                             f'x2="{data._arrow[1] * 2:.4f}" y2="0"/>')
            self._file.write('</reaction>')
        self._file.write('</MChemicalStruct></MDocument>\n')

    def __convert_structure(self, g):
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

        out = ['<atomArray>']
        for n, atom in g._atoms.items():
            x, y = gp[n]
            out.append(f'<atom id="a{n}" elementType="{atom.atomic_symbol}" '
                       f'x2="{x * 2:.4f}" y2="{y * 2:.4f}" mrvMap="{n}"')
            if gc[n]:
                out.append(f' formalCharge="{gc[n]}"')
            if gr[n]:
                out.append(' radical="monovalent"')
            if atom.isotope:
                out.append(f' isotope="{atom.isotope}"')
            out.append('/>')
        out.append('</atomArray>')
        out.extend(bonds)
        return ''.join(out)

    @classmethod
    def __convert_molecule(cls, g):
        bonds = g._bonds
        bond_map = cls.__bond_map
        wedge = defaultdict(set)
        out = ['<bondArray>']
        for n, (i, j, s) in enumerate(g._wedge_map, start=1):
            out.append(f'<bond id="b{n}" atomRefs2="a{i} a{j}" order="{bond_map[bonds[i][j].order]}">'
                       f'<bondStereo>{s == 1 and "W" or "H"}</bondStereo></bond>')
            wedge[i].add(j)
            wedge[j].add(i)
        for n, (i, j, bond) in enumerate(g.bonds(), start=len(out)):
            if j not in wedge[i]:
                out.append(f'<bond id="b{n}" atomRefs2="a{i} a{j}" order="{bond_map[bond.order]}"/>')
        out.append('</bondArray>')
        return out

    @classmethod
    def __convert_cgr(cls, g):
        bond_map = cls.__bond_map
        gpc = g._p_charges
        gpr = g._p_radicals

        out = ['<bondArray>']
        dyn = []
        for n, (i, j, bond) in enumerate(g.bonds(), start=1):
            if bond.order != bond.p_order:
                dyn.append((i, j, bond))
                out.append(f'<bond id="b{n}" atomRefs2="a{i} a{j}" order="1" queryType="Any"/>')
            else:
                out.append(f'<bond id="b{n}" atomRefs2="a{i} a{j}" order="{bond_map[bond.order]}"/>')
        out.append('</bondArray>')
        for n, (i, j, bond) in enumerate(dyn, start=1):
            out.append(f'<molecule id="sg{n}" role="DataSgroup" fieldName="dynbond" '
                       f'fieldData="{bond.order or 0}&gt;{bond.p_order or 0}" '
                       f'atomRefs="a{i} a{j}" x="0" y="{n / 3:.4f}"/>')
        n = len(dyn)
        for n, (m, c) in enumerate(g._charges.items(), start=n + 1):
            if c != gpc[m]:
                out.append(f'<molecule id="sg{n}" role="DataSgroup" fieldName="dynatom" '
                           f'fieldData="c{gpc[m] - c:+d}" atomRefs="a{m}" x="0" y="{n / 3:.4f}"/>')
        for n, (m, r) in enumerate(g._radicals.items(), start=n + 1):
            if r != gpr[m]:
                out.append(f'<molecule id="sg{n}" role="DataSgroup" fieldName="dynatom" '
                           f'fieldData="r1" atomRefs="a{m}" x="0" y="{n / 3:.4f}"/>')
        return out

    @classmethod
    def __convert_query(cls, g):
        out = cls.__convert_molecule(g)

        n = 0
        for n, (m, an) in enumerate(g._neighbors.items(), start=1):
            if an:
                an = ','.join(str(x) for x in an)
                out.append(f'<molecule id="sg{n}" role="DataSgroup" fieldName="neighbors" '
                           f'fieldData="{an}" atomRefs="a{m}" x="0" y="{n / 3:.4f}"/>')
        for n, (m, h) in enumerate(g._hybridizations.items(), start=n + 1):
            if h:
                h = ','.join(str(x) for x in h)
                out.append(f'<molecule id="sg{n}" role="DataSgroup" fieldName="hybridization" '
                           f'fieldData="{h}" atomRefs="a{m}" x="0" y="{n / 3:.4f}"/>')
        return out

    __bond_map = {8: '1" queryType="Any', 4: 'A', 1: '1', 2: '2', 3: '3', None: '0'}
    __finalized = False


class MRVread:
    def __init__(self, *args, **kwargs):
        warn('MRVread deprecated. Use MRVRead instead', DeprecationWarning)
        warning('MRVread deprecated. Use MRVRead instead')
        self.__obj = MRVRead(*args, **kwargs)

    def __getattr__(self, item):
        return getattr(self.__obj, item)

    def __iter__(self):
        return iter(self.__obj)

    def __next__(self):
        return next(self.__obj)

    def __enter__(self):
        return self.__obj.__enter__()

    def __exit__(self, _type, value, traceback):
        return self.__obj.__exit__(_type, value, traceback)


class MRVwrite:
    def __init__(self, *args, **kwargs):
        warn('MRVwrite deprecated. Use MRVWrite instead', DeprecationWarning)
        warning('MRVwrite deprecated. Use MRVWrite instead')
        self.__obj = MRVWrite(*args, **kwargs)

    def __getattr__(self, item):
        return getattr(self.__obj, item)

    def __enter__(self):
        return self.__obj.__enter__()

    def __exit__(self, _type, value, traceback):
        return self.__obj.__exit__(_type, value, traceback)


__all__ = ['MRVWrite', 'MRVwrite']

if find_spec('lxml'):
    from lxml.etree import iterparse, QName, tostring

    __all__.extend(['MRVRead', 'MRVread'])
else:
    warn('lxml library not installed', ImportWarning)
    del MRVRead, MRVread
