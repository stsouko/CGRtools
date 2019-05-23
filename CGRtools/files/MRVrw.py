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
from itertools import chain, count
from logging import warning
from traceback import format_exc
from warnings import warn
from ._CGRrw import CGRread, CGRwrite, WithMixin, cgr_keys
from ..containers.common import BaseContainer
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


class MRVread(CGRread, WithMixin):
    """
    ChemAxon MRV files reader. works similar to opened file object. support `with` context manager.
    on initialization accept opened in binary mode file, string path to file,
    pathlib.Path object or another binary buffered reader object
    """
    def __init__(self, file, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(CGRread, self).__init__(file, 'rb')
        self.__data = self.__reader()

    def read(self):
        """
        parse whole file

        :return: list of parsed molecules or reactions
        """
        return list(self.__data)

    def __iter__(self):
        return self.__data

    def __next__(self):
        return next(self.__data)

    def __reader(self):
        for _, element in iterparse(self._file, tag='{*}MChemicalStruct'):
            parsed = xml_dict(element)
            element.clear()
            if 'molecule' in parsed and isinstance(parsed['molecule'], dict):
                data = parsed['molecule']
                try:
                    record = self.__parse_molecule(data)
                except (KeyError, ValueError):
                    warning(f'record consist errors:\n{format_exc()}')
                else:
                    if 'propertyList' in data and 'property' in data['propertyList']:
                        record['meta'] = self.__parse_property(data['propertyList']['property'])
                    else:
                        record['meta'] = {}
                    try:
                        yield self._convert_structure(record)
                    except ValueError:
                        warning(f'record consist errors:\n{format_exc()}')
            elif 'reaction' in parsed and isinstance(parsed['reaction'], dict):
                data = parsed['reaction']
                try:
                    record = self.__parse_reaction(data)
                except (KeyError, ValueError):
                    warning(f'record consist errors:\n{format_exc()}')
                else:
                    if 'propertyList' in data and 'property' in data['propertyList']:
                        record['meta'] = self.__parse_property(data['propertyList']['property'])
                    else:
                        record['meta'] = {}
                    try:
                        yield self._convert_reaction(record)
                    except ValueError:
                        warning(f'record consist errors:\n{format_exc()}')
            else:
                warning('invalid MDocument')

    def __parse_reaction(self, data):
        reaction = dict(reactants=[], products=[], reagents=[])
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
        atoms, bonds, extra, cgr = [], [], [], []
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
                              'multiplicity': self.__radical_map[atom['@radical']] if '@radical' in atom else None,
                              'mapping': int(atom.get('@mrvMap', 0)),
                              'x': float(atom['@x3'] if '@x3' in atom else atom['@x2']),
                              'y': float(atom['@y3'] if '@y3' in atom else atom['@y2']),
                              'z': float(atom['@z3'] if '@z3' in atom else 0.)})
                if '@mrvQueryProps' in atom and atom['@mrvQueryProps'][0] == 'L':
                    al = atom['@mrvQueryProps']
                    _type = al[1]
                    al = al[2:-1]
                    if not al:
                        raise ValueError('invalid atomlist')
                    extra.append((n, 'atomlist' if _type == ',' else 'atomnotlist', al.split(_type)))
        else:
            atom = data['atomArray']
            for n, (_id, e, x, y) in enumerate(zip(atom['@atomID'].split(),
                                                   atom['@elementType'].split(),
                                                   (atom['@x3'] if '@x3' in atom else atom['@x2']).split(),
                                                   (atom['@y3'] if '@y3' in atom else atom['@y2']).split())):
                atom_map[_id] = n
                atoms.append({'element': e, 'charge': 0, 'x': float(x), 'y': float(y), 'z': 0., 'mapping': 0,
                              'isotope': None, 'multiplicity': None})
            if '@z3' in atom:
                for a, x in zip(atoms, atom['@z3'].split()):
                    a['z'] = float(x)
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
                        a['multiplicity'] = self.__radical_map[x]
            if '@mrvQueryProps' in atom:
                for n, x in enumerate(atom['@mrvQueryProps'].split()):
                    if x[0] == 'L':
                        _type = x[1]
                        x = x[2:-1]
                        if not x:
                            raise ValueError('invalid atomlist')
                        extra.append((n, 'atomlist' if _type == ',' else 'atomnotlist', x.split(_type)))
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
                        try:
                            stereo = self.__stereo_map[bond['bondStereo']['$']]
                        except KeyError:
                            warning('invalid or unsupported stereo')
                            stereo = None
                    else:
                        warning('incorrect bondStereo tag')
                        stereo = None
                else:
                    stereo = None
                bonds.append((atom_map[a1], atom_map[a2], order))

        if 'molecule' in data:
            dm = data['molecule']
            if isinstance(dm, dict):
                dm = (dm,)
            for cgr_dat in dm:
                if cgr_dat['@role'] == 'DataSgroup':
                    t = cgr_dat['@fieldName']
                    if t not in cgr_keys:
                        continue
                    a = tuple(atom_map[x] for x in cgr_dat['@atomRefs'].split())
                    v = cgr_dat['@fieldData'].replace('/', '').lower()
                    if len(a) != cgr_keys[t] or not v:
                        raise ValueError(f'CGR spec invalid: {cgr_dat}')
                    cgr.append((a, t, v))
        return {'atoms': atoms, 'bonds': bonds, 'extra': extra, 'cgr': cgr}

    __bond_map = {'Any': 8, 'any': 8, 'A': 4, 'a': 4, '1': 1, '2': 2, '3': 3}
    __radical_map = {'monovalent': 2, 'divalent': 1, 'divalent1': 1, 'divalent3': 3}
    __stereo_map = {'H': -1, 'W': 1}


class MRVwrite(CGRwrite, WithMixin):
    """
    ChemAxon MRV files writer. works similar to opened for writing file object. support `with` context manager.
    on initialization accept opened for writing in text mode file, string path to file,
    pathlib.Path object or another buffered writer object
    """
    def __init__(self, file, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(CGRwrite, self).__init__(file, 'w')

    def close(self, *args, **kwargs):
        """
        write close tag of MRV file and close opened file

        :param force: force closing of externally opened file or buffer
        """
        if not self.__finalized:
            self._file.write('</cml>')
            self.__finalized = True
        super().close(*args, **kwargs)

    def write(self, data):
        """
        write single molecule or reaction into file
        """
        self._file.write('<cml>')
        self.__write(data)
        self.write = self.__write

    def __write(self, data):
        self._file.write('<MDocument><MChemicalStruct>')

        if isinstance(data, BaseContainer):
            m = self._convert_structure(data)
            self._file.write('<molecule><propertyList>')
            for k, v in data.meta.items():
                if isinstance(v, str):
                    v = f'<![CDATA[{v}]]>'
                self._file.write(f'<property title="{k}"><scalar>{v}</scalar></property>')
            self._file.write('</propertyList>')
            self._file.write(self.__format_mol(*m))
            self._file.write('</molecule>')
        else:
            if not data._arrow:
                data.fix_positions()

            c = count(1)
            self._file.write('<reaction>')
            for i, j in ((data.reactants, 'reactantList'), (data.products, 'productList'),
                         (data.reagents, 'agentList')):
                self._file.write(f'<{j}>')
                for n, m in zip(c, i):
                    m = self._convert_structure(m)
                    self._file.write(f'<molecule molID="m{n}">')
                    self._file.write(self.__format_mol(*m))
                    self._file.write('</molecule>')
                self._file.write(f'</{j}>')

            self._file.write(f'<arrow type="DEFAULT" x1="{data._arrow[0]:.4f}" y1="1" x2="{data._arrow[1]:.4f}" '
                             f'y2="1"/><propertyList>')
            for k, v in data.meta.items():
                if isinstance(v, str):
                    v = f'<![CDATA[{v}]]>'
                self._file.write(f'<property title="{k}"><scalar>{v}</scalar></property>')
            self._file.write('</propertyList></reaction>')
        self._file.write('</MChemicalStruct></MDocument>')

    @staticmethod
    def __format_mol(atoms, bonds, cgr):
        return ''.join(chain(('<atomArray>',),
                             (f'<atom id="a{atom["id"]}" elementType="{atom["symbol"]}" x3="{atom["x"]:.4f}" '
                              f'y3="{atom["y"]:.4f}" z3="{atom["z"]:.4f}" mrvMap="{atom["mapping"]}"'
                              f'{atom["charge"]}{atom["multiplicity"]}{atom["isotope"]}/>'
                              for atom in atoms),
                             ('</atomArray><bondArray>',),
                             (f'<bond id="b{n}" atomRefs2="a{i} a{j}" order="{order}"{stereo}'
                              for n, (i, j, order, stereo) in enumerate(bonds, start=1)),
                             ('</bondArray>',),
                             (f'<molecule id="sg{n}" role="DataSgroup" fieldName="{t}" '
                              f'fieldData="{v.replace(">", "&gt;").replace("<", "&lt;")}" '
                              f'atomRefs="{" ".join(f"a{x}" for x in i)}" x="0" y="{n / 3}"/>'
                              for n, (i, t, v) in enumerate(cgr, start=1))))

    @staticmethod
    def _isotope_map(x, n):
        return x and f' isotope="{x}"' or ''

    @staticmethod
    def _multiplicity_map(x, n):
        if x == 3:
            return ' radical="divalent3"'
        elif x == 2:
            return ' radical="monovalent"'
        else:
            return ''

    _stereo_map = {-1: '><bondStereo>H</bondStereo></bond>', 1: '><bondStereo>W</bondStereo></bond>',
                   None: '/>'}.__getitem__
    _charge_map = {-3: ' formalCharge="-3"', -2: ' formalCharge="-2"', -1: ' formalCharge="-1"', 0: '',
                   1: ' formalCharge="1"', 2: ' formalCharge="2"', 3: ' formalCharge="3"'}.__getitem__
    _bond_map = {8: '1" queryType="Any', 4: 'A', 1: '1', 2: '2', 3: '3', 9: 's', None: '0'}.__getitem__
    __finalized = False


if find_spec('lxml'):
    from lxml.etree import iterparse, QName, tostring
    __all__ = ['MRVread', 'MRVwrite']
else:
    warn('lxml library not installed', ImportWarning)
    __all__ = ['MRVwrite']
    del MRVread
