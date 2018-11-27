# -*- coding: utf-8 -*-
#
#  Copyright 2017, 2018 Ramil Nugmanov <stsouko@live.ru>
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
from itertools import chain, count, repeat
from lxml.etree import iterparse, QName, tostring
from traceback import format_exc
from warnings import warn
from ._CGRrw import CGRread, CGRwrite, WithMixin, elements_set
from ..containers.common import BaseContainer


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
                text.append(tostring(element, encoding=str))
            else:
                elements_grouped[name].append(element)

            if element.tail:
                t = element.tail.strip()
                if t:
                    text.append(t)

        for element_tag, element_group in elements_grouped.items():
            if len(element_group) == 1:
                out[element_tag] = xml_dict(element_group[0])
            else:
                out[element_tag] = [xml_dict(x) for x in element_group]

    if parent_element.text:
        t = parent_element.text.strip()
        if t:
            text.insert(0, t)
    if text:
        out['$'] = ''.join(text)

    return out


class MRVread(CGRread, WithMixin):
    def __init__(self, file, *args, ignore=False, **kwargs):
        super().__init__(*args, ignore=ignore, **kwargs)
        super(CGRread, self).__init__(file, 'rb')
        self.__data = self.__reader()
        self.__ignore = ignore

    def read(self):
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
                try:
                    molecule = self.__parse_molecule(parsed['molecule'])
                except KeyError:
                    warn('molecule consist errors: %s' % format_exc(), ResourceWarning)
                else:
                    try:
                        yield self._get_molecule(molecule)
                    except Exception:
                        warn('record consist errors: %s' % format_exc(), ResourceWarning)
            elif 'reaction' in parsed and isinstance(parsed['reaction'], dict):
                try:
                    reaction = self.__parse_reaction(parsed['reaction'])
                except KeyError:
                    warn('reaction consist errors: %s' % format_exc(), ResourceWarning)
                else:
                    try:
                        yield self._get_reaction(reaction)
                    except Exception:
                        warn('record consist errors: %s' % format_exc(), ResourceWarning)
            else:
                warn('record invalid', ResourceWarning)

    @classmethod
    def __parse_reaction(cls, data):
        reaction = dict(reagents=[], products=[], reactants=[], meta={})
        if 'propertyList' in data and 'property' in data['propertyList']:
            reaction['meta'].update(cls.__parse_property(data['propertyList']))

        for tag, group in (('reactantList', 'reagents'), ('productList', 'products'), ('agentList', 'reactants')):
            if tag in data and 'molecule' in data[tag]:
                molecule = data[tag]['molecule']
                if isinstance(molecule, dict):
                    reaction[group].append(cls.__parse_molecule(molecule))
                else:
                    for m in molecule:
                        reaction[group].append(cls.__parse_molecule(m))
        return reaction

    @classmethod
    def __parse_property(cls, data):
        meta = defaultdict(list)
        dp = data['property']
        for x in (dp,) if isinstance(dp, dict) else dp:
            key = x['@title']
            val = x['scalar']['$']
            if key:
                meta[key].append(val)
        return meta

    @classmethod
    def __parse_molecule(cls, data):
        molecule = dict(atoms=[], bonds=[], CGR_DAT=[], meta={})

        if 'propertyList' in data and 'property' in data['propertyList']:
            molecule['meta'].update(cls.__parse_property(data['propertyList']))

        atom_map = {}
        if 'atom' in data['atomArray']:
            da = data['atomArray']['atom']
            for n, atom in (((1, da),) if isinstance(da, dict) else enumerate(da, start=1)):
                atom_map[atom['@id']] = n
                molecule['atoms'].append(dict(element=atom['@elementType'], isotope=0,
                                              charge=int(atom.get('@formalCharge', 0)),
                                              map=int(atom.get('@mrvMap', 0)), mark=atom.get('@ISIDAmark', '0'),
                                              x=float(atom['@x3'] if '@x3' in atom else atom['@x2']),
                                              y=float(atom['@y3'] if '@y3' in atom else atom['@y2']),
                                              z=float(atom['@z3'] if '@z3' in atom else atom.get('@z2', 0))))
                if '@isotope' in atom:
                    molecule['CGR_DAT'].append(dict(atoms=(n,), type='isotope', value=atom['@isotope']))
                if '@mrvQueryProps' in atom and atom['@mrvQueryProps'][0] == 'L':
                    _type = atom['@mrvQueryProps'][1]
                    molecule['CGR_DAT'].append(dict(atoms=(n,), type='atomlist' if _type == ',' else 'atomnotlist',
                                                    value=atom['@mrvQueryProps'][2:-1].split(_type)))
                if '@radical' in atom:
                    molecule['CGR_DAT'].append(dict(atoms=(n,), type='radical',
                                                    value=cls.__radical_map[atom['@radical']]))
        else:
            atom = data['atomArray']
            for n, (_id, el, iz, ch, mp, mk, al, rd, x, y, z) in \
                    enumerate(zip(atom['@atomID'].split(), atom['@elementType'].split(),
                                  atom['@isotope'].split() if '@isotope' in atom else repeat('0'),
                                  atom['@formalCharge'].split() if '@formalCharge' in atom else repeat(0),
                                  atom['@mrvMap'].split() if '@mrvMap' in atom else repeat(0),
                                  atom['@ISIDAmark'].split() if '@ISIDAmark' in atom else repeat('0'),
                                  atom['@mrvQueryProps'].split() if '@mrvQueryProps' in atom else repeat('0'),
                                  atom['@radical'].split() if '@radical' in atom else repeat('0'),
                                  (atom['@x3'] if '@x3' in atom else atom['@x2']).split(),
                                  (atom['@y3'] if '@y3' in atom else atom['@y2']).split(),
                                  (atom['@z3'].split() if '@z3' in atom else
                                   atom['@z2'].split() if '@z2' in atom else repeat(0))), start=1):
                atom_map[_id] = n
                molecule['atoms'].append(dict(element=el, isotope=0, charge=int(ch), map=int(mp), mark=mk,
                                              x=float(x), y=float(y), z=float(z)))
                if iz != '0':
                    molecule['CGR_DAT'].append(dict(atoms=(n,), type='isotope', value=iz))
                if al != '0' and al[0] == 'L':
                    _type = al[1]
                    molecule['CGR_DAT'].append(dict(atoms=(n,), type='atomlist' if _type == ',' else 'atomnotlist',
                                                    value=al[2:-1].split(_type)))
                if rd != '0':
                    molecule['CGR_DAT'].append(dict(atoms=(n,), type='radical', value=cls.__radical_map[rd]))

        if 'bond' in data['bondArray']:
            db = data['bondArray']['bond']
            for bond in ((db,) if isinstance(db, dict) else db):
                order = cls.__bond_map[bond['@queryType' if '@queryType' in bond else '@order']]
                a1, a2 = bond['@atomRefs2'].split()
                stereo = cls.__stereo_map.get(bond['bondStereo'].get('$'), 0) if 'bondStereo' in bond else 0
                molecule['bonds'].append((atom_map[a1], atom_map[a2], order, stereo))

        if 'molecule' in data:
            dm = data['molecule']
            for cgr_dat in ((dm,) if isinstance(dm, dict) else dm):
                if cgr_dat['@role'] == 'DataSgroup':
                    t = cgr_dat['@fieldName']
                    if t not in cls._cgr_keys:
                        continue

                    a = tuple(atom_map[x] for x in cgr_dat['@atomRefs'].split())
                    if len(a) == cls._cgr_keys[t]:
                        molecule['CGR_DAT'].append(dict(atoms=a, type=t,
                                                        value=cgr_dat['@fieldData'].replace('/', '').lower()))
        return molecule

    __bond_map = {'Any': 8, 'any': 8, 'A': 4, '1': 1, '2': 2, '3': 3}
    __radical_map = {'monovalent': '2', 'divalent': '1', 'divalent1': '1', 'divalent3': '3'}
    __stereo_map = {'H': -1, 'W': 1}


class MRVwrite(CGRwrite, WithMixin):
    def __init__(self, file, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(CGRwrite, self).__init__(file, 'w')

    def close(self):
        if not self.__finalized:
            self._file.write('</cml>')
            self.__finalized = True
        super().close()

    def write(self, data):
        self._file.write('<cml>')
        self.__write(data)
        self.write = self.__write

    def __write(self, data):
        self._file.write('<MDocument><MChemicalStruct>')

        if isinstance(data, BaseContainer):
            m = self._convert_structure(data)
            self._file.write('<molecule><propertyList>')
            for k, v in chain(m['colors'].items(), data.meta.items()):
                if '\n' in v:
                    v = f'<![CDATA[{v}]]>'
                self._file.write(f'<property title="{k}"><scalar>{v}</scalar></property>')
            self._file.write('</propertyList>')
            self._file.write(self.__format_mol(*m['structure']))
            self._file.write('</molecule>')
        else:
            colors = {}
            c = count(1)
            s = x = 0
            rl = len(data.reagents)
            self._file.write('<reaction>')
            for i, j in (('reagents', 'reactantList'), ('products', 'productList')):
                self._file.write(f'<{j}>')
                for cnext, m in zip(c, data[i]):
                    m = self._convert_structure(m, s)
                    if self._fix_position:
                        s = m['max_x'] + (3 if cnext == rl else 1)
                    if cnext == rl:  # get last reagent right atom position
                        x = m['max_x']
                    elif not rl and cnext == 1:
                        x = m['min_x'] - 3
                    self._file.write('<molecule>')
                    self._file.write(self.__format_mol(*m['structure']))
                    self._file.write('</molecule>')
                    colors.update({f'{k}.{cnext}': v for k, v in m['colors'].items()})
                self._file.write('</{j]>')

            self._file.write(f'<arrow type="DEFAULT" x1="{x + .5:.4f}" y1="0" x2="{x + 2.5:.4f}" y2="0"/>'
                             '<propertyList>')
            for k, v in chain(colors.items(), data.meta.items()):
                if '\n' in v:
                    v = f'<![CDATA[{v}]]>'
                    self._file.write(f'<property title="{k}"><scalar>{v}</scalar></property>')
            self._file.write('</propertyList></reaction>')
        self._file.write('</MChemicalStruct></MDocument>')

    @classmethod
    def __format_mol(cls, atoms, bonds, extra, cgr):
        isotope, atom_query, radical, isida, stereo = {}, {}, {}, {}, {}
        for i in extra:
            ia, it, iv = i
            if it == 'isotope':
                isotope[ia] = f' isotope="{iv}"'
            elif it == 'atomlist':
                atom_query[ia] = ' mrvQueryProps="L%s:"' % ''.join((f'!{x}' for x in elements_set.difference(iv))
                                                                   if len(iv) > cls._half_table else
                                                                   (',%s' % x for x in iv))
            elif it == 'radical':
                radical[ia] = f' radical="{iv}"'
        for n, atom in enumerate(atoms, start=1):
            if atom['mark']:
                isida[n] = f' ISIDAmark="{atom["mark"]}"'

        return ''.join(chain(('<atomArray>',),
                             (f'<atom id="a{n}" elementType="{atom["element"]}" x3="{atom["x"]:.4f}" '
                              f'y3="{atom["y"]:.4f}" z3="{atom["z"]:.4f}" mrvMap="{atom["map"]}" '
                              f'formalCharge="{atom["charge"]}"{radical.get(n, "")}{isotope.get(n, "")}'
                              f'{atom_query.get(n, "")}{isida.get(n, "")}/>'
                              for n, atom in enumerate(atoms, start=1)),
                             ('</atomArray><bondArray>',),
                             (f'<bond id="b{n}" atomRefs2="a{i} a{j}" order="{order}"{stereo}'
                              for n, (i, j, order, stereo) in enumerate(bonds, start=1)),
                             ('</bondArray>',),
                             (f'<molecule id="sg{n}" role="DataSgroup" fieldName="{t}" '
                              f'fieldData="{v.replace(">", "&gt;")}" '
                              f'atomRefs="{" ".join(f"a{x}" for x in i)}" x="0" y="{n / 3}"/>'
                              for n, (i, t, v) in enumerate(cgr, start=1))))

    _stereo_map = {-1: '><bondStereo>H</bondStereo></bond>', 1: '><bondStereo>W</bondStereo></bond>', None: '/>'}
    _charge_map = {-3: -3, -2: -2, -1: -1, 0: 0, 1: 1, 2: 2, 3: 3}
    _radical_map = {2: 'monovalent', 1: 'divalent1', 3: 'divalent3'}
    _bond_map = {8: '1" queryType="Any', 4: 'A', 1: '1', 2: '2', 3: '3', 9: 's'}
    __finalized = False


__all__ = ['MRVread', 'MRVwrite']
