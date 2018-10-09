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
from ..containers import MoleculeContainer


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
        reaction = dict(reagents=[], products=[], reactants=[], meta={}, colors={})
        if 'propertyList' in data and 'property' in data['propertyList']:
            meta, colors = cls.__parse_property(data['propertyList'], True)
            reaction['meta'].update(meta)
            reaction['colors'].update(colors)

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
    def __parse_property(cls, data, is_reaction=False):
        meta = defaultdict(list)
        colors = defaultdict(list)
        dp = data['property']
        for x in (dp,) if isinstance(dp, dict) else dp:
            key = x['@title']
            val = x['scalar']['$']
            col_key = key.split('.')[0] if is_reaction else key
            if col_key in ('PHTYP', 'FFTYP', 'PCTYP', 'EPTYP', 'HBONDCHG', 'CNECHG',
                           'dynPHTYP', 'dynFFTYP', 'dynPCTYP', 'dynEPTYP', 'dynHBONDCHG', 'dynCNECHG'):
                colors[key].append(val)
            elif key:
                meta[key].append(val)

        return meta, colors

    @classmethod
    def __parse_molecule(cls, data):
        molecule = dict(atoms=[], bonds=[], CGR_DAT=[], meta={}, colors={})

        if 'propertyList' in data and 'property' in data['propertyList']:
            meta, colors = cls.__parse_property(data['propertyList'])
            molecule['meta'].update(meta)
            molecule['colors'].update(colors)

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
        self.write = self.__init_write

    def close(self):
        if not self.__finalized:
            self._file.write('</cml>')
            self.__finalized = True
        super().close()

    def __init_write(self, data):
        self._file.write('<cml>')
        self.__write(data)
        self.write = self.__write

    def __write(self, data):
        self._file.write('<MDocument><MChemicalStruct>')

        if isinstance(data, MoleculeContainer):
            m = self.get_formatted_cgr(data)
            self._file.write('<molecule><propertyList>')
            for k, v in chain(m['colors'].items(), data.meta.items()):
                if '\n' in v:
                    v = '<![CDATA[%s]]>' % v
                self._file.write('<property title="%s"><scalar>%s</scalar></property>' % (k, v))

            self._file.write('</propertyList>')
            self._file.write(m['CGR'])
            self._file.write('</molecule>')
        else:
            colors = {}
            c = count(1)
            s = x = 0
            rl = len(data.reagents)
            self._file.write('<reaction>')
            for i, j in (('reagents', 'reactantList'), ('products', 'productList')):
                self._file.write('<%s>' % j)
                for cnext, m in zip(c, data[i]):
                    m = self.get_formatted_cgr(m, s)
                    if self._fix_position:
                        s = m['max_x'] + (3 if cnext == rl else 1)
                    if cnext == rl:  # get last reagent right atom position
                        x = m['max_x']
                    elif not rl and cnext == 1:
                        x = m['min_x'] - 3
                    self._file.write('<molecule>')
                    self._file.write(m['CGR'])
                    self._file.write('</molecule>')
                    colors.update({'%s.%d' % (k, cnext): v for k, v in m['colors'].items()})
                self._file.write('</%s>' % j)

            self._file.write('<arrow type="DEFAULT" x1="{:.4f}" y1="0" x2="{:.4f}" y2="0"/>'
                             '<propertyList>'.format(x + .5, x + 2.5))
            for k, v in chain(colors.items(), data.meta.items()):
                if '\n' in v:
                    v = '<![CDATA[%s]]>' % v
                    self._file.write('<property title="%s"><scalar>%s</scalar></property>' % (k, v))

            self._file.write('</propertyList></reaction>')

        self._file.write('</MChemicalStruct></MDocument>')

    @classmethod
    def _format_mol(cls, atoms, bonds, extended, cgr_dat):
        isotope, atom_query, radical = {}, {}, {}
        for i in extended:
            it, iv, ia = i['type'], i['value'], i['atom']
            if it == 'isotope':
                isotope[ia] = ' isotope="%d"' % iv
            elif it == 'atomlist':
                atom_query[ia] = ' mrvQueryProps="L%s:"' % ''.join(('!%s' % x for x in elements_set.difference(iv))
                                                                   if len(iv) > cls._half_table else
                                                                   (',%s' % x for x in iv))
            elif it == 'radical':
                radical[ia] = ' radical="%s"' % iv

        return ''.join(chain(('<atomArray>',),
                             ('<atom id="a{0}" elementType="{1[element]}" x3="{1[x]:.4f}" y3="{1[y]:.4f}" '
                              'z3="{1[z]:.4f}" mrvMap="{1[map]}" formalCharge="{1[charge]}"{2}{3}{4}{5}/>'
                              .format(i, j, radical.get(i, ''), isotope.get(i, ''), atom_query.get(i, ''),
                                      ' ISIDAmark="%s"' % j['mark'] if j['mark'] != '0' else '')
                              for i, j in enumerate(atoms, start=1)),
                             ('</atomArray><bondArray>',),
                             ('<bond id="b{0}" atomRefs2="a{1} a{2}" order="{3}"{4}'
                              .format(i, j, l, cls.__bond_map[k],
                                      '><bondStereo>%s</bondStereo></bond>' % s if s else '/>')
                              for i, (j, l, k, s) in enumerate(bonds, start=1)),
                             ('</bondArray>',),
                             ('<molecule id="sg{0}" role="DataSgroup" fieldName="{1}" fieldData="{2}" '
                              'atomRefs="{3}" x="{4[0]}" y="{4[1]}" '
                              '/>'.format(i, j['type'], j['value'].replace('>', '&gt;'),
                                          ' '.join('a%d' % x for x in j['atoms']),
                                          cls._get_position([atoms[i - 1] for i in j['atoms']]))
                              for i, j in enumerate(cgr_dat, start=1))))

    _stereo_map = {-1: 'H', 0: 0, 1: 'W', None: 0}
    _charge_map = {-3: -3, -2: -2, -1: -1, 0: 0, 1: 1, 2: 2, 3: 3}
    _radical_map = {2: 'monovalent', 1: 'divalent1', 3: 'divalent3'}
    __bond_map = {8: '1" queryType="Any', 4: 'A', 1: '1', 2: '2', 3: '3'}
    __finalized = False


__all__ = [MRVread.__name__, MRVwrite.__name__]
