# -*- coding: utf-8 -*-
#
#  Copyright 2017 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
from itertools import chain, count
from sys import stderr
from traceback import format_exc
from .CGRrw import CGRread, CGRwrite, fromMDL, EmptyMolecule, FinalizedFile, mendeleyset
from ..containers import MoleculeContainer, CGRContainer


class MRVread(CGRread):
    def __init__(self, file, remap=True):
        self.__file = file
        self.__data = self.__reader()
        CGRread.__init__(self, remap)
        raise Exception('NOT IMPLEMENTED')

    def read(self):
        return list(self.__data)

    def __iter__(self):
        return self.__data

    def __next__(self):
        return next(self.__data)

    def __reader(self):
        pass


class MRVwrite(CGRwrite):
    def __init__(self, file, extralabels=False, mark_to_map=False, xyz=False):
        CGRwrite.__init__(self, extralabels=extralabels, mark_to_map=mark_to_map, xyz=xyz)
        self.__file = file
        self.write = self.__init_write

    __finalized = False

    def close(self):
        if not self.__finalized:
            self.finalize()
        self.__file.close()

    def finalize(self):
        self.__file.write('</cml>')
        self.write = self.__write_adhoc
        self.__finalized = True

    @staticmethod
    def __write_adhoc(_):
        raise FinalizedFile('Writer closed')

    def __init_write(self, data):
        self.__file.write('<cml>')
        self.__write(data)
        self.write = self.__write

    def __write(self, data):
        self.__file.write('<MDocument><MChemicalStruct>')

        if isinstance(data, MoleculeContainer):
            m = self.get_formatted_cgr(data)
            self.__file.write('<molecule><propertyList>')
            for k, v in chain(m['colors'].items(), data.meta.items()):
                if '\n' in v:
                    v = '<![CDATA[%s]]>' % v
                self.__file.write('<property title="%s"><scalar>%s</scalar></property>' % (k, v))

            self.__file.write('</propertyList>')
            self.__file.write(m['CGR'])
            self.__file.write('</molecule>')
        else:
            colors = {}
            c = count(1)
            self.__file.write('<reaction>')
            for i, j in (('reagents', 'reactantList'), ('products', 'productList')):
                self.__file.write('<%s>' % j)
                for cnext, m in zip(c, data[i]):
                    m = self.get_formatted_cgr(m)
                    self.__file.write('<molecule>')
                    self.__file.write(m['CGR'])
                    self.__file.write('</molecule>')
                    colors.update({'%s.%d' % (k, cnext): v for k, v in m['colors'].items()})
                self.__file.write('</%s>' % j)

            self.__file.write('<propertyList>')
            for k, v in chain(colors.items(), data.meta.items()):
                if '\n' in v:
                    v = '<![CDATA[%s]]>' % v
                    self.__file.write('<property title="%s"><scalar>%s</scalar></property>' % (k, v))

            self.__file.write('</propertyList></reaction>')

        self.__file.write('</MChemicalStruct></MDocument>')

    @classmethod
    def _format_mol(cls, atoms, bonds, extended, cgr_dat):
        isotope, atom_query, radical = {}, {}, {}
        for i in extended:
            it, iv, ia = i['type'], i['value'], i['atom']
            if it == 'isotope':
                isotope[ia] = ' isotope="%d"' % iv
            elif it == 'atomlist':
                atom_query[ia] = ' mrvQueryProps="L%s:"' % ''.join(('!%s' % x for x in mendeleyset.difference(iv))
                                                                    if len(iv) > cls._half_table else iv)
            elif it == 'radical':
                radical[ia] = ' radical="%d"' % iv

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
                             ('<molecule id="sg{0}" role="DataSgroup" fieldName="{1[type]}" fieldData="{1[value]}" '
                              'atomRefs="{2}" x="{3[0]}" y="{3[1]}" '
                              '/>'.format(i, j, ' '.join('a%d' % x for x in j['atoms']),
                                          cls._get_position([atoms[i - 1] for i in j['atoms']]))
                              for i, j in enumerate(cgr_dat, start=1))))

    @staticmethod
    def _xyz_convert(x, y, z):
        return x * 2, y * 2, z * 2

    _stereo_map = {-1: 'H', 0: 0, 1: 'W', None: 0}
    _charge_map = {-3: -3, -2: -2, -1: -1, 0: 0, 1: 1, 2: 2, 3: 3}
    _radical_map = {2: 'monovalent', 1: 'divalent1', 3: 'divalent3'}
    __bond_map = {8: '1" queryType="Any', 4: 'A', 1: '1', 2: '2', 3: '3'}
