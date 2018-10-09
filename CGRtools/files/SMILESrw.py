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
from coho.smi import Parser
from re import split
from traceback import format_exc
from warnings import warn
from ._CGRrw import WithMixin, CGRread


class SMILESread(CGRread, WithMixin):
    """
    accept file wich consist smiles/smirks per lines.
    line should be start with smiles/smirks string and
    optionally continues with space/tab separated list of key:value [or key=value] data
    for reactions . [dot] in bonds should be used only for molecules separation
    example:
    C=C>>CC id:123 key=value\n
    """
    def __init__(self, file, *args, is_template=None, **kwargs):
        assert not is_template, 'is_tepmlate works unavailable for SMILES/SMIRKS'
        super().__init__(*args, is_template=False, **kwargs)
        super(CGRread, self).__init__(file)
        self.__parser = Parser()
        self.__data = self.__reader()

    def read(self):
        return list(self.__data)

    def __iter__(self):
        return self.__data

    def __next__(self):
        return next(self.__data)

    def __reader(self):
        for line in self._file:
            smi, *data = line.split()
            meta = {}
            for x in data:
                x = split('[=:]', x, 1)
                if len(x) == 2:
                    meta[x[0].strip()] = x[1].strip()

            if '>' in smi:
                reaction = dict(reagents=[], reactants=[], products=[], meta=meta, colors={})
                try:
                    reagents, reactants, products = smi.split('>')
                    if reagents:
                        reaction['reagents'].extend(self.__parse_molecule(x) for x in reagents.split('.'))
                    if products:
                        reaction['products'].extend(self.__parse_molecule(x) for x in products.split('.'))
                    if reactants:
                        reaction['reactants'].extend(self.__parse_molecule(x) for x in reactants.split('.'))
                except (ValueError, KeyError):
                    warn('line: %s\nconsist errors:\n%s' % (smi, format_exc()), ResourceWarning)
                else:
                    try:
                        yield self._get_reaction(reaction)
                    except Exception:
                        warn('record consist errors:\n%s' % format_exc(), ResourceWarning)
            else:
                try:
                    molecule = self.__parse_molecule(smi)
                except (ValueError, KeyError):
                    warn('line: %s\nconsist errors:\n%s' % (smi, format_exc()), ResourceWarning)
                else:
                    if meta:
                        molecule['meta'].update(meta)
                    try:
                        yield self._get_molecule(molecule)
                    except Exception:
                        warn('record consist errors:\n%s' % format_exc(), ResourceWarning)

    def __parse_molecule(self, smiles):
        self.__parser.parse(smiles)
        return {'atoms': [dict(element=a['symbol'], charge=a['charge'], map=a['aclass'] or 0,
                               x=0, y=0, z=0, isotope=0, mark='0') for a in self.__parser.atoms],
                'CGR_DAT': [dict(atoms=(n,), type='isotope', value=a['isotope'])
                            for n, a in enumerate(self.__parser.atoms, start=1) if a['isotope']],
                'bonds': [(b['a0'] + 1, b['a1'] + 1, self.__bond_map[b['order']], 0) for b in self.__parser.bonds],
                'meta': {}, 'colors': {}}

    __bond_map = {1: 1, 2: 2, 3: 3, 5: 4}
