# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
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
from importlib.util import find_spec
from logging import warning
from re import split
from traceback import format_exc
from warnings import warn
from ._CGRrw import WithMixin, CGRread
from ..periodictable import elements_list


class SMILESread(CGRread, WithMixin):
    """
    SMILES separated per lines files reader. works similar to opened file object. support `with` context manager.
    on initialization accept opened in text mode file, string path to file,
    pathlib.Path object or another buffered reader object.
    line should be start with SMILES string and
    optionally continues with space/tab separated list of key:value [or key=value] data.
    for reactions . [dot] in bonds should be used only for molecules separation.

    example:
    C=C>>CC id:123 key=value\n
    """
    def __init__(self, file, *args, **kwargs):
        super().__init__(*args, **kwargs)
        super(CGRread, self).__init__(file)
        self.__parser = Parser()
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
        for line in self._file:
            smi, *data = line.split()
            meta = {}
            for x in data:
                try:
                    k, v = split('[=:]', x, 1)
                    meta[k.strip()] = v.strip()
                except ValueError:
                    warning(f'invalid metadata entry: {x}')

            if '>' in smi:
                record = dict(reactants=[], reagents=[], products=[], meta=meta)
                try:
                    reactants, reagents, products = smi.split('>')
                except ValueError:
                    warning('invalid SMIRKS')
                    continue

                try:
                    if reactants:
                        for x in reactants.split('.'):
                            if not x and self._ignore:
                                warning('empty molecule ignored')
                            else:
                                record['reactants'].append(self.__parse_smiles(x))
                    if products:
                        for x in products.split('.'):
                            if not x and self._ignore:
                                warning('empty molecule ignored')
                            else:
                                record['products'].append(self.__parse_smiles(x))
                    if reagents:
                        for x in reagents.split('.'):
                            if not x and self._ignore:
                                warning('empty molecule ignored')
                            else:
                                record['reagents'].append(self.__parse_smiles(x))
                except ValueError:
                    warning(f'record consist errors:\n{format_exc()}')
                    continue

                try:
                    yield self._convert_reaction(record)
                except ValueError:
                    warning(f'record consist errors:\n{format_exc()}')
            else:
                try:
                    record = self.__parse_smiles(smi)
                except ValueError:
                    warning(f'line: {smi}\nconsist errors:\n{format_exc()}')
                    continue

                record['meta'] = meta
                try:
                    yield self._convert_structure(record)
                except ValueError:
                    warning(f'record consist errors:\n{format_exc()}')

    def __parse_smiles(self, smiles):
        self.__parser.parse(smiles)
        return {'atoms': [{'element': elements_list[a['atomic_number'] - 1], 'charge': a['charge'],
                           'mapping': a['atom_class'] or 0, 'x': 0., 'y': 0., 'z': 0., 'isotope': a['isotope'],
                           'multiplicity': None}
                          for a in self.__parser.atoms], 'extra': [], 'cgr': [],
                'bonds': [(b['atom0'], b['atom1'], self.__bond_map[b['order']]) for b in self.__parser.bonds]}

    __bond_map = {1: 1, 2: 2, 3: 3, 5: 4}


if find_spec('coho'):
    from coho.smiles import Parser
    __all__ = ['SMILESread']
else:
    warn('coho library not installed', ImportWarning)
    __all__ = []
    del SMILESread
