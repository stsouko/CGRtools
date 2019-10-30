# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
#  Copyright 2019 Artem Mukanov <nostro32@mail.ru>
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
from io import StringIO, TextIOWrapper
from logging import warning
from pathlib import Path
from re import split
from traceback import format_exc
from warnings import warn
from ._CGRrw import CGRRead, elements_list


main_atom_list = {'B': 11, 'C': 12, 'N': 14, 'O': 16, 'P': 31, 'S': 32, 'F': 19, 'Cl': 35, 'Br': 80, 'I': 127, 'c': 12}
# tokens structure:
# (type: int, value, optional)
# types:
# 0: atom
# 1: bond
# 2: open chain (
# 3: close chain )
# 4: dot bond .
# 5: in bracket raw data []
# 6: closure number
# 7: raw closure number

replace_dict = {'-': 1, '=': 2, '#': 3, ':': 4, '~': 8, '.': None, '(': 2, ')': 3}
dynamic_radical = {'^': False, '*': True}
upper_atom_symbols = {x[0] for x in elements_list}
lower_atom_symbols = {x[1] for x in elements_list if len(x) > 1}


class IncorrectSmiles(ValueError):
    pass


class SMILESRead(CGRRead):
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
        self.__data = self.__reader()

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

            if '>' in smi and smi[smi.index('>') + 1] not in ('+', '-', '.', '=', '#', ':', '~', '*', '^'):
                record = dict(reactants=[], reagents=[], products=[], meta=meta, title='')
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

    @classmethod
    def __description_parser(cls, atom_desc_list, dynamic_bond, token):
        token_list = [{'element': None, 'charge': 0, 'p_charge': None, 'isotope': None, 'mapping': 0, 'stereo': None, 'hydrogen': None,
                       'radical': None, 'p_radical': None}, None, None]
        actual_index = -1

        # atom description
        adl = atom_desc_list
        if any(set(elements_list).intersection(adl)):
            for n, i in enumerate(adl):
                if n <= actual_index:
                    continue
                #  atom
                elif i in elements_list:
                    if i != 'H':
                        if not token_list[0]['element']:
                            if adl[n - 1].isnumeric() or adl[n - 1] in bond_dict or n == 0:
                                token_list[0]['element'] = i
                            else:
                                raise IncorrectSmiles('incorrect description')
                        else:
                            raise IncorrectSmiles('two elements in description')
                    else:
                        if adl[n - 1].isnumeric() or adl[n - 1] in bond_dict or n == 0 or adl[n - 1] in elements_list:
                            if len(adl) > (n + 1) and adl[n + 1].isnumeric() and adl[n + 1] != '0':
                                token_list[0]['hydrogen'] = int(adl[n + 1])
                                actual_index = n + 1
                            else:
                                token_list[0]['hydrogen'] = 1
                #  charge
                elif '-' in i or '+' in i or i == '0':
                    if token_list[0]['element'] or token_list[0]['hydrogen']:
                        if len(adl) > (n + 1) and adl[n + 1].isnumeric() and adl[n + 1] != '0':
                            if len(i) == 1:
                                token_list[0]['charge'] = int(i + adl[n + 1])
                                actual_index = n + 1
                            else:
                                raise IncorrectSmiles('incorrect charge')
                        else:
                            if i != '0':
                                token_list[0]['charge'] = int(i[0] + str(len(i)))
                                actual_index = n
                            else:
                                actual_index = n
                        # dynamic
                        if len(adl) > (actual_index + 1) and adl[actual_index + 1] == '>':
                            if len(adl) > (actual_index + 2) and ('-' in adl[actual_index + 2] or '+' in adl[actual_index + 2] or adl[actual_index + 2] == '0'):
                                if adl[actual_index + 2] == '0':
                                    token_list[0]['p_charge'] = 0
                                    actual_index += 2
                                elif len(adl) > (actual_index + 3) and adl[actual_index + 3].isnumeric() and adl[actual_index + 3] != '0':
                                    if len(adl[actual_index + 2]) == 1:
                                        token_list[0]['p_charge'] = int(adl[actual_index + 2] + adl[actual_index + 3])
                                        actual_index += 3
                                    else:
                                        raise IncorrectSmiles('incorrect charge')
                                else:
                                    token_list[0]['p_charge'] = int(adl[actual_index + 2][0] + str(len(adl[actual_index + 2])))
                                    actual_index += 2
                            else:
                                raise IncorrectSmiles('incorrect dynamic charge')
                    else:
                        token_list[1] = bond_dict[i]
                    continue
                #  bond
                if i in bond_dict:
                    if n == 0:
                        token_list[1] = bond_dict[i]
                    elif i == ':' and adl[n + 1].isnumeric() and len(adl) == (n + 2):
                        token_list[0]['mapping'] = adl[n + 1]
                        actual_index = n + 2
                #  number
                elif i.isnumeric():
                    if len(adl) > (n + 1) and adl[n + 1] in elements_list:
                        token_list[0]['isotope'] = int(i)
                    else:
                        raise IncorrectSmiles('incorrect number in description')
                #  radical
                elif i == '*' or i == '^':
                    if token_list[0]['element']:
                        token_list[0]['radical'] = dynamic_radical[i]
                    else:
                        raise IncorrectSmiles('incorrect radical')
                    if len(adl) > (n + 1) and adl[n + 1] == '>':
                        if adl[n + 2] in dynamic_radical and dynamic_radical[adl[n + 2]] != token_list[0]['radical']:
                            token_list[0]['p_radical'] = dynamic_radical[adl[n + 2]]
                            actual_index = n + 2
                        else:
                            raise IncorrectSmiles('incorrect dynamic radical')
                    else:
                        token_list[0]['p_radical'] = dynamic_radical[i]
                #  chiral specification
                elif i == '@' or i == '@@':
                    token_list[0]['stereo'] = i
                else:
                    continue
            #  hydrogen
            if token_list[0]['element'] is None:
                if token_list[0]['hydrogen'] == 1:
                    token_list[0]['element'] = 'H'
                    token_list[0]['hydrogen'] = None
                else:
                    raise IncorrectSmiles('atom description contains only several H')
        #  dynamic_bond
        elif '>' in adl:
            if len(adl) == 3 and adl[0] in bond_dict and not dynamic_bond:
                if len(token) > 0:
                    dynamic_bond.append(bond_dict[adl[0]])
                    dynamic_bond.append(bond_dict[adl[-1]])
                else:
                    raise IncorrectSmiles('dynamic bond is the first element in smiles')
        else:
            raise IncorrectSmiles('atom description does not contain element or dynamic bond')
        if len(dynamic_bond) > 0 and token_list[0]['element'] is not None:
            if token_list[1] is None:
                token_list[1] = dynamic_bond[0]
                token_list[2] = dynamic_bond[1]
                dynamic_bond.clear()
            else:
                raise IncorrectSmiles('bond and dynamic bond')
        if not token_list[1]:
            token_list[1] = 1
        if token_list[0]['element']:
            return token_list

    @classmethod
    def __token_adder(cls, element, token, atom_desc_list, is_description, token_store, dynamic_bond):
        if not is_description:
            if element[-1] == ')' or element[-1] == '(':
                token.append(element)
            elif element.isnumeric():
                token.append([int(element), (bond_dict[token_store[0]] if token_store[0] in bond_dict else 1), None])
                if len(dynamic_bond) > 0:
                    if token_store[0] not in bond_dict:
                        token[-1][1] = dynamic_bond[0]
                        token[-1][2] = dynamic_bond[1]
                    else:
                        raise IncorrectSmiles('bond and dynamic bond')
                dynamic_bond.clear()
            elif not element.isnumeric() and (
                    (element >= 'a' and element <= 'z') or (element >= 'A' and element <= 'Z')):
                if element in main_atom_list:
                    token.append([{'element': element, 'charge': 0, 'p_charge': None, 'isotope': None, 'mapping': 0, 'stereo': None,
                                    'hydrogen': None, 'radical': None, 'p_radical': None},
                                  (1 if len(token_store) == 1 else bond_dict[token_store[0]]), None])
                    if len(dynamic_bond) > 0:
                        if token_store[0] not in bond_dict:
                            token[-1][1] = dynamic_bond[0]
                            token[-1][2] = dynamic_bond[1]
                        else:
                            raise IncorrectSmiles('bond and dynamic bond')
                    dynamic_bond.clear()
                else:
                    raise IncorrectSmiles('atom must be in brackets')
            elif element in bond_dict:
                atom_desc_list.append(element)
        else:
            if not element.isnumeric() and ((element >= 'a' and element <= 'z') or (element >= 'A' and element <= 'Z')):
                atom_desc_list.append(element)
            else:
                atom_desc_list.append(element)

    @classmethod
    def __make_token(cls, smiles):
        if smiles[0] in r'=#:-0123456789/\~.%':
            raise IncorrectSmiles('bond or cycle can not be the first element in smiles')
        elif smiles[-1] in r'=#:-/\~.%':
            raise IncorrectSmiles('bond or % can not be the last element in smiles')

        token_type = token = None
        tokens = []
        for n, s in enumerate(smiles):
            if s == '[':  # open complex token
                if token_type == 5:
                    raise IncorrectSmiles('[...[')
                elif token_type == 7:
                    raise IncorrectSmiles('invalid closure')
                elif token:
                    tokens.append((token_type, token))
                token = []
                token_type = 5
            elif s == ']':  # close complex token
                if token_type != 5:
                    raise IncorrectSmiles(']..]')
                elif not token:
                    raise IncorrectSmiles('empty [] brackets')
                tokens.append((5, token))
                token_type = token = None
            elif s in '()':
                if token_type == 5:
                    raise IncorrectSmiles('brackets in atom/bond token')
                elif token_type == 7:
                    raise IncorrectSmiles('invalid closure')
                elif token:
                    tokens.append((token_type, token))
                    token_type = token = None
                tokens.append((replace_dict[s], None))
            elif token_type == 5:  # grow token with brackets. skip validation
                token.append(s)
            elif s.isnumeric():  # closures
                if token_type == 7:  # % already found. collect cumber
                    if not token and s == '0':
                        raise IncorrectSmiles('number starts with 0')
                    token.append(s)
                elif s == '0':
                    raise IncorrectSmiles('number starts with 0')
                elif token:
                    tokens.append((token_type, token))
                    token_type = token = None
                tokens.append((6, int(s)))
            elif token_type == 7:
                raise IncorrectSmiles('invalid closure')
            elif s == '%':
                if token:
                    tokens.append((token_type, token))
                token_type = 7
                token = []
            elif s in '=#:-~':  # bonds found
                if token:
                    tokens.append((token_type, token))
                    token_type = token = None
                tokens.append((1, replace_dict[s]))
            elif s in r'\/':
                if token:
                    tokens.append((token_type, token))
                    token_type = token = None
                tokens.append((1, 1, s == '/'))  # Up is true
            elif s == '.':
                if token:
                    tokens.append((token_type, token))
                    token_type = token = None
                tokens.append((4, None))
            elif s in 'NOPSFI':  # organic atoms
                if token:
                    tokens.append((token_type, token))
                    token_type = token = None
                tokens.append((0, s))
            elif s in 'cnops':  # aromatic ring atom
                if token:
                    tokens.append((token_type, token))
                    token_type = token = None
                tokens.append((0, s.upper()))
            elif s in 'CB':  # flag possible Cl or Br
                if token:
                    tokens.append((token_type, token))
                token_type = 0
                token = s
            elif token_type == 0:
                if s == 'l':
                    if token == 'C':
                        tokens.append((0, 'Cl'))
                        token_type = token = None
                    else:
                        raise IncorrectSmiles('invalid element Bl')
                elif s == 'r':
                    if token == 'B':
                        tokens.append((0, 'Br'))
                        token_type = token = None
                    else:
                        raise IncorrectSmiles('invalid smiles for Cr')
                else:
                    raise IncorrectSmiles('invalid smiles')
            else:
                raise IncorrectSmiles('invalid smiles')

        if token_type == 5:
            raise IncorrectSmiles('atom description has not finished')
        return tokens

    @staticmethod
    def __check_token(token, strong_cycle=False):
        previous = ''
        dead_cycles = []
        cycles_check = {}
        bracket_numb = 0
        for n, i in enumerate(token):
            if isinstance(i[0], dict):  # atom
                previous = 'atom'
            elif i == '(':  # (((((((
                if previous == '(':
                    raise IncorrectSmiles('((')
                else:
                    previous = '('
                    bracket_numb += 1
            elif i == ')':  # ))))))
                if previous == '(':
                    raise IncorrectSmiles('()')
                else:
                    previous = ')'
                    bracket_numb -= 1
            else:  # cycles
                if previous == ')':  # (...)#1
                    raise IncorrectSmiles('cycle after brackets')
                elif i[0] in dead_cycles:  # finished cycles
                    raise IncorrectSmiles('this cycle has already been finished')
                elif i[0] not in cycles_check:
                    cycles_check[i[0]] = [n, i[1], i[2]]
                    previous = 'cycle'
                elif i[1] == cycles_check[i[0]][1]:  # close cycle
                    if i[2] != cycles_check[i[0]][2]:
                        if i[2] is None:
                            i[2] = cycles_check[i[0]][2]
                        elif cycles_check[i[0]][2] is None:
                            cycles_check[i[0]][2] = i[2]
                        else:
                            raise IncorrectSmiles('different dynamic bond')
                    dead_cycles.append(i[0])
                    del cycles_check[i[0]]
                    previous = 'cycle'
                elif i[1] == 1 and not i[2]:
                    if cycles_check[i[0]][2]:
                        i[1] = cycles_check[i[0]][1]
                        i[2] = cycles_check[i[0]][2]
                    else:
                        if not strong_cycle:
                            i[1] = cycles_check[i[0]][1]
                        else:
                            raise IncorrectSmiles('only one bond in cycle is defined')
                    dead_cycles.append(i[0])
                    del cycles_check[i[0]]
                    previous = 'cycle'
                elif cycles_check[i[0]][1] == 1 and not cycles_check[i[0]][2]:
                    if i[2]:
                        cycles_check[i[0]][1] = i[1]
                        cycles_check[i[0]][2] = i[2]
                    else:
                        if not strong_cycle:
                            cycles_check[i[0]][1] = i[1]
                        else:
                            raise IncorrectSmiles('only one bond in cycle is defined')
                    dead_cycles.append(i[0])
                    del cycles_check[i[0]]
                    previous = 'cycle'
                else:
                    raise IncorrectSmiles('different bonds in cycle')
        if any(cycles_check):
            raise IncorrectSmiles('cycle is not finished')
        if bracket_numb != 0:
            raise IncorrectSmiles('number of ( does not equal to number of )')
        return token

    @classmethod
    def __parse_smiles(cls, smiles):

        token = cls.__make_token(smiles)
        token = cls.__check_token(token)

        atoms = []
        bonds = []
        cgr = []
        atom_num = -1
        last_num = 0
        stack = []
        cycles = {}
        for i in token:
            if i == '(':  # ((((((
                stack.append(last_num)
            elif i == ')':  # ))))))
                last_num = stack.pop()
            elif isinstance(i[0], int):  # cycle
                c, b, f = i
                if c not in cycles:
                    cycles[c] = [atom_num, b, f]
                else:
                    bonds.append((atom_num, cycles[c][0], b))
                    if f:
                        cgr.append(((atom_num, cycles[c][0]), 'dynbond', (str(b) + '>' + str(f))))
                    del cycles[c]
            else:  # atom
                c, b, f = i
                atom_num += 1
                atoms.append({'element': c['element'],
                              'charge': (c['charge'][0] if isinstance(c['charge'], list) else c['charge']),
                              'mapping': c['mapping'], 'x': 0., 'y': 0., 'z': 0., 'isotope': c['isotope'],
                              'is_radical': (c['radical'] if c['radical'] is not None else False)})
                if atom_num > 0:
                    if b != 8 and not f:
                        bonds.append((atom_num, last_num, b))
                if c['p_charge'] is not None:
                    cgr.append(((atom_num,), 'dynatom', ('c' + str(int(c['p_charge']) - int(c['charge'])))))
                if f:
                    cgr.append(((atom_num, last_num), 'dynbond', (str(b) + '>' + str(f))))
                    bonds.append((atom_num, last_num, 8))
                if c['radical'] is False or c['p_radical'] is False:
                    cgr.append(((atom_num,), 'dynatom', 'r1'))
                last_num = atom_num

        return {'atoms': atoms,
                'bonds': bonds,
                'atoms_lists': {}, 'cgr': cgr, 'query': [], 'stereo': [], 'title': ''}


class SMILESread:
    def __init__(self, *args, **kwargs):
        warn('SMILESread deprecated. Use SMILESRead instead', DeprecationWarning)
        warning('SMILESread deprecated. Use SMILESRead instead')
        self.__obj = SMILESRead(*args, **kwargs)

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


__all__ = ['SMILESRead', 'SMILESread']
