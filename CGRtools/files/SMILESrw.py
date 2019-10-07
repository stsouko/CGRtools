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
from importlib.util import find_spec
from io import StringIO, TextIOWrapper
from logging import warning
from pathlib import Path
from re import split
from traceback import format_exc
from warnings import warn
from ._CGRrw import CGRRead, elements_list


main_atom_list = {'B': 11, 'C': 12, 'N': 14, 'O': 16, 'P': 31, 'S': 32, 'F': 19, 'Cl': 35, 'Br': 80, 'I': 127, 'c': 12}
bond_dict = {'-': 1, '=': 2, '#': 3, ':': 4, '.': 5}

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
        # self.__parser = Parser()
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

            if '>' in smi and smi[smi.index('>') + 1] not in {'+', '-', '.', '=', '#', ':', '~'}:
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

    @staticmethod
    def description_parser(atom_desc_list, token_store, dynamic_bond, token):
        token_list = [{'element': None, 'charge': 0, 'isotope': None, 'mapping': 0, 'stereo': None, 'hydrogen': None,
                       'mult': None}, None, []]
        dynamic_charge = []
        end_number = 0
        # dynamic charge/bond
        if '>' in atom_desc_list:
            # bond
            if len(atom_desc_list) == 3 and atom_desc_list[0] in bond_dict and not dynamic_bond:
                if len(token) > 0 and isinstance(token[-1][0], dict):
                    dynamic_bond.append(bond_dict[atom_desc_list[0]])
                    dynamic_bond.append(bond_dict[atom_desc_list[-1]])
                else:
                    raise IncorrectSmiles('dynamic bond is the first element in smiles')
            # charge
            elif any(set(elements_list).intersection(atom_desc_list)) or 'H' in atom_desc_list:
                for n, i in enumerate(atom_desc_list):
                    # element
                    if i in elements_list or i == 'H':
                        if not token_store:
                            token_store = [i]
                        elif token_store[-1] in bond_dict:
                            token_list[1] = bond_dict[token_store[-1]]
                            token_store = [i]
                        else:
                            raise IncorrectSmiles('incorrect dynamic charge')
                    # number
                    elif i.isnumeric() and i != '0':
                        if len(token_store) > 0:
                            if token_store[-1] == '-' or token_store[-1] == '+':
                                if not dynamic_charge:
                                    dynamic_charge = [token_store[-1] + i]
                                    token_store.clear()
                                else:
                                    dynamic_charge.append(token_store[-1] + i)
                                    token_store.clear()
                                    end_number = n
                                    break
                            else:
                                raise IncorrectSmiles('incorrect dynamic charge')
                        else:
                            raise IncorrectSmiles('incorrect dynamic charge')
                    # >
                    elif i == '>':
                        if not token_store:
                            token_store = [i]
                        else:
                            if '-' in token_store[-1] or '+' in token_store[-1]:
                                dynamic_charge = [token_store[-1][0] + str(len(token_store[-1]))]
                                token_store = [i]
                            elif token_store[-1] == '0':
                                dynamic_charge = ['0']
                                token_store = [i]
                            else:
                                raise IncorrectSmiles('incorrect dynamic charge')
                    # charge
                    elif '+' in i or ('-' in i and len(i) > 1):
                        if token_store[-1] in elements_list or token_store[-1] == 'H':
                            token_list[0]['element'] = token_store[-1]
                            token_store = [i]
                        elif token_store[-1] == '>':
                            if len(atom_desc_list) > (n + 1):
                                if atom_desc_list[n + 1].isnumeric() and len(i) == 1:
                                    dynamic_charge.append(i + atom_desc_list[n + 1])
                                    token_store.clear()
                                    end_number = n + 1
                                    break
                                else:
                                    raise IncorrectSmiles('incorrect dynamic charge')
                            else:
                                dynamic_charge.append(i[0] + str(len(i)))
                                token_store.clear()
                                end_number = n
                        else:
                            raise IncorrectSmiles('incorrect dynamic charge')
                    # bond
                    elif i in bond_dict:
                        # -1 charge
                        if i == '-' and n > 0:
                            if token_store[-1] in elements_list or token_store[-1] == 'H':
                                token_list[0]['element'] = token_store[-1]
                                token_store = [i]
                            elif token_store[-1] == '>':
                                if len(atom_desc_list) > (n + 1):
                                    if atom_desc_list[n + 1].isnumeric():
                                        dynamic_charge.append(i + atom_desc_list[n + 1])
                                        token_store.clear()
                                        end_number = n + 1
                                        break
                                else:
                                    dynamic_charge.append(i + str(len(i)))
                                    token_store.clear()
                                    end_number = n
                            else:
                                raise IncorrectSmiles('incorrect dynamic charge')
                        # bond
                        else:
                            token_store = [i]
                    # zero number
                    elif i == '0':
                        if any(token_store):
                            if token_store[-1] in elements_list or token_store[-1] == 'H':
                                token_list[0]['element'] = token_store[-1]
                                token_store = [i]
                            elif token_store[-1] == '>':
                                dynamic_charge.append(i)
                                token_store.clear()
                                end_number = n
                                break
                            else:
                                raise IncorrectSmiles('incorrect dynamic charge')
                        else:
                            raise IncorrectSmiles('incorrect dynamic charge')
                if len(dynamic_charge) == 2 and len(atom_desc_list) == (end_number + 1):
                    token_list[0]['charge'] = []
                    token_list[0]['charge'].append(dynamic_charge[0])
                    token_list[0]['charge'].append(dynamic_charge[1])
                    dynamic_charge.clear()
                else:
                    raise IncorrectSmiles('incorrect dynamic charge')
                # dynamic bond
                if any(dynamic_bond):
                    if token_list[1] is None:
                        token_list[2].append(dynamic_bond[0])
                        token_list[2].append(dynamic_bond[-1])
                        dynamic_bond.clear()
                    else:
                        raise IncorrectSmiles('bond and dynamic bond')
                return token_list
            else:
                raise IncorrectSmiles('incorrect dynamic bond/charge')
        # radical
        elif '*' in atom_desc_list and any(set(elements_list).intersection(atom_desc_list)) and \
                atom_desc_list[-1].isnumeric():
            for n, i in enumerate(atom_desc_list):
                if i in elements_list or i == 'H':
                    if not token_store:
                        token_store = [i]
                    elif token_store[-1] in bond_dict:
                        token_list[1] = bond_dict[token_store[-1]]
                        token_store = [i]
                    else:
                        raise IncorrectSmiles('incorrect element before atom in radical')
                elif i.isnumeric():
                    if len(token_store) > 0:
                        if token_store[-1] == '-' or token_store[-1] == '+':
                            token_store.append(i)
                        elif token_store[-1] == '*' and len(atom_desc_list) == (n + 1):
                            token_list[0]['mult'] = i
                        else:
                            raise IncorrectSmiles('incorrect element before number in radical')
                    else:
                        raise IncorrectSmiles('incorrect element before number in radical')
                elif ('-' in i and n > 0) or '+' in i:
                    if len(token_store) > 0:
                        if token_store[-1] in elements_list:
                            token_list[0]['element'] = token_store[-1]
                            token_store = [i]
                        else:
                            raise IncorrectSmiles('incorrect element before charge in radical')
                    else:
                        raise IncorrectSmiles('incorrect element before charge in radical')
                elif i in bond_dict and n == 0:
                    token_store = [i]
                elif i == '*':
                    if len(token_store) > 0:
                        if token_store[-1] in elements_list:
                            token_list[0]['element'] = token_store[-1]
                            token_store = [i]
                        elif token_store[-1].isnumeric() and len(token_store) == 2 and (
                                token_store[0] == '-' or token_store[0] == '+'):
                            token_list[0]['charge'] = token_store[0] + token_store[1]
                            token_store = [i]
                        elif '-' in token_store[-1] or '+' in token_store[-1]:
                            token_list[0]['charge'] = token_store[-1][0] + str(len(token_store[-1]))
                            token_store = [i]
                        else:
                            raise IncorrectSmiles('incorrect element before * in radical')
                    else:
                        raise IncorrectSmiles('incorrect element before * in radical')
                else:
                    raise IncorrectSmiles('error in radical')
            # dynamic bond
            if any(dynamic_bond):
                if token_list[1] is None:
                    token_list[2].append(dynamic_bond[0])
                    token_list[2].append(dynamic_bond[-1])
                    dynamic_bond.clear()
                else:
                    raise IncorrectSmiles('bond and dynamic bond')
            return token_list
        # atom description
        elif any(set(elements_list).intersection(atom_desc_list)) or 'H' in atom_desc_list:
            for n, i in enumerate(atom_desc_list):
                # number
                if i.isnumeric():
                    if len(token_store) == 0:
                        token_store = [i]
                    else:
                        # charge
                        if token_store[-1] == '+':
                            token_list[0]['charge'] = token_store[-1] + i
                            token_store.clear()
                        elif token_store[-1] in bond_dict:  # bond/mapping
                            if token_list[0]['element'] in elements_list:
                                if token_store[-1] == ':':
                                    token_list[0]['mapping'] = i
                                    token_store.clear()
                                else:
                                    token_list[0]['charge'] = '-' + i
                                    token_store.clear()
                            else:
                                token_list[1] = bond_dict[token_store[-1]]
                                token_store = [i]
                        elif token_store[-1] == 'H':
                            token_list[0]['hydrogen'] = int(i)
                            token_store.clear()
                        elif '-' in token_store[-1] and len(token_store[-1]) > 1:
                            raise IncorrectSmiles('incorrect charge')
                        elif token_store[-1] in elements_list:
                            raise IncorrectSmiles('incorrect atom description')
                        else:
                            raise IncorrectSmiles('incorrect atom description')
                # charge
                elif '+' in i or ('-' in i and len(i) > 1):
                    if len(token_store) > 0:
                        if token_store[-1] == 'H':
                            token_list[0]['hydrogen'] = 1
                            token_store = [i]
                        elif '@' in token_store[-1]:
                            token_list[0]['stereo'] = token_store[-1]
                            token_store = [i]
                        else:
                            token_list[0]['element'] = token_store[-1]
                            token_store = [i]
                    else:
                        token_store = [i]
                # bond
                elif i in bond_dict:
                    if len(token_store) > 0:  # ---------------
                        if token_store[-1] == 'H':
                            token_list[0]['hydrogen'] = 1
                            token_store = [i]
                        elif '@' in token_store[-1]:
                            token_list[0]['stereo'] = token_store[-1]
                            token_store = [i]
                        elif '-' in token_store[-1] or '+' in token_store[-1]:
                            token_list[0]['charge'] = token_store[-1][0] + str(len(token_store[-1]))
                            token_store = [i]
                        else:
                            token_list[0]['element'] = token_store[-1]
                            token_store = [i]
                    else:
                        token_store = [i]
                # chiral specification
                elif i == '@' or i == '@@':
                    token_list[0]['element'] = token_store[-1]
                    token_store = [i]
                # hydrogen
                elif i == 'H':
                    if len(token_store) == 0:
                        token_store = [i]
                    elif token_store[-1].isnumeric():
                        token_list[0]['isotope'] = token_store[-1]
                        token_store = [i]
                    elif '@' in token_store[-1]:
                        token_list[0]['stereo'] = token_store[-1]
                        token_store = [i]
                    elif token_store[-1] in bond_dict:
                        token_list[1] = bond_dict[token_store[-1]]
                        token_store = [i]
                    else:
                        token_list[0]['element'] = token_store[-1]
                        token_store = [i]
                # atom
                elif i in elements_list:
                    if len(token_store) == 0:
                        token_store = [i]
                    else:
                        if token_store[-1].isnumeric():
                            token_list[0]['isotope'] = token_store[-1]
                            token_store = [i]
                        elif token_store[-1] in bond_dict:
                            token_list[1] = bond_dict[token_store[-1]]
                            token_store = [i]
                        else:
                            token_list[0]['hydrogen'] = 1
                            token_store = [i]
                else:
                    raise IncorrectSmiles('impossible element in atom description')
            # left characters in description
            if len(token_store) > 0:
                if token_store[-1] == 'H':
                    token_list[0]['hydrogen'] = 1
                elif '@' in token_store[-1]:
                    token_list[0]['stereo'] = token_store[-1]
                elif '-' in token_store[-1] or '+' in token_store[-1]:
                    token_list[0]['charge'] = token_store[-1][0] + str(len(token_store[-1]))
                elif token_store[-1] in elements_list:
                    token_list[0]['element'] = token_store[-1]
                else:
                    raise IncorrectSmiles('atom description finished with incorrect element')
            if token_list[0]['element'] is None:  # Hydrogen
                if token_list[0]['hydrogen'] == 1:
                    token_list[0]['element'] = 'H'
                    token_list[0]['hydrogen'] = None
                else:
                    raise IncorrectSmiles('atom description contains only several H')
            if any(dynamic_bond) and token_list[0]['element'] is not None:  # dynamic bond
                if token_list[1] is None:
                    token_list[2].append(dynamic_bond[0])
                    token_list[2].append(dynamic_bond[-1])
                    dynamic_bond.clear()
                else:
                    raise IncorrectSmiles('bond and dynamic bond')
            return token_list
        else:
            raise IncorrectSmiles('atom description does not contain element or dynamic bond')

    @staticmethod
    def token_adder(element, token, atom_desc_list, is_description, token_store, dynamic_bond):
        if not is_description:
            if element[-1] == ')' or element[-1] == '(':
                token.append((element))
            elif element.isnumeric():
                token.append([int(element), (bond_dict[token_store[0]] if len(token_store) > 1 else None), []])
                if any(dynamic_bond):
                    token[-1][2].append(dynamic_bond[0])
                    token[-1][2].append(dynamic_bond[-1])
                dynamic_bond.clear()
            elif not element.isnumeric() and (
                    (element >= 'a' and element <= 'z') or (element >= 'A' and element <= 'Z')):
                if element in main_atom_list:
                    token.append([{'element': element, 'charge': 0, 'isotope': None, 'mapping': 0, 'stereo': None,
                                    'hydrogen': None, 'mult': None},
                                  (None if len(token_store) == 1 else bond_dict[token_store[0]]), []])
                    if len(dynamic_bond) > 0:
                        if token_store[0] not in bond_dict:
                            token[-1][2].append(dynamic_bond[0])
                            token[-1][2].append(dynamic_bond[-1])
                        else:
                            raise IncorrectSmiles('bond and dynamic bond')
                    dynamic_bond.clear()
                else:
                    raise IncorrectSmiles('atom must be in brackets')
            elif element in bond_dict:
                atom_desc_list.append(element)
        else:
            if not element.isnumeric() and ((element >= 'a' and element <= 'z') or (element >= 'A' and element <= 'Z')):
                if element in elements_list or element == 'H':
                    atom_desc_list.append(element)
                else:
                    raise IncorrectSmiles('element not in atom_list(token_add function/True')
            else:
                atom_desc_list.append(element)

    @staticmethod
    def make_token(smiles):

        is_description = False
        token_store = []
        atom_desc_list = []
        dynamic_bond = []
        token = []

        if smiles[0] in bond_dict or smiles[0].isnumeric():
            raise IncorrectSmiles('bond or cycle can not be the first element in smiles')
        elif smiles[-1] in bond_dict:
            raise IncorrectSmiles('bond can not be the last element in smiles')

        for n, i in enumerate(smiles):

            # atom_description
            if i == '[':  # [[[[[[[[[[[[[
                if not is_description:
                    if len(token_store) > 0:
                        SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                                               dynamic_bond)
                        token_store.clear()
                    is_description = True
                else:
                    raise IncorrectSmiles('[...[')
            elif i == ']':  # ]]]]]]]]]]]
                if is_description:
                    if len(token_store) > 0:
                        SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                                               dynamic_bond)
                        token_store.clear()
                    is_description = False
                    if len(atom_desc_list) > 0:
                        returned_token = SMILESRead.description_parser(atom_desc_list, token_store, dynamic_bond, token)
                        if isinstance(returned_token, list):
                            token.append(returned_token)
                        atom_desc_list.clear()
                    else:
                        raise IncorrectSmiles('empty [] brackets')
                else:
                    raise IncorrectSmiles(']..]')
            # atom
            elif not i.isnumeric() and ((i >= 'a' and i <= 'z') or (i >= 'A' and i <= 'Z')) and i != 'H':
                if not is_description:
                    if len(token_store) == 0:  # nothing in token_store
                        token_store.append(i)
                    elif len(token_store) == 1:  # one element in token_store
                        if token_store[-1] in bond_dict:
                            token_store.append(i)
                        elif not token_store[-1].isnumeric() and ((token_store[-1] >= 'a' and token_store[-1] <= 'z')
                                                                  or
                                                                  (token_store[-1] >= 'A' and token_store[-1] <= 'Z')):
                            if (token_store[-1] + i) in main_atom_list:
                                SMILESRead.token_adder((token_store[-1] + i), token, atom_desc_list, is_description,
                                                       token_store, dynamic_bond)
                                token_store.clear()
                            elif token_store[-1] in main_atom_list:
                                SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description,
                                                       token_store, dynamic_bond)
                                token_store.clear()
                                token_store = [i]
                            else:
                                raise IncorrectSmiles('invalid element in token_store1')
                        else:
                            SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                                                   dynamic_bond)
                            token_store = [i]
                    elif len(token_store) == 2:  # two elements in token_store
                        if any(dynamic_bond):
                            raise IncorrectSmiles('bond and dynamic bond')
                        elif (token_store[-1] + i) in main_atom_list:
                            SMILESRead.token_adder((token_store[-1] + i), token, atom_desc_list, is_description,
                                                   token_store, dynamic_bond)
                            token_store.clear()
                        elif token_store[-1] in main_atom_list:
                            SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                                                   dynamic_bond)
                            token_store = [i]
                        elif token_store[-1].isnumeric():
                            SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                                                   dynamic_bond)
                            token_store = [i]
                        else:
                            raise IncorrectSmiles('invalid element in token_store')
                    else:
                        raise IncorrectSmiles('token_store contains more than three elements')

                else:  # inside atom description
                    if any(set(elements_list).intersection(atom_desc_list)):  # atom_description contains an element
                        raise IncorrectSmiles('two atoms in one description')
                    elif len(token_store) > 0:
                        if not token_store[-1].isnumeric() and ((token_store[-1] >= 'a' and token_store[-1] <= 'z') or (
                                token_store[-1] >= 'A' and token_store[-1] <= 'Z')) and token_store[-1] != 'H':
                            if (token_store[-1] + i) in elements_list:
                                atom_desc_list.append(token_store[-1] + i)
                                token_store.clear()
                            else:
                                raise IncorrectSmiles('invalid element in atom_description')
                        elif token_store[-1].isnumeric() or token_store[-1] == 'H':
                            atom_desc_list.append(token_store[-1])
                            token_store = [i]
                    else:
                        token_store = [i]
            # Hydrogen
            elif i == 'H':
                if is_description:
                    if len(token_store) > 0:
                        if (token_store[-1] in elements_list and token_store[-1] != 'H') or '@' in token_store[-1] or (
                                token_store[-1].isnumeric() and not set(elements_list).intersection(atom_desc_list)):
                            SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                                                   dynamic_bond)
                            token_store = [i]
                        else:
                            raise IncorrectSmiles('impossible element before H')
                    else:
                        token_store = [i]
                else:
                    raise IncorrectSmiles('H outside atom description')
            # brackets
            elif i == '(' or i == ')':  # (((((((((((()))))))))))))))))
                if not is_description:
                    if len(token_store) == 0:
                        token_store.append(i)
                    else:
                        if token_store[-1] not in bond_dict:
                            SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                                                   dynamic_bond)
                            token_store = [i]
                        else:
                            raise IncorrectSmiles('bond before bracket')
                else:
                    raise IncorrectSmiles('brackets in atom description')
            # bond
            elif i in bond_dict and i != '-':
                if not is_description:
                    if len(token_store) == 0:
                        token_store.append(i)
                    else:
                        if token_store[-1] not in bond_dict:
                            SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                                                   dynamic_bond)
                            token_store = [i]
                        else:
                            raise IncorrectSmiles('bond after bond')
                else:
                    if len(token_store) > 0:
                        if token_store[-1] == '>' or (
                                i == ':' and smiles[n + 1].isnumeric()):
                            SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                                                   dynamic_bond)
                            token_store = [i]
                        else:
                            raise IncorrectSmiles('incorrect element before/after bond in description')
                    else:
                        if len(atom_desc_list) == 0 and smiles[n + 1] == '>':
                            token_store = [i]
                        elif i == ':' and any(set(elements_list).intersection(atom_desc_list)) and smiles[
                            n + 1].isnumeric():
                            token_store = [i]
                        else:
                            raise IncorrectSmiles('incorrect element before/after bond in description')
            # cycle/number
            elif i.isnumeric():
                if not is_description:
                    if len(token_store) == 0:
                        if i != '0':
                            token_store = [i]
                        else:
                            raise IncorrectSmiles('0 cycle')
                    else:
                        if token_store[-1] in bond_dict:
                            token_store.append(i)
                        elif token_store[-1].isnumeric():
                            SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                                                   dynamic_bond)
                            token_store = [i]
                        elif token_store[-1] == ')':
                            raise IncorrectSmiles('cycle after bracket')
                        else:
                            if i != '0':
                                SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description,
                                                       token_store, dynamic_bond)
                                token_store = [i]
                            else:
                                raise IncorrectSmiles('0 cycle')
                else:
                    if len(token_store) == 0:
                        token_store.append(i)
                    else:
                        if token_store[-1].isnumeric():
                            if token_store[-1] != '0':
                                token_store[-1] += i
                            else:
                                raise IncorrectSmiles('number starts with 0')
                        # elif token_store[-1] == '>':
                        #     if i == '0':
                        #         token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                        #                            dynamic_bond)
                        #         token_store = [i]
                        #     else:
                        #         raise IncorrectSmiles('0 after >')
                        elif '@' in token_store[-1]:
                            raise IncorrectSmiles('number after @')
                        else:
                            SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                                               dynamic_bond)
                            token_store = [i]
            # dynamic bond
            elif i == '>':
                if is_description and (
                        token_store[-1] == '0' or '-' in token_store[-1] or token_store[-1] in bond_dict or '+' in
                        token_store[-1] or (token_store[-1].isnumeric() and (
                        atom_desc_list[-1] == '-' or atom_desc_list[-1] == '+'))):
                    SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                                           dynamic_bond)
                    token_store = [i]
                else:
                    raise IncorrectSmiles('not bond > or outside atom description')
            # charge +
            elif i == '+':
                if is_description:
                    if len(token_store) == 0:
                        if len(atom_desc_list) > 0 and atom_desc_list[-1] in elements_list:
                            token_store = [i]
                        else:
                            raise IncorrectSmiles('impossible element before charge')
                    elif token_store[-1] in elements_list or token_store[-1] == 'H' or '@' in token_store[-1] or \
                            token_store[-1].isnumeric() or token_store[-1] == '>':
                        SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                                               dynamic_bond)
                        token_store = [i]
                    elif '+' in token_store[-1]:
                        token_store[-1] += i
                    else:
                        raise IncorrectSmiles('impossible element before charge')
                else:
                    raise IncorrectSmiles('charge is outside atom description')
            # bond and charge
            elif i == '-':
                if not is_description:
                    if len(token_store) == 0:
                        token_store.append(i)
                    else:
                        if token_store[-1] not in bond_dict:
                            SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                                                   dynamic_bond)
                            token_store = [i]
                        else:
                            raise IncorrectSmiles('bond before -')
                else:
                    if len(token_store) == 0:
                        if any(atom_desc_list) and atom_desc_list[-1] in elements_list:
                            token_store = [i]
                        elif not atom_desc_list and smiles[n + 1] == '>':
                            token_store = [i]
                        else:
                            raise IncorrectSmiles('impossible element before +-')
                    elif token_store[-1] in elements_list or token_store[-1] == 'H' or '@' in token_store[-1] or \
                            token_store[-1].isnumeric() or token_store[-1] == '>':
                        SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                                               dynamic_bond)
                        token_store = [i]
                    elif '-' in token_store[-1]:
                        token_store[-1] += i
                    else:
                        raise IncorrectSmiles('impossible element before +-')
            # chiral specification
            elif i == '@':
                if is_description:
                    if len(token_store) > 0:
                        if token_store[-1] in elements_list:
                            SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                                                   dynamic_bond)
                            token_store = [i]
                        elif token_store[-1] == '@':
                            token_store[-1] += i
                        else:
                            raise IncorrectSmiles('impossible element before @')
                    elif len(atom_desc_list) > 0 and atom_desc_list[-1] in elements_list:
                        token_store = [i]
                    else:
                        raise IncorrectSmiles('impossible element before @')
                else:
                    raise IncorrectSmiles('@ is outside atom description')
            # radical
            elif i == '*':
                if is_description:
                    if len(token_store) > 0:
                        if '-' in token_store[-1] or '+' in token_store[-1] or token_store[-1] in elements_list or \
                                token_store[-1].isnumeric():
                            SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store,
                                                   dynamic_bond)
                            token_store = [i]
                        else:
                            raise IncorrectSmiles('impossible element before *')
                    elif any(atom_desc_list):
                        if atom_desc_list[-1] in elements_list:
                            token_store = [i]
                        else:
                            raise IncorrectSmiles('impossible element before *')
                    else:
                        raise IncorrectSmiles('* can not be the first element in description')
                else:
                    raise IncorrectSmiles('* outside atom_description')
            else:
                raise IncorrectSmiles('impossible element')

        if len(token_store) > 0:
            SMILESRead.token_adder(token_store[-1], token, atom_desc_list, is_description, token_store, dynamic_bond)
        if len(dynamic_bond) > 0:
            raise IncorrectSmiles('dynamic bond left')
        if is_description:
            raise IncorrectSmiles('atom description has not finished')

        return token

    @staticmethod
    def check_token(token, strong_cycle=True):
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
                c, b, f = i
                if previous == ')':  # (...)#1
                    raise IncorrectSmiles('cycle after brackets')
                elif c in dead_cycles:  # finished cycles
                    raise IncorrectSmiles('this cycle has already been finished')
                elif c not in cycles_check:
                    cycles_check[c] = [n, b]
                    previous = 'cycle'
                elif b == cycles_check[c][1]:  # close cycle
                    dead_cycles.append(c)
                    del cycles_check[c]
                    previous = 'cycle'
                elif b is None or cycles_check[c][1] is None:  # one bond defined
                    if not strong_cycle:
                        dead_cycles.append(c)
                        previous = 'cycle'
                        if b is None:
                            i[1] = cycles_check[c][1]
                        elif cycles_check[c][1] is None:
                            token[cycles_check[c][0]][1] = b
                        del cycles_check[c]
                    else:  # strong smiles
                        raise IncorrectSmiles('only one bond in cycle is defined')
                else:  # incorrect smiles
                    raise IncorrectSmiles('different bonds in cycle')
        if any(cycles_check):
            raise IncorrectSmiles('cycle is not finished')
        if bracket_numb != 0:
            raise IncorrectSmiles('number of ( does not equal to number of )')

    @staticmethod
    def __parse_smiles(smiles):

        token = SMILESRead.make_token(smiles)
        SMILESRead.check_token(token)

        atoms = []
        adjacency = []
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
                    cycles[c] = [atom_num, b]
                else:
                    adjacency.append(())
                    adjacency[-1] += (atom_num, cycles[c][0], b)
                    del cycles[c]
            else:  # atom
                c, b, f = i
                atom_num += 1
                atoms.append({'element': c['element'],
                              'charge': (c['charge'][0] if isinstance(c['charge'], list) else c['charge']),
                              'p_charge': (c['charge'][-1] if isinstance(c['charge'], list) else None),
                              'mapping': c['mapping'], 'x': 0, 'y': 0, 'z': 0, 'isotope': c['isotope'],
                              'is_radical': False})
                if atom_num > 0:
                    adjacency.append((atom_num, last_num, (f if any(f) else b)))
                last_num = atom_num

        return {'atoms': atoms,
                'bonds': adjacency,
                'atoms_lists': {}, 'cgr': [], 'query': [], 'stereo': []}

        # self.__parser.parse(smiles)
        # return {'atoms': [{'element': elements_list[a['atomic_number'] - 1], 'charge': a['charge'],
        #                    'mapping': a['atom_class'] or 0, 'x': 0., 'y': 0., 'z': 0., 'isotope': a['isotope'],
        #                    'is_radical': False}
        #                   for a in self.__parser.atoms],
        #         'bonds': [(b['atom0'], b['atom1'], self.__bond_map[b['order']]) for b in self.__parser.bonds],
        #         'atoms_lists': {}, 'cgr': [], 'query': [], 'stereo': []}


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
