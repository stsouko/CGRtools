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
from itertools import takewhile, permutations
from io import StringIO, TextIOWrapper
from logging import warning
from pathlib import Path
from re import split
from traceback import format_exc
from warnings import warn
from ._CGRrw import CGRRead
from ..exceptions import IncorrectSmiles


# tokens structure:
# (type: int, value)
# types:
# 0: atom
# 1: bond
# 2: open chain (
# 3: close chain )
# 4: dot bond .
# 5: in bracket raw data []
# 6: closure number
# 7: raw closure number
# 8: aromatic atom
# 9: up down bond
# 10: dynamic bond
# 11: dynamic atom

replace_dict = {'-': 1, '=': 2, '#': 3, ':': 4, '~': 8, '.': None, '(': 2, ')': 3}
charge_dict = {'+': 1, '+1': 1, '++': 2, '+2': 2, '+3': 3, '+4': 4,
               '-': -1, '-1': -1, '--': -2, '-2': -2, '-3': -3, '-4': -4}
charge_dict = {tuple(k): v for k, v in charge_dict.items()}
charge_dict[()] = 0
dynamic_bonds = {'.>-': '0>1', '.>=': '0>2', '.>#': '0>3', '.>:': '0>4', '.>~': '0>8',
                 '->.': '1>0', '->=': '1>2', '->#': '1>3', '->:': '1>4', '->~': '1>8',
                 '=>.': '2>0', '=>-': '2>1', '=>#': '2>3', '=>:': '2>4', '=>~': '2>8',
                 '#>.': '3>0', '#>-': '3>1', '#>=': '3>2', '#>:': '3>4', '#>~': '3>8',
                 ':>.': '4>0', ':>-': '4>1', ':>=': '4>2', ':>#': '4>3', ':>~': '4>8',
                 '~>.': '8>0', '~>-': '8>1', '~>=': '8>2', '~>#': '8>3', '~>:': '8>4'}
dynamic_bonds = {tuple(k): v for k, v in dynamic_bonds.items()}

dyn_charge_dict = {'-4': -4, '-3': -3, '-2': -2, '-': -1, '0': 0, '+': 1, '++': 2, '+3': 3, '+4': 4}
tmp = {f'{i}>{j}': (x, f'c{y - x:+d}') for (i, x), (j, y) in permutations(dyn_charge_dict.items(), 2)}
dyn_charge_dict = {k: (v,) for k, v in dyn_charge_dict.items()}
dyn_charge_dict.update(tmp)

dyn_radical_dict = {(): (False,), ('*',): (True,), ('*', '>', '^'): (True, 'r1'), ('^', '>', '*'): (False, 'r1')}


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
        self._data = self.__reader()

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
        return list(iter(self))

    def __iter__(self):
        return (x for x in self._data if x is not None)

    def __next__(self):
        return next(iter(self))

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
                    yield None
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
                    yield None
                    continue

                try:
                    yield self._convert_reaction(record)
                except ValueError:
                    warning(f'record consist errors:\n{format_exc()}')
                    yield None
            else:
                try:
                    record = self.__parse_smiles(smi)
                except ValueError:
                    warning(f'line: {smi}\nconsist errors:\n{format_exc()}')
                    yield None
                    continue

                record['meta'] = meta
                try:
                    yield self._convert_structure(record)
                except ValueError:
                    warning(f'record consist errors:\n{format_exc()}')
                    yield None

    @staticmethod
    def __raw_tokenize(smiles):
        token_type = token = None
        tokens = []
        for s in smiles:
            if s == '[':  # open complex token
                if token_type == 5:  # two opened [
                    raise IncorrectSmiles('[...[')
                elif token:
                    tokens.append((token_type, token))
                elif token_type == 7:  # empty closure
                    raise IncorrectSmiles('invalid closure')
                token = []
                token_type = 5
            elif s == ']':  # close complex token
                if token_type != 5:
                    raise IncorrectSmiles(']..]')
                elif not token:
                    raise IncorrectSmiles('empty [] brackets')
                tokens.append((5, tuple(token)))
                token_type = token = None
            elif s in '()':
                if token_type == 5:
                    raise IncorrectSmiles('brackets in atom/bond token')
                elif token:
                    tokens.append((token_type, token))
                    token_type = token = None
                elif token_type == 7:  # empty closure
                    raise IncorrectSmiles('invalid closure')
                elif token_type == 2:  # barely opened
                    raise IncorrectSmiles('(( or ()')
                tokens.append((replace_dict[s], None))
            elif token_type == 5:  # grow token with brackets. skip validation
                token.append(s)
            elif s.isnumeric():  # closures
                if token_type == 7:  # % already found. collect cumber
                    if not token and s == '0':
                        raise IncorrectSmiles('number starts with 0')
                    token.append(s)
                else:
                    if s == '0':
                        raise IncorrectSmiles('number starts with 0')
                    elif token:
                        tokens.append((token_type, token))
                        token_type = token = None
                    tokens.append((6, int(s)))
            elif s == '%':
                if token:
                    tokens.append((token_type, token))
                elif token_type == 7:
                    raise IncorrectSmiles('%%')
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
                tokens.append((9, s == '/'))  # Up is true
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
                tokens.append((8, s.upper()))
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
        elif token:
            tokens.append((token_type, token))  # %closure or C or B
        return tokens

    @classmethod
    def __fix_tokens(cls, tokens):
        out = []
        for token_type, token in tokens:
            if token_type in (0, 8):  # simple atom
                out.append((token_type, {'element': token, 'charge': 0, 'isotope': None, 'is_radical': False,
                                         'mapping': 0, 'x': 0., 'y': 0., 'z': 0., 'hydrogen': 0, 'stereo': 0}))
            elif token_type == 5:
                if '>' in token:  # dynamic bond or atom
                    if len(token) == 3:  # bond only possible
                        try:
                            out.append((10, dynamic_bonds[token]))
                        except KeyError:
                            raise IncorrectSmiles('invalid dynamic bond token')
                    else:  # dynamic atom token
                        out.append((11, cls.__dynatom_parse(token)))
                elif '*' in token:  # CGR atom radical mark
                    out.append((11, cls.__dynatom_parse(token)))
                else:  # atom token
                    out.append(cls.__atom_parse(token))
            elif token_type == 7:  # composite closures folding
                out.append((6, int(''.join(token))))
            else:  # as is types: 1, 2, 3, 4, 6, 9
                out.append((token_type, token))
        return out

    @staticmethod
    def __atom_parse(token):
        isotope = ''.join(takewhile(lambda x: x.isnumeric(), token))
        if isotope:
            token = token[len(isotope):]
            if not token:
                raise IncorrectSmiles('atom token invalid')
            isotope = int(isotope)
        else:
            isotope = None

        element = token[0]
        if not element.isalpha():
            raise IncorrectSmiles('invalid atom token')
        elif len(token) == 1:
            hydrogen = charge = mapping = stereo = 0
        else:
            t1 = token[1]
            token = token[2:]
            if t1 == '@':  # only carbons stereo accepted
                if token:
                    t1 = token[0]
                    token = token[1:]
                    if t1 == '@':
                        if token:
                            t1 = token[0]
                            token = token[1:]
                            if t1 in ('H', ':') and element == 'C':
                                stereo = -1
                            else:
                                stereo = 0
                                warning('only neutral carbon stereo accepted')
                        else:
                            t1 = None
                            if element == 'C':
                                stereo = -1
                            else:
                                stereo = 0
                                warning('only neutral carbon stereo accepted')
                    elif t1 in ('H', ':') and element == 'C':
                        stereo = 1
                    else:
                        stereo = 0
                        warning('only neutral carbon stereo accepted')
                else:
                    t1 = None
                    if element == 'C':
                        stereo = 1
                    else:
                        stereo = 0
                        warning('only neutral carbon stereo accepted')
            else:
                stereo = 0

            if t1 == 'H':  # implicit H
                if token:
                    h = ''.join(takewhile(lambda x: x.isnumeric(), token))
                    if h:
                        token = token[len(h):]
                        hydrogen = int(h)
                    else:
                        hydrogen = 1
                else:
                    hydrogen = 1
            elif t1 is None:
                hydrogen = 0
            else:
                hydrogen = 0
                if t1.isalpha():
                    element += t1
                else:
                    token = (t1, *token)  # reset token

            if token:  # charge and mapping
                if ':' in token:
                    i = token.index(':')
                    try:
                        mapping = int(''.join(token[i + 1:]))
                    except ValueError:
                        raise IncorrectSmiles('invalid mapping token')
                    token = token[:i]
                else:
                    mapping = 0
                try:
                    charge = charge_dict[token]
                except KeyError:
                    raise IncorrectSmiles('charge token invalid')
            else:
                charge = mapping = 0

        if element in ('c', 'n', 'o', 'p', 's', 'as', 'se'):
            _type = 8
            element = element.capitatize()
        else:
            _type = 0
        return _type, {'element': element, 'charge': charge, 'isotope': isotope, 'is_radical': False,
                       'mapping': mapping, 'x': 0., 'y': 0., 'z': 0., 'hydrogen': hydrogen, 'stereo': stereo}

    @staticmethod
    def __dynatom_parse(token):
        # token contain * or > symbol allways
        isotope = ''.join(takewhile(lambda x: x.isnumeric(), token))
        if isotope:
            token = token[len(isotope):]
            isotope = int(isotope)
        else:
            isotope = None

        element = token[0]
        if not element.isalpha():
            raise IncorrectSmiles('invalid atom token')
        t1 = token[1]
        if t1.isalpha():  # extract element
            token = token[2:]
            element += t1
        else:
            token = token[1:]

        charge = ''.join(takewhile(lambda x: x in '+->01234', token))
        if charge:
            token = token[len(charge):]
            try:
                charge, *cgr = dyn_charge_dict[charge]
            except KeyError:
                raise IncorrectSmiles('invalid dynamic charge token')
        else:
            charge = 0
            cgr = []

        try:
            is_radical, *dyn = dyn_radical_dict[token]
        except KeyError:
            raise IncorrectSmiles('invalid dynamic radical token')
        else:
            cgr.extend(dyn)

        return {'element': element, 'charge': charge, 'isotope': isotope, 'is_radical': is_radical,
                'mapping': 0, 'x': 0., 'y': 0., 'z': 0., 'cgr': cgr}

    @staticmethod
    def __postprocess(tokens, strong_cycle=False):
        start = tokens[0][0]
        if start == 2:
            if len(tokens) == 1 or tokens[1][0] not in (0, 8, 11):
                raise IncorrectSmiles('smiles not started from (atom')
        elif start not in (0, 8, 11):
            raise IncorrectSmiles('smiles not started from atom')

        previous = None
        dead_cycles = []
        cycles_check = {}

        out = []
        for token_type, token in tokens:
            if token_type == 2:
                if previous:
                    if previous[0] != 4:
                        raise IncorrectSmiles('bond before side chain')
                    previous = None
                out.append((token_type, token))
            elif token_type == 3:
                if previous:
                    raise IncorrectSmiles('bond before closure')
                out.append((token_type, token))
            elif token_type in (1, 4, 9, 10):  # bonds. only keeping for atoms connecting
                previous = (token_type, token)

            elif token_type == 6:  # closure
                ...

            elif token_type in (0, 11):  # atoms
                if previous:
                    out.append((token_type, token, previous))  # atom with incoming bond
                    previous = None
                else:
                    out.append((token_type, token, (1, 1)))  # skipped single bond
            elif token_type == 8:  # aromatic
                if previous:
                    out.append((token_type, token, previous))  # atom with incoming bond
                    previous = None
                else:  # todo:
                    out.append((token_type, token, (1, 4)))  # skipped aromatic bond

        for n, i in enumerate(tokens):
            if isinstance(i[0], dict):  # atom
                previous = 'atom'
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

        if cycles_check:
            raise IncorrectSmiles('cycle is not finished')
        return tokens

    @classmethod
    def __parse_smiles(cls, smiles):
        if smiles[0] in r'=#:-0123456789/\~.%':
            raise IncorrectSmiles('bond or cycle can not be the first element in smiles')
        elif smiles[-1] in r'=#:-/\~.%':
            raise IncorrectSmiles('bond or % can not be the last element in smiles')

        tokens = cls.__raw_tokenize(smiles)
        tokens = cls.__fix_tokens(tokens)
        tokens = cls.__postprocess(tokens)

        atoms = []
        bonds = []
        cgr = []
        atom_num = -1
        last_num = 0
        stack = []
        cycles = {}
        previous = None
        for n, (token_type, token) in enumerate(tokens):
            if token_type == 2:  # ((((((
                stack.append(last_num)
            elif token_type == 3:  # ))))))
                last_num = stack.pop()
            elif token_type == 6:  # cycle
                if token not in cycles:
                    cycles[token] = atom_num
                else:
                    bonds.append((atom_num, cycles[token], b))
                    if f:
                        cgr.append(((atom_num, cycles[c][0]), 'dynbond', (str(b) + '>' + str(f))))
                    del cycles[c]
            else:  # atom
                c, b, f = token
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

        if stack:
            raise IncorrectSmiles('number of ( does not equal to number of )')
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
