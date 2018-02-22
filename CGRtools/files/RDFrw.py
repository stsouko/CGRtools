# -*- coding: utf-8 -*-
#
#  Copyright 2014-2018 Ramil Nugmanov <stsouko@live.ru>
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
from itertools import chain
from sys import stderr
from time import strftime
from traceback import format_exc
from ._CGRrw import fromMDL, WithMixin
from ._MDLrw import MOLwrite, MOLread
from ..containers import MoleculeContainer
from ..exceptions import EmptyMolecule, InvalidData


class RDFread(MOLread, WithMixin):
    def __init__(self, file, *args, ignore=False, **kwargs):
        WithMixin.__init__(self, file)
        MOLread.__init__(self, *args, **kwargs)
        self.__data = self.__reader()
        self.__ignore = ignore

    def read(self):
        return list(self.__data)

    def __iter__(self):
        return self.__data

    def __next__(self):
        return next(self.__data)

    def __reader(self):
        ir = im = atomcount = bondcount = n = reagents = products = spr = molcount = -1
        failkey = isreaction = True
        reaction = molecule = mkey = None
        for n, line in enumerate(self._file):
            if failkey and not line.startswith(("$RFMT", "$MFMT")):
                continue
            elif line.startswith("$RFMT"):
                if reaction:
                    try:
                        yield self._get_reaction(reaction) if isreaction else self._get_molecule(reaction)
                    except Exception:
                        print('line %d\n previous record consist errors: %s' % (n, format_exc()), file=stderr)
                isreaction = True
                ir = n + 5
                failkey = False
                reaction = mkey = None
            elif line.startswith("$MFMT"):
                if reaction:
                    try:
                        yield self._get_reaction(reaction) if isreaction else self._get_molecule(reaction)
                    except Exception:
                        print('line %d\n previous record consist errors: %s' % (n, format_exc()), file=stderr)
                reaction = dict(reagents=[], products=[], reactants=[], meta={}, colors={})
                reagents, products, spr, molcount = 1, 1, 1, 1
                mkey = None
                failkey = isreaction = False
                im = n + 4
                ir = -1
            elif n == ir:
                try:
                    reagents = int(line[:3])
                    products = int(line[3:6]) + reagents
                    spr = int(line[6:].rstrip() or 0) + products
                    reaction = dict(reagents=[], products=[], reactants=[], meta={}, colors={})
                    molcount = 0
                except ValueError:
                    failkey = True
                    reaction = molecule = None
                    print('line %d\n\n%s\n consist errors: %s' % (n, line, format_exc()), file=stderr)
            elif reaction:
                if line.startswith("$MOL"):
                    try:
                        if molcount == spr:
                            raise InvalidData('More then defined molecules')
                    except InvalidData:
                        failkey = True
                        reaction = molecule = None
                        print('line %d\n\n%s\n consist errors: %s' % (n, line, format_exc()), file=stderr)

                    molcount += 1
                    im = n + 4
                elif n == im:
                    try:
                        atoms = int(line[0:3])
                        if not atoms:
                            if self.__ignore:
                                continue
                            raise EmptyMolecule('Molecule without atoms')
                        atomcount = atoms + im
                        bondcount = int(line[3:6]) + atomcount
                        molecule = dict(atoms=[], bonds=[], CGR_DAT=[])
                    except (EmptyMolecule, ValueError):
                        failkey = True
                        reaction = molecule = None
                        print('line %d\n\n%s\n consist errors: %s' % (n, line, format_exc()), file=stderr)
                elif molecule:
                    if n <= atomcount:
                        try:
                            molecule['atoms'].append(dict(element=line[31:34].strip(), isotope=int(line[34:36]),
                                                          charge=fromMDL[int(line[38:39])], map=int(line[60:63]),
                                                          mark=line[54:57].strip(), x=float(line[:10]),
                                                          y=float(line[10:20]), z=float(line[20:30])))
                        except ValueError:
                            failkey = True
                            reaction = molecule = None
                            print('line %d\n\n%s\n consist errors: %s' % (n, line, format_exc()), file=stderr)
                    elif n <= bondcount:
                        try:
                            molecule['bonds'].append((int(line[:3]), int(line[3:6]), int(line[6:9]), int(line[9:12])))
                        except ValueError:
                            failkey = True
                            reaction = molecule = None
                            print('line %d\n\n%s\n consist errors: %s' % (n, line, format_exc()), file=stderr)
                    elif line.startswith("M  END"):
                        molecule['CGR_DAT'] = self._get_collected()
                        if molcount <= reagents:
                            reaction['reagents'].append(molecule)
                        elif molcount <= products:
                            reaction['products'].append(molecule)
                        else:
                            reaction['reactants'].append(molecule)
                        molecule = None
                    else:
                        try:
                            self._collect(line)
                        except ValueError:
                            self._flush_collected()
                            failkey = True
                            reaction = molecule = None
                            print('line %d\n\n%s\n consist errors: %s' % (n, line, format_exc()), file=stderr)

                elif line.startswith('$DTYPE'):
                    mkey = line[7:].strip()
                    if mkey.split('.')[0] in ('PHTYP', 'FFTYP', 'PCTYP', 'EPTYP', 'HBONDCHG', 'CNECHG', 'dynPHTYP',
                                              'dynFFTYP', 'dynPCTYP', 'dynEPTYP', 'dynHBONDCHG', 'dynCNECHG'):
                        target = 'colors'
                    elif mkey:
                        target = 'meta'
                    else:
                        continue
                    reaction[target][mkey] = []
                elif mkey:
                    data = line.lstrip("$DATUM").strip()
                    if data:
                        reaction[target][mkey].append(data)

        if reaction:
            try:
                yield self._get_reaction(reaction) if isreaction else self._get_molecule(reaction)
            except Exception:
                print('line %d\n previous record consist errors: %s' % (n, format_exc()), file=stderr)

    def _get_molecule(self, reaction):
        molecule = reaction['reagents'][0]
        molecule['meta'] = reaction['meta']
        molecule['colors'] = reaction['colors']
        return super()._get_molecule(molecule)


class RDFwrite(MOLwrite, WithMixin):
    def __init__(self, file, *args, **kwargs):
        WithMixin.__init__(self, file, 'w')
        MOLwrite.__init__(self, *args, **kwargs)
        self.write = self.__init_write

    def __init_write(self, data):
        self._file.write(strftime("$RDFILE 1\n$DATM    %m/%d/%y %H:%M\n"))
        self.__write(data)
        self.write = self.__write

    def __write(self, data):
        if isinstance(data, MoleculeContainer):
            m = self.get_formatted_cgr(data)
            self._file.write('$MFMT\n')
            self._file.write(m['CGR'])
            self._file.write("M  END\n")
            colors = m['colors']
        else:
            self._file.write('$RFMT\n$RXN\n\n  CGRtools. (c) Dr. Ramil I. Nugmanov\n\n%3d%3d\n' %
                             (len(data.reagents), len(data.products)))
            colors = {}
            for cnext, m in enumerate(chain(data.reagents + data.products), start=1):
                m = self.get_formatted_cgr(m)
                self._file.write('$MOL\n')
                self._file.write(m['CGR'])
                self._file.write("M  END\n")
                colors.update({'%s.%d' % (k, cnext): v for k, v in m['colors'].items()})

        for p in chain(colors.items(), data.meta.items()):
            self._file.write('$DTYPE %s\n$DATUM %s\n' % p)


__all__ = [RDFread.__name__, RDFwrite.__name__]
