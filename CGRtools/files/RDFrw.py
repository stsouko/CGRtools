# -*- coding: utf-8 -*-
#
#  Copyright 2014-2017 Ramil Nugmanov <stsouko@live.ru>
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
from .CGRrw import CGRread, CGRwrite, fromMDL, EmptyMolecule, FinalizedFile, InvalidData
from ..containers import MoleculeContainer


class RDFread(CGRread):
    def __init__(self, file, remap=True, ignore=False):
        self.__file = file
        self.__data = self.__reader()
        CGRread.__init__(self, remap, ignore=ignore)

    def read(self):
        return list(self.__data)

    def __iter__(self):
        return self.__data

    def __next__(self):
        return next(self.__data)

    def __reader(self):
        ir = im = atomcount = bondcount = n = substrats = products = spr = molcount = -1
        failkey = isreaction = True
        reaction = molecule = mkey = None
        for n, line in enumerate(self.__file):
            if failkey and not line.startswith(("$RFMT", "$MFMT")):
                continue
            elif line.startswith("$RFMT"):
                if reaction:
                    try:
                        yield self._get_reaction(reaction) if isreaction else self._get_molecule(reaction)
                    except:
                        print('line %d\n previous record consist errors: %s' % (n, format_exc()), file=stderr)
                isreaction = True
                ir = n + 5
                failkey = False
                reaction = mkey = None
            elif line.startswith("$MFMT"):
                if reaction:
                    try:
                        yield self._get_reaction(reaction) if isreaction else self._get_molecule(reaction)
                    except:
                        print('line %d\n previous record consist errors: %s' % (n, format_exc()), file=stderr)
                reaction = dict(substrats=[], products=[], reactants=[], meta={}, colors={})
                substrats, products, spr, molcount = 1, 1, 1, 1
                mkey = None
                failkey = isreaction = False
                im = n + 4
                ir = -1
            elif n == ir:
                try:
                    substrats = int(line[:3])
                    products = int(line[3:6]) + substrats
                    spr = int(line[6:].rstrip() or 0) + products
                    reaction = dict(substrats=[], products=[], reactants=[], meta={}, colors={})
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
                            raise EmptyMolecule('Molecule without atoms')
                        atomcount = atoms + im
                        bondcount = int(line[3:6]) + atomcount
                        molecule = dict(atoms=[], bonds=[], CGR_DAT={})
                    except (EmptyMolecule, ValueError):
                        failkey = True
                        reaction = molecule = None
                        print('line %d\n\n%s\n consist errors: %s' % (n, line, format_exc()), file=stderr)
                elif molecule:
                    if n <= atomcount:
                        try:
                            molecule['atoms'].append(dict(element=line[31:34].strip(), isotope=int(line[34:36]),
                                                          charge=fromMDL[int(line[38:39])],
                                                          map=int(line[60:63]), mark=line[54:57].strip(),
                                                          x=float(line[:10]), y=float(line[10:20]), z=float(line[20:30])))
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
                        if molcount <= substrats:
                            reaction['substrats'].append(molecule)
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
        else:
            if reaction:
                try:
                    yield self._get_reaction(reaction) if isreaction else self._get_molecule(reaction)
                except:
                    print('line %d\n previous record consist errors: %s' % (n, format_exc()), file=stderr)

    def _get_molecule(self, reaction):
        molecule = reaction['substrats'][0]
        molecule['meta'] = reaction['meta']
        molecule['colors'] = reaction['colors']
        return super()._get_molecule(molecule)


class RDFwrite(CGRwrite):
    def __init__(self, file, extralabels=False, mark_to_map=False, xyz=False):
        CGRwrite.__init__(self, extralabels=extralabels, mark_to_map=mark_to_map, xyz=xyz)
        self.__file = file
        self.write = self.__init_write

    def close(self):
        self.write = self.__write_adhoc
        self.__file.close()

    @staticmethod
    def __write_adhoc(_):
        raise FinalizedFile('Writer closed')

    def __init_write(self, data):
        self.__file.write(strftime("$RDFILE 1\n$DATM    %m/%d/%y %H:%M\n"))
        self.__write(data)
        self.write = self.__write

    def __write(self, data):
        if isinstance(data, MoleculeContainer):
            m = self.get_formatted_cgr(data)
            self.__file.write('$MFMT\n')
            self.__file.write(m['CGR'])
            self.__file.write("M  END\n")
            colors = m['colors']
        else:
            self.__file.write('$RFMT\n$RXN\n\n  CGRtools. (c) Dr. Ramil I. Nugmanov\n\n%3d%3d\n' %
                              (len(data.substrats), len(data.products)))
            colors = {}
            for cnext, m in enumerate(chain(data.substrats + data.products), start=1):
                m = self.get_formatted_cgr(m)
                self.__file.write('$MOL\n')
                self.__file.write(m['CGR'])
                self.__file.write("M  END\n")
                colors.update({'%s.%d' % (k, cnext): v for k, v in m['colors'].items()})

        for p in chain(colors.items(), data.meta.items()):
            self.__file.write('$DTYPE %s\n$DATUM %s\n' % p)
