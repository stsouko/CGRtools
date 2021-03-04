# -*- coding: utf-8 -*-
#
#  Copyright 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from .parser import CGRRead
from ...exceptions import NotChiral, IsChiral, ValenceError


class MDLStereo(CGRRead):
    def __init__(self, calc_cis_trans=False, ignore_stereo=False, **kwargs):
        super().__init__(**kwargs)
        self.__calc_cis_trans = calc_cis_trans
        self.__ignore_stereo = ignore_stereo

    def _convert_molecule(self, molecule, mapping):
        mol = super()._convert_molecule(molecule, mapping)
        if self.__ignore_stereo:
            return mol

        if self.__calc_cis_trans:
            mol.calculate_cis_trans_from_2d()

        stereo = [(mapping[n], mapping[m], s) for n, m, s in molecule['stereo']]
        while stereo:
            fail_stereo = []
            old_stereo = len(stereo)
            for n, m, s in stereo:
                try:
                    mol.add_wedge(n, m, s, clean_cache=False)
                except NotChiral:
                    fail_stereo.append((n, m, s))
                except IsChiral:
                    pass
                except ValenceError:
                    self._info('structure has errors, stereo data skipped')
                    mol.flush_cache()
                    break
            else:
                stereo = fail_stereo
                if len(stereo) == old_stereo:
                    break
                del mol.__dict__['_MoleculeStereo__chiral_centers']
                if self.__calc_cis_trans:
                    mol.calculate_cis_trans_from_2d(clean_cache=False)
                continue
            break
        return mol


__all__ = ['MDLStereo']
