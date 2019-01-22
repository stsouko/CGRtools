# -*- coding: utf-8 -*-
#
#  Copyright 2019 Ramil Nugmanov <stsouko@live.ru>
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
from rdkit.Chem import BondType
from ..containers import MoleculeContainer


def from_rdkit_molecule(data):
    """
    RDKit molecule object to MoleculeContainer converter
    """
    m = MoleculeContainer()
    atoms = []
    for a in data.GetAtoms():
        atom = {'element': a.GetSymbol(), 'charge': a.GetFormalCharge()}
        atoms.append(atom)

        mapping = a.GetAtomMapNum()
        if mapping:
            atom['mapping'] = mapping
        isotope = a.GetIsotope()
        if isotope:
            atom['isotope'] = isotope
        radical = a.GetNumRadicalElectrons()
        if radical:
            atom['multiplicity'] = 2 if radical == 1 else 3

    conformers = data.GetConformers()
    if conformers:
        for atom, (x, y, z) in zip(atoms, conformers[0].GetPositions()):
            atom['x'] = x
            atom['y'] = y
            atom['z'] = z

    for atom in atoms:
        m.add_atom(atom)

    for bond in data.GetBonds():
        m.add_bond(bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1, _rdkit_bond_map[bond.GetBondType()])

    return m


_rdkit_bond_map = {BondType.SINGLE: 1, BondType.DOUBLE: 2, BondType.TRIPLE: 3, BondType.AROMATIC: 4}


__all__ = ['from_rdkit_molecule']
