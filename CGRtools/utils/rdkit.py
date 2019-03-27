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
from rdkit.Chem import BondType, Atom, RWMol, SanitizeMol, Conformer
from ..containers import MoleculeContainer


def from_rdkit_molecule(data):
    """
    RDKit molecule object to MoleculeContainer converter
    """
    m = MoleculeContainer()
    atoms, mapping = [], []
    for a in data.GetAtoms():
        atom = {'element': a.GetSymbol(), 'charge': a.GetFormalCharge()}
        atoms.append(atom)
        mapping.append(a.GetAtomMapNum())

        isotope = a.GetIsotope()
        if isotope:
            atom['isotope'] = isotope
        radical = a.GetNumRadicalElectrons()
        if radical:
            atom['multiplicity'] = radical + 1

    conformers = data.GetConformers()
    if conformers:
        for atom, (x, y, z) in zip(atoms, conformers[0].GetPositions()):
            atom['x'] = x
            atom['y'] = y
            atom['z'] = z

    for atom, mapping in zip(atoms, mapping):
        a = m.add_atom(atom)
        if mapping:
            m.atom(a)._parsed_mapping = mapping

    for bond in data.GetBonds():
        m.add_bond(bond.GetBeginAtomIdx() + 1, bond.GetEndAtomIdx() + 1, _rdkit_bond_map[bond.GetBondType()])

    return m


def to_rdkit_molecule(data):
    """
    MoleculeContainer to RDKit molecule object converter
    """
    mol = RWMol()
    conf = Conformer()
    mapping = {}
    is_3d = False
    for n, a in data.atoms():
        ra = Atom(a.number)
        ra.SetAtomMapNum(n)
        if a.charge:
            ra.SetFormalCharge(a.charge)
        if a.isotope != a.common_isotope:
            ra.SetIsotope(a.isotope)
        if a.radical:
            ra.SetNumRadicalElectrons(a.radical)
        mapping[n] = m = mol.AddAtom(ra)
        conf.SetAtomPosition(m, (a.x, a.y, a.z))
        if a.z:
            is_3d = True
    if not is_3d:
        conf.Set3D(False)

    for n, m, b in data.bonds():
        mol.AddBond(mapping[n], mapping[m], _bond_map[b.order])

    mol.AddConformer(conf)
    SanitizeMol(mol)
    return mol


_rdkit_bond_map = {BondType.SINGLE: 1, BondType.DOUBLE: 2, BondType.TRIPLE: 3, BondType.AROMATIC: 4}
_bond_map = {1: BondType.SINGLE, 2: BondType.DOUBLE, 3: BondType.TRIPLE, 4: BondType.AROMATIC}

__all__ = ['from_rdkit_molecule', 'to_rdkit_molecule']
