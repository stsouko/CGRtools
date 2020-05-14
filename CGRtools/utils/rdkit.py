# -*- coding: utf-8 -*-
#
#  Copyright 2019 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from ..periodictable import Element


def from_rdkit_molecule(data):
    """
    RDKit molecule object to MoleculeContainer converter
    """
    mol = MoleculeContainer()
    parsed_mapping = mol._parsed_mapping
    mol_conformers = mol._conformers

    atoms, mapping = [], []
    for a in data.GetAtoms():
        e = Element.from_symbol(a.GetSymbol())
        isotope = a.GetIsotope()
        if isotope:
            e = e(isotope)
        else:
            e = e()
        atom = {'atom': e, 'charge': a.GetFormalCharge()}

        radical = a.GetNumRadicalElectrons()
        if radical:
            atom['is_radical'] = True

        atoms.append(atom)
        mapping.append(a.GetAtomMapNum())

    conformers = []
    c = data.GetConformers()
    if c:
        for atom, (x, y, _) in zip(atoms, c[0].GetPositions()):
            atom['xy'] = (x, y)
        for c in c:
            if c.Is3D():
                conformers.append(c.GetPositions())

    new_map = []
    for a, n in zip(atoms, mapping):
        a = mol.add_atom(**a)
        new_map.append(a)
        parsed_mapping[a] = n

    for bond in data.GetBonds():
        mol.add_bond(new_map[bond.GetBeginAtomIdx()], new_map[bond.GetEndAtomIdx()],
                     _rdkit_bond_map[bond.GetBondType()])

    for c in conformers:
        mol_conformers.append({k: tuple(v) for k, v in zip(new_map, c)})
    return mol


def to_rdkit_molecule(data):
    """
    MoleculeContainer to RDKit molecule object converter
    """
    mol = RWMol()
    mapping = {}

    for n, a in data.atoms():
        ra = Atom(a.atomic_number)
        ra.SetAtomMapNum(n)
        if a.charge:
            ra.SetFormalCharge(a.charge)
        if a.isotope:
            ra.SetIsotope(a.isotope)
        if a.is_radical:
            ra.SetNumRadicalElectrons(1)
        mapping[n] = mol.AddAtom(ra)

    for n, m, b in data.bonds():
        mol.AddBond(mapping[n], mapping[m], _bond_map[b.order])

    conf = Conformer()
    for n, a in data.atoms():
        conf.SetAtomPosition(mapping[n], (a.x, a.y, 0))
    conf.Set3D(False)
    mol.AddConformer(conf, assignId=True)

    for c in data._conformers:
        conf = Conformer()
        for n, xyz in c.items():
            conf.SetAtomPosition(mapping[n], xyz)
        mol.AddConformer(conf, assignId=True)

    SanitizeMol(mol)
    return mol


_rdkit_bond_map = {BondType.SINGLE: 1, BondType.DOUBLE: 2, BondType.TRIPLE: 3, BondType.AROMATIC: 4}
_bond_map = {1: BondType.SINGLE, 2: BondType.DOUBLE, 3: BondType.TRIPLE, 4: BondType.AROMATIC}

__all__ = ['from_rdkit_molecule', 'to_rdkit_molecule']
