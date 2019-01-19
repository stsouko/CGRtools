# -*- coding: utf-8 -*-
#
#  Copyright 2014-2019 Ramil Nugmanov <stsouko@live.ru>
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
from collections import defaultdict
from functools import reduce
from logging import warning
from operator import or_
from .containers import QueryContainer, QueryCGRContainer, MoleculeContainer, CGRContainer


class CGRreactor:
    def __init__(self, template, delete_atoms=False):
        pattern, atom_attrs, bond_attrs, conditional_element, is_cgr, to_delete = self.__prepare_template(template)
        self.__pattern = pattern
        self.__atom_attrs = atom_attrs
        self.__bond_attrs = bond_attrs
        self.__conditional_element = conditional_element
        self.__is_cgr = is_cgr
        self.__to_delete = delete_atoms and to_delete or set()
        self.__meta = template.meta.copy()

    def __call__(self, structure, limit=1, skip_intersection=True):
        if not isinstance(structure, (MoleculeContainer, CGRContainer)):
            raise TypeError('only Molecules and CGRs possible')

        mapping = self.__pattern.get_substructure_mapping(structure, limit)
        if limit == 1:
            if mapping:
                return self.__patcher(structure, mapping)
        elif limit > 1:
            if mapping:
                return [self.__patcher(structure, m) for m in mapping]
        elif skip_intersection:
            return self.__skip(structure, mapping)
        else:
            return (self.__patcher(structure, m) for m in mapping)

    def __skip(self, structure, mapping):
        found = set()
        for m in mapping:
            matched_atoms = set(m.values())
            if found.intersection(matched_atoms):
                continue
            found.update(matched_atoms)
            yield self.__patcher(structure, m)

    @staticmethod
    def __prepare_template(template):
        if not template.reagents or not template.products:
            raise ValueError('empty template')
        reagents = reduce(or_, template.reagents, QueryContainer())
        products = reduce(or_, template.products, QueryContainer())

        if isinstance(reagents, QueryCGRContainer):
            if isinstance(products, QueryContainer):
                products = QueryCGRContainer(products)
            is_cgr = True
        elif isinstance(products, QueryCGRContainer):
            if isinstance(reagents, QueryContainer):
                reagents = QueryCGRContainer(reagents)
            is_cgr = True
        else:
            is_cgr = False

        to_delete = set(reagents).difference(products)

        absolute_atom = defaultdict(dict)
        conditional_element = {}

        for n in set(products).difference(reagents):
            # if unique atoms in patch has variable properties exception is raised
            atom = products.atom(n)
            if atom.element is None or len(atom.element) > 1:
                raise ValueError('new atoms in patch should be static')
            elif atom.neighbors or atom.hybridization:
                warning('neighbors and hybridization for new atoms unusable')
            absolute_atom[n].update(element=atom.element[0], charge=atom.charge,
                                    isotope=atom.isotope, multiplicity=atom.multiplicity,
                                    x=atom.x, y=atom.y, z=atom.z, stereo=atom.stereo)
            if is_cgr:
                absolute_atom[n].update(p_charge=atom.p_charge, p_multiplicity=atom.p_multiplicity,
                                        p_x=atom.p_x, p_y=atom.p_y, p_z=atom.p_z, p_stereo=atom.p_stereo)

        for n in set(products).intersection(reagents):
            r_atom = reagents.atom(n)
            atom = products.atom(n)

            if atom.element:
                if len(r_atom.element) == 1:
                    absolute_atom[n]['element'] = atom.element[0]
                elif len(r_atom.element) != len(atom.element):
                    raise ValueError('element mapping invalid. possible A>A, E>A, A>E, [E]>A, [E]>[E]')
                else:
                    conditional_element[n] = dict(zip(r_atom.element, atom.element))
            if atom.isotope:
                absolute_atom[n]['isotope'] = atom.isotope

            absolute_atom[n]['charge'] = atom.charge
            absolute_atom[n]['stereo'] = atom.stereo
            if is_cgr:
                absolute_atom[n]['p_charge'] = atom.p_charge
                absolute_atom[n]['p_stereo'] = atom.p_stereo
                if atom.multiplicity or atom.p_multiplicity:
                    absolute_atom[n]['multiplicity'] = atom.multiplicity
                    absolute_atom[n]['p_multiplicity'] = atom.p_multiplicity
            elif atom.multiplicity:  # replace for fixed value
                absolute_atom[n]['multiplicity'] = atom.multiplicity

        bonds = []
        for n, m, bond in products.bonds():
            if is_cgr:
                bonds.append((n, m, {'order': bond.order, 'p_order': bond.p_order,
                                     'stereo': bond.stereo, 'p_stereo': bond.p_stereo}))
            else:
                bonds.append((n, m, {'order': bond.order, 'stereo': bond.stereo}))

        return reagents, dict(absolute_atom), bonds, conditional_element, is_cgr, to_delete

    def __patcher(self, structure, mapping):
        new = type(structure)()
        new.meta.update(self.__meta)
        to_delete = {mapping[x] for x in self.__to_delete}
        atoms = {}
        new_atoms = {}

        for n, atom in self.__atom_attrs.items():
            if n in mapping:
                n = mapping[n]
                atoms[n] = {**structure.atom(n), **atom}
            else:
                new_atoms[n] = atom
        for n, element in self.__conditional_element.items():
            n = mapping[n]
            atoms[n]['element'] = element[structure.atom(n).element]
        for n, atom in atoms.items():
            new.add_atom(atom, n)
        for n, atom in structure.atoms():  # add unmatched atoms
            if n not in atoms and n not in to_delete:
                new.add_atom(atom, n)
        for n, atom in new_atoms.items():
            mapping[n] = new.add_atom(atom)

        for n, m, bond in self.__bond_attrs:  # add patch bonds
            n = mapping[n]
            m = mapping[m]
            new.add_bond(n, m, bond)

        for n, m_bond in structure._adj.items():
            if n in to_delete:  # atoms for removing
                continue
            to_delete.add(n)
            for m, bond in m_bond.items():
                if m in to_delete or n in atoms and m in atoms:
                    continue
                new.add_bond(n, m, bond)

        # todo: calculate stereo mark based on new atom order
        return new


__all__ = ['CGRreactor']
