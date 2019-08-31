# -*- coding: utf-8 -*-
#
#  Copyright 2014-2019 Ramil Nugmanov <stsouko@live.ru>
#  Copyright 2019 Adelia Fatykhova <adelik21979@gmail.com>
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
from itertools import chain, count, islice, permutations, product
from logging import warning, info
from operator import or_
from typing import Union
from .containers import QueryContainer, QueryCGRContainer, MoleculeContainer, CGRContainer, ReactionContainer
from .periodictable import Element, DynamicElement


class CGRReactor:
    """
    CGR based editor for CGRs and molecules.
    generates transformation from input CGR/molecule
    using template (CGRtools ReactionContainer)
    -----------------------------------------------------
    input: transformation template, CGR/molecule to edit
    output: generator of edited CGRs/molecules
    -----------------------------------------------------
    template should contain one reactant and one product:
    CGRReactor allows only 1 -> 1 transformation

    CGRReactor init prepares the reactor container:

    >> reactor = CGRReactor(template, delete_atoms=True)

    CGRReactor calling transforms reactants to products and
    returns generator of all possible products.
    """
    def __init__(self, template: ReactionContainer, delete_atoms: bool = False):
        """
        :param template: CGRtools ReactionContainer
        :param delete_atoms: if True atoms exists in reactant but
                            not exists in product will be removed
        """
        pattern, elements, atom_attrs, bond_attrs, is_cgr, to_delete = self.__prepare_template(template)
        self.__pattern = pattern
        self.__elements = elements
        self.__atom_attrs = atom_attrs
        self.__bond_attrs = bond_attrs
        self.__is_cgr = is_cgr
        self.__to_delete = delete_atoms and to_delete or set()
        self.__meta = template.meta.copy()

    def __call__(self, structure: Union[MoleculeContainer, CGRContainer], automorphism_filter: bool = True):
        if not isinstance(structure, (MoleculeContainer, CGRContainer)):
            raise TypeError('only Molecules and CGRs possible')

        for mapping in self.__pattern.get_mapping(structure, automorphism_filter=automorphism_filter):
            yield self._patcher(structure, mapping)

    @staticmethod
    def __prepare_template(template):
        if not template.reactants or not template.products:
            raise ValueError('empty template')
        if any(isinstance(x, (CGRContainer, QueryCGRContainer)) for x in chain(template.reactants, template.products)):
            reactants = reduce(or_, template.reactants, QueryCGRContainer())
            products = reduce(or_, template.products, QueryCGRContainer())
            is_cgr = True
            e = DynamicElement
        else:
            reactants = reduce(or_, template.reactants, QueryContainer())
            products = reduce(or_, template.products, QueryContainer())
            is_cgr = False
            e = Element

        to_delete = set(reactants).difference(products)
        atoms = defaultdict(dict)
        elements = {}
        for n, atom in products.atoms():
            if atom.neighbors or atom.hybridization:
                warning('neighbors and hybridization for new atoms unusable')
            atoms[n].update(charge=atom.charge, is_radical=atom.is_radical)
            elements[n] = e.from_atomic_number(atom.atomic_number)(atom.isotope)
            if n not in reactants:
                atoms[n]['xy'] = atom.xy
            if is_cgr:
                atoms[n].update(p_is_radical=atom.p_is_radical, p_charge=atom.p_charge)

        bonds = list(products.bonds())
        return reactants, elements, dict(atoms), bonds, is_cgr, to_delete

    def _patcher(self, structure, mapping):
        mapping = mapping.copy()
        elements = self.__elements

        plane = structure._plane
        bonds = structure._bonds
        charges = structure._charges
        radicals = structure._radicals
        if self.__is_cgr:
            p_charges = structure._p_charges
            p_radicals = structure._p_radicals

        new = structure.__class__()
        new.meta.update(self.__meta)

        to_delete = {mapping[x] for x in self.__to_delete}
        if to_delete:
            # if deleted atoms have another path to remain fragment, the path is preserved
            remain = set(mapping.values()).difference(to_delete)
            delete, global_seen = set(), set()
            for x in to_delete:
                for n in bonds[x]:
                    if n in global_seen or n in remain:
                        continue
                    seen = {n}
                    global_seen.add(n)
                    stack = [x for x in bonds[n] if x not in global_seen]
                    while stack:
                        current = stack.pop()
                        if current in remain:
                            break
                        if current in to_delete:
                            continue
                        seen.add(current)
                        global_seen.add(current)
                        stack.extend([x for x in bonds[current] if x not in global_seen])
                    else:
                        delete.update(seen)

            to_delete.update(delete)

        new_atoms = {}
        for n, atom in self.__atom_attrs.items():
            if n in mapping:  # add matched atoms
                m = mapping[n]
                new.add_atom(elements[n].copy(), m, xy=plane[m], **atom)
            else:  # new atoms
                new_atoms[n] = {'atom': elements[n].copy(), **atom}

        old_atoms = set(new._atoms)
        if self.__is_cgr:
            for n, atom in structure.atoms():  # add unmatched atoms
                if n not in old_atoms and n not in to_delete:
                    new.add_atom(atom.copy(), n, charge=charges[n], is_radical=radicals[n], xy=plane[n],
                                 p_is_radical=p_radicals[n], p_charge=p_charges[n])
        else:
            for n, atom in structure.atoms():  # add unmatched atoms
                if n not in old_atoms and n not in to_delete:
                    new.add_atom(atom.copy(), n, charge=charges[n], is_radical=radicals[n],  xy=plane[n])

        for n, atom in new_atoms.items():
            mapping[n] = new.add_atom(**atom)

        for n, m, bond in self.__bond_attrs:  # add patch bonds
            n = mapping[n]
            m = mapping[m]
            new.add_bond(n, m, bond)

        for n, m_bond in bonds.items():
            if n in to_delete:  # atoms for removing
                continue
            to_delete.add(n)
            for m, bond in m_bond.items():
                if m in to_delete or n in old_atoms and m in old_atoms:
                    continue
                new.add_bond(n, m, bond)

        # todo: calculate stereo mark based on new atom order
        return new


class Reactor:
    """
    CGR based reactor for molecules/queries.
    generates reaction from input queries/molecules using
    transformation template (CGRtools ReactionContainer).
    -----------------------------------------------------
    input: transformation template, list of reactants
    output: reaction or list or generator of reactions
    -----------------------------------------------------
    reactor allows only this reaction transformations:
         ONE to ONE   # 1 -> 1
         ONE to MANY  # 1 -> 2
         MANY to ONE  # 3 -> 1
         MANY to MANY # 2 -> 2 (equal)

    reactor init prepares the reactor container:

    >> reactor = CGRreactor(template, delete_atoms=True)

    reactor calling transforms reactants to products and
    returns generator of reaction transformations with all
    possible products if limit=0, one reaction if limit=1,
    else limited to number list of reactions:

    >> reactions = reactor(structure, limit=0)  # generator
    >> reaction = reactor(structure, limit=1)   # one reaction
    >> reactions = reactor(structure, limit=5)  # list with 5 reactions

    """
    def __init__(self, template, delete_atoms=False):
        """
        :param template: CGRtools ReactionContainer
        :param delete_atoms: if True atoms exists in reactants but
                            not exists in products will be removed
        """
        reactants, products = template.reactants, template.products
        if not reactants or not products:
            raise ValueError('empty template')
        if any(not isinstance(structure, (MoleculeContainer, QueryContainer))
               for structure in chain(reactants, products)):
            raise TypeError('only Molecules and Queries possible')
        self.__split = False
        self.__single = False
        if len(reactants) == 1:
            if len(products) > 1:
                self.__split = True
            self.__single = True
        elif len(reactants) == len(products):
            self.__split = True
        elif len(products) != 1:
            raise ValueError('Only reactions with '
                             'ONE to ONE, '
                             'ONE to MANY, '
                             'MANY to ONE and '
                             'MANY to MANY (EQUAL) molecules allowed')
        self.__reactor = CGRreactor(template, delete_atoms)
        self.__patterns = [QueryContainer(r) for r in reactants]

    def __call__(self, structures, limit=1, skip_intersection=True):
        if any(not isinstance(structure, MoleculeContainer) for structure in structures):
            raise TypeError('only list of Molecules possible')
        if self.__single:
            patch = self.__reactor(structures[0], limit, skip_intersection)
            if limit == 1:
                if patch:
                    if self.__split:
                        return ReactionContainer(reactants=structures, products=patch.split())
                    return ReactionContainer(reactants=structures, products=[patch])
                else:
                    return
            if self.__split:
                g = (ReactionContainer(reactants=structures, products=x.split()) for x in patch)
            else:
                g = (ReactionContainer(reactants=structures, products=[x]) for x in patch)
            return list(g) if limit > 1 else g
        else:
            structures = self.__remap(structures)
            mapping = self.__get_mapping(structures)
            structure = reduce(or_, structures)
            if limit == 1:
                mapping = next(mapping, None)
                if mapping:
                    patch = self.__reactor.patcher(structure, mapping)
                    if self.__split:
                        return ReactionContainer(reactants=structures, products=patch.split())
                    return ReactionContainer(reactants=structures, products=[patch])
            else:
                if skip_intersection:
                    mapping = skip(mapping)
                g = (self.__reactor.patcher(structure, m) for m in mapping)

                if self.__split:
                    r = (ReactionContainer(reactants=structures, products=p.split()) for p in g)
                else:
                    r = (ReactionContainer(reactants=structures, products=[p]) for p in g)

                if limit > 1:
                    return list(islice(r, limit))
                else:
                    return r

    def __get_mapping(self, structures):
        """
        match each pattern to each molecule.
        if all patterns matches with all molecules
        return generator of all possible mapping.

        :param structures: disjoint molecules
        :return: mapping generator
        """
        for c in permutations(structures, len(self.__patterns)):
            for m in product(*(x.get_substructure_mapping(y, limit=0) for x, y in zip(self.__patterns, c))):
                mapping = {}
                for i in m:
                    mapping.update(i)
                if mapping:
                    yield mapping

    @staticmethod
    def __remap(structures):
        checked = []
        checked_atoms = set()
        for structure in structures:
            intersection = set(structure.atoms_numbers).intersection(checked_atoms)
            if intersection:
                mapping = {k: v for k, v in zip(intersection,
                                                count(max(max(checked_atoms), max(structure.atoms_numbers)) + 1))}
                structure = structure.remap(mapping, copy=True)
                structure.reset_query_marks()
                info("some atoms in input structures had the same numbers.\n"
                     f"atoms {list(mapping)} were remapped to {list(mapping.values())}")
            checked_atoms.update(structure.atoms_numbers)
            checked.append(structure)
        return checked


def skip(mapping):
    """
    :param mapping: generator
    :return: filtered generator
    """
    found = set()
    for m in mapping:
        matched_atoms = set(m.values())
        if found.intersection(matched_atoms):
            continue
        found.update(matched_atoms)
        yield m


__all__ = ['CGRReactor', 'Reactor']
