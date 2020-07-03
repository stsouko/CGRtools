# -*- coding: utf-8 -*-
#
#  Copyright 2014-2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from itertools import chain, count, permutations
from logging import info
from operator import or_
from typing import Union, Iterable
from ._functions import lazy_product
from .containers import QueryContainer, QueryCGRContainer, MoleculeContainer, CGRContainer, ReactionContainer
from .periodictable import Element, DynamicElement


class BaseReactor:
    def __init__(self, reactants, products, delete_atoms):
        if isinstance(reactants, QueryContainer):
            self.__is_cgr = is_cgr = False
            e = Element
        else:
            self.__is_cgr = is_cgr = True
            e = DynamicElement

        self.__to_delete = set(reactants).difference(products) if delete_atoms else set()
        self.__elements = elements = {}
        atoms = defaultdict(dict)
        for n, atom in products.atoms():
            if atom.neighbors or atom.hybridization:
                info('neighbors and hybridization for new atoms unusable')
            atoms[n].update(charge=atom.charge, is_radical=atom.is_radical)
            elements[n] = e.from_atomic_number(atom.atomic_number)(atom.isotope)
            if n not in reactants:
                atoms[n]['xy'] = atom.xy
            if is_cgr:
                atoms[n].update(p_is_radical=atom.p_is_radical, p_charge=atom.p_charge)

        self.__atom_attrs = dict(atoms)
        self.__bond_attrs = list(products.bonds())

    def _patcher(self, structure, mapping):
        elements = self.__elements

        plane = structure._plane
        bonds = structure._bonds
        charges = structure._charges
        radicals = structure._radicals
        if self.__is_cgr:
            p_charges = structure._p_charges
            p_radicals = structure._p_radicals

        new = structure.__class__()

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

        max_atom = max(charges) + 1
        for n, atom in self.__atom_attrs.items():
            if n in mapping:  # add matched atoms
                m = mapping[n]
                new.add_atom(elements[n].copy(), m, xy=plane[m], **atom)
            else:  # new atoms
                mapping[n] = new.add_atom(elements[n].copy(), max_atom, **atom)
                max_atom += 1

        old_atoms = set(new._atoms)
        if self.__is_cgr:
            for n, atom in structure.atoms():  # add unmatched atoms
                if n not in old_atoms and n not in to_delete:
                    new.add_atom(atom.copy(), n, charge=charges[n], is_radical=radicals[n], xy=plane[n],
                                 p_is_radical=p_radicals[n], p_charge=p_charges[n])
        else:
            for n, atom in structure.atoms():  # add unmatched atoms
                if n not in old_atoms and n not in to_delete:
                    new.add_atom(atom.copy(), n, charge=charges[n], is_radical=radicals[n], xy=plane[n])

        for n, m, bond in self.__bond_attrs:  # add patch bonds
            n = mapping[n]
            m = mapping[m]
            new.add_bond(n, m, bond.copy())

        for n, m_bond in bonds.items():
            if n in to_delete:  # atoms for removing
                continue
            to_delete.add(n)
            for m, bond in m_bond.items():
                if m in to_delete or n in old_atoms and m in old_atoms:
                    continue
                new.add_bond(n, m, bond.copy())

        # todo: calculate stereo mark based on new atom order
        return new

    def __getstate__(self):
        return {'elements': self.__elements, 'atom_attrs': self.__atom_attrs, 'bond_attrs': self.__bond_attrs,
                'is_cgr': self.__is_cgr, 'to_delete': self.__to_delete}

    def __setstate__(self, state):
        self.__elements = state['elements']
        self.__atom_attrs = state['atom_attrs']
        self.__bond_attrs = state['bond_attrs']
        self.__is_cgr = state['is_cgr']
        self.__to_delete = state['to_delete']


class CGRReactor(BaseReactor):
    """
    CGR based editor for CGRs and molecules.
    generates transformation from input CGR/molecule
    using template (CGRtools ReactionContainer).
    Template should contain one reactant and one product:

    CGRReactor calling transforms reactants to products and
    returns generator of all possible products.
    """
    def __init__(self, template: ReactionContainer, delete_atoms: bool = False):
        """
        :param template: CGRtools ReactionContainer
        :param delete_atoms: if True atoms exists in reactant but
                            not exists in product will be removed
        """
        reactants, products = template.reactants, template.products
        if not reactants or not products:
            raise ValueError('empty template')
        if any(isinstance(x, (CGRContainer, QueryCGRContainer)) for x in chain(template.reactants, template.products)):
            reactants = reduce(or_, reactants, QueryCGRContainer())
            products = reduce(or_, products, QueryCGRContainer())
        else:
            reactants = reduce(or_, reactants, QueryContainer())
            products = reduce(or_, products, QueryContainer())

        self.__pattern = reactants
        self.__meta = template.meta.copy()
        super().__init__(reactants, products, delete_atoms)

    def __call__(self, structure: Union[MoleculeContainer, CGRContainer], automorphism_filter: bool = True):
        if not isinstance(structure, (MoleculeContainer, CGRContainer)):
            raise TypeError('only Molecules and CGRs possible')

        for mapping in self.__pattern.get_mapping(structure, automorphism_filter=automorphism_filter):
            new = self._patcher(structure, mapping)
            new.meta.update(self.__meta)
            yield new

    def __getstate__(self):
        return {'pattern': self.__pattern, 'meta': self.__meta, **super().__getstate__()}

    def __setstate__(self, state):
        self.__pattern = state['pattern']
        self.__meta = state['meta']
        super().__setstate__(state)


class Reactor(BaseReactor):
    """
    CGR based reactor for molecules/queries.
    generates reaction from input queries/molecules using
    transformation template (CGRtools ReactionContainer).

    Reactor calling transforms reactants to products and
    returns generator of reaction transformations with all
    possible products
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
        if not all(isinstance(x, (QueryContainer, MoleculeContainer)) for x in chain(products, reactants)):
            raise TypeError('only Molecules and Queries possible')

        self.__patterns = reactants = tuple(QueryContainer() | x for x in reactants)
        self.__split = len(products)

        products = reduce(or_, products, QueryContainer())
        reactants = reduce(or_, reactants)
        self.__meta = template.meta.copy()
        super().__init__(reactants, products, delete_atoms)

    def __call__(self, structures: Iterable[MoleculeContainer], automorphism_filter: bool = True):
        if any(not isinstance(structure, MoleculeContainer) for structure in structures):
            raise TypeError('only list of Molecules possible')

        structures = self.__remap(structures)
        s_nums = set(range(len(structures)))
        for chosen in permutations(s_nums, len(self.__patterns)):
            ignored = [structures[x] for x in s_nums.difference(chosen)]
            ignored_numbers = {x for x in ignored for x in x}
            chosen = [structures[x] for x in chosen]
            united_chosen = reduce(or_, chosen)
            for match in lazy_product(*(x.get_mapping(y, automorphism_filter=automorphism_filter)
                                        for x, y in zip(self.__patterns, chosen))):
                mapping = match[0]
                for m in match[1:]:
                    mapping.update(m)
                new = self._patcher(united_chosen, mapping)
                collision = set(new).intersection(ignored_numbers)
                if collision:
                    new.remap(dict(zip(collision, count(max(max(ignored_numbers), max(new.atoms_numbers)) + 1))))
                if self.__split > 1:
                    new = new.split()
                    if len(new) != self.__split:
                        info(f'expected {self.__split} molecules in reaction products, but {len(new)} formed.\n'
                             'input molecules has disconnected components')
                else:
                    new = [new]
                yield ReactionContainer(structures, new + ignored, meta=self.__meta)

    @staticmethod
    def __remap(structures):
        checked = []
        checked_atoms = set()
        for structure in structures:
            intersection = set(structure).intersection(checked_atoms)
            if intersection:
                mapping = dict(zip(intersection, count(max(max(checked_atoms), max(structure.atoms_numbers)) + 1)))
                structure = structure.remap(mapping, copy=True)
                info('some atoms in input structures had the same numbers.\n'
                     f'atoms {list(mapping)} were remapped to {list(mapping.values())}')
            checked_atoms.update(structure)
            checked.append(structure)
        return checked

    def __getstate__(self):
        return {'patterns': self.__patterns, 'meta': self.__meta, 'split': self.__split, **super().__getstate__()}

    def __setstate__(self, state):
        self.__patterns = state['patterns']
        self.__meta = state['meta']
        self.__split = state['split']
        super().__setstate__(state)


__all__ = ['CGRReactor', 'Reactor']
