# -*- coding: utf-8 -*-
#
#  Copyright 2017-2019 Ramil Nugmanov <stsouko@live.ru>
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
from collections.abc import MutableSequence
from functools import reduce
from operator import or_
from .cgr import CGRContainer
from .molecule import MoleculeContainer
from .query import QueryCGRContainer
from ..algorithms import HashableSmiles, DepictReaction
from ..cache import cached_method


class MindfulList(MutableSequence):
    """list with self-checks of modification. need for control of ReactionContainer caches actuality"""
    def __init__(self, data=None):
        self.__data = [] if data is None else list(data)
        self.__check = True

    def get_state(self):
        """return True if structure data list changed from previous checking time"""
        tmp = self.__check
        self.__check = False
        return tmp

    def insert(self, index, obj):
        self.__check = True
        self.__data.insert(index, obj)

    def __delitem__(self, index):
        self.__check = True
        del self.__data[index]

    def __add__(self, other):
        return self.__class__(self.__data + list(other))

    def __getitem__(self, index):
        return self.__data[index]

    def __len__(self):
        return len(self.__data)

    def __setitem__(self, key, value):
        self.__check = True
        self.__data[key] = value

    def __str__(self):
        return '[%s]' % ', '.join(str(x) for x in self.__data)


class ReactionContainer(DepictReaction, HashableSmiles):
    """
    reaction storage. contains reactants, products and reagents lists.

    reaction storages hashable and comparable. based on reaction unique signature (SMIRKS).
    for reactions with query containers hash and comparison may give errors due to non-uniqueness.
    query containers itself not support hashing and comparison.
    """
    __slots__ = ('__reactants', '__products', '__reagents', '__meta', '_arrow')

    def __init__(self, reactants=None, products=None, reagents=None, meta=None):
        """
        new empty or filled reaction object creation

        :param reactants: list of MoleculeContainers [or other Structure Containers] in left side of reaction
        :param products: right side of reaction. see reactants
        :param reagents: middle side of reaction: solvents, catalysts, etc. see reactants
        :param meta: dictionary of metadata. like DTYPE-DATUM in RDF

        """
        self.__reactants = MindfulList(reactants)
        self.__products = MindfulList(products)
        self.__reagents = MindfulList(reagents)
        if meta is None:
            self.__meta = {}
        else:
            self.__meta = dict(meta)
        self._arrow = None

    def __getitem__(self, item):
        if item in ('reactants', 0):
            return self.__reactants
        elif item in ('products', 1):
            return self.__products
        elif item in ('reagents', 2):
            return self.__reagents
        elif item == 'meta':
            return self.__meta
        raise AttributeError('invalid attribute')

    def __getstate__(self):
        return dict(reactants=list(self.__reactants), products=list(self.__products), reagents=list(self.__reagents),
                    meta=self.meta)

    def __setstate__(self, state):
        if next(iter(state)) == 'reagents':  # 3.0 compatibility
            state['reagents'], state['reactants'] = state['reactants'], state['reagents']
        self.__init__(**state)

    @property
    def reactants(self):
        """reactants list. see products"""
        return self.__reactants

    @property
    def reagents(self):
        """reagents list. see products"""
        return self.__reagents

    @property
    def products(self):
        """list of CGRs or/and Molecules in products side"""
        return self.__products

    @property
    def meta(self):
        """dictionary of metadata. like DTYPE-DATUM in RDF"""
        return self.__meta

    def copy(self):
        """
        get copy of object

        :return: ReactionContainer
        """
        return type(self)(reagents=[x.copy() for x in self.__reagents], meta=self.__meta.copy(),
                          products=[x.copy() for x in self.__products],
                          reactants=[x.copy() for x in self.__reactants])

    def implicify_hydrogens(self):
        """
        remove explicit hydrogens if possible

        :return: number of removed hydrogens
        """
        total = 0
        for ml in (self.__reagents, self.__reactants, self.__products):
            for m in ml:
                if hasattr(m, 'implicify_hydrogens'):
                    total += m.implicify_hydrogens()
        if total:
            self.flush_cache()
        return total

    def explicify_hydrogens(self):
        """
        add explicit hydrogens to atoms

        :return: number of added atoms
        """
        total = 0
        for ml in (self.__reagents, self.__reactants, self.__products):
            for m in ml:
                if hasattr(m, 'explicify_hydrogens'):
                    total += m.explicify_hydrogens()
        if total:
            self.flush_cache()
        return total

    def aromatize(self):
        """
        convert structures to aromatic form. works only for Molecules

        :return: number of processed molecules
        """
        total = 0
        for ml in (self.__reagents, self.__reactants, self.__products):
            for m in ml:
                if hasattr(m, 'aromatize'):
                    if m.aromatize():
                        total += 1
        if total:
            self.flush_cache()
        return total

    def standardize(self):
        """
        standardize functional groups and convert structures to aromatic form. works only for Molecules

        :return: number of processed molecules
        """
        total = 0
        for ml in (self.__reagents, self.__reactants, self.__products):
            for m in ml:
                if hasattr(m, 'standardize'):
                    if m.standardize():
                        total += 1
        if total:
            self.flush_cache()
        return total

    def reset_query_marks(self):
        """
        set or reset hyb and neighbors marks to atoms.
        """
        for ml in (self.__reagents, self.__reactants, self.__products):
            for m in ml:
                if hasattr(m, 'reset_query_marks'):
                    m.reset_query_marks()
        self.flush_cache()

    @cached_method
    def compose(self):
        """
        get CGR of reaction

        reagents will be presented as unchanged molecules
        :return: CGRContainer
        """
        rr = self.__reagents + self.__reactants
        if rr:
            if not all(isinstance(x, (MoleculeContainer, CGRContainer)) for x in rr):
                raise TypeError('Queries not composable')
            r = reduce(or_, rr)
        else:
            r = MoleculeContainer()
        if self.__products:
            if not all(isinstance(x, (MoleculeContainer, CGRContainer)) for x in self.__products):
                raise TypeError('Queries not composable')
            p = reduce(or_, self.__products)
        else:
            p = MoleculeContainer()
        return r ^ p

    def __invert__(self):
        """
        get CGR of reaction
        """
        return self.compose()

    def calculate2d(self, force=True):
        """
        recalculate 2d coordinates. currently rings can be calculated badly.

        :param force: ignore existing coordinates of atoms
        """
        for ml in (self.__reagents, self.__reactants, self.__products):
            for m in ml:
                m.calculate2d(force)
        self.fix_positions()

    def fix_positions(self):
        """
        fix coordinates of molecules in reaction
        """
        shift_x = 0
        for m in self.__reactants:
            max_x = self.__fix_positions(m, shift_x, 0)
            shift_x = max_x + 1
        arrow_min = shift_x

        if self.__reagents:
            for m in self.__reagents:
                max_x = self.__fix_positions(m, shift_x, 1.5)
                shift_x = max_x + 1
        else:
            shift_x += 3
        arrow_max = shift_x - 1

        for m in self.__products:
            max_x = self.__fix_positions(m, shift_x, 0)
            shift_x = max_x + 1
        self._arrow = (arrow_min, arrow_max)
        self.flush_cache()

    @staticmethod
    def __fix_positions(molecule, shift_x, shift_y):
        min_x = min(atom.x for atom in molecule._node.values()) - shift_x
        max_x = max(atom.x for atom in molecule._node.values()) - min_x
        min_y = min(atom.y for atom in molecule._node.values()) - shift_y
        for atom in molecule._node.values():
            atom.x = atom.x - min_x
            atom.y = atom.y - min_y
        return max_x

    @cached_method
    def __str__(self):
        """
        SMIRKS of reaction. query containers in reaction {surrounded by curly braces}
        """
        sig = []
        for ml in (self.__reactants, self.__reagents, self.__products):
            sig.append(self.__get_smiles(ml) if ml else '')
        return '>'.join(sig)

    @staticmethod
    def __get_smiles(molecules):
        smiles = []
        union = MoleculeContainer()  # need for whole atoms ordering
        queries = []
        atoms = {}
        for m in molecules:
            if isinstance(m, (MoleculeContainer, CGRContainer)):
                union._node.update(m._node)
                union._adj.update(m._adj)
                atoms.update(dict.fromkeys(m, m))
            elif isinstance(m, QueryCGRContainer):  # queries added as is without ordering
                queries.append('{%s}' % m)
            else:
                queries.append(str(m))

        order_atoms = set(atoms)
        order = union.atoms_order  # whole atoms order
        while order_atoms:
            next_molecule = atoms[min(order_atoms, key=order.__getitem__)]  # get molecule with smallest atom
            order_atoms.difference_update(next_molecule)
            smiles.append('{%s}' % next_molecule if isinstance(next_molecule, CGRContainer) else str(next_molecule))
        smiles.extend(queries)  # queries always in the end of list
        return '.'.join(smiles)

    def flush_cache(self):
        self.__dict__.clear()


__all__ = ['ReactionContainer']
