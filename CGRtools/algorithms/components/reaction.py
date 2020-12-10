# -*- coding: utf-8 -*-
#
#  Copyright 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2020 Ravil Mukhametgaleev <sonic-mc@mail.ru>
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
from CachedMethods import cached_property
from collections import ChainMap
from itertools import chain, product
from typing import Tuple, Iterator, TYPE_CHECKING
from ...containers import molecule  # cyclic imports resolve


if TYPE_CHECKING:
    from CGRtools import ReactionContainer


class ReactionComponents:
    __slots__ = ()

    @cached_property
    def centers_list(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Reaction centers with leaving and coming groups.
        """
        if not self.reactants or not self.products:
            return ()  # no rc
        elif not isinstance(self.reactants[0], molecule.MoleculeContainer):
            raise TypeError('Only Molecules supported')

        cgr = self.compose()
        bonds = cgr._bonds
        all_groups = {x for x in self.reactants for x in x} ^ {x for x in self.products for x in x}
        all_groups = cgr._connected_components({n: bonds[n].keys() & all_groups for n in all_groups})
        centers_list = list(cgr.centers_list)

        for x in all_groups:
            intersection = []
            for i, y in enumerate(centers_list):
                if not x.isdisjoint(y):
                    intersection.append(i)
            if intersection:
                for i in reversed(intersection):
                    x.update(centers_list.pop(i))
                centers_list.append(x)
        return tuple(tuple(x) for x in centers_list)

    @cached_property
    def extended_centers_list(self) -> Tuple[Tuple[int, ...], ...]:
        """
        Additionally to `centers_list` include:
        * First environment of dynamic atoms.
        * Whole formed cycles. For condensed cycles smallest is taken.
        * Whole aromatic cycle with at least one dynamic atom.
        * Whole small (3, 4) cycle with at least one dynamic atom.
        * Double or triple bonds connected to previous atoms.

        Note for multiple RCs intersection possible. Use `enumerate_centers` to prevent unobvious RCs.
        """
        cgr = self.compose()
        bonds = cgr._bonds
        center_atoms = set(cgr.center_atoms)

        formed_rings = {}
        small_aromatic_rings = set()
        for r in cgr.sssr:
            if len(r) < 5 and not center_atoms.isdisjoint(r):
                small_aromatic_rings.add(frozenset(r))

            n, m = r[0], r[-1]
            if bonds[n][m].order is None:
                nm = frozenset((n, m))
                if nm not in formed_rings or len(formed_rings[nm]) > len(r):
                    formed_rings[nm] = r
            for n, m in zip(r, r[1:]):
                if bonds[n][m].order is None:
                    nm = frozenset((n, m))
                    if nm not in formed_rings or len(formed_rings[nm]) > len(r):
                        formed_rings[nm] = r

        for r in cgr.aromatic_rings:
            if not center_atoms.isdisjoint(r):
                small_aromatic_rings.add(frozenset(r))

        out = []
        for rc in self.centers_list:
            c = center_atoms.intersection(rc)
            fe = {m for n in c for m in bonds[n]}  # add first environment
            # add double or triple bond to first env
            fe |= {m for n in fe for m, b in bonds[n].items() if m not in fe and b.order in (2, 3)}
            fe.update(rc)  # add leaving|coming groups and alone center atoms

            for rk, r in formed_rings.items():  # add formed rings to RC
                if c.issuperset(rk):
                    fe.update(r)
            for r in small_aromatic_rings:  # add small or aromatic rings with dyn atoms
                if not c.isdisjoint(r):
                    fe.update(r)
            out.append(tuple(fe))
        return tuple(out)

    def enumerate_centers(self) -> Iterator['ReactionContainer']:
        """
        Get all possible single stage reactions from multistage.
        Note multicomponent molecules (salts etc) can be treated incorrectly.
        """
        if len(self.centers_list) > 1:
            centers_list = self.centers_list

            charges = ChainMap(*(x._charges for x in self.reactants))
            radicals = ChainMap(*(x._radicals for x in self.reactants))
            bonds = ChainMap(*(x._bonds for x in self.reactants))
            atoms = ChainMap(*(x._atoms for x in self.reactants))
            p_charges = ChainMap(*(x._charges for x in self.products))
            p_radicals = ChainMap(*(x._radicals for x in self.products))
            p_bonds = ChainMap(*(x._bonds for x in self.products))
            p_atoms = ChainMap(*(x._atoms for x in self.products))

            centers = {x for x in centers_list for x in x}
            common = {x for x in chain(self.reactants, self.products) for x in x if x not in centers}
            reactants = {x for x in self.reactants for x in x}
            products = {x for x in self.products for x in x}

            common_molecule = molecule.MoleculeContainer()
            for n in common:
                common_molecule.add_atom(atoms[n].copy(), n, charge=charges[n], is_radical=radicals[n])
            seen = set()
            for n in common:
                seen.add(n)
                for m, b in bonds[n].items():
                    if m not in seen and m in common:
                        common_molecule.add_bond(n, m, b.copy())

            products_bonds = {}
            reactants_bonds = {}
            common_bonds = []
            seen = set()
            p_seen = set()
            for c in centers_list:
                not_rc = centers.difference(c)
                reactants_bonds[c] = (c_bonds, c_atoms) = [], reactants.intersection(c)
                for n in c_atoms:
                    seen.add(n)
                    for m, b in bonds[n].items():
                        if m not in seen and m in reactants:
                            if m in not_rc:
                                common_bonds.append((n, m, b))
                            else:
                                c_bonds.append((n, m, b))
                products_bonds[c] = (c_bonds, c_atoms) = [], products.intersection(c)
                for n in c_atoms:
                    p_seen.add(n)
                    for m, b in p_bonds[n].items():
                        if m not in p_seen and m in products and m not in not_rc:
                            c_bonds.append((n, m, b))

            for rc in range(len(centers_list)):
                not_rc = centers_list[:rc] + centers_list[rc + 1:]
                rc = centers_list[rc]
                for combo in list(product((False, True), repeat=len(not_rc))):
                    r = common_molecule.copy()
                    p = common_molecule.copy()

                    for n in reactants_bonds[rc][1]:
                        r.add_atom(atoms[n].copy(), n, charge=charges[n], is_radical=radicals[n])
                    for n in products_bonds[rc][1]:
                        p.add_atom(p_atoms[n].copy(), n, charge=p_charges[n], is_radical=p_radicals[n])

                    for is_p, center in zip(combo, not_rc):
                        if is_p:
                            c_bonds, c_atoms = products_bonds[center]
                            for n in c_atoms:
                                r.add_atom(p_atoms[n].copy(), n, charge=p_charges[n], is_radical=p_radicals[n])
                                p.add_atom(p_atoms[n].copy(), n, charge=p_charges[n], is_radical=p_radicals[n])
                        else:
                            c_bonds, c_atoms = reactants_bonds[center]
                            for n in c_atoms:
                                r.add_atom(atoms[n].copy(), n, charge=charges[n], is_radical=radicals[n])
                                p.add_atom(atoms[n].copy(), n, charge=charges[n], is_radical=radicals[n])
                        for n, m, b in c_bonds:
                            r.add_bond(n, m, b.copy())
                            p.add_bond(n, m, b.copy())

                    for n, m, b in products_bonds[rc][0]:
                        p.add_bond(n, m, b.copy())
                    for n, m, b in reactants_bonds[rc][0]:
                        r.add_bond(n, m, b.copy())
                    for n, m, b in common_bonds:
                        r.add_bond(n, m, b.copy())
                        p.add_bond(n, m, b.copy())
                    yield self.__class__(r.split(), p.split(), [x.copy() for x in self.reagents])
        else:
            cp = self.copy()
            cp.meta.clear()
            yield cp


__all__ = ['ReactionComponents']
