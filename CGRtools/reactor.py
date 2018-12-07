# -*- coding: utf-8 -*-
#
#  Copyright 2014-2018 Ramil Nugmanov <stsouko@live.ru>
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
from itertools import repeat
from logging import warning
from operator import or_
from .containers import QueryContainer, QueryCGRContainer, MoleculeContainer, CGRContainer


class CGRreactor:
    def __init__(self, template):
        pattern, absolute_atom, absolute_bond, conditional_atom, conditional_bond, is_cgr = \
            self.__prepare_template(template)
        self.__pattern = pattern
        self.__absolute_atom = absolute_atom
        self.__absolute_bond = absolute_bond
        self.__conditional_atom = conditional_atom
        self.__conditional_bond = conditional_bond
        self.__is_cgr = is_cgr
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

        absolute_atom = defaultdict(dict)
        conditional_atom = defaultdict(dict)

        for n in set(products).difference(reagents):
            # if unique atoms in patch has variable properties exception is raised
            atom = products._node[n]
            if len(atom.element) > 1 or atom.element == ('A',) or any(len(atom[x]) > 1 for x in
                                                                      ('charge', 'multiplicity', 'isotope')):
                raise ValueError('new atoms in patch should be static')
            elif atom.neighbors or atom.hybridization:
                warning('neighbors and hybridization for new atoms unusable')
            absolute_atom[n].update(element=atom.element[0], charge=atom.charge[0],
                                    isotope=atom.isotope and atom.isotope[0] or None,
                                    multiplicity=atom.multiplicity and atom.multiplicity[0] or None,
                                    x=atom.x, y=atom.y, z=atom.z, stereo=atom.stereo)
            if is_cgr:
                absolute_atom[n].update(p_charge=atom.p_charge[0],
                                        p_multiplicity=atom.p_multiplicity and atom.p_multiplicity[0] or None,
                                        p_x=atom.p_x, p_y=atom.p_y, p_z=atom.p_z, p_stereo=atom.p_stereo)

        for n in set(products).intersection(reagents):
            r_atom = reagents._node[n]
            atom = products._node[n]

            absolute_atom[n]['stereo'] = atom.stereo
            if is_cgr:
                absolute_atom[n]['p_stereo'] = atom.p_stereo

            if atom.element != ('A',):
                if len(r_atom.element) == 1:
                    absolute_atom[n]['element'] = atom.element[0]
                elif len(r_atom.element) != len(atom.element):
                    raise ValueError('element mapping invalid. possible A>A, E>A, A>E, [E]>A, [E]>[E]')
                else:
                    conditional_atom[n]['element'] = dict(zip(r_atom.element, atom.element))
            if len(atom.isotope) == 1:
                absolute_atom[n]['isotope'] = atom.isotope[0]
            elif len(r_atom.isotope) != len(atom.isotope) != 0:
                raise ValueError(f'isotope mapping invalid. possible n>1, n>n, n>0')
            elif atom.isotope:  # mapped replace
                conditional_atom[n]['isotope'] = dict(zip(r_atom.isotope, atom.isotope))

            if len(atom.charge) == 1:  # replace for fixed value
                absolute_atom[n]['charge'] = atom.charge[0]
                if is_cgr:
                    absolute_atom[n]['p_charge'] = atom.p_charge[0]
            elif len(r_atom.charge) != len(atom.charge):
                raise ValueError(f'charge mapping invalid. possible n>1, n>n')
            else:
                if is_cgr:
                    conditional_atom[n]['charge'] = dict(zip(zip(r_atom.charge, r_atom.p_charge),
                                                        zip(atom.charge, atom.p_charge)))
                else:
                    conditional_atom[n]['charge'] = dict(zip(r_atom.charge, atom.charge))

            if len(atom.multiplicity) == 1:  # replace for fixed value
                absolute_atom[n]['multiplicity'] = atom.multiplicity[0]
                if is_cgr:
                    absolute_atom[n]['p_multiplicity'] = atom.p_multiplicity[0]
            elif len(r_atom.multiplicity) != len(atom.multiplicity) != 0:
                raise ValueError(f'multiplicity mapping invalid. possible n>1, n>n, n>0')
            elif atom.multiplicity:  # mapped replace
                if is_cgr:
                    conditional_atom[n]['multiplicity'] = dict(zip(zip(r_atom.multiplicity, r_atom.p_multiplicity),
                                                              zip(atom.multiplicity, atom.p_multiplicity)))
                else:
                    conditional_atom[n]['multiplicity'] = dict(zip(r_atom.multiplicity, atom.multiplicity))
            # save found state if empty atom[multiplicity or charge]

        absolute_bond = defaultdict(dict)
        conditional_bond = defaultdict(dict)

        seen = set()
        for n, m_bond in products._adj.items():
            seen.add(n)
            for m, bond in m_bond.items():
                if m in seen:
                    continue
                elif n not in reagents or m not in reagents or n not in products._adj[m]:
                    if is_cgr:
                        if bond.order and len(bond.order) == 1 or bond.p_order and len(bond.p_order) == 1:
                            # not sequential are more common
                            absolute_bond[n][m] = absolute_bond[m][n] = \
                                {'order': bond.order and bond.order[0], 'p_order': bond.p_order and bond.p_order[0],
                                 'stereo': bond.stereo, 'p_stereo': bond.p_stereo}
                        else:
                            raise ValueError('new bonds in patch should be static')
                    elif len(bond.order) == 1:
                        absolute_bond[n][m] = absolute_bond[m][n] = {'order': bond.order[0], 'stereo': bond.stereo}
                    else:
                        raise ValueError('new bonds in patch should be static')
                else:
                    r_bond = reagents._adj[n][m]
                    if is_cgr:
                        if bond.order and len(bond.order) == 1 or bond.p_order and len(bond.p_order) == 1:
                            absolute_bond[n][m] = absolute_bond[m][n] = \
                                {'order': bond.order and bond.order[0], 'p_order': bond.p_order and bond.p_order[0],
                                 'stereo': bond.stereo, 'p_stereo': bond.p_stereo}
                        elif not bond.order:
                            if not r_bond.p_order:
                                raise ValueError('None bond not mappable')
                            elif len(bond.p_order) != len(r_bond.p_order):
                                raise ValueError('bond order mapping invalid. possible n>n, n>1')
                            conditional_bond[n][m] = conditional_bond[m][n] = \
                                {'order': dict(zip(zip(repeat(None), r_bond.p_order),
                                                   zip(repeat(None), bond.p_order))),
                                 'stereo': bond.stereo, 'p_stereo': bond.p_stereo}
                        elif not bond.p_order:
                            if not r_bond.order:
                                raise ValueError('None bond not mappable')
                            elif len(bond.order) != len(r_bond.order):
                                raise ValueError('bond order mapping invalid. possible n>n, n>1')
                            conditional_bond[n][m] = conditional_bond[m][n] = \
                                {'order': dict(zip(zip(r_bond.order, repeat(None)),
                                                   zip(bond.order, repeat(None)))),
                                 'stereo': bond.stereo, 'p_stereo': bond.p_stereo}
                        elif len(bond.order) != len(r_bond.order):
                            raise ValueError('bond order mapping invalid. possible n>n, n>1')
                        else:
                            conditional_bond[n][m] = conditional_bond[m][n] = \
                                {'order': dict(zip(zip(r_bond.order, r_bond.p_order),
                                                   zip(bond.order, bond.p_order))),
                                 'stereo': bond.stereo, 'p_stereo': bond.p_stereo}
                    elif len(bond.order) == 1:
                        absolute_bond[n][m] = absolute_bond[m][n] = {'order': bond.order[0], 'stereo': bond.stereo}
                    elif len(bond.order) != len(r_bond.order):
                        raise ValueError(f'bond order mapping invalid. possible n>n, n>1')
                    else:
                        conditional_bond[n][m] = conditional_bond[m][n] = \
                            {'order': dict(zip(r_bond.order, bond.order)), 'stereo': bond.stereo}

        return reagents, dict(absolute_atom), dict(absolute_bond), dict(conditional_atom), dict(conditional_bond), \
            is_cgr

    def __patcher(self, structure, mapping):
        new = type(structure)()
        new.meta.update(self.__meta)

        new_atoms = []
        for n, atom in self.__absolute_atom.items():
            if n in mapping:
                n = mapping[n]
                new.add_atom({**structure._node[n], **atom}, n)
            else:
                new_atoms.append(n)
        for n, atom in self.__conditional_atom:
            n = mapping[n]
            attr = {}
            r_atom = structure._node[n]
            for k, replace in atom.items():
                if k in ('element', 'isotope'):
                    attr[k] = replace[r_atom[k]]
                elif self.__is_cgr:
                    p_k = f'p_{k}'
                    attr[k], attr[p_k] = replace[(r_atom[k], r_atom[p_k])]
                else:
                    attr[k] = replace[r_atom[k]]
            new.add_atom(attr, n)
        for n in structure._node.keys() - new._node.keys():  # add unmatched atoms
            new.add_atom(structure._node[n], n)
        for n in new_atoms:
            mapping[n] = new.add_atom(self.__absolute_atom[n])

        seen_a = set()
        for n, m_bond in self.__absolute_bond.items():
            n = mapping[n]
            seen_a.add(n)
            for m, bond in m_bond.items():
                m = mapping[m]
                if m in seen_a:
                    continue
                new.add_bond(n, m, bond)
        seen_c = set()
        for n, m_bond in self.__conditional_bond.items():
            n = mapping[n]
            seen_c.add(n)
            for m, bond in m_bond.items():
                m = mapping[m]
                if m in seen_c:
                    continue
                r_bond = structure._adj[n][m]
                if self.__is_cgr:
                    order, p_order = bond['order'][(r_bond['order'], r_bond['p_order'])]
                    new.add_bond(n, m, {'order': order, 'p_order': p_order, 'stereo': bond['stereo']})
                else:
                    new.add_bond(n, m, {'order': bond['order'][r_bond['order']],
                                        'stereo': bond['stereo']})
        seen_matched = seen_a | seen_c
        seen = set()
        for n, m_bond in structure._adj.items():
            seen.add(n)
            for m, bond in m_bond.items():
                if m in seen or n in seen_matched and m in seen_matched:
                    continue
                new.add_bond(n, m, bond)

        return new


__all__ = ['CGRreactor']
