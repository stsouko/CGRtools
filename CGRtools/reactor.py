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
from functools import reduce


class CGRreactor:
    def get_template_searcher(self, templates):
        def searcher(g, skip_intersection=True):
            if skip_intersection:
                found = set()

            for i in templates:
                gm = self.get_cgr_matcher(g, i.pattern)
                for j in gm.subgraph_isomorphisms_iter():
                    matched_atoms = set(j)
                    if skip_intersection:
                        if found.intersection(matched_atoms):
                            continue
                        else:
                            found.update(matched_atoms)

                    yield MatchContainer(j, self.__remap_group(i.patch, g, {y: x for x, y in j.items()})[0], i.meta)

        return searcher

    @staticmethod
    def __remap_group(g, h, mapping):
        newmap = mapping.copy()
        newmap.update({x: y for x, y in zip(set(g).difference(newmap), set(range(1, 1000)).difference(h))})
        return g.remap(newmap, copy=True), newmap

    @staticmethod
    def prepare_templates(raw_templates):
        a_marks = {'element', 's_charge', 's_hyb', 's_neighbors', 'p_charge', 'p_hyb', 'p_neighbors',
                   's_stereo', 'p_stereo', 's_radical', 'p_radical'}
        b_marks = {'s_bond', 'p_bond', 's_stereo', 'p_stereo'}
        x_marks = ('s_x', 's_y', 's_z', 'p_x', 'p_y', 'p_z')
        templates = []
        for template in raw_templates:
            if isinstance(template, CGRTemplate):
                templates.append(template)
                continue

            products = reduce(CGRcore.union, template.products).copy()
            reagents = reduce(CGRcore.union, template.reagents).copy()
            if not (isinstance(reagents, QueryContainer) and isinstance(products, QueryContainer)):
                raise InvalidTemplate('templates should be QueryContainer')

            new_atoms = set(products).difference(reagents)
            for n in new_atoms:  # if unique atoms in patch has variable properties exception is raised
                pnn = products.nodes[n]
                for j in a_marks.intersection(pnn):
                    if isinstance(pnn[j], list):
                        raise InvalidTemplate("new atoms can't be variable")

            common = set(products).intersection(reagents)
            for n in common:
                pnn = products.nodes[n]
                rnn = reagents.nodes[n]
                for j in a_marks.intersection(pnn):
                    if isinstance(pnn[j], list):
                        pnn[j] = dict(zip(rnn[j], pnn[j]))
                for j in x_marks:
                    pnn.pop(j)

            for m, n, a in products.edges(data=True):
                if reagents.has_edge(m, n):
                    rmn = reagents[m][n]
                    for j in b_marks.intersection(a):
                        if isinstance(a[j], list):
                            a[j] = dict(zip(rmn[j], a[j]))
                else:
                    for j in b_marks.intersection(a):
                        if isinstance(a[j], list):
                            raise InvalidTemplate("new bonds can't be variable")

            reagents.remap({x: x + 1000 for x in reagents})
            products.remap({x: x + 1000 for x in products})

            templates.append(CGRTemplate(reagents, products, template.meta.copy()))
        return templates

    @classmethod
    def patcher(cls, structure, patch):
        """
        remove edges bw common nodes. add edges from template and replace nodes data

        :param structure: MoleculeContainer or CGRContainer
        :param patch: MoleculeContainer or CGRContainer with replacement data
        """
        if isinstance(structure, CGRContainer):
            node_marks = cls.__cgr_node_marks
            full_node_marks = cls.__cgr_full_node_marks
            bond_marks = cls.__cgr_bond_marks
        elif isinstance(structure, MoleculeContainer):
            node_marks = cls.__node_marks
            full_node_marks = cls.__full_node_marks
            bond_marks = cls.__bond_marks
        else:
            raise ValueError('invalid data type of g: %s' % type(structure))

        p = structure.fresh_copy()
        p.meta.update(structure.meta)
        p.add_nodes_from(structure.nodes(data=True))
        p.add_edges_from((m, n, attr) for m, n, attr in structure.edges(data=True)
                         if m not in patch or n not in patch)

        for i, attr in patch.nodes(data=True):
            if i not in structure:
                p.add_node(i, **{x: y for x, y in attr.items() if x in full_node_marks})
            else:
                p.add_node(i, **{x: y[structure.nodes[i][x]] if isinstance(y, dict) else y
                                 for x, y in attr.items() if x in node_marks})

        for m, n, attr in patch.edges(data=True):
            p.add_edge(m, n, **{x: y[structure[m][n][x]] if isinstance(y, dict) else y
                                for x, y in attr.items() if x in bond_marks})
        return p


__all__ = ['CGRreactor']
