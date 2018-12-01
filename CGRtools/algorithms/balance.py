# -*- coding: utf-8 -*-
#
#  Copyright 2017, 2018 Ramil Nugmanov <stsouko@live.ru>
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
from itertools import cycle, product, combinations
from networkx import compose, has_path
from ..periodictable import radical_unmap


class Balance:
    def _balance(self, ):
        """ calc unbalanced charges and radicals for skin atoms
        """
        meta = h.meta
        for n in (skin_reagent.keys() | skin_product.keys()):
            lost = skin_reagent[n]
            cycle_lost = cycle(lost)
            new = skin_product[n]
            cycle_new = cycle(new)
            atom = h._node[n]
            dr = atom.p_radical - atom.radical
            # radical balancing
            if dr > 0:  # radical added or increased.
                for _, m in zip(range(dr), cycle_lost):  # homolysis
                    s_atom = h._node[m]
                    s_atom.p_multiplicity = radical_unmap[s_atom.p_radical + 1]
                    meta.setdefault('rule #14. atom lost. common atom radical added or increased. '
                                    'lost atom radical added', []).append((m, n))
                for m in lost[dr:]:
                    meta.setdefault('rule #15. atom lost. common atom radical added or increased. '
                                    'lost atom radical unchanged', []).append((m, n))
            elif dr < 0:  # radical removed or decreased.
                if n in skin_product:
                    for m in lost:
                        meta.setdefault('rule #20. atom lost. common atom radical removed or decreased. '
                                        'lost atom radical unchanged', []).append((m, n))
                else:
                    for _, m in zip(range(-dr), cycle_lost):  # radical elimination
                        s_atom = h._node[m]
                        s_atom.p_multiplicity = radical_unmap[s_atom.p_radical + 1]
                        meta.setdefault('rule #21. atom lost. common atom radical removed or decreased. '
                                        'lost atom radical added', []).append((m, n))
                    for m in lost[-dr:]:
                        meta.setdefault('rule #20. atom lost. common atom radical removed or decreased. '
                                        'lost atom radical unchanged', []).append((m, n))
            else:
                env = h.environment(n)
                sv = atom.get_valence([(b.reagent, a.reagent) for b, a in env if b.order])
                pv = atom.p_get_valence([(b.product, a.product) for b, a in env if b.p_order])
                sh, ph = h.atom_total_h(n)

                dv = pv - sv
                dh = ph - sh
                dc = atom.p_charge - atom.charge

                if not (dv or dh or dc):  # common atom unchanged. Substitution, Elimination
                    for m in skins:
                        meta.setdefault('rule #1. atom lost. common atom unchanged. '
                                        'substitution, elimination, addition', []).append((m, n))
                elif dv == dh == dc < 0:  # explicit hydrogen removing
                    for m in skins:
                        h._node[m].p_charge = 1
                        meta.setdefault('rule #4. atom lost. common atom deprotonation', []).append((m, n))
                else:
                    for m in skins:
                        meta.setdefault('rule #5. atom lost. common atom changed. '
                                        'convert to reduction or oxidation', []).append((m, n))

                    pth = ph + sum(h.atom_total_h(x)[1] for x in skins)
                    if n in skin_product:
                        sth = sh + sum(h.atom_total_h(x)[0] for x in skin_product[n])
                    else:
                        sth = sh
                    dth = pth - sth

        for n, skins in skin_product.items():
            cycle_skins = cycle(skins)
            atom = h._node[n]
            dr = atom.p_radical - atom.radical
            # radical balancing
            if dr > 0:  # radical added or increased.
                if n in skin_reagent:
                    for m in skins:
                        meta.setdefault('rule #16. atom new. common atom radical added or increased. '
                                        'new atom radical unchanged', []).append((m, n))
                else:
                    for _, m in zip(range(dr), cycle_skins):  # radical addition
                        s_atom = h._node[m]
                        s_atom.multiplicity = radical_unmap[s_atom.radical + 1]
                        meta.setdefault('rule #17. atom new. common atom radical added or increased. '
                                        'new atom radical added', []).append((m, n))
                    for m in skins[dr:]:
                        meta.setdefault('rule #16. atom new. common atom radical added or increased. '
                                        'new atom radical unchanged', []).append((m, n))
            elif dr < 0:  # radical removed or decreased.
                for _, m in zip(range(-dr), cycle_skins):  # recombination
                    s_atom = h._node[m]
                    s_atom.multiplicity = radical_unmap[s_atom.radical + 1]
                    meta.setdefault('rule #18. atom new. common atom radical removed or decreased. '
                                    'new atom radical added', []).append((m, n))
                for m in skins[-dr:]:
                    meta.setdefault('rule #19. atom new. common atom radical removed or decreased. '
                                    'new atom radical unchanged', []).append((m, n))
            else:
                env = h.environment(n)
                sv = atom.get_valence([(b.reagent, a.reagent) for b, a in env if b.order])
                pv = atom.p_get_valence([(b.product, a.product) for b, a in env if b.p_order])
                sh, ph = h.atom_total_h(n)

                dv = pv - sv
                dh = ph - sh
                dc = atom.p_charge - atom.charge

                if not (dv or dh or dc):  # common atom unchanged. Substitution, Addition
                    for m in skins:
                        meta.setdefault('rule #2. atom new. common atom unchanged. '
                                        'substitution, elimination, addition', []).append((m, n))
                elif dv == dh == dc > 0:  # explicit hydrogen addition
                    for m in skins:
                        h._node[m].charge = 1
                        h.meta.setdefault('rule #3. atom new. common atom protonation', []).append((m, n))
                else:
                    for m in skins:
                        meta.setdefault('rule #6. atom new. common atom changed. '
                                        'convert to reduction or oxidation', []).append((m, n))

                    sth = sh + sum(h.atom_total_h(x)[0] for x in skins)
                    if n in skin_reagent:
                        pth = ph + sum(h.atom_total_h(x)[1] for x in skin_reagent[n])
                    else:
                        pth = ph
                    dth = pth - sth

        for n, sp in reverse_ext.items():

            # charge neutralization
            if dc > 0:
                for _ in range(dc):
                    h.meta.setdefault('rule #7. charge neutralization. hydroxide radical added',
                                      []).append(h.add_atom(O(multiplicity=2), O(charge=-1)))
            elif dc < 0:
                for _ in range(-dc):
                    h.meta.setdefault('rule #8. charge neutralization. hydrogen radical added',
                                      []).append(h.add_atom(H(multiplicity=2), H(charge=1)))

            # hydrogen balancing
            if dth > 0:
                red_e = 0
                for m in sp['products']:
                    if h.nodes[m]['element'] == 'H':  # set reduction H if explicit H count increased
                        h.nodes[m]['s_radical'] = 2
                        red_e += 1
                        h.meta.setdefault('rule #11. protonation. new explicit hydrogen radical added',
                                          []).append(m)

                red = []
                for _ in range(dth - red_e):  # add reduction agents
                    m = h.add_atom(H(multiplicity=2), H())
                    red.append(m)
                    h.meta.setdefault('rule #10. protonation. hydrogen radical added', []).append(m)
                red = iter(red)

                dih = sub(*h.atom_implicit_h(n))
                if dih < 0:  # attach reduction H to central atom if implicit H atoms count increased
                    for _ in range(-dih):
                        m = next(red)
                        h.add_bond(m, n, None)
                        h.meta.setdefault('rule #12. protonation. new implicit hydrogen radical added',
                                          []).append(m)

                for m in sp['reagents']:  # attach reduction H if detached group implicit H count increased
                    dih = sub(*h.atom_implicit_h(m))
                    if dih < 0:
                        for _ in range(-dih):
                            o = next(red)
                            h.add_bond(o, m, None)
            elif dth < 0:
                oxo = []
                for _ in range(-dth):
                    m = h.add_atom(O(multiplicity=2), O())
                    oxo.append(m)
                    h.meta.setdefault('rule #9. deprotonation. hydroxide radical added', []).append(m)
                oxo = iter(oxo)

                for m in sp['reagents']:
                    if h.nodes[m]['element'] == 'H':
                        o = next(oxo)
                        h.add_bond(o, m, None)
                        h.meta.setdefault('rule #13. hydrogen accepting by hydroxide radical added',
                                          []).append(m)

        return h

    def clone_subgraphs(self, g):
        if not isinstance(g, CGRContainer):
            raise InvalidData('only CGRContainer acceptable')

        r_group = []
        x_group = {}
        r_group_clones = []
        newcomponents = []

        ''' search bond breaks and creations
        '''
        components, lost_bonds, term_atoms = self.__split_graph(g)
        lost_map = {x: y for x, y in lost_bonds}
        ''' extract subgraphs and sort by group type (R or X)
        '''
        x_terminals = set(lost_map.values())
        r_terminals = set(lost_map)

        for i in components:
            x_terminal_atom = x_terminals.intersection(i)
            if x_terminal_atom:
                x_group[x_terminal_atom.pop()] = i
                continue

            r_terminal_atom = r_terminals.intersection(i)
            if r_terminal_atom:
                r_group.append([r_terminal_atom, i])
                continue

            newcomponents.append(i)
        ''' search similar R groups and patch.
        '''
        tmp = g
        for i in newcomponents:
            for k, j in r_group:
                gm = GraphMatcher(j, i, node_match=self.__node_match_products,
                                  edge_match=self.__edge_match_products)
                ''' search for similar R-groups started from bond breaks.
                '''
                mapping = next((x for x in gm.subgraph_isomorphisms_iter() if k.issubset(x) and
                                all(x[y] in term_atoms for y in k)), None)
                if mapping:
                    r_group_clones.append([k, mapping])
                    tmp = compose(tmp, self.__remap_group(j, tmp, mapping)[0])
                    break

        ''' add lose X groups to R groups
        '''
        for i, j in r_group_clones:
            for k in i:
                remappedgroup, mapping = self.__remap_group(x_group[lost_map[k]], tmp, {})
                tmp = CGRcore.union(tmp, remappedgroup)
                tmp.add_edge(j[k], mapping[lost_map[k]], s_bond=1, sp_bond=(1, None))

        if r_group_clones:
            tmp.meta.update(g.meta)
            return tmp

        return tmp.copy()

    @classmethod
    def __split_graph(cls, g):
        g = g.copy()
        lost_bonds = []
        term_atoms = []

        for l, n, m in list(cls.__get_substitution_paths(g)):
            nl = (n, l)
            if nl in lost_bonds:
                continue

            lost_bonds.append(nl)
            g.remove_edge(n, l)
            g.remove_edge(n, m)

        for n, m in list(cls.__get_broken_paths(g)):
            if not any(has_path(g, *x) for x in product((y for x in lost_bonds for y in x), (n, m))):
                g.remove_edge(n, m)
                term_atoms.append(n)
                term_atoms.append(m)

        return CGRcore.split(g), lost_bonds, term_atoms

    @staticmethod
    def __get_substitution_paths(g):
        """
        get atoms paths from detached atom to attached

        :param g: CGRContainer
        :return: tuple of atoms numbers
        """
        for n, nbrdict in g.adjacency():
            for m, l in combinations(nbrdict, 2):
                nms = nbrdict[m]['sp_bond']
                nls = nbrdict[l]['sp_bond']
                if nms == (1, None) and nls == (None, 1):
                    yield m, n, l
                elif nms == (None, 1) and nls == (1, None):
                    yield l, n, m

    @staticmethod
    def __get_broken_paths(g):
        for m, n, attr in g.edges(data=True):
            if attr['sp_bond'] == (None, 1):
                yield n, m
