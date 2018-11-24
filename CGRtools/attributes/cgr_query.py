# -*- coding: utf-8 -*-
#
#  Copyright 2018 Ramil Nugmanov <stsouko@live.ru>
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
from .cgr import DynAtom, DynBond, DynamicContainer
from .molecule import Atom, Bond
from .query import QueryAtom, QueryBond
from ..periodictable import elements_numbers


class DynQueryAtom(DynAtom):
    def __setattr__(self, key, value):
        if key in self._p_static:
            raise AttributeError(f'{key} is invalid')
        elif key.startswith('p_'):
            key = key[2:]
            value = getattr(self._atom_factory, f'_{key}_check')(value)
            if key != 'stereo' and len(value) != len(getattr(self._reagent, key)):
                raise ValueError(f'{key} lists in reagent and product should be equal')
            setattr(self._product, key, value)
        else:
            value = getattr(self._atom_factory, f'_{key}_check')(value)
            if key in self._static:
                setattr(self._product, key, value)
            elif key != 'stereo' and len(value) != len(getattr(self._product, key)):
                raise ValueError(f'{key} lists in reagent and product should be equal')
            setattr(self._reagent, key, value)

    def _update(self, value, kwargs):
        if isinstance(value, DynQueryAtom):
            r, p = self._split_check_kwargs(kwargs)
            p_value = value._product
            value = value._reagent
        elif isinstance(value, QueryAtom):
            r, p = self._split_check_kwargs(kwargs)
            p_value = value
        else:
            return super()._update(value, kwargs)

        self._reagent._update(value, r)
        self._product._update(p_value, p)

    def stringify(self, atom=True, isotope=True, stereo=True, hybridization=True, neighbors=True):
        rmi, pmi = [], []
        if stereo:
            if self.stereo:
                rmi.append(self._stereo_str[self.stereo])
            if self.p_stereo:
                pmi.append(self._stereo_str[self.p_stereo])
        if hybridization:
            if len(self.hybridization) > 1:
                r, p = zip(*sorted(zip(self.hybridization, self.p_hybridization)))
                rmi.append('<%s>' % ''.join(self._hybridization_str[x] for x in r))
                pmi.append('<%s>' % ''.join(self._hybridization_str[x] for x in p))
            else:
                if self.hybridization:
                    rmi.append(self._hybridization_str[self.hybridization[0]])
                if self.p_hybridization:
                    pmi.append(self._hybridization_str[self.p_hybridization[0]])
        if neighbors:
            if len(self.neighbors) > 1:
                r, p = zip(*sorted(zip(self.neighbors, self.p_neighbors)))
                rmi.append('<%s>' % ''.join(str(x) for x in r))
                pmi.append('<%s>' % ''.join(str(x) for x in p))
            else:
                if self.neighbors:
                    rmi.append(str(self.neighbors[0]))
                if self.p_neighbors:
                    pmi.append(str(self.p_neighbors[0]))
        if rmi:
            rmi.append(';')
            rmi.insert(0, ';')
        if pmi:
            pmi.append(';')
            pmi.insert(0, ';')

        if atom:
            if self.element == ('A',):
                atom = False
                rmi.insert(0, '*')
                pmi.insert(0, '*')
            elif len(self.element) > 1:
                atom = False
                tmp = ','.join(sorted(self.element, key=elements_numbers.get))
                rmi.insert(0, tmp)
                pmi.insert(0, tmp)
            else:
                if self.element[0] not in ('C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B'):
                    atom = False
                rmi.insert(0, self.element[0])
                pmi.insert(0, self.element[0])

            if len(self.charge) > 1:
                r, p = zip(*sorted(zip(self.charge, self.p_charge)))
                rmi.append('<%s>' % ''.join(self._charge_str[x] for x in r))
                pmi.append('<%s>' % ''.join(self._charge_str[x] for x in p))
            else:
                if self.charge != (0,):
                    rmi.append(self._charge_str[self.charge[0]])
                if self.p_charge != (0,):
                    pmi.append(self._charge_str[self.p_charge[0]])
            if len(self.multiplicity) > 1:
                r, p = zip(*sorted(zip(self.multiplicity, self.p_multiplicity)))
                rmi.append('<%s>' % ''.join(self._multiplicity_str[x] for x in r))
                pmi.append('<%s>' % ''.join(self._multiplicity_str[x] for x in p))
            else:
                if self.multiplicity:
                    rmi.append(self._multiplicity_str[self.multiplicity[0]])
                if self.p_multiplicity:
                    pmi.append(self._multiplicity_str[self.p_multiplicity[0]])
            if isotope:
                if len(self.isotope) > 1:
                    tmp = '<%s>' % ''.join(str(x) for x in sorted(self.isotope))
                    rmi.insert(0, tmp)
                    pmi.insert(0, tmp)
                elif self.isotope:
                    tmp = str(self.isotope[0])
                    rmi.insert(0, tmp)
                    pmi.insert(0, tmp)
        else:
            rmi.insert(0, '*')
            pmi.insert(0, '*')

        if len(rmi) != 1 or not atom:
            rmi.insert(0, '[')
            rmi.append(']')
        if len(pmi) != 1 or not atom:
            pmi.insert(0, '[')
            pmi.append(']')
        return ''.join(rmi), ''.join(pmi)

    def weight(self, atom=True, isotope=False, stereo=False, hybridization=False, neighbors=False):
        r_weight, p_weight = [], []
        if atom:
            tmp = tuple(sorted(elements_numbers[x] for x in self.element))
            r_weight.append(tmp)
            p_weight.append(tmp)
            if isotope:
                tmp = tuple(sorted(self.isotope))
                r_weight.append(tmp)
                p_weight.append(tmp)
            r, p = zip(*sorted(zip(self.charge, self.p_charge)))
            r_weight.append(tuple(r))
            p_weight.append(tuple(p))
            if self.multiplicity:
                r, p = zip(*sorted(zip(self.multiplicity, self.p_multiplicity)))
                r_weight.append(tuple(r))
                p_weight.append(tuple(p))
            else:
                r_weight.append(())
                p_weight.append(())
        if stereo:
            r_weight.append(self.stereo or 0)
            p_weight.append(self.p_stereo or 0)
        if hybridization:
            if self.hybridization:
                r, p = zip(*sorted(zip(self.hybridization, self.p_hybridization)))
                r_weight.append(tuple(r))
                p_weight.append(tuple(p))
            else:
                r_weight.append(())
                p_weight.append(())
        if neighbors:
            if self.neighbors:
                r, p = zip(*sorted(zip(self.neighbors, self.p_neighbors)))
                r_weight.append(tuple(r))
                p_weight.append(tuple(p))
            else:
                r_weight.append(())
                p_weight.append(())
        return tuple(r_weight), tuple(p_weight)

    def __eq__(self, other):
        if isinstance(other, DynQueryAtom):
            return all(sorted(self._reagent[x]) == sorted(other._reagent[x]) for x in self._static) and \
                       all(sorted(zip(self._reagent[x], self._product[x])) ==
                           sorted(zip(other._reagent[x], other._product[x]))
                           for x in ('charge', 'multiplicity', 'hybridization', 'neighbors'))
        elif isinstance(other, DynAtom):
            return all(other._reagent[x] in self._reagent[x] for x in self._static) and \
                all((other._reagent[x], other._product[x]) in zip(self._reagent[x], self._product[x])
                    for x in ('charge', 'multiplicity', 'hybridization', 'neighbors'))
        elif isinstance(other, (QueryAtom, Atom)):
            return self._reagent == self._product == other
        return False

    def __ne__(self, other):
        """
        != equality checks with stereo
        """
        if self == other:
            if self.stereo is self.p_stereo is None:
                return True
            elif isinstance(other, DynAtom):
                return self.stereo == other.stereo and self.p_stereo == other.p_stereo
            return self.stereo == self.p_stereo == other.stereo
        return False

    @classmethod
    def _split_check_kwargs(cls, kwargs):
        r, p = super()._split_check_kwargs(kwargs)
        if not all(len(r.get(x, ())) == len(p.get(x, ())) for x in
                   ('charge', 'multiplicity', 'neighbors', 'hybridization')):
            raise ValueError('charge, multiplicity, neighbors, hybridization should be presented in both states '
                             'with same number of values')
        return r, p

    _hybridization_str = Atom._hybridization_str
    _stereo_str = Atom._stereo_str
    _multiplicity_str = Atom._multiplicity_str
    _charge_str = Atom._charge_str
    _static = {'element', 'isotope'}
    _p_static = {f'p_{x}' for x in _static}
    _atom_factory = QueryAtom


class DynQueryBond(DynBond):
    def _update(self, value, kwargs):
        if isinstance(value, (DynQueryBond, QueryBond)):
            super()._update(value, kwargs)
        elif isinstance(value, (DynBond, Bond)):
            r, p = self._split_check_kwargs(kwargs)
            if isinstance(value, DynBond):
                r = {'stereo': None,
                     **{k: getattr(self._bond_factory, f'_{k}_check')(v) for k, v in value._reagent.items()}, **r}
                p = {'stereo': None,
                     **{k: getattr(self._bond_factory, f'_{k}_check')(v) for k, v in value._product.items()}, **p}
            else:
                value = {k: getattr(self._bond_factory, f'_{k}_check')(v) for k, v in value.items()}
                r = {'stereo': None, **value, **r}
                p = {'stereo': None, **value, **p}
            if not (r['order'] or p['order']):
                raise ValueError('empty bond not allowed')
            if r['stereo'] and not r['order'] or p['stereo'] and not p['order']:
                raise ValueError('stereo bond not nullable')

            self._reagent._update(r, {})
            self._product._update(p, {})
        else:
            super()._update(value, kwargs)

    def stringify(self, stereo=False):
        if self._reagent.order is None:
            return DynamicContainer('.', self._product.stringify(stereo=stereo))
        elif self.p_order is None:
            return DynamicContainer(self._reagent.stringify(stereo=stereo), '.')

        r, p = zip(*sorted(zip(self.order, self.p_order)))
        r_order = '<%s>' % ''.join(sorted(self._order_str[x] for x in r))
        p_order = '<%s>' % ''.join(sorted(self._order_str[x] for x in p))
        if stereo:
            if self.stereo:
                r_order += self._stereo_str[self.stereo]
            if self.p_stereo:
                p_order += self._stereo_str[self.p_stereo]

        return DynamicContainer(r_order, p_order)

    def weight(self, stereo=False):
        if self.order is None:
            r = ()
            p = tuple(sorted(self.p_order))
        elif self.p_order is None:
            r = tuple(sorted(self.order))
            p = ()
        else:
            r, p = zip(*sorted(zip(self.order, self.p_order)))

        if stereo:
            return r, p, self.stereo or 0, self.p_stereo or 0
        return r, p

    def __eq__(self, other):
        if isinstance(other, DynQueryBond):
            if self.order is None:
                if other.order is not None:
                    return False
                return sorted(self.p_order) == sorted(other.p_order)
            elif other.order is None:
                return False
            elif self.p_order is None:
                if other.p_order is not None:
                    return False
                return sorted(self.order) == sorted(other.order)
            elif other.p_order is None:
                return False
            return sorted(zip(self.order, self.p_order)) == sorted(zip(other.order, other.p_order))

        elif isinstance(other, DynBond):
            if self.order is None:
                return other.order is None and other.p_order in self.p_order
            elif self.p_order is None:
                return other.p_order is None and other.order in self.order
            return (other.order, other.p_order) in list(zip(self.order, self.p_order))

        elif isinstance(other, QueryBond):
            return self.order == self.p_order and sorted(self.order) == sorted(other.order)

        elif isinstance(other, Bond):
            # before check for not None
            return self.order == self.p_order and other.order in self.order
        return False

    def __ne__(self, other):
        """
        != equality checks with stereo
        """
        if self == other:
            if self.stereo is self.p_stereo is None:
                return True
            if isinstance(other, DynBond):
                return self.stereo == other.stereo and self.p_stereo == other.p_stereo
            return self.stereo == self.p_stereo == other.stereo
        return False

    _bond_factory = QueryBond
    _order_str = Bond._order_str
    _stereo_str = Bond._stereo_str


__all__ = ['DynQueryAtom', 'DynQueryBond']
