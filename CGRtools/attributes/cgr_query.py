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
        ratom = self._reagent._atom
        patom = self._product._atom

        if stereo:
            if ratom['stereo']:
                rmi.append(self._stereo_str[ratom['stereo']])
            if patom['stereo']:
                pmi.append(self._stereo_str[patom['stereo']])
        if hybridization:
            if len(ratom['hybridization']) > 1:
                r, p = zip(*sorted(zip(ratom['hybridization'], patom['hybridization'])))
                rmi.append('<%s>' % ''.join(self._hybridization_str[x] for x in r))
                pmi.append('<%s>' % ''.join(self._hybridization_str[x] for x in p))
            else:
                if ratom['hybridization']:
                    rmi.append(self._hybridization_str[ratom['hybridization'][0]])
                if patom['hybridization']:
                    pmi.append(self._hybridization_str[patom['hybridization'][0]])
        if neighbors:
            if len(ratom['neighbors']) > 1:
                r, p = zip(*sorted(zip(ratom['neighbors'], patom['neighbors'])))
                rmi.append('<%s>' % ''.join(str(x) for x in r))
                pmi.append('<%s>' % ''.join(str(x) for x in p))
            else:
                if ratom['neighbors']:
                    rmi.append(str(ratom['neighbors'][0]))
                if patom['neighbors']:
                    pmi.append(str(patom['neighbors'][0]))
        if rmi:
            rmi.append(';')
            rmi.insert(0, ';')
        if pmi:
            pmi.append(';')
            pmi.insert(0, ';')

        if atom:
            if ratom['element'] == ('A',):
                atom = False
                rmi.insert(0, '*')
                pmi.insert(0, '*')
            elif len(ratom['element']) > 1:
                atom = False
                tmp = ','.join(sorted(ratom['element'], key=elements_numbers.get))
                rmi.insert(0, tmp)
                pmi.insert(0, tmp)
            else:
                if ratom['element'][0] not in ('C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B'):
                    atom = False
                rmi.insert(0, ratom['element'][0])
                pmi.insert(0, ratom['element'][0])

            if len(ratom['charge']) > 1:
                r, p = zip(*sorted(zip(ratom['charge'], patom['charge'])))
                rmi.append('<%s>' % ''.join(self._charge_str[x] for x in r))
                pmi.append('<%s>' % ''.join(self._charge_str[x] for x in p))
            else:
                if ratom['charge'] != (0,):
                    rmi.append(self._charge_str[ratom['charge'][0]])
                if patom['charge'] != (0,):
                    pmi.append(self._charge_str[patom['charge'][0]])
            if len(ratom['multiplicity']) > 1:
                r, p = zip(*sorted(zip(ratom['multiplicity'], patom['multiplicity'])))
                rmi.append('<%s>' % ''.join(self._multiplicity_str[x] for x in r))
                pmi.append('<%s>' % ''.join(self._multiplicity_str[x] for x in p))
            else:
                if ratom['multiplicity']:
                    rmi.append(self._multiplicity_str[ratom['multiplicity'][0]])
                if patom['multiplicity']:
                    pmi.append(self._multiplicity_str[patom['multiplicity'][0]])
            if isotope:
                if len(ratom['isotope']) > 1:
                    tmp = '<%s>' % ''.join(str(x) for x in sorted(ratom['isotope']))
                    rmi.insert(0, tmp)
                    pmi.insert(0, tmp)
                elif ratom['isotope']:
                    tmp = str(ratom['isotope'][0])
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
        ratom = self._reagent._atom
        patom = self._product._atom
        r_weight, p_weight = [], []
        if atom:
            tmp = tuple(sorted(elements_numbers[x] for x in ratom['element']))
            r_weight.append(tmp)
            p_weight.append(tmp)
            if isotope:
                tmp = tuple(sorted(ratom['isotope']))
                r_weight.append(tmp)
                p_weight.append(tmp)
            r, p = zip(*sorted(zip(ratom['charge'], patom['charge'])))
            r_weight.append(tuple(r))
            p_weight.append(tuple(p))
            if ratom['multiplicity']:
                r, p = zip(*sorted(zip(ratom['multiplicity'], patom['multiplicity'])))
                r_weight.append(tuple(r))
                p_weight.append(tuple(p))
            else:
                r_weight.append(())
                p_weight.append(())
        if stereo:
            r_weight.append(ratom['stereo'] or 0)
            p_weight.append(patom['stereo'] or 0)
        if hybridization:
            if ratom['hybridization']:
                r, p = zip(*sorted(zip(ratom['hybridization'], patom['hybridization'])))
                r_weight.append(tuple(r))
                p_weight.append(tuple(p))
            else:
                r_weight.append(())
                p_weight.append(())
        if neighbors:
            if ratom['neighbors']:
                r, p = zip(*sorted(zip(ratom['neighbors'], patom['neighbors'])))
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
                r = {**self._bond_factory._defaults,
                     **{k: getattr(self._bond_factory, f'_{k}_check')(v) for k, v in value._reagent.items()}, **r}
                p = {**self._bond_factory._defaults,
                     **{k: getattr(self._bond_factory, f'_{k}_check')(v) for k, v in value._product.items()}, **p}
            else:
                value = {k: getattr(self._bond_factory, f'_{k}_check')(v) for k, v in value.items()}
                r = {**self._bond_factory._defaults, **value, **r}
                p = {**self._bond_factory._defaults, **value, **p}
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
        elif self._product.order is None:
            return DynamicContainer(self._reagent.stringify(stereo=stereo), '.')

        r, p = zip(*sorted(zip(self._reagent.order, self._product.order)))
        r_order = '<%s>' % ''.join(sorted(self._order_str[x] for x in r))
        p_order = '<%s>' % ''.join(sorted(self._order_str[x] for x in p))
        if stereo:
            if self._reagent.stereo:
                r_order += self._stereo_str[self._reagent.stereo]
            if self._product.stereo:
                p_order += self._stereo_str[self._product.stereo]

        return DynamicContainer(r_order, p_order)

    def weight(self, stereo=False):
        if self._reagent.order is None:
            r = ()
            p = tuple(sorted(self._product.order))
        elif self._product.order is None:
            r = tuple(sorted(self._reagent.order))
            p = ()
        else:
            r, p = zip(*sorted(zip(self._reagent.order, self._product.order)))

        if stereo:
            return r, p, self._reagent.stereo or 0, self._product.stereo or 0
        return r, p

    def __eq__(self, other):
        if isinstance(other, DynQueryBond):
            if self._reagent.order is None:
                if other._reagent.order is not None:
                    return False
                return sorted(self._product.order) == sorted(other._product.order)
            elif other._reagent.order is None:
                return False
            elif self._product.order is None:
                if other._product.order is not None:
                    return False
                return sorted(self._reagent.order) == sorted(other._reagent.order)
            elif other._product.order is None:
                return False
            return sorted(zip(self._reagent.order, self._product.order)) == \
                sorted(zip(other._reagent.order, other._product.order))

        elif isinstance(other, DynBond):
            if self._reagent.order is None:
                return other._reagent.order is None and other._product.order in self._product.order
            elif self._product.order is None:
                return other._product.order is None and other._reagent.order in self._reagent.order
            return (other._reagent.order, other._product.order) in list(zip(self._reagent.order, self._product.order))

        elif isinstance(other, QueryBond):
            return self._reagent.order == self._product.order and sorted(self._reagent.order) == sorted(other.order)

        elif isinstance(other, Bond):
            # before check for not None
            return self._reagent.order == self._product.order and other.order in self._reagent.order
        return False

    def __ne__(self, other):
        """
        != equality checks with stereo
        """
        if self == other:
            if self._reagent.stereo is self._product.stereo is None:
                return True
            if isinstance(other, DynBond):
                return self._reagent.stereo == other._reagent.stereo and self._product.stereo == other._product.stereo
            return self._reagent.stereo == self._product.stereo == other.stereo
        return False

    _bond_factory = QueryBond
    _order_str = Bond._order_str
    _stereo_str = Bond._stereo_str


__all__ = ['DynQueryAtom', 'DynQueryBond']
