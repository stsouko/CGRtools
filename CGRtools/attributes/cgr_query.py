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
from .cgr import DynAtom, DynBond, DynAtomAttribute, DynBondAttribute
from .molecule import Atom, Bond
from .query import QueryAtom, QueryBond
from ..periodictable import Element


class DynQueryAtom(DynAtomAttribute):
    def __setattr__(self, key, value):
        if key in self._p_static:
            raise AttributeError(f'{key} is invalid')
        elif key.startswith('p_'):
            key = key[2:]
            value = getattr(self._reagent, f'_{key}_check')(value)
            if key not in ('stereo', 'x', 'y', 'z'):
                if len(value) != len(self[key]):
                    raise ValueError(f'{key} lists in reagent and product should be equal size')
                value = sorted(zip(self[key], value))
                if len(value) != len(set(value)):
                    raise ValueError('duplicates found')
                self._reagent[key], value = zip(*value)
            self._product[key] = value
        else:
            value = getattr(self._reagent, f'_{key}_check')(value)
            if key in self._static:
                self._product[key] = value
            elif key not in ('stereo', 'x', 'y', 'z'):
                if len(value) != len(self._product[key]):
                    raise ValueError(f'{key} lists in reagent and product should be equal size')
                value = sorted(zip(value, self._product[key]))
                if len(value) != len(set(value)):
                    raise ValueError('duplicates found')
                value, self._product[key] = zip(*value)
            self._reagent[key] = value

    def _update(self, value, kwargs):
        if isinstance(value, (DynQueryAtom, DynAtom)):
            r, p = self._split_check_kwargs(kwargs)
            p_value = value._product
            value = value._reagent
        elif isinstance(value, (QueryAtom, Atom, Element)):
            r, p = self._split_check_kwargs(kwargs)
            p_value = value
        elif isinstance(value, type):
            if not issubclass(value, Element):
                ValueError('only CGRtools.periodictable.Element subclasses allowed')
            r, p = self._split_check_kwargs(kwargs)
            p_value = value
        else:
            if not isinstance(value, dict):
                try:
                    value = dict(value)
                except (TypeError, ValueError):
                    raise TypeError('invalid attrs sequence')

            value.update(kwargs)
            if not value:  # ad-hoc for add_nodes_from method
                return
            value, p_value = self._split_check_kwargs(value)
            r = p = {}

        self._reagent._update(value, r)
        self._product._update(p_value, p)

    def __eq__(self, other):
        if isinstance(other, DynAtom):
            return (self.element == ('A',) or other.element in self.element) and \
                   (self.isotope == () or other.isotope in self.isotope) and \
                   all((other[x], other[y]) in zip(self[x], self[y])
                       for x, y in (('charge', 'p_charge'),
                                    ('multiplicity', 'p_multiplicity'),
                                    ('hybridization', 'p_hybridization'),
                                    ('neighbors', 'p_neighbors')) if self[x])
        elif isinstance(other, DynQueryAtom):
            return (self.element == ('A',) or set(self.element).issuperset(other.element)) and \
                   (self.isotope == () or set(self.isotope).issuperset(other.isotope)) and \
                   all(set(zip(self[x], self[y])).issuperset(zip(other[x], other[y]))
                       for x, y in (('charge', 'p_charge'),
                                    ('multiplicity', 'p_multiplicity'),
                                    ('hybridization', 'p_hybridization'),
                                    ('neighbors', 'p_neighbors')) if self[x])
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

    def _split_check_kwargs(self, kwargs):
        r, p = super()._split_check_kwargs(kwargs)
        for k in ('charge', 'multiplicity', 'neighbors', 'hybridization'):
            if k in r:
                if k not in p:
                    raise ValueError(f'{k} attribute should be presented with p_{k}')
                elif len(r[k]) != len(p[k]):
                    raise ValueError(f'{k} lists in reagent and product should be equal size')
                value = sorted(zip(r[k], p[k]))
                if value:
                    if len(value) != len(set(value)):
                        raise ValueError('duplicates found')
                    r[k], p[k] = zip(*value)
        return r, p

    _factory = QueryAtom
    _static = {'element', 'isotope'}
    _p_static = {'p_element', 'p_isotope'}


class DynQueryBond(DynBondAttribute):
    def __setattr__(self, key, value):
        if key == 'p_order':
            if value is None:
                if not self.order:
                    raise ValueError('empty bond not allowed')
                elif self.p_stereo:
                    raise ValueError('stereo bond not nullable')
                self._reagent.order = tuple(sorted(self.order))
            else:
                value = self._reagent._order_check(value)
                if self.order:
                    if len(value) != len(self.order):
                        raise ValueError('order lists in reagent and product should be equal size')
                    value = sorted(zip(self.order, value))
                    if len(value) != len(set(value)):
                        raise ValueError('duplicates found')
                    self._reagent.order, value = zip(*value)
                elif len(value) != len(set(value)):
                    raise ValueError('duplicates found')
                else:
                    value = tuple(sorted(value))
            self._product.order = value
        elif key == 'order':
            if value is None:
                if not self.p_order:
                    raise ValueError('empty bond not allowed')
                elif self.stereo:
                    raise ValueError('stereo bond not nullable')
                self._product.order = tuple(sorted(self.p_order))
            else:
                value = self._reagent._order_check(value)
                if self.p_order:
                    if len(value) != len(self.p_order):
                        raise ValueError('order lists in reagent and product should be equal size')
                    value = sorted(zip(value, self.p_order))
                    if len(value) != len(set(value)):
                        raise ValueError('duplicates found')
                    value, self._product.order = zip(*value)
                elif len(value) != len(set(value)):
                    raise ValueError('duplicates found')
                else:
                    value = tuple(sorted(value))
            self._reagent.order = value
        elif key == 'p_stereo':
            if value is None:
                if not self.p_order:
                    raise ValueError('null-bond stereo impossible')
            else:
                value = self._reagent._stereo_check(value)
            self._product.stereo = value
        elif key == 'stereo':
            if value is None:
                if not self.order:
                    raise ValueError('null-bond stereo impossible')
            else:
                value = self._reagent._stereo_check(value)
            self._reagent.stereo = value
        else:
            raise AttributeError('invalid bond attribute')

    def _update(self, value, kwargs):
        if isinstance(value, (DynBond, Bond)):
            r, p = self._split_check_kwargs(kwargs)
            if r or p:
                if isinstance(value, DynQueryBond):
                    r = {'stereo': value.stereo, 'order': value.order, **r}
                    p = {'stereo': value.p_stereo, 'order': value.p_order, **p}
                elif isinstance(value, QueryBond):
                    r = {'stereo': value.stereo, 'order': value.order, **r}
                    p = {'stereo': value.stereo, 'order': value.order, **p}
                elif isinstance(value, DynBond):
                    r = {'stereo': value.stereo, 'order': value.order and (value.order,), **r}
                    p = {'stereo': value.p_stereo, 'order': value.p_order and (value.p_order,), **p}
                else:
                    r = {'stereo': value.stereo, 'order': (value.order,), **r}
                    p = {'stereo': value.stereo, 'order': (value.order,), **p}

                if r['order'] and p['order']:
                    if len(r['order']) != len(p['order']):
                        raise ValueError('order lists in reagent and product should be equal size')
                    value = sorted(zip(r['order'], p['order']))
                    if len(value) != len(set(value)):
                        raise ValueError('duplicates found')
                    r['order'], p['order'] = zip(*value)
                elif r['order']:
                    if len(r['order']) != len(set(r['order'])):
                        raise ValueError('duplicates found')
                    r['order'] = tuple(sorted(r['order']))
                elif p['order']:
                    if len(p['order']) != len(set(p['order'])):
                        raise ValueError('duplicates found')
                    p['order'] = tuple(sorted(p['order']))
                else:
                    raise ValueError('empty bond not allowed')

                if r['stereo'] and not r['order'] or p['stereo'] and not p['order']:
                    raise ValueError('stereo bond not nullable')
            elif isinstance(value, DynQueryBond):
                r, p = value._reagent, value._product
            elif isinstance(value, QueryBond):
                r = p = value
            elif isinstance(value, DynBond):
                r = {'stereo': value.stereo, 'order': value.order and (value.order,)}
                p = {'stereo': value.p_stereo, 'order': value.p_order and (value.p_order,)}
            else:
                r = {'stereo': value.stereo, 'order': (value.order,)}
                p = {'stereo': value.stereo, 'order': (value.order,)}
        else:
            if not isinstance(value, dict):
                try:
                    value = dict(value)
                except (TypeError, ValueError):
                    raise TypeError('invalid attrs sequence')

            value.update(kwargs)
            if not value:
                return  # ad-hoc for add_edges_from method

            r, p = self._split_check_kwargs(value)
            if r.get('order') and p.get('order'):
                if len(r['order']) != len(p['order']):
                    raise ValueError('order lists in reagent and product should be equal size')
                value = sorted(zip(r['order'], p['order']))
                if len(value) != len(set(value)):
                    raise ValueError('duplicates found')
                r['order'], p['order'] = zip(*value)
            elif r.get('order'):
                if 'order' in p or not self.p_order:
                    if len(r['order']) != len(set(r['order'])):
                        raise ValueError('duplicates found')
                    r['order'] = tuple(sorted(r['order']))
                else:
                    if len(r['order']) != len(self.p_order):
                        raise ValueError('order lists in reagent and product should be equal size')
                    value = sorted(zip(r['order'], self.p_order))
                    if len(value) != len(set(value)):
                        raise ValueError('duplicates found')
                    r['order'], p['order'] = zip(*value)
            elif p.get('order'):
                if 'order' in r or not self.order:
                    if len(p['order']) != len(set(p['order'])):
                        raise ValueError('duplicates found')
                    p['order'] = tuple(sorted(p['order']))
                else:
                    if len(p['order']) != len(self.order):
                        raise ValueError('order lists in reagent and product should be equal size')
                    value = sorted(zip(self.order, p['order']))
                    if len(value) != len(set(value)):
                        raise ValueError('duplicates found')
                    r['order'], p['order'] = zip(*value)
            elif 'order' in r and 'order' in p:
                raise ValueError('empty bond not allowed')

            if r.get('stereo') and (not r.get('order', True) or not self.order) or \
                    p.get('stereo') and (not p.get('order', True) or not self.p_order):
                raise ValueError('stereo bond not nullable')

        self._reagent._update(r, {})
        self._product._update(p, {})

    def weight(self, stereo=False):
        if stereo:
            return self.order or (), self.p_order or (), self.stereo or 0, self.p_stereo or 0
        return self.order or (), self.p_order or ()

    def __eq__(self, other):
        if isinstance(other, DynBond):
            if self.order is None:
                return other.order is None and other.p_order in self.p_order
            elif self.p_order is None:
                return other.p_order is None and other.order in self.order
            return (other.order, other.p_order) in zip(self.order, self.p_order)

        elif isinstance(other, DynQueryBond):
            if self.order is None:
                if other.order is not None:
                    return False
                return set(self.p_order).issuperset(other.p_order)
            elif other.order is None:
                return False
            elif self.p_order is None:
                if other.p_order is not None:
                    return False
                return set(self.order).issuperset(other.order)
            elif other.p_order is None:
                return False
            return set(zip(self.order, self.p_order)).issuperset(zip(other.order, other.p_order))

        elif isinstance(other, QueryBond):
            return self.order == self.p_order and set(self.order).issuperset(other.order)

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

    _factory = QueryBond


__all__ = ['DynQueryAtom', 'DynQueryBond']
