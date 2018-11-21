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
from .cgr import DynamicContainer, DynBond
from .molecule import Bond
from .query import QueryBond


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


__all__ = ['DynQueryBond']
