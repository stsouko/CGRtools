# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
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
from collections.abc import MutableMapping
from itertools import chain
from .molecule import Bond, Atom
from ..periodictable import Element


class DynAttribute(MutableMapping):
    __slots__ = ('_reagent', '_product')

    def __init__(self):
        super().__setattr__('_reagent', self._factory(skip_checks=True))
        super().__setattr__('_product', self._factory(skip_checks=True))

    def __init_copy__(self, reagent, product):
        super().__setattr__('_reagent', reagent)
        super().__setattr__('_product', product)

    def update(self, *args, **kwargs):
        if len(args) > 1:
            raise TypeError('update expected at most 1 arguments')
        elif args:
            value = args[0]
        else:
            value = kwargs
            kwargs = ()
        self._update(value, kwargs)
        self.__dict__.clear()

    def __len__(self):
        return sum(1 for _ in self)

    def __setitem__(self, key, value):
        try:
            setattr(self, key, value)
        except AttributeError as e:
            raise KeyError from e

    def __delitem__(self, key):
        raise TypeError('attribute deletion impossible')

    def __getstate__(self):
        return {'reagent': self._reagent, 'product': self._product}

    def __setstate__(self, state):
        super().__setattr__('_reagent', state['reagent'])
        super().__setattr__('_product', state['product'])

    def copy(self):
        cls = type(self)
        copy = cls.__new__(cls)
        copy.__init_copy__(self._reagent.copy(), self._product.copy())
        return copy


class DynAtomAttribute(DynAttribute):
    def __getattr__(self, key):
        if key.startswith('p_'):
            if key in self._p_static:
                raise AttributeError(f'{key} is invalid')
            value = getattr(self._product, key[2:])
        else:
            value = getattr(self._reagent, key)
        self.__dict__[key] = value
        return value

    def __getitem__(self, key):
        if key.startswith('p_'):
            if key in self._p_static:
                raise KeyError(f'{key} is invalid')
            return self._product[key[2:]]
        return self._reagent[key]

    def __contains__(self, key):
        if key in self._p_static:
            return False
        if key.startswith('p_'):
            return key[2:] in self._product
        return key in self._reagent

    def __iter__(self):
        return chain(self._reagent, (f'p_{x}' for x in self._product if x not in self._static))

    def _split_check_kwargs(self, kwargs):
        r, p = {}, {}
        for k, v in kwargs.items():
            if k in self._p_static:
                raise KeyError(f'{k} is invalid')
            elif k in self._static:
                r[k] = p[k] = getattr(self._reagent, f'_{k}_check')(v)
            elif k.startswith('p_'):
                k = k[2:]
                p[k] = getattr(self._reagent, f'_{k}_check')(v)
            else:
                r[k] = getattr(self._reagent, f'_{k}_check')(v)
        return r, p


class DynAtom(DynAtomAttribute):
    def __setattr__(self, key, value):
        if key in self._p_static:
            raise AttributeError(f'{key} is invalid')
        elif key.startswith('p_'):
            key = key[2:]
            value = getattr(self._reagent, f'_{key}_check')(value)
            self._product[key] = value
        else:
            value = getattr(self._reagent, f'_{key}_check')(value)
            self._reagent[key] = value
            if key in self._static:
                self._product[key] = value
        self.__dict__.clear()

    def _update(self, value, kwargs):
        if isinstance(value, DynAtom):
            r, p = self._split_check_kwargs(kwargs)
            p_value = value._product
            value = value._reagent
        elif isinstance(value, (Atom, Element)):
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

    def __int__(self):
        return int(self._reagent) << 21 | int(self._product)

    def __eq__(self, other):
        if isinstance(other, DynAtom):
            return self._reagent == other._reagent and self._product == other._product
        elif isinstance(other, Atom):
            return self._reagent == other and self._product == other
        return False

    _factory = Atom
    _static = {'element', 'isotope', 'mapping'}
    _p_static = {f'p_{x}' for x in _static}


class DynBond(DynAttribute):
    def __setattr__(self, key, value):
        if key == 'p_order':
            if value is None:
                if not self.order:
                    raise ValueError('empty bond not allowed')
                elif self.p_stereo:
                    raise ValueError('stereo bond not nullable')
            else:
                value = self._reagent._order_check(value)
            self._product.order = value
        elif key == 'order':
            if value is None:
                if not self.p_order:
                    raise ValueError('empty bond not allowed')
                elif self.stereo:
                    raise ValueError('stereo bond not nullable')
            else:
                value = self._reagent._order_check(value)
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
        self.__dict__.clear()

    def _update(self, value, kwargs):
        if isinstance(value, (DynBond, Bond)):
            r, p = self._split_check_kwargs(kwargs)
            if r or p:
                if isinstance(value, DynBond):
                    r = {'stereo': value.stereo, 'order': value.order, **r}
                    p = {'stereo': value.p_stereo, 'order': value.p_order, **p}
                else:
                    r = {'stereo': value.stereo, 'order': value.order, **r}
                    p = {'stereo': value.stereo, 'order': value.order, **p}

                if not (r['order'] or p['order']):
                    raise ValueError('empty bond not allowed')
                if r['stereo'] and not r['order'] or p['stereo'] and not p['order']:
                    raise ValueError('stereo bond not nullable')
            elif isinstance(value, DynBond):
                r, p = value._reagent, value._product
            else:
                r = p = value
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
            if not ((r.get('order', True) or p.get('order') or self.p_order) and
                    (p.get('order', True) or r.get('order') or self.order)):
                raise ValueError('empty bond not allowed')
            if r.get('stereo') and (not r.get('order', True) or not self.order) or \
                    p.get('stereo') and (not p.get('order', True) or not self.p_order):
                raise ValueError('stereo bond not nullable')

        self._reagent._update(r, {})
        self._product._update(p, {})

    def __getattr__(self, key):
        value = getattr(self._product, key[2:]) if key.startswith('p_') else getattr(self._reagent, key)
        self.__dict__[key] = value
        return value

    def __getitem__(self, key):
        if key.startswith('p_'):
            return self._product[key[2:]]
        return self._reagent[key]

    def __contains__(self, key):
        if key.startswith('p_'):
            return key[2:] in self._product
        return key in self._reagent

    def __iter__(self):
        return chain(self._reagent, (f'p_{x}' for x in self._product))

    def __int__(self):
        return (self.order or 0) << 3 | (self.p_order or 0)

    def __eq__(self, other):
        """
        equality checks without stereo
        """
        if isinstance(other, DynBond):
            return self._reagent == other._reagent and self._product == other._product
        return self._reagent == self._product == other

    def _split_check_kwargs(self, kwargs):
        r, p = {}, {}
        for k, v in kwargs.items():
            if k == 'p_order':
                p['order'] = None if v is None else self._reagent._order_check(v)
            elif k == 'order':
                r['order'] = None if v is None else self._reagent._order_check(v)
            elif k == 'p_stereo':
                p['stereo'] = self._reagent._stereo_check(v)
            elif k == 'stereo':
                r['stereo'] = self._reagent._stereo_check(v)
            else:
                raise AttributeError('invalid bond attribute')
        return r, p

    _factory = Bond


__all__ = ['DynAtom', 'DynBond']
