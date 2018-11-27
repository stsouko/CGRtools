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
from collections import namedtuple
from collections.abc import MutableMapping
from itertools import chain
from .molecule import Bond, Atom
from .query import QueryAtom
from ..periodictable import Element


DynamicContainer = namedtuple('DynamicContainer', ['reagent', 'product'])


class DynAtom(MutableMapping):
    __slots__ = ('_reagent', '_product')

    def __init__(self):
        super().__setattr__('_reagent', self._atom_factory(skip_checks=True))
        super().__setattr__('_product', self._atom_factory(skip_checks=True))

    def __init_copy__(self, reagent, product):
        super().__setattr__('_reagent', reagent)
        super().__setattr__('_product', product)

    def __setattr__(self, key, value):
        if key in self._p_static:
            raise AttributeError(f'{key} is invalid')
        elif key.startswith('p_'):
            key = key[2:]
            value = getattr(self._atom_factory, f'_{key}_check')(value)
            setattr(self._product, key, value)
        else:
            value = getattr(self._atom_factory, f'_{key}_check')(value)
            setattr(self._reagent, key, value)
            if key in self._static:
                setattr(self._product, key, value)

    def update(self, *args, **kwargs):
        """
        update atom

        :param args: tuple with 1 or 0 elements. element can be dict of atom attrs or atom object or atom class
        or DynAtom object.
        :param kwargs: atom attrs. has precedence other args[0]
        """
        if len(args) > 1:
            raise TypeError('update expected at most 1 arguments')
        elif args:
            value = args[0]
        else:
            value = kwargs
            kwargs = ()
        self._update(value, kwargs)

    def _update(self, value, kwargs):
        if isinstance(value, DynAtom):
            r, p = self._split_check_kwargs(kwargs)
            p_value = value._product
            value = value._reagent
            if isinstance(value, QueryAtom):
                raise TypeError('QueryAtom not supported')
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

    def stringify(self, *args, **kwargs):
        return DynamicContainer(self._reagent.stringify(*args, **kwargs), self._product.stringify(*args, **kwargs))

    def weight(self, *args, **kwargs):
        return DynamicContainer(self._reagent.weight(*args, **kwargs), self._product.weight(*args, **kwargs))

    def __eq__(self, other):
        if isinstance(other, DynAtom):
            return self._reagent == other._reagent and self._product == other._product
        elif isinstance(other, Atom):
            return self._reagent == other == self._product
        return False

    def __ne__(self, other):
        """
        != equality checks with stereo
        """
        if isinstance(other, DynAtom):
            return self._reagent != other._reagent and self._product != other._product
        return self._reagent != self._product != other

    def __getitem__(self, key):
        if key.startswith('p_'):
            if key in self._p_static:
                raise KeyError(f'{key} is invalid')
            return self._product[key[2:]]
        return self._reagent[key]

    def __setitem__(self, key, value):
        try:
            setattr(self, key, value)
        except AttributeError as e:
            raise KeyError from e

    def __getattr__(self, key):
        if key.startswith('p_'):
            if key in self._p_static:
                raise AttributeError(f'{key} is invalid')
            return getattr(self._product, key[2:])
        return getattr(self._reagent, key)

    def __len__(self):
        return len(self._reagent) + len(self._product)

    def __iter__(self):
        return chain(self._reagent, (f'p_{x}' for x in self._product))

    def __contains__(self, key):
        if key in self._p_static:
            return False
        if key.startswith('p_'):
            return key[2:] in self._product
        return key in self._reagent

    def __delitem__(self, key):
        raise TypeError('attribute deletion impossible')

    def __str__(self):
        return '%s>>%s' % self.stringify(stereo=True)

    def copy(self):
        cls = type(self)
        copy = cls.__new__(cls)
        copy.__init_copy__(self._reagent.copy(), self._product.copy())
        return copy

    @classmethod
    def _split_check_kwargs(cls, kwargs):
        r, p = {}, {}
        for k, v in kwargs.items():
            if k in cls._p_static:
                raise KeyError(f'{k} is invalid')
            elif k in cls._static:
                r[k] = p[k] = getattr(cls._atom_factory, f'_{k}_check')(v)
            elif k.startswith('p_'):
                k = k[2:]
                p[k] = getattr(cls._atom_factory, f'_{k}_check')(v)
            else:
                r[k] = getattr(cls._atom_factory, f'_{k}_check')(v)
        return r, p

    _static = {'element', 'isotope', 'mark', 'mapping'}
    _p_static = {f'p_{x}' for x in _static}
    _atom_factory = Atom


class DynBond(MutableMapping):
    __slots__ = ('_reagent', '_product')

    def __init__(self):
        super().__setattr__('_reagent', self._bond_factory(skip_checks=True))
        super().__setattr__('_product', self._bond_factory(skip_checks=True))

    def __init_copy__(self, reagent, product):
        super().__setattr__('_reagent', reagent)
        super().__setattr__('_product', product)

    def __setattr__(self, key, value):
        if key == 'p_order':
            if value is None:
                if not self.order:
                    raise ValueError('empty bond not allowed')
                elif self.p_stereo:
                    raise ValueError('stereo bond not nullable')
            else:
                value = self._bond_factory._order_check(value)
            self._product.order = value
        elif key == 'order':
            if value is None:
                if not self.p_order:
                    raise ValueError('empty bond not allowed')
                elif self.stereo:
                    raise ValueError('stereo bond not nullable')
            else:
                value = self._bond_factory._order_check(value)
            self._reagent.order = value
        else:
            if key.startswith('p_'):
                key = key[2:]
                container = self._product
            else:
                container = self._reagent
            setattr(container, key, getattr(self._bond_factory, f'_{key}_check')(value))

    def update(self, *args, **kwargs):
        """
        update bond

        :param args: tuple with 1 or 0 elements. element can be Bond|DynBond object or dict of bond attrs.
        :param kwargs: bond attrs. has precedence other args[0]
        """
        if len(args) > 1:
            raise TypeError('update expected at most 1 arguments')
        elif args:
            value = args[0]
        else:
            value = kwargs
            kwargs = ()
        self._update(value, kwargs)

    def _update(self, value, kwargs):
        if isinstance(value, (DynBond, Bond)):
            r, p = self._split_check_kwargs(kwargs)
            if r or p:
                if isinstance(value, DynBond):
                    r = {'stereo': None, **value._reagent, **r}
                    p = {'stereo': None, **value._product, **p}
                else:
                    r = {'stereo': None, **value, **r}
                    p = {'stereo': None, **value, **p}
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
            if r.get('stereo') and (not r.get('order', True) or not self.stereo) or \
                    p.get('stereo') and (not p.get('order', True) or not self.p_stereo):
                raise ValueError('stereo bond not nullable')

        self._reagent._update(r, {})
        self._product._update(p, {})

    def stringify(self, stereo=False):
        if not self.order:
            return DynamicContainer('.', self._product.stringify(stereo=stereo))
        elif not self.p_order:
            return DynamicContainer(self._reagent.stringify(stereo=stereo), '.')
        return DynamicContainer(self._reagent.stringify(stereo=stereo), self._product.stringify(stereo=stereo))

    def weight(self, stereo=False):
        if stereo:
            return self.order or 0, self.p_order or 0, self.stereo or 0, self.p_stereo or 0
        return self.order or 0, self.p_order or 0

    def __eq__(self, other):
        """
        == equality checks without stereo
        """
        if isinstance(other, DynBond):
            return self._reagent == other._reagent and self._product == other._product
        return self._reagent == self._product == other

    def __ne__(self, other):
        """
        != equality checks with stereo
        """
        if isinstance(other, DynBond):
            return self._reagent != other._reagent and self._product != other._product
        return self._reagent != self._product != other

    def __getitem__(self, key):
        if key.startswith('p_'):
            return self._product[key[2:]]
        return self._reagent[key]

    def __setitem__(self, key, value):
        try:
            setattr(self, key, value)
        except AttributeError as e:
            raise KeyError from e

    def __getattr__(self, key):
        if key.startswith('p_'):
            return getattr(self._product, key[2:])
        return getattr(self._reagent, key)

    def __len__(self):
        return len(self._reagent) + len(self._product)

    def __iter__(self):
        return chain(self._reagent, (f'p_{x}' for x in self._product))

    def __contains__(self, key):
        if key.startswith('p_'):
            return key[2:] in self._product
        return key in self._reagent

    def __delitem__(self, key):
        raise TypeError('attribute deletion impossible')

    def __str__(self):
        return '%s>>%s' % self.stringify(stereo=True)

    @classmethod
    def _split_check_kwargs(cls, kwargs):
        r, p = {}, {}
        for k, v in kwargs.items():
            if k.startswith('p_'):
                k = k[2:]
                d = p
            else:
                d = r
            d[k] = getattr(cls._bond_factory, f'_{k}_check')(v)
        return r, p

    def copy(self):
        cls = type(self)
        copy = cls.__new__(cls)
        copy.__init_copy__(self._reagent.copy(), self._product.copy())
        return copy

    _bond_factory = Bond


__all__ = ['DynamicContainer', 'DynAtom', 'DynBond']
