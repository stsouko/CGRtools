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
from logging import warning
from .molecule import Bond, Atom
from ..periodictable import Element


DynamicContainer = namedtuple('DynamicContainer', ['reagent', 'product'])


class DynAtom(MutableMapping):
    __slots__ = ('reagent', 'product')

    def __init__(self, atom=None, p_atom=None):
        if atom is None or p_atom is None:
            if atom != p_atom:
                warning('only one atom state passed. ignored')
            atom = Atom()
            p_atom = Atom()
        else:
            if isinstance(atom, DynAtom):
                atom = Atom(atom.reagent)
            elif isinstance(atom, (Atom, Element)):
                atom = Atom(atom)
            else:
                raise TypeError('invalid atom passed')
            if isinstance(p_atom, DynAtom):
                p_atom = Atom(p_atom.product)
            elif isinstance(p_atom, (Atom, Element)):
                p_atom = Atom(p_atom)
            else:
                raise TypeError('invalid p_atom passed')
            if atom.symbol != p_atom.symbol or atom.isotope != p_atom.isotope:
                raise ValueError('different atoms impossible')
        super().__setattr__('reagent', atom)
        super().__setattr__('product', p_atom)

    def __getitem__(self, key):
        if key.startswith('p_'):
            if key in self.__p_static:
                raise KeyError(f'{key} is invalid')
            return self.product[key[2:]]
        return self.reagent[key]

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __getattr__(self, key):
        if key.startswith('p_'):
            if key in self.__p_static:
                raise KeyError(f'{key} is invalid')
            return getattr(self.product, key[2:])
        return getattr(self.reagent, key)

    def __setattr__(self, key, value):
        if key.startswith('p_'):
            if key in self.__p_static:
                raise KeyError(f'{key} is invalid')
            setattr(self.product, key[2:], value)
        else:
            setattr(self.reagent, key, value)
            if key in self.__static:
                setattr(self.product, key, value)

    def __len__(self):
        return len(self.reagent) + len(self.product)

    def __iter__(self):
        return chain(self.reagent, (f'p_{x}' for x in self.product))

    def __contains__(self, key):
        if key.startswith('p_'):
            return key[2:] in self.product
        return key in self.reagent

    def __delitem__(self, key):
        raise TypeError('attribute deletion impossible')

    def __eq__(self, other):
        if isinstance(other, DynAtom):
            return self.reagent == other.reagent and self.product == other.product
        elif isinstance(other, Atom):
            return self.reagent == other == self.product
        return False

    def __str__(self):
        return f'{self.reagent}>>{self.product}'

    def stringify(self, *args, **kwargs):
        return DynamicContainer(self.reagent.format(*args, **kwargs), self.product.format(*args, **kwargs))

    def weight(self, *args, **kwargs):
        return DynamicContainer(self.reagent.weight(*args, **kwargs), self.product.weight(*args, **kwargs))

    @property
    def atom(self):
        """
        reagent and product state atoms
        """
        return DynamicContainer(self.reagent.atom, self.product.atom)

    def copy(self):
        return type(self)(self.reagent, self.product)

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

        if isinstance(value, DynAtom):
            p_value = value.product
            value = value.reagent
        else:
            p_value = value

        if kwargs:
            if kwargs.keys() & self.__p_static:
                raise KeyError(f'{self.__static} keys is static')
            self.reagent.update(value, **{k: v for k, v in kwargs.items() if not k.startswith('p_')})
            self.product.update(p_value, **{k[2:]: v for k, v in kwargs.items() if k.startswith('p_')},
                                **{k: v for k, v in kwargs.items() if k in self.__static})
        else:
            self.reagent.update(value)
            self.product.update(p_value)

    __static = {'color', 'element', 'isotope', 'mark', 'mapping'}
    __p_static = {f'p_{x}' for x in __static}


class DynBond(MutableMapping):
    __slots__ = ('_reagent', '_product')

    def __init__(self):
        super().__setattr__('_reagent', self._bond_factory(skip_checks=True))
        super().__setattr__('_product', self._bond_factory(skip_checks=True))

    def __init_copy__(self, parent):
        super().__setattr__('_reagent', parent._reagent.copy())
        super().__setattr__('_product', parent._product.copy())

    def __setattr__(self, key, value):
        if key == 'p_order':
            if value is None:
                if not self._reagent.order:
                    raise ValueError('empty bond not allowed')
                elif self._product.stereo:
                    raise ValueError('stereo bond not nullable')
            else:
                value = self._bond_factory._order_check(value)
            self._product.order = value
        elif key == 'order':
            if value is None:
                if not self._product.order:
                    raise ValueError('empty bond not allowed')
                elif self._reagent.stereo:
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
                    r = {**self._bond_factory._defaults, **value._reagent, **r}
                    p = {**self._bond_factory._defaults, **value._product, **p}
                else:
                    r = {**self._bond_factory._defaults, **value, **r}
                    p = {**self._bond_factory._defaults, **value, **p}
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
            if not ((r.get('order', True) or p.get('order') or self._product.order) and
                    (p.get('order', True) or r.get('order') or self._reagent.order)):
                raise ValueError('empty bond not allowed')
            if r.get('stereo') and (not r.get('order', True) or not self._reagent.stereo) or \
                    p.get('stereo') and (not p.get('order', True) or not self._product.stereo):
                raise ValueError('stereo bond not nullable')

        self._reagent._update(r, {})
        self._product._update(p, {})

    def stringify(self, stereo=False):
        if not self._reagent.order:
            return DynamicContainer('.', self._product.stringify(stereo=stereo))
        elif not self._product.order:
            return DynamicContainer(self._reagent.stringify(stereo=stereo), '.')
        return DynamicContainer(self._reagent.stringify(stereo=stereo), self._product.stringify(stereo=stereo))

    def weight(self, stereo=False):
        if stereo:
            return self._reagent.order or 0, self._product.order or 0,\
                   self._reagent.stereo or 0, self._product.stereo or 0
        return self._reagent.order or 0, self._product.order or 0

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
        copy.__init_copy__(self)
        return copy

    _bond_factory = Bond


__all__ = ['DynamicContainer', 'DynAtom', 'DynBond']
