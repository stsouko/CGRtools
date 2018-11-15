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

    def format(self, *args, **kwargs):
        return DynamicContainer(self.reagent.format(*args, **kwargs), self.product.format(*args, **kwargs))

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
    __slots__ = ('reagent', 'product')

    def __init__(self, bond=None, p_bond=None):
        if bond is None or p_bond is None:
            if bond != p_bond:
                warning('only one atom state passed. ignored')
            bond = Bond(allow_none=True)
            p_bond = Bond(allow_none=True)
        else:
            if isinstance(bond, DynBond):
                bond = bond.reagent.copy()
            elif isinstance(bond, Bond):
                bond = Bond(**bond, allow_none=True)
            else:
                raise TypeError('invalid bond passed')
            if isinstance(p_bond, DynBond):
                p_bond = p_bond.product.copy()
            elif isinstance(p_bond, Bond):
                p_bond = Bond(**p_bond, allow_none=True)
            else:
                raise TypeError('invalid p_bond passed')
            if bond.order == p_bond.order is None:
                raise ValueError('empty bond not allowed')
        super().__setattr__('reagent', bond)
        super().__setattr__('product', p_bond)

    def __getitem__(self, key):
        if key.startswith('p_'):
            return self.product[key[2:]]
        return self.reagent[key]

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __getattr__(self, key):
        if key.startswith('p_'):
            return getattr(self.product, key[2:])
        return getattr(self.reagent, key)

    def __setattr__(self, key, value):
        if key.startswith('p_'):
            if key == 'p_order' and value == self.reagent.order is None:
                raise ValueError('empty bond not allowed')
            setattr(self.product, key[2:], value)
        else:
            if key == 'order' and value == self.product.order is None:
                raise ValueError('empty bond not allowed')
            setattr(self.reagent, key, value)

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
        if isinstance(other, DynBond):
            return self.reagent == other.reagent and self.product == other.product
        elif isinstance(other, Bond):
            return self.reagent == other == self.product
        return False

    def __str__(self):
        return f'{self.reagent}>>{self.product}'

    def format(self, *args, **kwargs):
        return DynamicContainer(self.reagent.format(*args, **kwargs), self.product.format(*args, **kwargs))

    def __repr__(self):
        return f'{self.reagent}>>{self.product}'

    def __hash__(self):
        return hash((self.reagent, self.product))

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

        if isinstance(value, DynBond):
            p_value = value.product
            value = value.reagent
        elif isinstance(value, Bond):
            if value.order is None:
                raise ValueError('empty bond not allowed')
            p_value = value
        else:
            p_value = value = dict(value)
            if value.get('order', True) is None:
                raise ValueError('empty bond not allowed')

        if kwargs:
            if kwargs.get('p_order', True) == kwargs.get('order', True) is None:
                raise ValueError('empty bond not allowed')
            self.reagent.update(value, **{k: v for k, v in kwargs.items() if not k.startswith('p_')})
            self.product.update(p_value, **{k[2:]: v for k, v in kwargs.items() if k.startswith('p_')})
        else:
            self.reagent.update(value)
            self.product.update(p_value)

    def copy(self):
        return type(self)(self.reagent, self.product)


__all__ = ['DynamicContainer', 'DynAtom', 'DynBond']
