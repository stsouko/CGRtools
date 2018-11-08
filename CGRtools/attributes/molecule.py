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
from collections.abc import MutableMapping
from ..periodictable import Element, elements_classes


class Atom(MutableMapping):
    __slots__ = '__atom'

    def __init__(self):
        super().__setattr__('_Atom__atom', None)

    def __getitem__(self, key):
        """
        dict like access to atom's attrs
        """
        if key == 'element':
            return self.__atom.symbol
        elif key not in self.__acceptable:
            raise KeyError('unknown atom attribute')
        return getattr(self.__atom, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __getattr__(self, key):
        if key == 'element':
            key = 'symbol'
        return getattr(self.__atom, key)

    def __setattr__(self, key, value):
        if key in ('element', 'symbol'):
            if value not in elements_classes or value == 'A':
                raise ValueError('invalid atom symbol')
            if self.__atom is None:
                super().__setattr__('_Atom__atom', elements_classes[value]())
            elif value != self.__atom.symbol:
                attrs = {k: getattr(self.__atom, k) for k in self.__acceptable}
                del attrs['isotope']
                super().__setattr__('_Atom__atom', elements_classes[value](**attrs))
        elif self.__atom is None:
            raise TypeError('any atom not allowed')
        elif key in self.__unmutable:
            if not self.__unmutable[key](value):
                raise ValueError('invalid value of unmutable attribute')
            attrs = {k: getattr(self.__atom, k) for k in self.__acceptable}
            attrs[key] = value
            super().__setattr__('_Atom__atom', elements_classes[self.__atom.symbol](**attrs))
        elif key in self.__mutable:
            if not self.__mutable[key](value):
                raise ValueError('invalid value of mutable attribute')
            setattr(self.__atom, key, value)
        else:
            raise KeyError('unknown atom attributes not allowed')

    def __len__(self):
        return sum(1 for _ in self)

    def __iter__(self):
        """
        iterate other non-default atom's init attrs
        """
        if self.__atom is None:
            raise StopIteration
        yield 'element'
        if self.__atom.isotope != self.__atom.common_isotope:
            yield 'isotope'
        for k, d in self.__defaults.items():
            if d != getattr(self, k):
                yield k

    def __contains__(self, key):
        return key in self.__possible

    def __delitem__(self, key):
        raise NotImplemented('attribute deletion impossible')

    def __eq__(self, other):
        return self.__atom == other

    def __gt__(self, other):
        return self.__atom > other

    def __ge__(self, other):
        return self.__atom >= other

    def __lt__(self, other):
        return self.__atom < other

    def __le__(self, other):
        return self.__atom <= other

    def __repr__(self):
        return repr(self.__atom)

    @property
    def atom(self):
        """
        wrapped atom
        """
        return self.__atom

    def copy(self):
        atom = self.__atom
        return type(atom)(charge=atom.charge, multiplicity=atom.multiplicity, isotope=atom.isotope,
                          mapping=atom.mapping, mark=atom.mark, x=atom.x, y=atom.y, z=atom.z, stereo=atom.stereo)

    def update(self, *args, **kwargs):
        """
        update atom

        :param args: tuple with 1 or 0 elements. element can be dict of atom attrs or atom object or atom class.
        :param kwargs: atom attrs. has precedence other args[0]
        """
        if len(args) > 1:
            raise TypeError('update expected at most 1 arguments')
        elif args:
            value = args[0]
        else:
            value = kwargs
            kwargs = ()

        if isinstance(value, type):
            if not issubclass(value, Element):
                ValueError('only CGRtools.periodictable.Element subclasses allowed')
            if kwargs.keys() - self.__mutable.keys():  # not need to check values
                raise KeyError('unmutable atom attributes not allowed')

            if self.__atom is None:
                super().__setattr__('_Atom__atom', value(**kwargs))
            else:
                attrs = {k: getattr(self.__atom, k) for k in self.__mutable}
                attrs.update(kwargs)
                super().__setattr__('_Atom__atom', value(**attrs))

        elif isinstance(value, Element):
            if kwargs.keys() - self.__mutable.keys():
                raise KeyError('unmutable atom attributes not allowed')

            attrs = {k: getattr(value, k) for k in self.__mutable}
            attrs.update(kwargs)
            if self.__atom is None or self.__atom != value:
                super().__setattr__('_Atom__atom', type(value)(charge=value.charge, isotope=value.isotope,
                                                               multiplicity=value.multiplicity, **attrs))
            else:
                if not all(self.__mutable[k](v) for k, v in kwargs.items()):
                    raise ValueError('invalid value of mutable attribute')
                for k, v in attrs.items():
                    setattr(self.__atom, k, v)
        else:
            if not isinstance(value, dict):
                try:
                    value = dict(value)
                except (TypeError, ValueError):
                    raise TypeError('invalid attrs sequence')

            value.update(kwargs)
            if not value:  # ad-hoc for add_nodes_from method
                return
            if value.keys() - self.__possible:
                raise ValueError('unknown atom attributes not allowed')

            if self.__atom is None:
                if 'element' not in value or value['element'] not in elements_classes or value['element'] == 'A':
                    raise ValueError('invalid atom symbol')

                super().__setattr__('_Atom__atom', elements_classes[value['element']](**{k: v for k, v in value.items()
                                                                                         if k != 'element'}))
            elif value.keys() - self.__mutable.keys():
                attrs = {k: getattr(self.__atom, k) for k in self.__acceptable}
                attrs.update(value)

                if 'element' in value:
                    if value['element'] not in elements_classes or value['element'] == 'A':
                        raise ValueError('invalid atom symbol')

                    del attrs['element']
                    if 'isotope' not in value:
                        del attrs['isotope']
                    super().__setattr__('_Atom__atom', elements_classes[value['element']](**attrs))
                else:
                    super().__setattr__('_Atom__atom', elements_classes[self.__atom.symbol](**attrs))
            else:
                if not all(self.__mutable[k](v) for k, v in value.items()):
                    raise ValueError('invalid value of mutable attribute')
                for k, v in value.items():
                    setattr(self.__atom, k, v)

    __defaults = {'mapping': None, 'mark': '0', 'x': 0, 'y': 0, 'z': 0, 'stereo': None, 'charge': 0,
                  'multiplicity': None}
    __mutable = {'mapping': lambda x: x is None or isinstance(x, int), 'mark': lambda x: isinstance(x, str),
                 'x': lambda x: isinstance(x, (float, int)), 'y': lambda x: isinstance(x, (float, int)),
                 'z': lambda x: isinstance(x, (float, int)), 'stereo': lambda x: x in {None, -1, 1, 0},
                 'color': lambda x: x is None or isinstance(x, dict) and all(isinstance(y, int) for y in x) and
                                        all(isinstance(y, str) for y in x.values())}
    __unmutable = {'isotope': lambda x: isinstance(x, int), 'charge': lambda x: isinstance(x, int),
                   'multiplicity': lambda x: x is None or isinstance(x, int)}
    __acceptable = {*__mutable, *__unmutable}
    __possible = {'element', *__acceptable}


class Bond(MutableMapping):
    __slots__ = ('order', 'stereo', '__allow_none')

    def __init__(self, order=1, stereo=None, allow_none=False):
        self.order = order
        self.stereo = stereo
        super().__setattr__('_Bond__allow_none', allow_none)

    def __getitem__(self, key):
        if key not in self.__acceptable:
            raise KeyError('unknown bond attribute')
        return getattr(self, key)

    def __setitem__(self, key, value):
        setattr(self, key, value)

    def __setattr__(self, key, value):
        if key not in self.__acceptable:
            raise AttributeError('unknown bond attribute')
        if value not in self.__acceptable[key] and value is None and not self.__allow_none:
            raise ValueError(f"attribute '{key}' value should be from acceptable list: {self.__acceptable[key]}")
        super().__setattr__(key, value)

    def __len__(self):
        return sum(1 for _ in self)

    def __iter__(self):
        """
        iterate other non-default bonds attrs

        need for nx.readwrite.json_graph.node_link_data
        """
        if self.order != self.__defaults['order']:
            yield 'order'
        if self.stereo != self.__defaults['stereo']:
            yield 'stereo'

    def __contains__(self, key):
        return key in self.__acceptable

    def __delitem__(self, key):
        raise NotImplemented('attribute deletion impossible')

    def __eq__(self, other):
        if isinstance(other, Bond):
            return self.order == other.order
        return False

    def __repr__(self):
        s = '' if self.stereo is None else f', stereo={self.stereo}'
        return f'{type(self).__name__}({self.order}{s})'

    def __hash__(self):
        return hash((self.order, self.stereo))

    def update(self, *args, **kwargs):
        """
        update bond

        :param args: tuple with 1 or 0 elements. element can be Bond object or dict of bond attrs.
        :param kwargs: bond attrs. has precedence other args[0]
        """
        if len(args) > 1:
            raise TypeError('update expected at most 1 arguments')
        elif args:
            value = args[0]
        else:
            value = kwargs
            kwargs = ()

        if isinstance(value, Bond):
            try:
                if not all(v in self.__acceptable[k] or v is None and self.__allow_none for k, v in kwargs.items()):
                    raise ValueError('invalid attribute value')
            except KeyError:
                raise KeyError('unknown bond attributes not allowed')

            super().__setattr__('order', value.order)
            super().__setattr__('stereo', value.stereo)
            for k, v in kwargs.items():
                super().__setattr__(k, v)

        else:
            if not isinstance(value, dict):
                try:
                    value = dict(value)
                except (TypeError, ValueError):
                    raise TypeError('invalid attrs sequence')

            value.update(kwargs)
            if not value:
                return  # ad-hoc for add_edges_from method

            try:
                if not all(v in self.__acceptable[k] or v is None and self.__allow_none for k, v in value.items()):
                    raise ValueError('invalid attribute value')
            except KeyError:
                raise KeyError('unknown bond attributes not allowed')

            for k, v in value.items():
                super().__setattr__(k, v)

    def copy(self):
        return type(self)(self.order, self.stereo)

    __acceptable = {'order': {1, 2, 3, 4, 9}, 'stereo': {None, -1, 1, 0}}
    __defaults = {'order': 1, 'stereo': None}


__all__ = ['Atom', 'Bond']
