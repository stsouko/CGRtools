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
from .cgr import DynAtom
from .molecule import Bond, Atom
from ..periodictable import Element, elements_classes, elements_numbers


class QueryAtom(MutableMapping):
    __slots__ = '_atom'

    def __init__(self, atom=None):
        if atom is None:
            super().__setattr__('_atom', self.__defaults.copy())
        elif isinstance(atom, Atom):
            super().__setattr__('_atom', {'element': (atom.symbol,),  'charge': (atom.charge,), 'stereo': atom.stereo,
                                          'multiplicity': (atom.multiplicity,) if atom.multiplicity else (),
                                          'neighbors': (atom.neighbors,) if atom.neighbors else (),
                                          'hybridization': (atom.hybridization,) if atom.hybridization else (),
                                          'isotope': (atom.isotope,) if atom.isotope != atom.common_isotope else ()})
        elif isinstance(atom, Element):
            super().__setattr__('_atom', {'element': (atom.symbol,), 'charge': (atom.charge,),
                                          'multiplicity': (atom.multiplicity,) if atom.multiplicity else (),
                                          'stereo': None, 'neighbors': (), 'hybridization': (),
                                          'isotope': (atom.isotope,) if atom.isotope != atom.common_isotope else ()})
        elif isinstance(atom, dict):
            if not all(self.__possible[k](v) for k, v in atom.items()):
                raise ValueError('invalid query atom data')
            super().__setattr__('_atom', {**self.__defaults, **atom})
        else:
            raise TypeError('invalid atom passed')

    def __getitem__(self, key):
        """
        dict like access to atom's attrs
        """
        if key == 'symbol':
            return self._atom['element']
        return self._atom[key]

    def __setitem__(self, key, value):
        if key == 'symbol':
            key = 'element'
        elif key not in self.__possible:
            raise KeyError('unknown atom attributes not allowed')
        if not self.__possible[key](value):
            raise ValueError('invalid atom attribute')
        self._atom[key] = value

    def __getattr__(self, key):
        if key == '__dict__':
            raise AttributeError()
        elif key == 'symbol':
            return self._atom['element']
        return self._atom[key]

    def __setattr__(self, key, value):
        self.__setitem__(key, value)

    def __len__(self):
        return sum(1 for _ in self)

    def __iter__(self):
        """
        iterate other non-default atom's init attrs
        """
        for k, d in self.__defaults.items():
            if d != self._atom[k]:
                yield k

    def __delitem__(self, key):
        raise TypeError('attribute deletion impossible')

    def __eq__(self, other):
        if isinstance(other, QueryAtom):
            return all(sorted(self._atom[attr], key=lambda x: x or 0) == sorted(other[attr], key=lambda x: x or 0)
                       for attr in ('element', 'isotope', 'charge', 'multiplicity', 'hybridization', 'neighbors')) \
                   and self._atom['stereo'] == other['stereo']
        elif isinstance(other, Atom):
            if (other.element in self._atom['element'] or self._atom['element'][0] == 'A') and \
                    all(getattr(other, attr) in self._atom[attr] for attr in
                        ('isotope', 'charge', 'multiplicity', 'hybridization', 'neighbors') if self._atom[attr]):
                if self._atom['stereo']:
                    return self._atom['stereo'] == other.stereo
                return True
        return False

    def __str__(self):
        return self.stringify()

    def stringify(self, atom=True, isotope=True, stereo=True, hybridization=True, neighbors=True):
        smi = []
        satom = self._atom
        if stereo and satom['stereo']:
            smi.append(Atom._stereo_str[satom['stereo']])
        if hybridization:
            if len(satom['hybridization']) > 1:
                smi.append('<%s>' % ''.join(Atom._hybridization_str[x] for x in sorted(satom['hybridization'])))
            elif satom['hybridization']:
                smi.append(Atom._hybridization_str[satom['hybridization'][0]])
        if neighbors:
            if len(satom['neighbors']) > 1:
                smi.append('<%s>' % ''.join(str(x) for x in sorted(satom['neighbors'])))
            elif satom['neighbors']:
                smi.append(str(satom['neighbors'][0]))
        if smi:
            smi.append(';')
            smi.insert(0, ';')

        if atom:
            if satom['element'] == ('A',):
                atom = False
                smi.insert(0, '*')
            elif len(satom['element']) > 1:
                atom = False
                smi.insert(0, ','.join(sorted(satom['element'], key=elements_numbers)))
            else:
                if satom['element'][0] not in ('C', 'N', 'O', 'P', 'S', 'F', 'Cl', 'Br', 'I', 'B'):
                    atom = False
                smi.insert(0, satom['element'][0])

            if len(satom['charge']) > 1:
                smi.append('<%s>' % ''.join(Atom._charge_str[x] for x in sorted(satom['charge'])))
            elif satom['charge'] != (0,):
                smi.append(Atom._charge_str[satom['charge'][0]])

            if len(satom['multiplicity']) > 1:
                smi.append('<%s>' % ''.join(Atom._multiplicity_str[x] for x in sorted(satom['multiplicity'])))
            elif satom['multiplicity']:
                smi.append(Atom._multiplicity_str[satom['multiplicity'][0]])

            if isotope:
                if len(satom['isotope']) > 1:
                    smi.insert(0, '<%s>' % ''.join(str(x) for x in sorted(satom['isotope'])))
                elif satom['isotope']:
                    smi.insert(0, str(satom['isotope'][0]))
        else:
            smi.insert(0, '*')

        if len(smi) != 1 or not atom:
            smi.insert(0, '[')
            smi.append(']')

        return ''.join(smi)

    def weight(self, atom=True, isotope=False, stereo=False, hybridization=False, neighbors=False):
        satom = self._atom
        weight = []
        if atom:
            weight.append(tuple(sorted(elements_numbers[x] for x in satom['element'])))
            if isotope:
                weight.append(tuple(sorted(satom['isotope'])))
            weight.append(tuple(sorted(satom['charge'])))
            weight.append(satom['multiplicity'] and tuple(sorted(satom['multiplicity'])))
        if stereo:
            weight.append(satom['stereo'] or 0)
        if hybridization:
            weight.append(satom['hybridization'] and tuple(sorted(satom['hybridization'])))
        if neighbors:
            weight.append(satom['neighbors'] and tuple(sorted(satom['neighbors'])))
        return tuple(weight)

    def copy(self):
        raise type(self)(self._atom)

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
            if not all(self.__possible[k](v) for k, v in kwargs.items()):
                raise ValueError('invalid query atom data')
            self._atom['element'] = (value.symbol,)
            self._atom.update(kwargs)

        elif isinstance(value, Element):
            if not all(self.__possible[k](v) for k, v in kwargs.items()):
                raise ValueError('invalid query atom data')
            self._atom.update(element=(value.symbol,), charge=(value.charge,),
                              multiplicity=(value.multiplicity,) if value.multiplicity else (),
                              isotope=(value.isotope,) if value.isotope != value.common_isotope else ())
            self._atom.update(kwargs)
        elif isinstance(value, Atom):
            if not all(self.__possible[k](v) for k, v in kwargs.items()):
                raise ValueError('invalid query atom data')
            self._atom.update(element=(value.element,), charge=(value.charge,), stereo=value.stereo,
                              multiplicity=(value.multiplicity,) if value.multiplicity else (),
                              isotope=(value.isotope,) if value.isotope != value.common_isotope else (),
                              neighbors=(value.neighbors,) if value.neighbors else (),
                              hybridization=(value.hybridization,) if value.hybridization else ())
            self._atom.update(kwargs)
        elif isinstance(value, DynAtom):
            raise TypeError('CGR unsupported')
        else:
            if not isinstance(value, dict):
                try:
                    value = dict(value)
                except (TypeError, ValueError):
                    raise TypeError('invalid attrs sequence')

            value.update(kwargs)
            if not value:  # ad-hoc for add_nodes_from method
                return

            if not all(self.__possible[k](v) for k, v in value.items()):
                raise ValueError('invalid query atom data')
            self._atom.update(value)

    __defaults = {'element': ('A',), 'isotope': (), 'charge': (0,), 'multiplicity': (),
                  'neighbors': (), 'hybridization': (), 'stereo': None}
    __possible = {'element': lambda x: x == ('A',) or isinstance(x, tuple) and
                                       all(x != 'A' and x in elements_classes for x in x) and len(set(x)) == len(x),
                  'stereo': lambda x: x in (None, -1, 1),
                  'isotope': lambda x: isinstance(x, tuple) and all(isinstance(x, int) and x > 0 for x in x) and
                                       len(set(x)) == len(x),
                  'charge': lambda x: isinstance(x, tuple) and all(isinstance(x, int) and -3 <= x <= 3 for x in x) and
                                      len(set(x)) == len(x) > 0,
                  'multiplicity': lambda x: isinstance(x, tuple) and all(isinstance(x, int) and 1 <= x <= 3 for x in x)
                                            and len(set(x)) == len(x),
                  'neighbors': lambda x: isinstance(x, tuple) and all(isinstance(x, int) and 0 <= x < 999 for x in x)
                                         and len(set(x)) == len(x),
                  'hybridization': lambda x: isinstance(x, tuple) and all(isinstance(x, int) and 1 <= x <= 4 for x in x)
                                             and len(set(x)) == len(x)}


class QueryBond(Bond):
    def __init__(self, order=None, stereo=None):
        super().__init__(order, stereo, True)

    def __eq__(self, other):
        if isinstance(other, QueryBond):
            return sorted(self.order, key=lambda x: x or 0) == sorted(other.order, key=lambda x: x or 0) and \
                   self.stereo == other.stereo
        elif isinstance(other, Bond):
            if other.order in self.order:
                if self.stereo:
                    return self.stereo == other.stereo
                return True
        return False

    def stringify(self, stereo=True):
        order = '<%s>' % ''.join(sorted(self._order_str[x] for x in sorted(self.order)))
        if stereo and self.stereo:
            return order + self._stereo_str[self.stereo]
        return order

    def __hash__(self):
        return hash((tuple(sorted(self.order)), self.stereo))

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
                if not all(self._acceptable[k](v) for k, v in kwargs.items()):
                    raise ValueError('invalid attribute value')
            except KeyError:
                raise KeyError('unknown bond attributes not allowed')
            if isinstance(value, QueryBond):
                super().__setattr__('order', value.order)
                super().__setattr__('stereo', value.stereo)
            else:
                super().__setattr__('order', [value.order])
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
                if not all(self._acceptable[k](v) for k, v in value.items()):
                    raise ValueError('invalid attribute value')
            except KeyError:
                raise KeyError('unknown bond attributes not allowed')

            for k, v in value.items():
                super().__setattr__(k, v)

    _acceptable = {'order': lambda x, _=None: isinstance(x, tuple) and all(x in (1, 2, 3, 4, 9) for x in x),
                   'stereo': lambda x, _=None: x in (None, -1, 1)}
    _defaults = {'order': (1,), 'stereo': None}


__all__ = ['QueryAtom', 'QueryBond']
