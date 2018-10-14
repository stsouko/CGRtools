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
"""
contains periodic table of elements classes
"""
from operator import ge, le, gt, lt
from .data import *
from .shells import *
from .tables import *


class PeriodicMeta(type):
    """
    metaclass for creation of all classes of periodic table
    """
    def __new__(mcs, defined_name, bases, attrs, **extra):
        name = extra.get('name', defined_name)
        return mcs.classes.get(name) or mcs.classes.setdefault(name, super().__new__(mcs, name, bases, attrs))

    def __repr__(cls):
        return f"<class '{__package__}.{cls.__name__}'>"

    classes = {}


class ElementMeta(PeriodicMeta):
    def __init__(cls, defined_name, bases, attrs, **extra):
        cls.__number = extra.get('number', 0)
        cls.__symbol = extra.get('name', defined_name)
        super().__init__(defined_name, bases, attrs)

    @property
    def number(cls):
        return cls.__number

    @property
    def symbol(cls):
        return cls.__symbol

    def __eq__(cls, other):
        if isinstance(other, str):
            return cls.__symbol == other
        elif isinstance(other, int):
            return cls.__number == other
        elif isinstance(other, type):
            if issubclass(other, Element):
                return cls.__number == other.number
        elif isinstance(other, Element):
            return cls.__number == other.number
        return False

    def __gt__(cls, other):
        return cls.__ops(other, gt)

    def __ge__(self, other):
        return self.__ops(other, ge)

    def __lt__(self, other):
        return self.__ops(other, lt)

    def __le__(self, other):
        return self.__ops(other, le)

    def __ops(cls, other, op):
        if isinstance(other, str):
            if other in elements_list:
                return op(cls.__number, elements_numbers[other])
        elif isinstance(other, int):
            return op(cls.__number, other)
        elif isinstance(other, type):
            if issubclass(other, Element):
                return op(cls.__number, other.number)
            raise TypeError(f'unorderable types {cls} and {other}')
        elif isinstance(other, Element):
            return op(cls.__number, other.number)
        raise TypeError(f'unorderable types {cls} and {type(other)}')


class Periodic(metaclass=PeriodicMeta):
    """
    Base class of elements periodic classes
    """
    def __repr__(self):
        return f'{type(self).__name__}()'


class Element(Periodic, metaclass=ElementMeta):
    """
    Base class for all elements
    """


def arab2roman(number):
    roma = []
    for arabic, roman in ((1000, "M"), (900, "CM"), (500, "D"), (400, "CD"), (100, "C"), (90, "XC"), (50, "L"),
                          (40, "XL"), (10, "X"), (9, "IX"), (5, "V"), (4, "IV"), (1, "I")):
        factor, number = divmod(number, arabic)
        roma.append(roman * factor)
    return ''.join(roma)


def get_element(symbol, number):
    common_isotope = common_isotopes[symbol]
    electrons = valence_electrons[symbol]
    configuration = ' '.join('%d%s%d' % (n, orbitals_names[l], e) for (n, l), e in
                             electrons_configuration[symbol].items())

    class Period(Periodic, name=f'Period{arab2roman(periods[symbol])}'):
        pass

    class Group(Periodic, name=f'Group{arab2roman(groups[symbol])}'):
        pass

    class Type(Periodic, name=classes[symbol]):
        pass

    class ElementClass(Element, Period, Group, Type, name=symbol, number=number):
        __slots__ = ('_ElementClass__charge', '_ElementClass__multiplicity', '_ElementClass__isotope',
                     'mapping', 'mark', 'x', 'y', 'z')

        def __init__(self, charge: int=0, multiplicity: int=None, isotope: int=common_isotope,
                     mapping: int=None, mark: str='0', x: float=0, y: float=0, z: float=0):
            self.__charge = charge
            self.__multiplicity = multiplicity
            self.__isotope = isotope
            self.mapping = mapping
            self.mark = mark
            self.x = x
            self.y = y
            self.z = z

        @property
        def charge(self):
            return self.__charge

        @property
        def multiplicity(self):
            return self.__multiplicity

        @property
        def isotope(self):
            return self.__isotope

        @property
        def radical(self):
            return radical_map[self.__multiplicity]

        @property
        def number(self):
            return number

        @property
        def symbol(self):
            return symbol

        @property
        def electron_configuration(self):
            return configuration

        @property
        def electrons(self):
            return electrons

        def check_atom(self):
            """
            check possibility of current charge and radical state of atom
            """
            return (symbol, self.__charge, self.radical) in atom_charge_radical

        def check_valence(self, neighbors):
            """
            check possibility of atom with current charge and radical state to have given bonded neighbors

            :param neighbors: list of pairs of (bond, symbol or Element object)
            """
            return self.get_valence(neighbors) is not None

        def get_valence(self, neighbors):
            """
            get total valence (include implicit hydrogens) of atom with current charge and radical state
            and given bonded neighbors

            :param neighbors: list of pairs of (bond, symbol or Element object)
            :return: valence number or None if valence impossible
            """
            if not self.check_atom():
                return None
            bonds = [x for x, _ in neighbors]
            res = atom_valences.get((symbol, self.__charge, self.radical, self.__bonds_sum(bonds)))
            if res is None:
                key = atom_valences_exceptions.get((symbol, self.__charge, self.radical, len(bonds)))
                if key:
                    neighbors = [x if isinstance(x, str) else x.symbol for _, x in neighbors]
                    res = key.get(tuple(sorted(zip(neighbors, bonds))))
                    if not res and not all(isinstance(x, Element) or isinstance(x, type) and issubclass(x, Element)
                                           or x in elements_list for x in neighbors):
                        raise TypeError('invalid neighbors type')
            return res

        def get_implicit_h(self, bonds):
            """
            get implicit hydrogens count of atom  with current charge and radical state and given bonds.
            Note: if atom valence invalid return 0. use check_valence for validation

            :param bonds: list of bonds
            """
            return atom_implicit_h.get((symbol, self.__charge, self.radical, self.__bonds_sum(bonds)), 0)

        @staticmethod
        def __bonds_sum(bonds):
            return int(sum(bonds_map[x] for x in bonds))

        def __eq__(self, other):
            if isinstance(other, str):
                return symbol == other
            elif isinstance(other, int):
                return number == other
            elif isinstance(other, type):
                if issubclass(other, Element):
                    return number == other.number
            elif isinstance(other, Element):
                return number == other.number and self.__isotope == other.isotope and \
                       self.__charge == other.charge and self.__multiplicity == other.multiplicity
            return False

        def __gt__(self, other):
            return self.__ops(other, gt)

        def __ge__(self, other):
            return self.__ops(other, ge)

        def __lt__(self, other):
            return self.__ops(other, lt)

        def __le__(self, other):
            return self.__ops(other, le)

        def __ops(self, other, op):
            if isinstance(other, str):
                if other in elements_list:
                    return op(number, elements_numbers[other])
            elif isinstance(other, int):
                return op(number, other)
            elif isinstance(other, type):
                if issubclass(other, Element):
                    return op(number, other.number)
                raise TypeError(f'unorderable types {type(self)} and {other}')
            elif isinstance(other, Element):
                return op((number, self.__isotope, self.__charge, self.__multiplicity),
                          (other.number, other.isotope, other.charge, other.multiplicity))
            raise TypeError(f'unorderable types {type(self)} and {type(other)}')

        def __repr__(self):
            r = []
            if self.__charge:
                r.append(str(self.__charge))
            if self.__multiplicity:
                r.append(f'multiplicity={self.__multiplicity}')
            if self.__isotope != common_isotope:
                r.append(f'isotope={self.__isotope}')
            if self.mapping:
                r.append(f'mapping={self.mapping}')
            if self.mark != '0':
                r.append(f'mark={self.mark}')
            if self.x or self.y or self.z:
                r.append(f'x={self.x}, y={self.y}, x={self.x}')

            r = ', '.join(r)
            return f'{type(self).__name__}({r})'

        def __hash__(self):
            return hash((number, self.__charge, self.__isotope, self.__multiplicity))

    return ElementClass


elements_set = set(elements)
elements_list = elements
elements_numbers = {s: n for n, s in enumerate(elements, start=1)}
groups = {elements[z - 1]: x for x, y in enumerate(groups, start=1) for z in y}
periods = {elements[z - 1]: x for x, y in enumerate(periods, start=1) for z in y}
classes = {elements[z - 1]: x.capitalize() for x, y in classes.items() for z in y}
radical_map = {1: 2, 2: 1, 3: 2, None: 0}
radical_unmap = {None: None, 0: None, 1: 2, 2: 3}
common_isotopes = dict(zip(elements_list, isotopes))
aromatic = ('B', 'C', 'N', 'P', 'O', 'S')
elements_classes = {s: get_element(s, n) for s, n in elements_numbers.items()}

del elements

locals().update(PeriodicMeta.classes)
__all__ = [Element.__name__] + list(PeriodicMeta.classes)
