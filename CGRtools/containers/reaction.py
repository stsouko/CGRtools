# -*- coding: utf-8 -*-
#
#  Copyright 2017, 2018 Ramil Nugmanov <stsouko@live.ru>
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
from warnings import warn
from .cgr import CGRContainer
from .query import QueryContainer
from .molecule import MoleculeContainer
from ..algorithms import hash_cgr_string


class MindfulList:
    """list with self-checks of modification. need for control of ReactionContainer caches actuality"""
    def __init__(self, data=None):
        self.__data = [] if data is None else list(data)
        self.__check = True

    def get_state(self):
        """return True if structure data list changed from previous checking time"""
        tmp = self.__check
        self.__check = False
        return tmp

    def append(self, obj):
        self.__check = True
        self.__data.append(obj)

    def clear(self):
        self.__check = True
        self.__data.clear()

    def extend(self, iterable):
        self.__check = True
        self.__data.extend(iterable)

    def insert(self, index, obj):
        self.__check = True
        self.__data.insert(index, obj)

    def pop(self, index=None):
        self.__check = True
        return self.__data.pop(index)

    def __delitem__(self, index):
        self.__check = True
        del self.__data[index]

    def __iadd__(self, obj):
        self.__check = True
        self.__data.append(obj)
        return self

    def __add__(self, other):
        return self.__class__(self.__data + list(other))

    def __getitem__(self, index):
        return self.__data[index]

    def __iter__(self):
        return iter(self.__data)

    def __len__(self):
        return len(self.__data)

    def __setitem__(self, key, value):
        self.__check = True
        self.__data[key] = value

    def __repr__(self):
        return '%s([%s])' % (self.__class__.__name__, ', '.join(repr(x) for x in self.__data))

    def __str__(self):
        return '[%s]' % ', '.join(str(x) for x in self.__data)


class ReactionContainer:
    """reaction storage. contains reagents, products and reactants lists"""
    __slots__ = ('__reagents', '__products', '__reactants', '__meta', '__signatures', '__pickle')

    def __init__(self, reagents=None, products=None, reactants=None, meta=None, substrats=None):
        """
        new empty or filled reaction object creation

        :param reagents: list of MoleculeContainers [or CGRContainers] in left side of reaction
        :param products: right side of reaction. see reagents
        :param reactants: middle side of reaction: solvents, catalysts, etc. see reagents
        :param meta: dictionary of metadata. like DTYPE-DATUM in RDF

        """
        if substrats:
            warn('deprecated key. use reagents instead', DeprecationWarning)

        self.__reagents = MindfulList(reagents if substrats is None else substrats)
        self.__products = MindfulList(products)
        self.__reactants = MindfulList(reactants)
        self.__meta = meta or {}
        self.__signatures = {}
        self.__pickle = None

    def __getitem__(self, item):
        if item == 'substrats':
            warn('deprecated key. use reagents instead', DeprecationWarning)
            return self.__reagents
        elif item == 'reagents':
            return self.__reagents
        elif item == 'products':
            return self.__products
        elif item == 'reactants':
            return self.__reactants
        elif item == 'meta':
            return self.__meta
        else:
            raise Exception('invalid key: %s' % item)

    def __getstate__(self):
        return dict(reagents=list(self.__reagents), meta=self.meta,
                    products=list(self.__products), reactants=list(self.__reactants))

    def __setstate__(self, state):
        self.__init__(**state)

    def pickle(self):
        """return json serializable reaction"""
        return dict(reagents=[x.pickle() for x in self.__reagents], meta=self.meta,
                    products=[x.pickle() for x in self.__products], reactants=[x.pickle() for x in self.__reactants])

    @staticmethod
    def unpickle(data):
        """convert json serializable reaction into ReactionContainer object instance"""
        return ReactionContainer(reagents=[CGRContainer.unpickle(x) for x in
                                           (data['reagents'] if 'reagents' in data else data['substrats'])],
                                 products=[CGRContainer.unpickle(x) for x in data['products']], meta=data['meta'],
                                 reactants=[CGRContainer.unpickle(x) for x in data.get('reactants', [])])

    @property
    def substrats(self):  # reverse compatibility
        warn('deprecated key. use reagents instead', DeprecationWarning)
        return self.__reagents

    @property
    def reagents(self):
        """reagents list. see products"""
        return self.__reagents

    @property
    def products(self):
        """list of CGRs or/and Molecules in products side"""
        return self.__products

    @property
    def reactants(self):
        """reactants list. see products"""
        return self.__reactants

    @property
    def meta(self):
        """dictionary of metadata. like DTYPE-DATUM in RDF"""
        return self.__meta

    def copy(self):
        """
        get copy of object

        :return: ReactionContainer
        """
        return self.__class__(reagents=[x.copy() for x in self.__reagents], meta=self.__meta.copy(),
                              products=[x.copy() for x in self.__products],
                              reactants=[x.copy() for x in self.__reactants])

    def get_signature_hash(self, *args, **kwargs):
        """
        get 40bytes hash of signature string. see get_signature

        :return: bytes
        """
        return hash_cgr_string(self.get_signature(*args, **kwargs))

    def get_signature(self, isotope=False, stereo=False, hybridization=False, neighbors=False, element=True,
                      flush_cache=False, hyb=False):
        """
        return string representation of reaction with molecules
        in order same as in lists of reagents, reactants, products.
        CAUTION: if reaction contains CGRs. signature will be unobvious

        :param isotope: set isotope marks
        :param stereo: set stereo marks
        :param hybridization: set hybridization mark of atom
        :param neighbors: set neighbors count mark of atom
        :param element: set elements marks and charges of atoms
        :param flush_cache: recalculate signature if True
        :param hyb: deprecated. see hybridization arg
        """
        if hyb:
            warn('attr hyb is deprecated, use hybridization instead', DeprecationWarning)
            hybridization = hyb

        if flush_cache or self.__signatures is None or any(x.get_state() for x in
                                                           (self.__reagents, self.__reactants, self.__products)):
            self.__signatures = {}

        k = (isotope, element, stereo, hybridization, neighbors)
        out = self.__signatures.get(k)
        if not out:
            sig = []
            for ml in (self.__reagents, self.__reactants, self.__products):
                ms = []
                for m in ml:
                    mol = m.get_signature(isotope=isotope, stereo=stereo, hybridization=hybridization,
                                          neighbors=neighbors, element=element)
                    ms.append('{%s}' % mol if isinstance(m, QueryContainer) else mol)
                sig.append('.'.join(ms))
            self.__signatures[k] = out = '>'.join(sig)
        return out

    def implicify_hydrogens(self):
        """
        remove explicit hydrogens if possible

        :return: number of removed hydrogens
        """
        total = 0
        for ml in (self.__reagents, self.__reactants, self.__products):
            for m in ml:
                if isinstance(m, MoleculeContainer):
                    total += m.implicify_hydrogens()
        if total:
            self.flush_cache()
        return total

    def explicify_hydrogens(self):
        """
        add explicit hydrogens to atoms

        :return: number of added atoms
        """
        total = 0
        for ml in (self.__reagents, self.__reactants, self.__products):
            for m in ml:
                if isinstance(m, MoleculeContainer):
                    total += m.explicify_hydrogens()
        if total:
            self.flush_cache()
        return total

    def aromatize(self):
        """
        convert structures to aromatic form. works only for Molecules and CGRs

        :return: number of processed molecules
        """
        total = 0
        for ml in (self.__reagents, self.__reactants, self.__products):
            for m in ml:
                if isinstance(m, MoleculeContainer):
                    res = m.aromatize()
                    if isinstance(m, CGRContainer):
                        if res[0] or res[1]:
                            total += 1
                    elif res:
                        total += 1
        if total:
            self.flush_cache()
        return total

    def flush_cache(self):
        """clear cached signatures and representation strings. use if structures objects in reaction object changed"""
        self.__pickle = self.__signatures = None

    def __str__(self):
        return self.get_signature(isotope=True, stereo=True, hybridization=True, neighbors=True)

    def __repr__(self):
        if self.__pickle is None or any(x.get_state() for x in (self.__reagents, self.__reactants, self.__products)):
            self.__pickle = '%s.unpickle(%s)' % (self.__class__.__name__, self.pickle())
        return self.__pickle


class MergedReaction:
    """represent reactions as single disjointed reagents and single disjointed products graphs"""
    __slots__ = ('__reagents', '__products', '__meta', '__signatures', '__pickle')

    def __init__(self, reagents=None, products=None, meta=None):
        self.__reagents = reagents
        self.__products = products
        self.__meta = meta or {}
        self.__signatures = {}
        self.__pickle = None

    def __getstate__(self):
        return dict(reagents=self.__reagents, products=self.__products, meta=self.__meta)

    def __setstate__(self, state):
        self.__init__(**state)

    @property
    def reagents(self):
        """disjointed reagents graph"""
        return self.__reagents

    @property
    def products(self):
        """disjointed products graph"""
        return self.__products

    @property
    def meta(self):
        """dictionary of metadata. like DTYPE-DATUM in RDF"""
        return self.__meta

    def copy(self):
        """
        get copy of object

        :return: MergedReaction
        """
        return self.__class__(self.reagents.copy(), self.products.copy(), self.meta.copy())

    def get_signature_hash(self, *args, **kwargs):
        """
        get 40bytes hash of signature string. see get_signature

        :return: bytes
        """
        return hash_cgr_string(self.get_signature(*args, **kwargs))

    def get_signature(self, isotope=False, stereo=False, hybridization=False, neighbors=False, element=True,
                      flush_cache=False, hyb=False):
        """
        return string representation of reaction with unique atoms and molecules order
        CAUTION: if reaction contains CGRs. signature will be unobvious

        :param isotope: set isotope marks to string
        :param stereo: set stereo marks
        :param hybridization: set hybridization mark of atom
        :param neighbors: set neighbors count mark of atom
        :param element: set elements marks and charges of atoms
        :param flush_cache: recalculate signature if True
        :param hyb: deprecated. see hybridization arg
        """
        if hyb:
            warn('attr hyb is deprecated, use hybridization instead', DeprecationWarning)
            hybridization = hyb

        if flush_cache or self.__signatures is None:
            self.__signatures = {}

        k = (isotope, element, stereo, hybridization, neighbors)
        out = self.__signatures.get(k)
        if not out:
            r = self.reagents.get_signature(isotope=isotope, stereo=stereo, hybridization=hybridization,
                                            neighbors=neighbors, element=element)
            p = self.products.get_signature(isotope=isotope, stereo=stereo, hybridization=hybridization,
                                            neighbors=neighbors, element=element)
            self.__signatures[k] = out = '%s>>%s' % ('{%s}' % r if isinstance(self.reagents, QueryContainer) else r,
                                                     '{%s}' % p if isinstance(self.products, QueryContainer) else p)
        return out

    def flush_cache(self):
        """clear cached signatures and representation strings. use if structure changed"""
        self.__pickle = self.__signatures = None

    def __str__(self):
        return self.get_signature(isotope=True, stereo=True, hybridization=True, neighbors=True)

    def __repr__(self):
        if self.__pickle is None:
            self.__pickle = '%s(%s, %s)' % (self.__class__.__name__, repr(self.reagents), repr(self.products))
        return self.__pickle


__all__ = [ReactionContainer.__name__, MergedReaction.__name__]
