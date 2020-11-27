# -*- coding: utf-8 -*-
#
#  Copyright 2019, 2020 Ramil Nugmanov <nougmanoff@protonmail.com>
#  Copyright 2019 Tagir Akhmetshin <tagirshin@gmail.com>
#  Copyright 2019 Dayana Bashirova <dayana.bashirova@yandex.ru>
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
from typing import Union
from .core import *
from .element import *
from .query import *
from .dynamic import DynamicElement
from .dynamic_query import *


AnyAtom = Union[Element, DynamicElement, QueryElement, DynamicQueryElement, AnyElement, DynamicAnyElement, ListElement]


__all__ = ['Core', 'Element', 'DynamicElement', 'QueryElement', 'DynamicQueryElement',
           'AnyElement', 'DynamicAnyElement', 'AnyAtom', 'ListElement']
