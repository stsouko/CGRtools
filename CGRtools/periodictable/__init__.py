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
"""
periodic table of elements classes
"""
from .element import Element
from .groups import *
from .periods import *
from .I import *
from .II import *
from .III import *
from .IV import *
from .V import *
from .VI import *
from .VIII import *
from .IX import *
from .X import *
from .XI import *
from .XII import *
from .XIII import *
from .XIV import *
from .XV import *
from .XVI import *
from .XVII import *
from .XVIII import *


__all__ = ['Element']  # ordered wildcard
__all__.extend(k for k in locals() if k.startswith('Group'))
__all__.extend(k for k in locals() if k.startswith('Period'))
__all__.extend(k for k, v in locals().items() if isinstance(v, Element))
