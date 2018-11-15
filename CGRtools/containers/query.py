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
from .common import BaseContainer
from ..algorithms import CGRstring
from ..attributes import QueryAtom, QueryBond


class QueryContainer(BaseContainer):
    def add_stereo(self, *args, **kwargs):
        pass

    node_attr_dict_factory = QueryAtom
    edge_attr_dict_factory = QueryBond

    @property
    def _morgan_init(self):
        return

    def _prepare_stereo(self):
        return {}

    def _signature_generator(self, *args, **kwargs):
        return CGRstring(*args, is_cgr=True, **kwargs)
