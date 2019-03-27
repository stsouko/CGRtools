# -*- coding: utf-8 -*-
#
#  Copyright 2018, 2019 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU Lesser General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, see <https://www.gnu.org/licenses/>.
#
from abc import abstractmethod
from itertools import islice


class Isomorphism:
    def __lt__(self, other):
        if len(self) >= len(other):
            return False
        return self.is_substructure(other)

    def __le__(self, other):
        return self.is_equal(other) if len(self) == len(other) else self.is_substructure(other)

    def __gt__(self, other):
        if len(self) <= len(other):
            return False
        return other.is_substructure(self)

    def __ge__(self, other):
        return other.is_equal(self) if len(self) == len(other) else other.is_substructure(self)

    @abstractmethod
    def _matcher(self, other):
        pass

    def is_substructure(self, other):
        """
        test self is substructure of other
        """
        return self._matcher(other).subgraph_is_isomorphic()

    def is_equal(self, other):
        """
        test self is structure of other
        """
        return self._matcher(other).is_isomorphic()

    def get_mapping(self, other):
        """
        get self to other mapping
        """
        m = next(self._matcher(other).isomorphisms_iter(), None)
        if m:
            return {v: k for k, v in m.items()}

    def get_substructure_mapping(self, other, limit=1):
        """
        get self to other substructure mapping

        :param limit: number of matches. if 0 return iterator for all possible; if 1 return dict or None;
            if > 1 return list of dicts
        """
        i = self._matcher(other).subgraph_isomorphisms_iter()
        if limit == 1:
            m = next(i, None)
            if m:
                return {v: k for k, v in m.items()}
            return
        elif limit == 0:
            return ({v: k for k, v in m.items()} for m in i)
        return [{v: k for k, v in m.items()} for m in islice(i, limit)]


__all__ = ['Isomorphism']
