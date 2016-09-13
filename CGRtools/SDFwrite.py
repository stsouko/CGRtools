#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
#  Copyright 2014-2016 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of cgrtools.
#
#  cgrtools is free software; you can redistribute it and/or modify
#  it under the terms of the GNU Affero General Public License as published by
#  the Free Software Foundation; either version 3 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU Affero General Public License for more details.
#
#  You should have received a copy of the GNU Affero General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
#  MA 02110-1301, USA.
#
from itertools import chain
from CGRtools.CGRrw import CGRWrite


class SDFwrite(CGRWrite):
    def __init__(self, output, extralabels=False, flushmap=False):
        CGRWrite.__init__(self, extralabels=extralabels, flushmap=flushmap)
        self.__file = output

    def close(self):
        self.__file.close()

    def writedata(self, data):
        m = self.getformattedcgr(data)
        self.__file.write(m['CGR'])
        self.__file.write("M  END\n")

        for i in chain(m['colors'].items(), m['meta'].items()):
            self.__file.write(">  <%s>\n%s\n" % i)
        self.__file.write("$$$$\n")
