#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# Copyright 2014-2016 Ramil Nugmanov <stsouko@live.ru>
# This file is part of FEAR (Find Errors in Automapped Reactions).
#
# FEAR is free software; you can redistribute it and/or modify
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
import time
from itertools import chain
from CGRtools.CGRrw import CGRWrite


class RDFwrite(CGRWrite):
    def __init__(self, file, extralabels=False, flushmap=False):
        CGRWrite.__init__(self, extralabels=extralabels, flushmap=flushmap)
        self.__file = file
        self.writedata = self.__initwrite

    def close(self):
        self.__file.close()

    def __initwrite(self, data):
        self.__file.write(time.strftime("$RDFILE 1\n$DATM    %m/%d/%y %H:%M\n"))
        self.__writedata(data)
        self.writedata = self.__writedata

    def __writedata(self, data):
        self.__file.write('$RFMT\n$RXN\n\n  CGRtools. (c) Dr. Ramil I. Nugmanov\n\n%3d%3d\n' %
                          (len(data['substrats']), len(data['products'])))
        colors = {}
        for cnext, m in enumerate(data['substrats'] + data['products'], start=1):
            m = self.getformattedcgr(m)
            self.__file.write('$MOL\n')
            self.__file.write(m['CGR'])
            self.__file.write("M  END\n")
            colors.update({'%s.%d' % (k, cnext): v for k, v in m['colors'].items()})

        for p in chain(colors.items(), data['meta'].items()):
            self.__file.write('$DTYPE %s\n$DATUM %s\n' % p)
