#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# func.py
#
# Copyright 2014 Ramil Nugmanov <stsouko@live.ru>
# This file is part of condenser.
#
#  condenser is free software; you can redistribute it and/or modify
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
import argparse
import re
import subprocess as sp
from condenserpkg.version import version


def main():
    rawopts = argparse.ArgumentParser(description="CGR generator", epilog="created by stsouko")
    rawopts.add_argument("--version", "-v", action="version", version=version(), default=False)
    rawopts.add_argument("--input", "-i", type=str, default="input.rdf", help="RDF inputfile")
    rawopts.add_argument("--output", "-o", type=str, default="output.rdf", help="RDF outputfile")
    options = vars(rawopts.parse_args())

    repdict = dict(TRANS='E', CIS='Z', UNKNOWN='U')

    with open(options['input']) as f, open(options['output'], 'w') as o:
        buff = []
        for line in f:
            if "$MOL" in line[0:4]:
                buff.append('')
                o.write(line)
            elif "M  END" in line:
                buff.append(line)
                try:
                    t = open('tmp.mol', 'w')
                    t.write(''.join(buff))
                except IOError:
                    print('error write to tmp file')
                    return 0
                finally:
                    t.close()
                    proc = sp.check_output(['cxcalc', 'stereoanalysis', 'tmp.mol'])
                    res = [x for x in proc.split('\n') if x and x.split(' ')[0] in 'TETRAHEDRAL CISTRANS']

                    if res:
                        buff.pop(-1)
                        for j in xrange(1000):
                            sty = res[j * 8:j * 8 + 8]
                            if sty:
                                stydat = ' '.join(['%3d DAT' % (x + 1) for x in range(len(sty))])
                                buff.append('M  STY  %d %s\n' % (len(sty), stydat))
                            else:
                                break
                        for i, j in enumerate(res):
                            j = re.search('\[\d+(,\s\d+)?\]', j).group().strip('[]').split(', ') + [
                                repdict.get(j.split(': ')[1], 'U'), 'stereo']
                        #print j, j[:-2]
                            buff.append(
                                'M  SAL %3d%3d %s\n' % (i + 1, len(j) - 2, ' '.join(['%3d' % (int(x) + 1) for x in j[:-2]])))
                            buff.append('M  SDT %3d %s\n' % (i + 1, j[-1]))
                            buff.append('M  SED %3d %s\n' % (i + 1, j[-2]))
                        buff.append('M  END\n')
                o.write(''.join(buff))


                buff = []
            elif buff:
                buff.append(line)
            else:
                o.write(line)


if __name__ == '__main__':
    main()