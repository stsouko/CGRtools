#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# preparerdf.py
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
import copy
import re
import subprocess as sp
from condenserpkg.version import version


def rings(connections, trace, inter):
    '''
    поиск всех атомов колец
    '''
    if len(trace) < 8:
        iterlist = set(connections[inter]).difference(trace[-2:])
        for i in iterlist:
            if i == trace[0]:
                # детект цикла
                return True
            if rings(connections, copy.copy(trace + [i]), i):
                return True
    return False


def main():
    rawopts = argparse.ArgumentParser(description="CGR generator", epilog="created by stsouko")
    rawopts.add_argument("--version", "-v", action="version", version=version(), default=False)
    rawopts.add_argument("--input", "-i", type=str, default="input.rdf", help="RDF inputfile")
    rawopts.add_argument("--output", "-o", type=str, default="output.rdf", help="RDF outputfile")
    options = vars(rawopts.parse_args())

    repdict = dict(TRANS='E', CIS='Z', UNKNOWN='U', CYCLE='C', SYMM='X', ODD='R', EVEN='S')

    with open(options['input']) as f, open(options['output'], 'w') as o:
        buff = []
        proc = sp.check_output(['cxcalc', 'stereoanalysis', options['input']])
        rescalc = (x for x in proc.split('\n\n'))

        for line in f:
            if "$RXN" in line[:4]:
                stereodata = [x for x in rescalc.next().split('\n') if x.split(' ')[0] in 'TETRAHEDRAL CISTRANS']
                keygen = 0
                o.write(line)
            elif "$MOL" in line[0:4]:
                buff.append('')
                o.write(line)
            elif "M  END" in line:
                res = []
                atcount = int(buff[4][:3])

                listed = []
                forscan = []

                for j in stereodata:
                    i = re.search('\[\d+(,\s\d+)?\]', j).group().strip('[]').split(', ')
                    if keygen <= int(i[0]) < atcount + keygen:
                        #print keygen, i
                        corrected = j.split(' ')
                        corrected[1] = str(sorted([int(x) - keygen for x in i]))
                        #TETRAHEDRAL [15] - [11, 14, 16] : UNKNOWN
                        res.append(' '.join(corrected))
                        if len(i) == 2:
                            listed.append(sorted([int(x) + 1 - keygen for x in i]))

                keygen += atcount
                #some magic on rings and symmetry with double bond
                a = atcount + 5
                b = a + int(buff[4][3:6])
                connections = {x: {} for x in range(1, a - 4)}
                for i in range(a, b):
                    x, y = int(buff[i][:3]), int(buff[i][3:6])
                    connections[x][y] = connections[y][x] = '*'
                    if buff[i][6:9] == '  2' and sorted([int(buff[i][:3]), int(buff[i][3:6])]) not in listed:
                        forscan.append((int(buff[i][:3]), int(buff[i][3:6])))

                for i, j in forscan:
                    if rings(connections, [i], i):
                        res.append('CISTRANS [%d, %d] : CYCLE' % tuple(sorted([i - 1, j - 1])))
                    else:
                        res.append('CISTRANS [%d, %d] : SYMM' % tuple(sorted([i - 1, j - 1])))

                if res:
                    for j in xrange(1000):
                        sty = res[j * 8:j * 8 + 8]
                        if sty:
                            stydat = ' '.join(['%3d DAT' % (x + 1 + j * 8) for x in range(len(sty))])
                            buff.append('M  STY  %d %s\n' % (len(sty), stydat))
                        else:
                            break
                    for i, j in enumerate(res):
                        j = re.search('\[\d+(,\s\d+)?\]', j).group().strip('[]').split(', ') + [
                            repdict.get(j.split(': ')[1], 'U'), 'stereo']

                        buff.append(
                            'M  SAL %3d%3d %s\n' % (i + 1, len(j) - 2, ' '.join(['%3d' % (int(x) + 1) for x in j[:-2]])))
                        buff.append('M  SDT %3d %s\n' % (i + 1, j[-1]))
                        buff.append('M  SED %3d %s\n' % (i + 1, j[-2]))
                buff.append(line)


                o.write(''.join(buff))
                #print connections, forscan
                buff = []
            elif buff:
                buff.append(line)
            else:
                o.write(line)


if __name__ == '__main__':
    main()