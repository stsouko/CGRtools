#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2014-2020 Ramil Nugmanov <nougmanoff@protonmail.com>
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
from distutils.util import get_platform
from pathlib import Path
from setuptools import setup
from wheel.bdist_wheel import bdist_wheel


platform = get_platform()
if platform == 'win-amd64':
    libinchi = ['INCHI/libinchi.dll']
elif platform == 'linux-x86_64':
    libinchi = ['INCHI/libinchi.so']
else:
    libinchi = []


class _bdist_wheel(bdist_wheel):
    def finalize_options(self):
        super().finalize_options()
        self.root_is_pure = False


setup(
    name='CGRtools',
    version='4.0.18',
    packages=['CGRtools', 'CGRtools.algorithms', 'CGRtools.containers', 'CGRtools.files',
              'CGRtools.periodictable', 'CGRtools.utils', 'CGRtools.attributes'],
    url='https://github.com/cimm-kzn/CGRtools',
    license='LGPLv3',
    author='Dr. Ramil Nugmanov',
    author_email='nougmanoff@protonmail.com',
    python_requires='>=3.6.1',
    cmdclass={'bdist_wheel': _bdist_wheel},
    install_requires=['CachedMethods>=0.1.4,<0.2'],
    extras_require={'mrv': ['lxml>=4.1,<4.5'], 'clean2d': ['numpy>=1.18'],
                    'clean2djit': ['numpy>=1.18', 'numba>=0.46']},
    data_files=[('lib', libinchi)],
    zip_safe=False,
    long_description=(Path(__file__).parent / 'README.md').read_text(),
    classifiers=['Environment :: Plugins',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 3 :: Only',
                 'Programming Language :: Python :: 3.6',
                 'Programming Language :: Python :: 3.7',
                 'Topic :: Scientific/Engineering',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Scientific/Engineering :: Information Analysis',
                 'Topic :: Software Development',
                 'Topic :: Software Development :: Libraries',
                 'Topic :: Software Development :: Libraries :: Python Modules'],
    command_options={'build_sphinx': {'source_dir': ('setup.py', 'doc'),
                                      'build_dir':  ('setup.py', 'build/doc'),
                                      'all_files': ('setup.py', True)}}
)
