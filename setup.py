#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2014-2019 Ramil Nugmanov <stsouko@live.ru>
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
from pathlib import Path
from setuptools import setup


setup(
    name='CGRtools',
    version='3.1.6',
    packages=['CGRtools', 'CGRtools.algorithms', 'CGRtools.attributes', 'CGRtools.containers', 'CGRtools.files',
              'CGRtools.files.dll', 'CGRtools.periodictable', 'CGRtools.utils'],
    url='https://github.com/cimm-kzn/CGRtools',
    license='LGPLv3',
    author='Dr. Ramil Nugmanov',
    author_email='stsouko@live.ru',
    python_requires='>=3.6.1',
    install_requires=['networkx>=2.3,<2.4'],
    extras_require={'smiles': ['coho>=0.4,<0.5'], 'mrv': ['lxml>=4.1,<4.4']},
    package_data={'CGRtools.files.dll': ['LICENCE', 'readme.txt', 'libinchi.so', 'libinchi.dll']},
    zip_safe=False,
    long_description=(Path(__file__).parent / 'README.md').open().read(),
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
