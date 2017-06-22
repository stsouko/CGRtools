#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2014-2017 Ramil Nugmanov <stsouko@live.ru>
#  This file is part of CGRtools.
#
#  CGRtools is free software; you can redistribute it and/or modify
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
from CGRtools.version import version
from pathlib import Path
from setuptools import setup, find_packages

setup(
    name='CGRtools',
    version=version(),
    packages=find_packages(),
    url='https://github.com/stsouko/CGRtools',
    license='AGPLv3',
    author='Dr. Ramil Nugmanov',
    author_email='stsouko@live.ru',
    description='CGR tools',
    entry_points={'console_scripts': ['cgrtools=CGRtools.CLI:launcher']},
    package_data={'CGRtools.utils': ['aromatize.rdf']},
    install_requires=['networkx>=2.0.dev', 'periodictable'],
    extras_require={'autocomplete': ['argcomplete']},
    dependency_links=['git+https://github.com/networkx/networkx.git@master#egg=networkx-2.0.dev'],
    long_description=(Path(__file__).parent / 'README.md').open().read(),
    keywords="tools cgr cli",
    classifiers=['Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'Intended Audience :: Developers',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Software Development :: Libraries :: Python Modules',
                 'License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 3',
                 'Programming Language :: Python :: 3.4',
                 'Programming Language :: Python :: 3.5',
                 ]
)
