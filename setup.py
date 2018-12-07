#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#
#  Copyright 2014-2018 Ramil Nugmanov <stsouko@live.ru>
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
from CGRtools.version import version
from pathlib import Path
from setuptools import setup, find_packages


setup(
    name='CGRtools',
    version=version(),
    packages=find_packages(),
    url='https://github.com/stsouko/CGRtools',
    license='LGPLv3',
    author='Dr. Ramil Nugmanov',
    author_email='stsouko@live.ru',
    description='CGR tools',
    entry_points={'console_scripts': ['cgrtools=CGRtools.CLI:launcher']},
    install_requires=['networkx>=2.3,<2.4', 'lxml>=4.1.1,<4.3', 'coho>=0.3,<0.4'],
    extras_require={'autocomplete': ['argcomplete'], 'sphinx': ['sphinx>=1.6']},
    dependency_links=['git+https://github.com/networkx/networkx.git@master#egg=networkx-2.3'],
    zip_safe=True,
    long_description=(Path(__file__).parent / 'README.md').open().read(),
    classifiers=['Environment :: Console',
                 'Intended Audience :: Science/Research',
                 'License :: OSI Approved :: GNU Lesser General Public License v3 or later (LGPLv3+)',
                 'Operating System :: OS Independent',
                 'Programming Language :: Python',
                 'Programming Language :: Python :: 3 :: Only',
                 'Programming Language :: Python :: 3.7',
                 'Topic :: Scientific/Engineering',
                 'Topic :: Scientific/Engineering :: Chemistry',
                 'Topic :: Scientific/Engineering :: Information Analysis',
                 'Topic :: Software Development',
                 'Topic :: Software Development :: Libraries',
                 'Topic :: Software Development :: Libraries :: Python Modules'],
    command_options={'build_sphinx': {'project': ('setup.py', 'CGRtools'),
                                      'version': ('setup.py', version()), 'source_dir': ('setup.py', 'doc'),
                                      'build_dir':  ('setup.py', 'build/doc'),
                                      'all_files': ('setup.py', True),
                                      'copyright': ('setup.py', 'Dr. Ramil Nugmanov <stsouko@live.ru>')},
                     'easy_install': {'allow_hosts': ('setup.py', 'github.com, pypi.python.org')}},
)
