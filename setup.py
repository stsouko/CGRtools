#!/usr/bin/env python3
from CGRtools.version import version
from setuptools import setup, find_packages


setup(name='cgrtools',
      version=version(),
      packages=find_packages(),
      url='https://github.com/stsouko/CGRtools',
      license='AGPLv3',
      author='Dr. Ramil Nugmanov',
      author_email='stsouko@live.ru',
      description='CGR tools',
      scripts=['cgrtools.py'],
      requires=['networkx', 'periodictable'],
      long_description='CGR tools distributive',

      keywords="tools cgr cli",
      classifiers=['Environment :: Console',
                   'Intended Audience :: End Users/Desktop',
                   'Intended Audience :: Developers',
                   ('License :: OSI Approved :: GNU Affero General Public License'
                    ' v3 or later (AGPLv3+)'),
                   'Operating System :: OS Independent',
                   'Programming Language :: Python',
                   'Programming Language :: Python :: 3',
                   'Programming Language :: Python :: 3.2',
                   'Programming Language :: Python :: 3.3',
                   'Programming Language :: Python :: 3.4',
                   'Programming Language :: Python :: 3.5',
                   ]
      )
