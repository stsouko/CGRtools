CGRTools
========
Tools for processing of reactions based on Condensed Graph of Reaction (CGR) approach.

Basic operations:
   - Read/write/convert formats: MDL .RDF (RXN) and .SDF (MOL), .MRV, SMILES, INCHI (Linux and Windows), .XYZ, .PDB
   - Standardize molecules and reactions and valid structures checker.
   - Duplicate searching.
   - Tetrahedron, Allene and CIS-TRANS stereo checking.
   - Produce CGRs.
   - Perform subgraph search.
   - Build/edit molecules and reactions.
   - Produce template based reactions and molecules.
   - Atom-to-atom mapping checker and rule-based fixer.
   - Perform MCS search.
   - 2d depiction.

INSTALL
=======

Highly recommended to use python 3.8+. Python 3.6 and 3.7 deprecated.


Linux Debian based
------------------
* Install python3.8, virtualenv and git::

    sudo apt install python3.8 python3.8-dev git python3-virtualenv
    
* Create new environment and activate it::

    virtualenv -p python3.8 venv
    source venv/bin/activate

Mac
---
* Install python3.8 and git using [brew](<https://brew.sh>)::

    brew install git
    brew install python3

* Install virtualenv::

    pip install virtualenv

* Create new environment and activate it::

    virtualenv -p python3.8 venv
    source venv/bin/activate
    
Windows
-------
* Install python3.8 and git using [Chocolatey](<https://chocolatey.org/>)::

    choco install git
    choco install python3
    
* Install virtualenv::

    pip install virtualenv

* Create new environment and activate it::

    virtualenv venv
    venv\Scripts\activate

General part
------------

* **stable version available through PyPI**::

    pip install CGRTools

* Install CGRtools with MRV files parsing support::

    pip install CGRTools[mrv]

* Install CGRtools with structures `clean2d` support (optimized version) \[numba and numpy used\]::

    pip install CGRtools[clean2djit]

* Install CGRtools with structures `clean2d` support slow version \[without numba\]::

    pip install CGRtools[clean2d]

* Install CGRtools library DEV version for features that are not well tested::

    pip install -U git+https://github.com/cimm-kzn/CGRtools.git@master#egg=CGRtools


**If you still have questions, please open issue within github.**

PACKAGING
=========

For wheel generation just type next command in source root::

    python setup.py bdist_wheel

On Linux additionally do repairing of package::

    pip install auditwheel
    auditwheel repair dist/CGRtools-<version>-<python_version>-linux_x86_64.whl

COPYRIGHT
=========

* 2014-2020 Ramil Nugmanov <nougmanoff@protonmail.com> main developer
* 2014-2019 Timur Madzhidov <tmadzhidov@gmail.com> features and API discussion
* 2014-2019 Alexandre Varnek <varnek@unistra.fr> base idea of CGR approach

CONTRIBUTORS
============

* Dinar Batyrshin <batyrshin-dinar@mail.ru>
* Timur Gimadiev <timur.gimadiev@gmail.com>
* Adelia Fatykhova <adelik21979@gmail.com>
* Tagir Akhmetshin <tagirshin@gmail.com>
* Ravil Mukhametgaleev <sonic-mc@mail.ru>

CITE THIS
=========

CGRtools: Python Library for Molecule, Reaction, and Condensed Graph of Reaction Processing.
Journal of Chemical Information and Modeling 2019 59 (6), 2516-2521.
DOI: 10.1021/acs.jcim.9b00102 
