CGRTools
========
Tools for processing of reactions based on Condensed Graph of Reaction (CGR) approach.

Basic operations:
   - Read/write/convert formats: MDL .RDF and .SDF, .MRV, SMILES, INCHI (Linux and Windows), .XYZ, .PDB
   - Standardize molecules and reactions and valid structures checker.
   - Produce CGRs.
   - Perform subgraph search.
   - Build/edit molecules and reactions.
   - Produce template based reactions and molecules.
   - Atom-to-atom mapping checker and rule-based fixer.
   - Perform MCS search.
   - 2d depiction.

INSTALL
=======

Linux Debian based
------------------

* Install python3.7, virtualenv and git

    ```
    sudo apt install python3.7 python3.7-dev git python3-virtualenv
    ```
    
* Create new environment and activate it.

    ```
    virtualenv -p python3.7 venv
    source venv/bin/activate
    ```

Mac
---
* Install python3.7 and git using [brew](<https://brew.sh>)

    ```
    brew install git
    brew install python3
    ```
    
* Install virtualenv.

    ```
    pip install virtualenv
    ```

* Create new environment and activate it.

    ```
    virtualenv -p python3.7 venv
    source venv/bin/activate
    ```
    
Windows
-------

* Install python3.7 and git using [Chocolatey](<https://chocolatey.org/>)

    ```
    choco install git
    choco install python3
    ```
    
* Install virtualenv.

    ```
    pip install virtualenv
    ```

* Create new environment and activate it.

    ```
    virtualenv venv
    venv\Scripts\activate
    ```

General part
------------

* **stable version will be available through PyPI**

    ```
    pip install CGRTools
    ```    
    
* Install CGRtools library DEV version for features that are not well tested

    ```
    pip install -U git+https://github.com/cimm-kzn/CGRtools.git@master#egg=CGRtools
    ```

* Install CGRtools with MRV files parsing support

    ```
    pip install CGRTools[mrv]
    ```

* Install CGRtools with structures `clean2d` support \[only DEV\] optimized version \[numba and numpy used\]

    ```
    pip install -U git+https://github.com/cimm-kzn/CGRtools.git@master#egg=CGRTools[clean2djit]
    ```

* Install CGRtools with structures `clean2d` support \[only DEV\] slow version \[numpy used\]

    ```
    pip install -U git+https://github.com/cimm-kzn/CGRtools.git@master#egg=CGRtools[clean2d]
    ```

* Jupyter integration:

    ```
    pip install jupyter
    jupyter notebook
    ```
    
* Download tutorial files to same directory

   <https://github.com/cimm-kzn/CGRtools/tree/master/tutorial>

* Open .ipynb files from tutorial directory in Jupyter browser

**If you still have questions, please open issue within github.**

PACKAGING
=========

For wheel generation just type next command in source root

    python setup.py bdist_wheel

On Linux additionally do repairing of package

    pip install auditwheel
    auditwheel repair dist/CGRtools-<version>-<python_version>-linux_x86_64.whl

COPYRIGHT
=========

2014-2020 Ramil Nugmanov <nougmanoff@protonmail.com> main developer  
2014-2019 Timur Madzhidov <tmadzhidov@gmail.com> features and API discussion  
2014-2019 Alexandre Varnek <varnek@unistra.fr> base idea of CGR approach

CONTRIBUTORS
============

* Dinar Batyrshin <batyrshin-dinar@mail.ru>
* Timur Gimadiev <timur.gimadiev@gmail.com>
* Adelia Fatykhova <adelik21979@gmail.com>
* Tagir Akhmetshin <tagirshin@gmail.com>
* Ravil Mukhametgaleev <sonic-mc@mail.ru>
