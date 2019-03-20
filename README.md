CGRTools
========
Tools for processing of reactions based on Condensed Graph of Reaction (CGR) approach.

Basic opertions:
   - Read /write /convert formats MDL .RDF and .SDF, SMILES, .MRV
   - Standardize reactions and valid structures checker.
   - Produce CGRs.
   - Perfrom subgraph search.
   - Build /correct molecules and reactions.
   - Produce template based reactions.
    
INSTALL
=======

Linux Debian based
==================

* Install python3.7 and git

    ```
    sudo apt install python3.7 python3.7-dev git python3-virtualenv
    ```
    
* Install virtualenv.

    ```
    sudo apt install virtualenv
    ```

* Create new environment and activate it.

    ```
    virtualenv -p python3.7 venv
    source venv/bin/activate
    ```

Mac
===
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
=======

* Install python3.7 and git using [Chocolately](<https://chocolatey.org/>)

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
============
* Into activated environment install networkx library (Currently DEV version, due to special 
parts of code that was contributed to networkx and will appear in next release)

    ```
    pip install -U git+https://github.com/networkx/networkx.git@master#egg=networkx
    ```
    
* **stable version will be available through PyPI (The same as DEV for now)**

    ```
    pip install CGRTools
    ```    
    
* Install CGRtools library DEV version for features that are not well tested (Currently DEV version and stable version is the same as 
bugs fixing is going on)

    ```
    pip install -U git+https://github.com/cimm-kzn/CGRtools.git@master#egg=CGRtools
    ```

* Jupyter integration:

    ```
    pip install jupyter
    jupyter notebook
    ```
    
* Download tutorial files

   <https://github.com/cimm-kzn/CGRtools/tree/master/tutorial>

* Open .ipynb file in jupyter browser


**If you still have questions, please open issue within github.**

COPYRIGHT
=========

2014-2019 Ramil Nugmanov <stsouko@live.ru> main developer  
2014-2019 Timur Madzhidov <tmadzhidov@gmail.com> atom ordering algorithm and API discussion  
2014-2019 Alexandre Varnek <varnek@unistra.fr> base idea of CGR approach

CONTRIBUTORS
============

* Timur Gimadiev <timur.gimadiev@gmail.com>
* Ravil Mukhametgaleev <sonic-mc@mail.ru>
* Tagir Akhmetshin <tagirshin@gmail.com>
* Adelia Fatykhova <adelik21979@gmail.com>
