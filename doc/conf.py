# -*- coding: utf-8 -*-
from os.path import abspath
from sys import path
parent = abspath('..')
if parent not in path:
    path.insert(0, parent)
from CGRtools.periodictable import C, DynamicC, QueryC, DynamicQueryC, ListElement, AnyElement, DynamicAnyElement

author = 'Dr. Ramil Nugmanov'
copyright = '2014-2020, Dr. Ramil Nugmanov <nougmanoff@protonmail.com>'
version = '4.0'
project = 'CGRtools'

needs_sphinx = '1.8'
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.autosummary', 'nbsphinx']

exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store', '**.ipynb_checkpoints']
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'

language = None
pygments_style = 'sphinx'
todo_include_todos = False
autoclass_content = 'both'

html_theme = 'sphinx_rtd_theme'
html_logo = 'logo.jpg'
html_favicon = 'favicon.ico'
html_theme_options = {'github_user': 'stsouko', 'github_repo': 'CGRtools', 'show_related': True}
html_show_copyright = True
html_show_sourcelink = False
html_sidebars = {
    '**': [
        'about.html',
        'navigation.html',
        'relations.html',  # needs 'show_related': True theme option to display
        'searchbox.html',
    ]
}

nbsphinx_execute = 'never'
