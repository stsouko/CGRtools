# -*- coding: utf-8 -*-
from os.path import abspath
from sys import path
parent = abspath('..')
if parent not in path:
    path.insert(0, parent)

author = 'Dr. Ramil Nugmanov'
copyright = '2014-2019, Dr. Ramil Nugmanov <stsouko@live.ru>'
version = '4.0'
project = 'CGRtools'

needs_sphinx = '1.6'
extensions = ['sphinx.ext.autodoc', 'sphinx.ext.autosummary', 'sphinx.ext.coverage']


exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
templates_path = ['_templates']
source_suffix = '.rst'
master_doc = 'index'

language = None
pygments_style = 'sphinx'
todo_include_todos = False
autoclass_content = 'both'

html_theme = 'alabaster'
html_theme_options = {'github_user': 'cimm-kzn', 'github_repo': 'CGRtools', 'show_related': True}
html_static_path = ['_static']
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
