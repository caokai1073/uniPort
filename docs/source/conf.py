# Configuration file for the Sphinx documentation builder.

# -- Project information

import inspect
import uniport

import sphinx_autodoc_typehints

project = 'uniPort'
copyright = u'2022, Kai Cao'
author = 'Kai Cao'

release = uniport.__version__
version = uniport.__version__

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',
    'sphinx_markdown_tables',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx_copybutton',
    'nbsphinx',
    'myst_parser',
]

# Generate the API documentation when building
autosummary_generate = True
autodoc_member_order = 'bysource'

napoleon_google_docstring = False
napoleon_numpy_docstring = True
napoleon_include_init_with_doc = False
napoleon_use_rtype = True  # having a separate entry generally helps readability
napoleon_use_param = True
napoleon_custom_sections = [('Params', 'Parameters')]
todo_include_todos = False


intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']
html_static_path = ['_static']
html_logo = '_static/uniPort.png'

# -- Options for HTML output

html_theme = 'sphinx_rtd_theme'

# -- Options for EPUB output
epub_show_urls = 'footnote'

# latex_engine = 'xelatex'
# latex_use_xindy = False
# latex_elements = {
#     'preamble': '\\usepackage[UTF8]{ctex}\n',
# }

html_theme_options = dict(navigation_depth=4, logo_only=True)  # Only show the logo

from recommonmark.parser import CommonMarkParser
source_parsers = {
    '.md': CommonMarkParser,
}
source_suffix = ['.rst', '.md']