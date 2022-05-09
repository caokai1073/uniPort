# Configuration file for the Sphinx documentation builder.

# -- Project information

import inspect

import sphinx_autodoc_typehints

project = 'uniPort'
copyright = u'2022, Kai Cao'
author = 'Kai Cao'

release = '1.0.5'
version = '1.0.5'

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


intersphinx_mapping = {
    'python': ('https://docs.python.org/3/', None),
    'sphinx': ('https://www.sphinx-doc.org/en/master/', None),
}
intersphinx_disabled_domains = ['std']

templates_path = ['_templates']

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