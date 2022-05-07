# Configuration file for the Sphinx documentation builder.

# -- Project information

import inspect

import sphinx_autodoc_typehints

project = 'uniPort'
author = 'Kai Cao'

release = '1.0'
version = '1.0.0'

# -- General configuration

extensions = [
    'sphinx.ext.duration',
    'sphinx.ext.doctest',
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.intersphinx',
    'sphinx.ext.viewcode',
    'recommonmark',
    'sphinx_markdown_tables',
    'sphinx.ext.coverage',
    'sphinx.ext.mathjax',
    'sphinx.ext.napoleon',
    'sphinx_copybutton',
    'nbsphinx'
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

latex_engine = 'xelatex'
latex_use_xindy = False
latex_elements = {
    'preamble': '\\usepackage[UTF8]{ctex}\n',
}

html_theme_options = dict(navigation_depth=4, logo_only=True)  # Only show the logo
html_context = dict(
    display_github=True,  # Integrate GitHub
    github_user='jsxlei',  # Username
    github_repo='SCALEX',  # Repo name
    github_version='main',  # Version
    conf_py_path='/docs/',  # Path in the checkout to the docs root
)
html_static_path = ['_static']
html_show_sphinx = False

fa_orig = sphinx_autodoc_typehints.format_annotation
def format_annotation(annotation, fully_qualified=True):  # pylint: disable=unused-argument
    r"""
    Adapted from https://github.com/agronholm/sphinx-autodoc-typehints/issues/38#issuecomment-448517805
    """
    if inspect.isclass(annotation):
        full_name = f'{annotation.__module__}.{annotation.__qualname__}'
        override = qualname_overrides.get(full_name)
        if override is not None:
            return f':py:class:`~{override}`'
    return fa_orig(annotation)
sphinx_autodoc_typehints.format_annotation = format_annotation