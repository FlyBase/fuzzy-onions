# -- General configuration ------------------------------------------------

source_suffix = { '.md': 'markdown' }
master_doc = 'scrnaseq-sop'
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
pygments_style = 'sphinx'
extensions = ['myst_parser']

project = 'Curation of scRNAseq datasets'
copyright = u'2022 Damien Goutte-Gattat'
author = u'Damien Goutte-Gattat <dpg44@cam.ac.uk>'
language = 'en'


# -- Options for HTML output ----------------------------------------------

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
