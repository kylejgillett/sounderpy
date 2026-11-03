# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'SounderPy'
copyright = '2026, Kyle J Gillett'
author = 'Kyle J Gillett'
release = 'v3.3.0'
html_logo = 'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/images/sounderpy_logo.png'
html_favicon = 'https://raw.githubusercontent.com/kylejgillett/sounderpy/main/images/sounderpy_logo.png'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = ['nbsphinx', 'sphinx_reredirects']

exclude_patterns = []

nbsphinx_execute = "never"

pygments_style = "zenburn"

# Long function signatures list one parameter per line.
maximum_signature_line_length = 70

# Pages moved into section folders keep their old URLs working.
redirects = {
    "getting_started": "user_guide/getting_started.html",
    "about": "user_guide/about.html",
    "troubleshooting": "user_guide/troubleshooting.html",
    "community": "user_guide/community.html",
    "additional_resources": "user_guide/additional_resources.html",
    "gettingdata": "reference/gettingdata.html",
    "customdatasources": "reference/customdatasources.html",
    "clean_data_scheme": "reference/clean_data_scheme.html",
    "plottingdata": "reference/plottingdata.html",
    "helpertools": "reference/helpertools.html",
    "cli": "reference/cli.html",
    "examplegallery": "gallery/examplegallery.html",
    "examplescripts": "gallery/examplescripts.html",
    "beginners_sounderpy_cookbook": "tutorials/beginners_sounderpy_cookbook.html",
    "acars_data_example": "notebooks/acars_data_example.html",
    "bufkit_data_example": "notebooks/bufkit_data_example.html",
    "composite_example": "notebooks/composite_example.html",
    "custom_data_example": "notebooks/custom_data_example.html",
    "hodograph_example": "notebooks/hodograph_example.html",
    "raob_data_example": "notebooks/raob_data_example.html",
    "reanalysis_data_example": "notebooks/reanalysis_data_example.html",
    "sounding_example": "notebooks/sounding_example.html",
    "hodograph-example": "notebooks/hodograph_example.html",
    "notebook_gallery": "gallery/examplescripts.html",
    "index_old": "index.html",
}

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']
html_theme_options = {
    'analytics_anonymize_ip': False,
    'logo_only': False,
    'prev_next_buttons_location': 'bottom',
    'style_external_links': False,
    'vcs_pageview_mode': '',
    # Toc options
    'collapse_navigation': True,
    'sticky_navigation': True,
    'navigation_depth': 4,
    'includehidden': True,
    'titles_only': False,
}

html_css_files = [
    "css/homepage.css",
    "css/gallery.css",
    "css/case_studies.css",
    "css/sounderpy-theme.css",
]

html_js_files = ["js/param-layout.js"]

html_extra_path = ["_extra"]