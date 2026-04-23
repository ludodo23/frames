## docs/conf.py — Sphinx configuration for frames

project   = "frames"
copyright = "2026, Ludovic Andrieux"
author    = "Ludovic Andrieux"
release   = "0.1.0"

extensions = [
    "breathe",
    "myst_parser",
    "sphinx.ext.autodoc",
    "sphinx.ext.viewcode",
    "sphinx_rtd_theme",
]

myst_enable_extensions = ["colon_fence", "deflist"]
source_suffix = {".rst": "restructuredtext", ".md": "markdown"}

breathe_projects        = {"frames": "../doxygen/xml"}
breathe_default_project = "frames"

html_theme = "sphinx_rtd_theme"
html_theme_options = {
    "navigation_depth": 4,
    "collapse_navigation": False,
    "titles_only": False,
}
html_static_path = ["_static"]

exclude_patterns = ["_build", "Thumbs.db", ".DS_Store"]