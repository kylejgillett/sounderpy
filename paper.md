---
title: "SounderPy: A sounding visualization tool for severe-weather analysis and forecasting"
tags:
  - Python
  - meteorology
  - atmospheric-sciences
  - data-analysis
  - weather-data
authors:
  - name: Kyle J. Gillett
    orcid: 0009-0001-5528-5037
    affiliation: 1
    ror: 04a5szx83
affiliations:
  - name: College of Aerospace Sciences, University of North Dakota, United States
    index: 1
date: 9 December 2024
bibliography: paper.bib
---


# Summary

SounderPy is a simple, open-source Python package for retrieving and plotting 
vertical profile (sounding) data. Built for simplicity and reliability for all 
uses and users, this project’s goal is to provide a uniform method for sounding 
analysis across multiple data types. Severe weather analysis and forecasting 
requires a sound comprehension of thermodynamic and kinematic properties of the 
environment. SounderPy makes this possible with robust access to data and custom 
visualizations. The tool creates complex yet effective sounding and hodograph 
plots with high readability which are designed specifically for severe weather 
analysis and forecasting. SounderPy is capable of retrieving and plotting model 
forecast data, observed radiosonde data, Aircraft Communications Addressing and 
Reporting System (ACARS) observation data, and model reanalysis data. All of 
this functionality can be completed in three simple lines of code or less, 
making SounderPy an accessible tool for both Python experts and novices. A number
of scientific Python libraries build the base of SounderPy’s efficient and 
durable functionality, such as NumPy, Matplotlib, xarray, Metpy, and SHARPpy. 
SounderPy is available through GitHub and PyPi and is distributed under an 
MIT license.

# Statement of need

`Gala` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for `Gala` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. `Gala` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the `Astropy` package [@astropy] (`astropy.units` and
`astropy.coordinates`).

`Gala` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in `Gala` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$

You can also use plain \LaTeX for equations
\begin{equation}\label{eq:fourier}
\hat f(\omega) = \int_{-\infty}^{\infty} f(x) e^{i\omega x} dx
\end{equation}
and refer to \autoref{eq:fourier} from text.

# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

If you want to cite a software repository URL (e.g. something on GitHub without a preferred
citation) then you can do it with the example BibTeX entry below for @fidgit.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this:
![Caption for example figure.\label{fig:example}](figure.png)
and referenced from text using \autoref{fig:example}.

Figure sizes can be customized by adding an optional second parameter:
![Caption for example figure.](figure.png){ width=20% }

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
