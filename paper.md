---
title: "SounderPy: An atmospheric sounding visualization and analysis tool for Python"
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

SounderPy is a simple, open-source Python package for retrieving, processing, & 
plotting atmospheric vertical profile (sounding) data. Built for simplicity and
reliability for all users and use cases, this project’s goal is to provide a 
uniform method for sounding analysis across all data types. Severe weather analysis
and forecasting requires a sound comprehension of thermodynamic and kinematic properties of the 
environment. SounderPy makes this possible with robust access to data and custom 
visualizations. The tool creates complex yet effective sounding and hodograph 
plots with high readability which are designed specifically for severe weather 
analysis and forecasting. All of this functionality can be completed in three simple lines of code or less, 
making SounderPy an accessible tool for both Python experts and novices. A number
of scientific Python libraries build the base of SounderPy’s efficient and 
durable functionality, such as NumPy, Matplotlib, xarray, Metpy, and SHARPpy. 
SounderPy is available through GitHub and PyPi and is distributed under an 
MIT license.

# Statement of need

Meteorological data from varying sources may be stored in a variety of file types 
with a wide range of data structures. Such diversity in data formats & availability 
can make consistent and effectient processing of complex meteorological data difficult.
Additionally, thorough yet effectient meteorological analysis of atmospheric properties is 
vital in describing the past, current, and future state of the atmosphere. Such statements 
have been made since the dawn of "free-air" observations (soundings) by kites, balloons and 
aircraft in the early 20th century [@Byers:1934]. Consistent calculations and displays of 
meteorological data normalizes data analysis such that meaningful comparisons can be 
drawn from different data types and sources. Reliable statisics and analogs can be developed
from normalized data analysis which can aid forecasters and researchers with pattern recognition 
and context. An example of such would be the comparison of numerical weather predicition output 
to observations of past events. A forecaster may determine the similarities between 
past and future events and factor those similarities into their forecast.

SounderPy allows for simple access to multiple data sources, such as National Weather Service 
observed radiosonde observations, Aircraft Communications Addressing and Reporting System 
observations, numerical weather prediction forecast data, and numerical weather prediction 
reanalysis data. Each data type has it's own source, file type, and data structure. 


# Basics of the API





