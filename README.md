
# Overview
This repository contains python files and jupiter notebooks for the analysis of images from sliced and masked CT and MRI images from organs found in the human abdomen. 

* [Introduction](#introduction)
* [Installation](#installation)
* [Usage](#usage)
* [Contributing](#contributing)
* [License](#license)
* [Contact](#contact)
* [Acknowledgements](#acknowledgements)

# Introduction

## `3D Viewer`
A simple 3D object viewer using plotly. The input is a singlefolder containing a stack of `.tif` images. The output is a 3D object viewer that can be rotated and zoomed in and out. 

## `surrace_classifier`
The Surface Classifier tool analyzes 3D image stacks to quantify surface characteristics, providing surface statistics and exporting a 3D mesh (`.obj` file) for further insights or 3D printing.

## `2D_distance`
The 2D Distance tool calculates the distance between two points in a 2D image stack, providing insights the size of features or objects in the organ.

## `Location_finder_EIT`
The Location Finder tool uses a trained nn-Unet segmentation model to predict where in the abdomen your organ can be. 

## `volume_compare` 
Plots your .obj file from surface classifier on the same axis as common objects and against a to scale human body outline. Helpful for gauging relative size.

## `surface_analysis`
Computes texture features using the gray-level co-occurrence matrix (GLCM) method (extracting metrics like contrast, correlation, energy, and homogeneity), saves these results to a CSV file, and finally displays each image with its corresponding texture analysis results overlaid.

## `shape_classifier`
Analyzes organ geometry and intensity to compute centroid, min/max radius, area, and intensity variation. Also estimates the organ's total mass using pixel intensity, slice thickness, and uncertainty propagation, saving results to a CSV file for further analysis.

# Installation
Each tool is designed to be run in its own jupyter notebook. Run all cells in the notebook in order to use the tool.

# Usage

# Contributing
Feel free to submit a PR on GitHub

# License
Currently not distributed under any license.

# Contact
* @agadin on github

# Acknowledgements
