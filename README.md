### SHARP

This repository contains code and examples for the Spatial Heterogeneity Analysis by Recursive Partitioning (SHARP) algorithm. The algorithm recursively segments a tumor on 2D or 3D images and outputs the amount of heterogeneity at various distance scales.

Created by Michael Gensheimer and Andrew Trister at the University of Washington.
Direct questions to michael.gensheimer@gmail.com

### How to install

The code was tested in MATLAB R2012b in Linux. It requires the Image Processing Toolbox and Sameer Agarwal's free Spectral Clustering Toolbox, available here: http://homes.cs.washington.edu/~sagarwal/code.html
To install, simply add the SHARP and Spectral Clustering Toolbox directories to the MATLAB path.

### How to use

Run mri.m for a demonstration using one slice from MRI of a breast tumor.
Run synthetic.m for a demonstration using synthetic images.

