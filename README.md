# SigmaTransform-Matlab
Implements the Sigma Transform in Matlab/Octave.

## Contents
- [General Information](#general-information)
- [Installation](#installation)
- [Usage](#usage)

# General information
This repository shows an exemplary implementation of the "SigmaTransform", as defined in the thesis _"Quantum Frames and Uncertainty Principles arising from Symplectomorphisms"_, in MATLAB/Octave.

Note that this library is not intended to show maximal performance, but shows the usability of the universal interface of the "Sigma Transform" to perform well-known signal processing transforms / algorithms like the Short-Time Fourier Transform, Wavelet Transform,
Curvelet Transform, Shearlet Transform, etc., differing only by single paramater - the "spectral diffeomorphism".

# Installation
The code, written and tested on _MATLAB R2016a_, with the Signal Processing toolbox installed, should be self-contained and, thus, cloning the repository or copying _all_ the files, found in this repository, to a local folder and changing the Matlab-path to this folder should be enough to use and test the code provided.

The code should also run on GNU Octave, tested with _GNU Octave 4.4.1_, with the _signal and image packages_ installed and _loaded_.

# Usage
The two examples

    Examples_SigmaTransform1D.m
    Examples_SigmaTransform2D.m

show how to use the implementation.
