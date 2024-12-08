# Overview
This package is an extension of the module geometricVOF of OpenFOAM.
It presents a new interface reconstruction scheme called plicRDFFast faster than both the existing isoAlpha and plicRDF. The scheme can be used with the interIsoFoam solver employing the isoAdvector advection algorithm.

# Version
A pre-installed OpenFOAM version is required. The implementation and testing of this extension was done with the v2312 OpenFOAM version.

# Installation

``wmake -j`` in ``openfoam-plicrdffast/`` compiles the library into $(FOAM_LIBBIN)/libgeometricVoF.so. Careful! This overrides the originally installed libgeometricVOF.so.

After cloning the repository, one can use the new PLIC scheme by setting the parameter "reconstructionScheme" to "plicRDFFast" in the dictionary of the VOF variable (alpha.\*) in "solvers" in the project file "system/fvSolution". 

# Benchmark cases
This extension was tested on various benchmark cases like the quasi-2D disc in constant and reversed vortex flow, 3D sphere in constant and reversed vortex flow, as well as both a 2D and 3D Dam Break test case.

# License
This work is released as an extension of OpenFOAM: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version. See the file COPYING in this directory or http://www.gnu.org/licenses/, for a description of the GNU General Public License terms under which you may redistribute files.
