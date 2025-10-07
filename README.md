# Cubical Ridge Surface Finder

This repositrory contains the source code used to calculate a cubibal ridge surface, described in the paper [A Local Iterative Approach for the Extraction of 2D Manifolds from Strongly Curved and Folded Thin-Layer Structures](https://arxiv.org/abs/2308.07070). A ridge surface can be seen as a generalized version of an isosurface.Given a scalar field, each point inside the ridge surface is a local maxima point in the direction of the ridge surface. Cubical means that the surface is voxelized i.e. the surface only consists of faces of the voxels of the scalar field (see images).

## Installation

To compile, one should use cmake with features of c++23 or newer.

## Folder and File Structure

- src: contains the source code library
- includes: extern libraries used in this project and included with the project for ease of installation and use
- main.cxx: entry point for the terminal program
- test.cpp: entry point for the test program, created with catch2

## License

This project is licensed under a Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License.
You should have received a copy of the license along with this work.  If not, see [CC BY-NC-SA](http://creativecommons.org/licenses/by-nc-sa/3.0/).