# Exact Voxelisation of SPH Particles

This repository contains a C++ implementation of several methods for voxelising smoothed particle hydrodynamics (SPH) particle distributions onto regular Cartesian grids.

The main method implemented here is a finite-volume voxelisation scheme based on precomputed cumulative smoothing kernels. The method accompanies the manuscript

> *Mass-conservative voxelisation of smoothed particles using cumulative kernels* (in preparation).

## Overview

An SPH particle with position \((x,y,z)\), smoothing length \(h\), and mass \(m\) represents a continuous density distribution rather than a point mass. The code provides several methods for converting a collection of such particles into a three-dimensional Cartesian density grid:

- nearest-grid-point deposition (`ngp`)
- sampling at voxel centres (`centre`)
- regular subsampling with \(3^3\) or \(5^3\) samples per voxel (`subsample3`, `subsample5`)
- Monte Carlo sampling with 10, 100, or 1000 samples per voxel (`montecarlo10`, `montecarlo100`, `montecarlo1000`)
- cumulative-kernel finite-volume voxelisation (`cumkernel`)

The cumulative-kernel method computes the mass contained inside each voxel using a precomputed three-dimensional cumulative representation of the adopted smoothing kernel.

The implementation currently includes three smoothing kernels:

- cubic spline
- quartic spline
- scaled Gaussian

The smoothing kernels are spherically symmetric, have finite support for normalized radius \(u=r/h \leq 1\), and are normalized to unit total mass.

## Directory structure

The code assumes the following directory structure:

```text
.
├── Code/
│   ├── VoxeliseSimulation.cpp
│   ├── CalculateCumulativeKernel.cpp
│   ├── Voxelise.cpp
│   ├── Voxelise.hpp
│   ├── SmoothingKernel.cpp
│   ├── SmoothingKernel.hpp
│   ├── CubicSplineSmoothingKernel.cpp
│   ├── CubicSplineSmoothingKernel.hpp
│   ├── QuarticSplineSmoothingKernel.cpp
│   ├── QuarticSplineSmoothingKernel.hpp
│   ├── ScaledGaussianSmoothingKernel.cpp
│   ├── ScaledGaussianSmoothingKernel.hpp
│   └── Makefile
│
├── Kernels/
│   └── ...
│
├── Data/
│   └── ...
│
└── Results/
    └── ...
```

The programs use relative paths and are therefore intended to be run from the `Code` directory.

- `Code/` contains the C++ source code and Makefile.
- `Kernels/` contains the cumulative-kernel lookup tables used by the `cumkernel` method.
- `Data/` contains input particle distributions.
- `Results/` contains the generated FITS density cubes.

## Dependencies

The code requires:

- a C++11-compatible compiler
- CFITSIO
- CCfits

The supplied Makefile currently assumes that CCfits and CFITSIO are installed under `/usr/local`.

If they are installed elsewhere, modify `INCLUDES` and `LIBS` in the Makefile accordingly.

## Compilation

From the `Code` directory, run

```bash
make
```

This builds two executables:

```text
VoxeliseSimulation
CalculateCumulativeKernel
```

To build only one executable, use

```bash
make VoxeliseSimulation
```

or

```bash
make CalculateCumulativeKernel
```

To remove compiled object files and executables, use

```bash
make clean
```

## Input particle format

The voxelisation code uses particles with five quantities:

```text
x   y   z   h   m
```

where:

- `x`, `y`, `z` are the Cartesian particle coordinates
- `h` is the smoothing length
- `m` is the particle mass

The current benchmark program reads these quantities from an ASCII text file with one particle per line.

Lines beginning with `#` are ignored. Malformed lines are skipped with a warning.

In the current example setup, positions and smoothing lengths are expressed in kpc and masses in solar masses. The corresponding output density is therefore expressed in `Msun/kpc3`.

## Voxelisation interface

The main interface is

```cpp
void voxelise(
    const std::string& method,
    const std::vector<Particle>& particles,
    const SmoothingKernel& kernel,
    std::valarray<double>& densities,
    int numx,
    int numy,
    int numz,
    double voxelsize);
```

The supported method strings are:

```text
ngp
centre
subsample3
subsample5
montecarlo10
montecarlo100
montecarlo1000
cumkernel
```

An unknown method name results in an exception.

The output density array is flattened according to

```text
i + numx*j + numx*numy*k
```

with the \(x\) index varying fastest.

## Smoothing kernels

All smoothing kernels derive from the abstract `SmoothingKernel` class and implement

```cpp
double density(double u) const;
```

where

```text
u = r / h
```

is the normalized radius.

The currently implemented classes are:

```text
CubicSplineSmoothingKernel
QuarticSplineSmoothingKernel
ScaledGaussianSmoothingKernel
```

For cumulative-kernel voxelisation, a kernel object is constructed with the desired lookup-table resolution, for example

```cpp
CubicSplineSmoothingKernel kernel(129);
```

This loads the corresponding cumulative-kernel FITS file from the `Kernels` directory.

## Cumulative-kernel lookup tables

The cumulative method uses a precomputed function

```text
Phi(X,Y,Z)
```

representing the integral of the normalized smoothing-kernel density over the rectangular region extending from the origin to \((X,Y,Z)\).

Only the positive octant needs to be stored because the smoothing kernels are spherically symmetric. Values for signed coordinates are obtained using the corresponding symmetry of the cumulative integral.

The lookup tables are stored as three-dimensional FITS cubes and are read by `SmoothingKernel`.

The expected filename convention is

```text
<KernelName>_N<N>.fits
```

where `<N>` is the number of lookup-table grid points along each coordinate axis.

For example:

```text
CubicSplineSmoothingKernel_N129.fits
QuarticSplineSmoothingKernel_N129.fits
ScaledGaussianSmoothingKernel_N129.fits
```

The required lookup-table resolution must match the value used when constructing the smoothing-kernel object. For example,

```cpp
CubicSplineSmoothingKernel kernel(129);
```

requires

```text
../Kernels/CubicSplineSmoothingKernel_N129.fits
```

If the required file is not present, the program terminates with an error rather than generating the table automatically.

Between tabulated points, `Phi` is evaluated using trilinear interpolation.

## Generating cumulative-kernel tables

Lookup tables are generated with the separate

```text
CalculateCumulativeKernel
```

executable.

The kernel type and lookup-table resolutions are currently selected directly in `CalculateCumulativeKernel.cpp`. For example,

```cpp
const int num = 129;
QuarticSplineSmoothingKernel kernel(num, false);
```

The second constructor argument, `false`, prevents the kernel object from attempting to read an existing cumulative-kernel table. This is required because the purpose of this program is to generate that table.

The cumulative integral is evaluated using three-dimensional Gauss-Legendre quadrature and written to a FITS cube in the `Kernels` directory.

Run the generator from the `Code` directory:

```bash
./CalculateCumulativeKernel
```

The generated lookup tables need to be calculated only once for each combination of smoothing kernel and lookup-table resolution.

Table generation, particularly at high resolution, can be computationally expensive.

## Benchmark program

`VoxeliseSimulation.cpp` is currently a benchmark and experiment driver rather than a general command-line application.

Parameters such as

- the input particle distribution
- field of view
- grid resolution
- smoothing kernel
- cumulative-kernel resolution
- voxelisation method

are currently specified directly in `VoxeliseSimulation.cpp`.

For example, the current code constructs a cubic spline kernel with a cumulative-table resolution of 129:

```cpp
int num = 129;
CubicSplineSmoothingKernel kernel(num);
```

The benchmark reads particle data from

```text
../Data/
```

and writes resulting FITS cubes to

```text
../Results/
```

The output filenames have the form

```text
<dataset>_<method>_n<grid-size>.fits
```

The FITS files contain the three-dimensional voxel density together with basic spatial metadata.

## Running the benchmark

After compiling, run from the `Code` directory:

```bash
./VoxeliseSimulation
```

Before running the cumulative-kernel method, make sure that the appropriate lookup table exists in `../Kernels/`.

For example, if `VoxeliseSimulation.cpp` contains

```cpp
CubicSplineSmoothingKernel kernel(129);
```

then the following file must exist:

```text
../Kernels/CubicSplineSmoothingKernel_N129.fits
```

## Adding another smoothing kernel

A new smoothing kernel can be added by deriving a class from `SmoothingKernel` and implementing

```cpp
double density(double u) const override;
```

The kernel should be spherically symmetric, normalized to unit total mass, and have finite support within \(u \leq 1\).

The new class should assign a unique kernel name and support construction without reading an existing cumulative table, following the pattern of the existing smoothing-kernel classes.

A cumulative lookup table can then be generated with `CalculateCumulativeKernel` and subsequently used by the `cumkernel` voxelisation method.

## Notes

This repository is primarily a research implementation accompanying the development and validation of the cumulative-kernel voxelisation method. The emphasis is on numerical experiments and comparison between voxelisation algorithms rather than on providing a general-purpose command-line interface.

The cumulative-kernel lookup tables are intentionally generated separately from the main voxelisation program. They depend only on the smoothing-kernel profile and lookup-table resolution and can therefore be reused for different particle distributions and voxel grids.

## License

This software is distributed under the BSD 3-Clause License. See the LICENSE file for details.

If you use this software in scientific work, please cite the accompanying publication.
