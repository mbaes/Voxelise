//////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <valarray>
#include <vector>

#include "SmoothingKernel.hpp"

//////////////////////////////////////////////////////////////////////

struct Particle
{
    double x;
    double y;
    double z;
    double h;
    double m;
};

//////////////////////////////////////////////////////////////////////

void voxelise(const std::string& method, const std::vector<Particle>& particles, const SmoothingKernel& kernel, std::valarray<double>& densityv, int numx, int numy, int numz, double voxelsize);

void voxelise_ngp(const std::vector<Particle>& particles, std::valarray<double>& densityv, int numx, int numy, int numz, double voxelsize);

void voxelise_centre(const std::vector<Particle>& particles, const SmoothingKernel& kernel, std::valarray<double>& densityv, int numx, int numy, int numz, double voxelsize);

void voxelise_subsample(const std::vector<Particle>& particles, const SmoothingKernel& kernel, std::valarray<double>& densityv, int numx, int numy, int numz, double voxelsize, int numss);

void voxelise_montecarlo(const std::vector<Particle>& particles, const SmoothingKernel& kernel, std::valarray<double>& densityv, int numx, int numy, int numz, double voxelsize, int nummc);

void voxelise_cumkernel(const std::vector<Particle>& particles, const SmoothingKernel& kernel, std::valarray<double>& densityv, int numx, int numy, int numz, double voxelsize);

//////////////////////////////////////////////////////////////////////
