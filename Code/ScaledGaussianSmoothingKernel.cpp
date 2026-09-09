//////////////////////////////////////////////////////////////////////

#include <cmath>
#include "ScaledGaussianSmoothingKernel.hpp"

//////////////////////////////////////////////////////////////////////

ScaledGaussianSmoothingKernel::ScaledGaussianSmoothingKernel(int num, bool readkernel) :
    SmoothingKernel(num)
{
    _name = "ScaledGaussianSmoothingKernel";
    if (readkernel)
        readCumulativeKernelFile();
}

//////////////////////////////////////////////////////////////////////

double ScaledGaussianSmoothingKernel::density(double u) const
{
    if (u < 0.0 || u >= 1.0) return 0.;
    return 2.56810060330949540082 * std::exp(-5.85836755024609305208 * u * u);
}

//////////////////////////////////////////////////////////////////////
