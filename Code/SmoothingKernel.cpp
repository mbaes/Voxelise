//////////////////////////////////////////////////////////////////////

#include "SmoothingKernel.hpp"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <memory>
#include <stdexcept>
#include <CCfits>
#include <CCfits/CCfits>

using namespace CCfits;

//////////////////////////////////////////////////////////////////////

SmoothingKernel::SmoothingKernel(int num):
    _num(num)
{
    if (num < 2)
    {
        throw std::invalid_argument("Error in SmoothingKernel::SmoothingKernel: lookup table must contain at least 2 grid points per axis.");
    }
}

//////////////////////////////////////////////////////////////////////

const std::string& SmoothingKernel::name() const
{
    return _name;
}

//////////////////////////////////////////////////////////////////////

int SmoothingKernel::numGridPoints() const
{
    return _num;
}

//////////////////////////////////////////////////////////////////////

void SmoothingKernel::readCumulativeKernelFile()
{
    const std::string filename = "../Kernels/" + _name + "_N" + std::to_string(_num) + ".fits";
    std::ifstream file(filename);
    if (!file.good())
    {
        throw std::runtime_error("Error in SmoothingKernel::readCumulativeKernelFile: FITS file " + filename + " does not exist.");
    }
    std::unique_ptr<FITS> pFits(new FITS(filename, Read, true));
    PHDU& image = pFits->pHDU();
    if (image.axes() != 3)
    {
        throw std::runtime_error("Error in SmoothingKernel::readCumulativeKernelFile: FITS file " + filename + " is not a 3D cube.");
    }
            
    long num1 = image.axis(0);
    long num2 = image.axis(1);
    long num3 = image.axis(2);
    if (num1 != num2 || num2 != num3 || num3 != _num)
    {
        throw std::runtime_error("Error in SmoothingKernel::readCumulativeKernelFile: FITS file " + filename + " does not contain a data cube with the desired dimensions.");
    }
    
    image.read(_cumkernel);
}

//////////////////////////////////////////////////////////////////////

double SmoothingKernel::PhiSigned(double X, double Y, double Z) const
{
    double s = 1.0;
    if (X < 0.0) s *= -1.0;
    if (Y < 0.0) s *= -1.0;
    if (Z < 0.0) s *= -1.0;
    return s * Phi(std::abs(X), std::abs(Y), std::abs(Z));
}

//////////////////////////////////////////////////////////////////////

double SmoothingKernel::Phi(double X, double Y, double Z) const
{
    if (_cumkernel.size() == 0)
    {
        throw std::runtime_error("Error in SmoothingKernel::Phi: cumulative kernel has not been loaded.");
    }
    if (X < 0.0 || Y < 0.0 || Z < 0.0)
    {
        throw std::runtime_error("Error in SmoothingKernel::Phi: invalid arguments.");
    }
    if (X >= 1.0 && Y >= 1.0 && Z >= 1.0) return 0.125;
    X = std::min(X, 1.0);
    Y = std::min(Y, 1.0);
    Z = std::min(Z, 1.0);
    double scale = (_num - 1.0);
    X *= scale;
    Y *= scale;
    Z *= scale;
    int i0 = static_cast<int>(std::floor(X));
    int j0 = static_cast<int>(std::floor(Y));
    int k0 = static_cast<int>(std::floor(Z));
    if (i0 >= _num - 1) i0 = _num - 2;
    if (j0 >= _num - 1) j0 = _num - 2;
    if (k0 >= _num - 1) k0 = _num - 2;
    int i1 = i0 + 1;
    int j1 = j0 + 1;
    int k1 = k0 + 1;
    double tx = X - i0;
    double ty = Y - j0;
    double tz = Z - k0;
    
    int num = _num;
    auto idx = [num](int i, int j, int k)
    {
        return i + num*j + num*num*k;
    };
    double v000 = _cumkernel[idx(i0,j0,k0)];
    double v100 = _cumkernel[idx(i1,j0,k0)];
    double v010 = _cumkernel[idx(i0,j1,k0)];
    double v110 = _cumkernel[idx(i1,j1,k0)];
    double v001 = _cumkernel[idx(i0,j0,k1)];
    double v101 = _cumkernel[idx(i1,j0,k1)];
    double v011 = _cumkernel[idx(i0,j1,k1)];
    double v111 = _cumkernel[idx(i1,j1,k1)];
    double v00 = (1.0 - tx)*v000 + tx*v100;
    double v10 = (1.0 - tx)*v010 + tx*v110;
    double v01 = (1.0 - tx)*v001 + tx*v101;
    double v11 = (1.0 - tx)*v011 + tx*v111;
    double v0 = (1.0 - ty)*v00 + ty*v10;
    double v1 = (1.0 - ty)*v01 + ty*v11;
    return (1.0 - tz)*v0 + tz*v1;
}

//////////////////////////////////////////////////////////////////////
