//////////////////////////////////////////////////////////////////////

#include <cmath>
#include "CubicSplineSmoothingKernel.hpp"

//////////////////////////////////////////////////////////////////////

CubicSplineSmoothingKernel::CubicSplineSmoothingKernel(int num, bool readkernel) :
    SmoothingKernel(num)
{
    _name = "CubicSplineSmoothingKernel";
    if (readkernel)
        readCumulativeKernelFile();
}

//////////////////////////////////////////////////////////////////////

double CubicSplineSmoothingKernel::density(double u) const
{
    if (u < 0.0 || u >= 1.0)
        return 0.0;
    else if (u < 0.5)
        return 2.5464790894703253723 * (1.0 - 6.0 * u * u * (1.0 - u));
    else
        return 5.0929581789406507446 * (1.0 - u) * (1.0 - u) * (1.0 - u);
}

//////////////////////////////////////////////////////////////////////
