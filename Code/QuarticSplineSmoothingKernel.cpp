//////////////////////////////////////////////////////////////////////

#include <cmath>
#include "QuarticSplineSmoothingKernel.hpp"

//////////////////////////////////////////////////////////////////////

QuarticSplineSmoothingKernel::QuarticSplineSmoothingKernel(int num, bool readkernel) : SmoothingKernel(num)
{
    _name = "QuarticSplineSmoothingKernel";
    if (readkernel)
        readCumulativeKernelFile();
}

//////////////////////////////////////////////////////////////////////

double QuarticSplineSmoothingKernel::density(double u) const
{
    if (u < 0. || u >= 1.) return 0.;

    double u2 = u * u;
    double u3 = u2 * u;
    double u4 = u2 * u2;

    if (u < 0.2) return 58.284280917442139563 * u4 - 23.313712366976855825 * u2 + 3.5747692296031178931;

    if (u < 0.6)
        return -38.856187278294759708 * u4 + 77.712374556589519415 * u3 - 46.627424733953711649 * u2
               + 3.1084949822635807767 * u + 3.4193444804899388544;

    return 9.7140468195736899271 * u4 - 38.856187278294759708 * u3 + 58.284280917442139563 * u2
           - 38.856187278294759708 * u + 9.7140468195736899271;
}

//////////////////////////////////////////////////////////////////////
