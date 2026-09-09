//////////////////////////////////////////////////////////////////////

#pragma once

#include <string>
#include <valarray>

//////////////////////////////////////////////////////////////////////

class SmoothingKernel
{
public:
    SmoothingKernel(int num);
    virtual ~SmoothingKernel() = default;
    const std::string& name() const;
    int numGridPoints() const;
    virtual double density(double u) const = 0;
    double PhiSigned(double X, double Y, double Z) const;
    double Phi(double X, double Y, double Z) const;

protected:
    void readCumulativeKernelFile();

protected:
    int _num;
    std::string _name;
    std::valarray<double> _cumkernel;
};

//////////////////////////////////////////////////////////////////////
