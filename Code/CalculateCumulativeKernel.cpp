//////////////////////////////////////////////////////////////////////

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <string>
#include <valarray>
#include <vector>
#include <CCfits>
#include <CCfits/CCfits>

#include "SmoothingKernel.hpp"
#include "CubicSplineSmoothingKernel.hpp"
#include "QuarticSplineSmoothingKernel.hpp"
#include "ScaledGaussianSmoothingKernel.hpp"

using namespace CCfits;

//////////////////////////////////////////////////////////////////////

int main(void)
{
    try
    {
        const int num = 129;
        CubicSplineSmoothingKernel kernel(num, false);

        std::string filename = "../Kernels/" + kernel.name() + "_N" + std::to_string(num) + ".fits";
        const size_t num3 = static_cast<std::size_t>(num) * num * num;
        std::valarray<double> cumkernel(num3);
        
        // Definition of the Gauss-Legendre nodes and weights.
        
        int numg = 64;
        std::vector<double> thetav = {0.00034747913211389999, 0.001829941614022, 0.0044933142616279996, 0.0083318730576869993, 0.01333658610504, 0.01949560017397, 0.026794312570800002, 0.03521541393403, 0.044738931460750003, 0.055342277002440002, 0.067000300922950007, 0.079685351873709995, 0.093367342438599996, 0.10801382052833, 0.1235900463697, 0.1400590749142, 0.15738184347289999, 0.17551726437269999, 0.1944223224138, 0.21405217689869999, 0.23436026799010001, 0.25529842714649997, 0.27681699137330001, 0.298864921018, 0.32138992083119999, 0.34433856400489998, 0.36765641889560002, 0.39128817813, 0.41517778978800002, 0.43926859035190002, 0.46350343910609998, 0.48782485366829997, 0.51217514633170003, 0.53649656089389997, 0.56073140964809998, 0.58482221021199998, 0.60871182186999995, 0.63234358110440003, 0.65566143599509996, 0.67861007916880001, 0.701135078982, 0.72318300862670004, 0.74470157285350003, 0.76563973200989999, 0.78594782310129996, 0.80557767758620002, 0.82448273562730001, 0.84261815652710004, 0.8599409250858, 0.87640995363030005, 0.8919861794717, 0.90663265756139999, 0.92031464812630004, 0.93299969907699998, 0.94465772299759998, 0.95526106853930004, 0.96478458606600004, 0.97320568742919999, 0.98050439982600002, 0.98666341389500001, 0.9916681269423, 0.99550668573839995, 0.99817005838600004, 0.99965252086790002};
        std::vector<double> weightv = {0.00089164036084821597, 0.002073516630281234, 0.003252228984489181, 0.0044233799131819743, 0.0055840697300655641, 0.0067315239483593196, 0.0078630152380123608, 0.0089758578878486699, 0.010067411576765099, 0.01113508690419163, 0.01217635128435544, 0.013188734857527329, 0.01416983630712974, 0.015117328536201241, 0.016028964177425779, 0.0169025809185708, 0.017736106628441189, 0.01852756427012002, 0.01927507658930781, 0.019976870566360171, 0.02063128162131176, 0.021236757561826788, 0.021791862264661729, 0.02229527908187828, 0.022745813963709071, 0.023142398290657212, 0.02348409140810501, 0.023770082857415151, 0.023999694298229152, 0.024172381117401481, 0.024287733720751711, 0.0243454785045699, 0.0243454785045699, 0.024287733720751711, 0.024172381117401481, 0.023999694298229152, 0.023770082857415151, 0.02348409140810501, 0.023142398290657212, 0.022745813963709071, 0.02229527908187828, 0.021791862264661729, 0.021236757561826788, 0.02063128162131176, 0.019976870566360171, 0.01927507658930781, 0.01852756427012002, 0.017736106628441189, 0.0169025809185708, 0.016028964177425779, 0.015117328536201241, 0.01416983630712974, 0.013188734857527329, 0.01217635128435544, 0.01113508690419163, 0.010067411576765099, 0.0089758578878486699, 0.0078630152380123608, 0.0067315239483593196, 0.0055840697300655641, 0.0044233799131819743, 0.003252228984489181, 0.002073516630281234, 0.00089164036084821597};
        
        // Calculation of the function Phi(X,Y,Z) in the grid points.
        
        std::vector<double> Xv(num,0.0);
        std::vector<double> Yv(num,0.0);
        std::vector<double> Zv(num,0.0);
        for (int i=0; i<num; i++) Xv[i] = i/(num-1.0);
        for (int j=0; j<num; j++) Yv[j] = j/(num-1.0);
        for (int k=0; k<num; k++) Zv[k] = k/(num-1.0);
        
        for (int i=0; i<num; i++)
        {
            std::cout << "  i = " << i+1 << " / " << num << std::endl;
            double X = Xv[i];
            for (int j=0; j<num; j++)
            {
                double Y = Yv[j];
                for (int k=0; k<num; k++)
                {
                    double Z = Zv[k];
                    int l = i + num*j + num*num*k;
                    
                    if (X == 0.0 || Y == 0.0 || Z == 0.0)
                    {
                        cumkernel[l] = 0.0;
                        continue;
                    }
                    
                    if (X >= 1.0 && Y >= 1.0 && Z >= 1.0)
                    {
                        cumkernel[l] = 0.125;
                        continue;
                    }
                    
                    double sum = 0.0;
                    for (int im=0; im<numg; im++)
                    {
                        double Xm = thetav[im] * X;
                        double weightX = weightv[im];
                        for (int jm=0; jm<numg; jm++)
                        {
                            double Ym = thetav[jm] * Y;
                            double weightY = weightv[jm];
                            for (int km=0; km<numg; km++)
                            {
                                double Zm = thetav[km] * Z;
                                double weightZ = weightv[km];
                                sum += weightX * weightY * weightZ * kernel.density(sqrt(Xm*Xm+Ym*Ym+Zm*Zm));
                            }
                        }
                    }
                    cumkernel[l] = sum * X * Y * Z;
                }
            }
        }
        
        // Write the precalculated values to a FITS file
        
        long naxis = 3;
        long naxes[3];
        naxes[0] = num;
        naxes[1] = num;
        naxes[2] = num;
        long numelements = static_cast<long>(num3);
        std::unique_ptr<FITS> pFits;

        pFits.reset(new FITS("!" + filename, DOUBLE_IMG, naxis, naxes));
        pFits->pHDU().write(1, numelements, cumkernel);
        pFits->pHDU().addKey("ORIGIN", "Generated with Voxelise", "");
        pFits->pHDU().addKey("BUNIT", "dimensionless", "");
        pFits->pHDU().addKey("KERNEL", kernel.name(), "");
        pFits->pHDU().addKey("NGRID", num, "Number of grid points per axis");
        
        pFits->pHDU().addKey("CRPIX1", 1.0, "");
        pFits->pHDU().addKey("CRVAL1", 0.0, "");
        pFits->pHDU().addKey("CDELT1", 1.0 / (num - 1.0), "");
        pFits->pHDU().addKey("CTYPE1", "X", "");
        pFits->pHDU().addKey("CRPIX2", 1.0, "");
        pFits->pHDU().addKey("CRVAL2", 0.0, "");
        pFits->pHDU().addKey("CDELT2", 1.0 / (num - 1.0), "");
        pFits->pHDU().addKey("CTYPE2", "Y", "");
        pFits->pHDU().addKey("CRPIX3", 1.0, "");
        pFits->pHDU().addKey("CRVAL3", 0.0, "");
        pFits->pHDU().addKey("CDELT3", 1.0 / (num - 1.0), "");
        pFits->pHDU().addKey("CTYPE3", "Z", "");
        
        std::cout << "Lookup table written to " << filename << "." << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
        
//////////////////////////////////////////////////////////////////////
