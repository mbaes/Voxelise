//////////////////////////////////////////////////////////////////////

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <string>
#include <valarray>
#include <vector>
#include <CCfits>
#include <CCfits/CCfits>

#include "Voxelise.hpp"
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
        std::string galaxy = "NGC628";
        double fov = 80;
        
        // The voxelisation method
        
        const std::string method = "cumkernel";
        
        // The number of voxels (in one direction)
        
        int n = 64;
        
        // The smoothing kernel
        
        int num = 129;
        CubicSplineSmoothingKernel kernel(num);
        
        // Read the data
        
        std::string inputfilename = "../Data/" + galaxy + ".txt";
        std::ifstream infile(inputfilename.c_str());
        if (!infile.good())
        {
            std::cerr << "Input file " << inputfilename << " does not exist.\n";
            return EXIT_FAILURE;
        }
        std::string line;
        Particle p;
        std::vector<Particle> particles;
        while (getline(infile, line))
        {
            if (line.size() != 0 && line[0] != '#')
            {
                std::istringstream iss(line);
                if (!(iss >> p.x >> p.y >> p.z >> p.h >> p.m))
                {
                    std::cerr << "Skipping malformed line: " << line << '\n';
                    continue;
                }
                particles.push_back(p);
            }
        }
        
        // Perform the voxelisation
        
        double voxelsize = fov/n;
        
        std::cout << "Starting voxelisation with method " << method << " and n = " << n << ".\n";
        std::valarray<double> densities;
        auto start = std::chrono::steady_clock::now();
        voxelise(method, particles, kernel, densities, n, n, n, voxelsize);
        auto stop = std::chrono::steady_clock::now();
        std::chrono::duration<double> elapsed = stop - start;
        std::cout << "Elapsed time: " << elapsed.count() << " s\n";
        std::cout << "Finished voxelisation.\n";
        
        // Write the results to a FITS file
        
        std::string outputfilename = "../Results/" + galaxy + "_" + method + "_n" + std::to_string(n) + ".fits";
        long naxis = 3;
        long naxes[3];
        const long n3 = static_cast<long>(n) * n * n;
        naxes[0] = n;
        naxes[1] = n;
        naxes[2] = n;
        std::unique_ptr<FITS> pFits;
        pFits.reset( new FITS("!" + outputfilename , DOUBLE_IMG , naxis , naxes ) );
        pFits->pHDU().write(1, n3, densities);
        pFits->pHDU().addKey("ORIGIN", "Generated with Voxelise", "");
        pFits->pHDU().addKey("BUNIT", "Msun/kpc3", "");
        pFits->pHDU().addKey("CRPIX1", 0.5*(n+1), "");
        pFits->pHDU().addKey("CRVAL1", 0.0, "");
        pFits->pHDU().addKey("CDELT1", voxelsize, "");
        pFits->pHDU().addKey("CUNIT1", "kpc", "");
        pFits->pHDU().addKey("CRPIX2", 0.5*(n+1), "");
        pFits->pHDU().addKey("CRVAL2", 0.0, "");
        pFits->pHDU().addKey("CDELT2", voxelsize, "");
        pFits->pHDU().addKey("CUNIT2", "kpc", "");
        pFits->pHDU().addKey("CRPIX3", 0.5*(n+1), "");
        pFits->pHDU().addKey("CRVAL3", 0.0, "");
        pFits->pHDU().addKey("CDELT3", voxelsize, "");
        pFits->pHDU().addKey("CUNIT3", "kpc", "");
        std::cout << "File " << outputfilename << " generated." << std::endl;
    }
    catch (const std::exception& e)
    {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }
    return EXIT_SUCCESS;
}
        
//////////////////////////////////////////////////////////////////////
