//////////////////////////////////////////////////////////////////////

#include "SmoothingKernel.hpp"
#include "Voxelise.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>

//////////////////////////////////////////////////////////////////////

void voxelise(const std::string& method, const std::vector<Particle>& particles, const SmoothingKernel& kernel, std::valarray<double>& densities, int numx, int numy, int numz, double voxelsize)
{
    if (method=="ngp")
        voxelise_ngp(particles, densities, numx, numy, numz, voxelsize);
    else if (method=="centre")
        voxelise_centre(particles, kernel, densities, numx, numy, numz, voxelsize);
    else if (method=="subsample3")
        voxelise_subsample(particles, kernel, densities, numx, numy, numz, voxelsize, 3);
    else if (method=="subsample5")
        voxelise_subsample(particles, kernel, densities, numx, numy, numz, voxelsize, 5);
    else if (method=="montecarlo10")
        voxelise_montecarlo(particles, kernel, densities, numx, numy, numz, voxelsize, 10);
    else if (method=="montecarlo100")
        voxelise_montecarlo(particles, kernel, densities, numx, numy, numz, voxelsize, 100);
    else if (method=="montecarlo1000")
        voxelise_montecarlo(particles, kernel, densities, numx, numy, numz, voxelsize, 1000);
    else if (method=="cumkernel")
        voxelise_cumkernel(particles, kernel, densities, numx, numy, numz, voxelsize);
    else
    {
        throw std::invalid_argument("Unknown voxelisation method: " + method);
    }
    return;
}
    
//////////////////////////////////////////////////////////////////////

void voxelise_ngp(const std::vector<Particle>& particles, std::valarray<double>& densities, int numx, int numy, int numz, double voxelsize)
{
    if (numx <= 0 || numy <= 0 || numz <= 0)
        throw std::invalid_argument("Grid dimensions must be positive");
    if (voxelsize <= 0.0)
        throw std::invalid_argument("Voxel size must be positive");
    const double xfov = numx * voxelsize;
    const double yfov = numy * voxelsize;
    const double zfov = numz * voxelsize;
    const std::size_t numCells = static_cast<std::size_t>(numx) * numy * numz;
    densities = std::valarray<double>(0.0, numCells);

    for (const auto& particle : particles)
    {
        const double x = particle.x;
        const double y = particle.y;
        const double z = particle.z;
        const double h = particle.h;
        const double m = particle.m;
        if (h <= 0.0 || m <= 0.0) continue;
        const int i = static_cast<int>(std::floor((x + 0.5*xfov) / voxelsize));
        const int j = static_cast<int>(std::floor((y + 0.5*yfov) / voxelsize));
        const int k = static_cast<int>(std::floor((z + 0.5*zfov) / voxelsize));
        if (i < 0 || i >= numx || j < 0 || j >= numy || k < 0 || k >= numz) continue;
        const int l = i + numx*j + numx*numy*k;
        densities[l] += m / (voxelsize * voxelsize * voxelsize);
    }
    return;
}

//////////////////////////////////////////////////////////////////////

void voxelise_centre(const std::vector<Particle>& particles, const SmoothingKernel& kernel, std::valarray<double>& densities, int numx, int numy, int numz, double voxelsize)
{
    if (numx <= 0 || numy <= 0 || numz <= 0)
        throw std::invalid_argument("Grid dimensions must be positive");
    if (voxelsize <= 0.0)
        throw std::invalid_argument("Voxel size must be positive");
    const double xfov = numx * voxelsize;
    const double yfov = numy * voxelsize;
    const double zfov = numz * voxelsize;
    const std::size_t numCells = static_cast<std::size_t>(numx) * numy * numz;
    densities = std::valarray<double>(0.0, numCells);

    std::vector<double> xv(numx);
    std::vector<double> yv(numy);
    std::vector<double> zv(numz);
    for (int i=0; i<numx; i++) xv[i] = -0.5*xfov + (i+0.5)*voxelsize;
    for (int j=0; j<numy; j++) yv[j] = -0.5*yfov + (j+0.5)*voxelsize;
    for (int k=0; k<numz; k++) zv[k] = -0.5*zfov + (k+0.5)*voxelsize;

    for (const auto& particle : particles)
    {
        const double x = particle.x;
        const double y = particle.y;
        const double z = particle.z;
        const double h = particle.h;
        const double m = particle.m;
        if (h <= 0.0 || m <= 0.0) continue;

        const int imin = std::max(0, static_cast<int>(std::floor((x-h+0.5*xfov) / voxelsize)));
        const int imax = std::min(numx-1, static_cast<int>(std::floor((x+h+0.5*xfov) / voxelsize)));
        const int jmin = std::max(0, static_cast<int>(std::floor((y-h+0.5*yfov) / voxelsize)));
        const int jmax = std::min(numy-1, static_cast<int>(std::floor((y+h+0.5*yfov) / voxelsize)));
        const int kmin = std::max(0, static_cast<int>(std::floor((z-h+0.5*zfov) / voxelsize)));
        const int kmax = std::min(numz-1, static_cast<int>(std::floor((z+h+0.5*zfov) / voxelsize)));
        if (imin > imax || jmin > jmax || kmin > kmax)
            continue;
        
        for (int i=imin; i<=imax; i++)
        {
            const double xi = xv[i];
            for (int j=jmin; j<=jmax; j++)
            {
                const double yj = yv[j];
                for (int k=kmin; k<=kmax; k++)
                {
                    const double zk = zv[k];
                    const int l = i + numx*j + numx*numy*k;
                    const double u2 = ((x-xi)*(x-xi) + (y-yj)*(y-yj) + (z-zk)*(z-zk)) / (h*h);
                    const double u = std::sqrt(u2);
                    densities[l] += m / (h*h*h) * kernel.density(u);
                }
            }
        }
    }
    return;
}

//////////////////////////////////////////////////////////////////////

void voxelise_subsample(const std::vector<Particle>& particles, const SmoothingKernel& kernel, std::valarray<double>& densities, int numx, int numy, int numz, double voxelsize, int numss)
{
    if (numx <= 0 || numy <= 0 || numz <= 0)
        throw std::invalid_argument("Grid dimensions must be positive");
    if (numss <= 0)
        throw std::invalid_argument("Subsampling rate must be positive");
    if (voxelsize <= 0.0)
        throw std::invalid_argument("Voxel size must be positive");
    const double xfov = numx * voxelsize;
    const double yfov = numy * voxelsize;
    const double zfov = numz * voxelsize;
    const std::size_t numCells = static_cast<std::size_t>(numx) * numy * numz;
    densities = std::valarray<double>(0.0, numCells);
    const double invnumss = 1.0 / numss;
    const double norm = 1.0 / (numss*numss*numss);

    std::vector<double> xv(numx);
    std::vector<double> yv(numy);
    std::vector<double> zv(numz);
    for (int i=0; i<numx; i++) xv[i] = -0.5*xfov + (i+0.5)*voxelsize;
    for (int j=0; j<numy; j++) yv[j] = -0.5*yfov + (j+0.5)*voxelsize;
    for (int k=0; k<numz; k++) zv[k] = -0.5*zfov + (k+0.5)*voxelsize;

    for (const auto& particle : particles)
    {
        const double x = particle.x;
        const double y = particle.y;
        const double z = particle.z;
        const double h = particle.h;
        const double m = particle.m;
        if (h <= 0.0 || m <= 0.0) continue;

        const int imin = std::max(0, static_cast<int>(std::floor((x-h+0.5*xfov) / voxelsize)));
        const int imax = std::min(numx-1, static_cast<int>(std::floor((x+h+0.5*xfov) / voxelsize)));
        const int jmin = std::max(0, static_cast<int>(std::floor((y-h+0.5*yfov) / voxelsize)));
        const int jmax = std::min(numy-1, static_cast<int>(std::floor((y+h+0.5*yfov) / voxelsize)));
        const int kmin = std::max(0, static_cast<int>(std::floor((z-h+0.5*zfov) / voxelsize)));
        const int kmax = std::min(numz-1, static_cast<int>(std::floor((z+h+0.5*zfov) / voxelsize)));
        if (imin > imax || jmin > jmax || kmin > kmax)
            continue;
        
        for (int i=imin; i<=imax; i++)
        {
            const double xi = xv[i];
            for (int j=jmin; j<=jmax; j++)
            {
                const double yj = yv[j];
                for (int k=kmin; k<=kmax; k++)
                {
                    const double zk = zv[k];
                    const int l = i + numx*j + numx*numy*k;
                    for (int iss=0; iss<numss; iss++)
                    {
                        const double xr = xi + ((iss + 0.5)*invnumss - 0.5) * voxelsize;
                        for (int jss=0; jss<numss; jss++)
                        {
                            const double yr = yj + ((jss + 0.5)*invnumss - 0.5) * voxelsize;
                            for (int kss=0; kss<numss; kss++)
                            {
                                const double zr = zk + ((kss + 0.5)*invnumss - 0.5) * voxelsize;
                                const double u2 = ((x-xr)*(x-xr) + (y-yr)*(y-yr) + (z-zr)*(z-zr)) / (h*h);
                                const double u = std::sqrt(u2);
                                densities[l] += m / (h*h*h) * kernel.density(u) * norm;
                            }
                        }
                    }
                }
            }
        }
    }
    return;
}

//////////////////////////////////////////////////////////////////////

void voxelise_montecarlo(const std::vector<Particle>& particles, const SmoothingKernel& kernel, std::valarray<double>& densities, int numx, int numy, int numz, double voxelsize, int nummc)
{
    if (numx <= 0 || numy <= 0 || numz <= 0)
        throw std::invalid_argument("Grid dimensions must be positive");
    if (nummc <= 0)
        throw std::invalid_argument("Monte Carlo subsampling rate must be positive");
    if (voxelsize <= 0.0)
        throw std::invalid_argument("Voxel size must be positive");
    const double xfov = numx * voxelsize;
    const double yfov = numy * voxelsize;
    const double zfov = numz * voxelsize;
    const std::size_t numCells = static_cast<std::size_t>(numx) * numy * numz;
    densities = std::valarray<double>(0.0, numCells);
    std::mt19937 gen(123456);
    std::uniform_real_distribution<double> dis(-0.5, 0.5);
    
    std::vector<double> xv(numx);
    std::vector<double> yv(numy);
    std::vector<double> zv(numz);
    for (int i=0; i<numx; i++) xv[i] = -0.5*xfov + (i+0.5)*voxelsize;
    for (int j=0; j<numy; j++) yv[j] = -0.5*yfov + (j+0.5)*voxelsize;
    for (int k=0; k<numz; k++) zv[k] = -0.5*zfov + (k+0.5)*voxelsize;

    for (const auto& particle : particles)
    {
        const double x = particle.x;
        const double y = particle.y;
        const double z = particle.z;
        const double h = particle.h;
        const double m = particle.m;
        if (h <= 0.0 || m <= 0.0) continue;

        const int imin = std::max(0, static_cast<int>(std::floor((x-h+0.5*xfov) / voxelsize)));
        const int imax = std::min(numx-1, static_cast<int>(std::floor((x+h+0.5*xfov) / voxelsize)));
        const int jmin = std::max(0, static_cast<int>(std::floor((y-h+0.5*yfov) / voxelsize)));
        const int jmax = std::min(numy-1, static_cast<int>(std::floor((y+h+0.5*yfov) / voxelsize)));
        const int kmin = std::max(0, static_cast<int>(std::floor((z-h+0.5*zfov) / voxelsize)));
        const int kmax = std::min(numz-1, static_cast<int>(std::floor((z+h+0.5*zfov) / voxelsize)));
        if (imin > imax || jmin > jmax || kmin > kmax)
            continue;
        
        for (int i=imin; i<=imax; i++)
        {
            const double xi = xv[i];
            for (int j=jmin; j<=jmax; j++)
            {
                const double yj = yv[j];
                for (int k=kmin; k<=kmax; k++)
                {
                    const double zk = zv[k];
                    const int l = i + numx*j + numx*numy*k;
                    for (int q=0; q<nummc; q++)
                    {
                        const double xr = xi + dis(gen) * voxelsize;
                        const double yr = yj + dis(gen) * voxelsize;
                        const double zr = zk + dis(gen) * voxelsize;
                        const double u2 = ((x-xr)*(x-xr) + (y-yr)*(y-yr) + (z-zr)*(z-zr)) / (h*h);
                        const double u = std::sqrt(u2);
                        densities[l] += m / (h*h*h) * kernel.density(u) / nummc;
                    }
                }
            }
        }
    }
    return;
}

//////////////////////////////////////////////////////////////////////

void voxelise_cumkernel(const std::vector<Particle>& particles, const SmoothingKernel& kernel, std::valarray<double>& densities, int numx, int numy, int numz, double voxelsize)
{
    if (numx <= 0 || numy <= 0 || numz <= 0)
        throw std::invalid_argument("Grid dimensions must be positive");
    if (voxelsize <= 0.0)
        throw std::invalid_argument("Voxel size must be positive");
    const double xfov = numx * voxelsize;
    const double yfov = numy * voxelsize;
    const double zfov = numz * voxelsize;
    const std::size_t numCells = static_cast<std::size_t>(numx) * numy * numz;
    densities = std::valarray<double>(0.0, numCells);
    const double invVoxelVolume = 1.0 / (voxelsize * voxelsize * voxelsize);
    
    std::vector<double> xedge(numx + 1);
    std::vector<double> yedge(numy + 1);
    std::vector<double> zedge(numz + 1);
    for (int i=0; i <= numx; ++i)
        xedge[i] = -0.5 * xfov + i * voxelsize;
    for (int j=0; j <= numy; ++j)
        yedge[j] = -0.5 * yfov + j * voxelsize;
    for (int k=0; k <= numz; ++k)
        zedge[k] = -0.5 * zfov + k * voxelsize;
    
    std::vector<double> cumkernel;
    std::vector<double> Xs;
    std::vector<double> Ys;
    std::vector<double> Zs;
    
    for (const auto& particle : particles)
    {
        const double x = particle.x;
        const double y = particle.y;
        const double z = particle.z;
        const double h = particle.h;
        const double m = particle.m;
        if (h <= 0.0 || m <= 0.0) continue;
        const double invh = 1.0 / h;
        
        const int imin = std::max(0, static_cast<int>(std::floor((x-h+0.5*xfov) / voxelsize)));
        const int imax = std::min(numx-1, static_cast<int>(std::floor((x+h+0.5*xfov) / voxelsize)));
        const int jmin = std::max(0, static_cast<int>(std::floor((y-h+0.5*yfov) / voxelsize)));
        const int jmax = std::min(numy-1, static_cast<int>(std::floor((y+h+0.5*yfov) / voxelsize)));
        const int kmin = std::max(0, static_cast<int>(std::floor((z-h+0.5*zfov) / voxelsize)));
        const int kmax = std::min(numz-1, static_cast<int>(std::floor((z+h+0.5*zfov) / voxelsize)));
        if (imin > imax || jmin > jmax || kmin > kmax)
            continue;
        
        const int numX = imax - imin + 2;
        const int numY = jmax - jmin + 2;
        const int numZ = kmax - kmin + 2;
        const std::size_t numCorners = static_cast<std::size_t>(numX) * static_cast<std::size_t>(numY) * static_cast<std::size_t>(numZ);
        cumkernel.resize(numCorners);
        Xs.resize(numX);
        Ys.resize(numY);
        Zs.resize(numZ);
        for (int i=imin; i<=imax+1; i++)
            Xs[i-imin] = std::max(-1.0, std::min(1.0, (xedge[i] - x) * invh));
        for (int j=jmin; j<=jmax+1; j++)
            Ys[j-jmin] = std::max(-1.0, std::min(1.0, (yedge[j] - y) * invh));
        for (int k=kmin; k<=kmax+1; k++)
            Zs[k-kmin] = std::max(-1.0, std::min(1.0, (zedge[k] - z) * invh));
        for (int k = 0; k < numZ; ++k)
            for (int j = 0; j < numY; ++j)
                for (int i = 0; i < numX; ++i)
                    cumkernel[i + numX*j + numX*numY*k] = kernel.PhiSigned(Xs[i], Ys[j], Zs[k]);
        
        for (int k = kmin; k <= kmax; ++k)
            for (int j = jmin; j <= jmax; ++j)
                for (int i = imin; i <= imax; ++i)
                {
                    const int ii = i - imin;
                    const int jj = j - jmin;
                    const int kk = k - kmin;
                    const double weight =
                    cumkernel[(ii+1) + numX*(jj+1) + numX*numY*(kk+1)]
                    - cumkernel[ii     + numX*(jj+1) + numX*numY*(kk+1)]
                    - cumkernel[(ii+1) + numX*jj     + numX*numY*(kk+1)]
                    - cumkernel[(ii+1) + numX*(jj+1) + numX*numY*kk]
                    + cumkernel[ii     + numX*jj     + numX*numY*(kk+1)]
                    + cumkernel[ii     + numX*(jj+1) + numX*numY*kk]
                    + cumkernel[(ii+1) + numX*jj     + numX*numY*kk]
                    - cumkernel[ii     + numX*jj     + numX*numY*kk];
                    const int l = i + numx*j + numx*numy*k;
                    densities[l] += m * weight * invVoxelVolume;
                }
    }
    return;
}

//////////////////////////////////////////////////////////////////////
