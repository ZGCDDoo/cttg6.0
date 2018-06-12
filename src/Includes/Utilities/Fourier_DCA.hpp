#pragma once

#include "Integrator.hpp"
#include "Utilities.hpp"
#include "MPIUtilities.hpp"
#include "GreenMat.hpp"
#include "ABC_SelfConsistency.hpp"

namespace FourierDCA
{
using DataK_t = ClusterMatrixCD_t; //first index (row) = K, second index (col) =iwn

DataK_t RtoK(const ClusterCubeCD_t &greenR, const ClusterSites_t &RSites, const ClusterSites_t &KWaveVectors)
{
    assert(RSites.size() == KWaveVectors.size());
    DataK_t greenK(KWaveVectors.size(), greenR.n_slices);

    //for each matsubara freq.
    const cd_t im = cd_t(0.0, 1.0);
    for (size_t nn = 0; nn < greenR.n_slices; nn++)
    {
        for (size_t RIndex = 0; RIndex < RSites.size(); RIndex++)
        {
            for (size_t KIndex = 0; KIndex < KWaveVectors.size(); KIndex++)
            {
                const double Rx = RSites.at(RIndex)(0);
                const double Ry = RSites.at(RIndex)(1);
                greenK(KIndex, nn) += std::exp(-im * dot(KWaveVectors.at(KIndex), RSites.at(RIndex))) * greenR(Rx, Ry, nn);
            }
        }
    }
    return greenK;
}

ClusterCubeCD_t KtoR(const DataK_t &greenK, const ClusterSites_t &RSites, const ClusterSites_t &KWaveVectors)
{
    assert(RSites.size() == KWaveVectors.size());
    ClusterCubeCD_t greenR(KWaveVectors.size(), RSites.size(), greenK.n_cols);
    greenR.zeros();

    //for each matsubara freq.
    const cd_t im = cd_t(0.0, 1.0);
    for (size_t nn = 0; nn < greenR.n_slices; nn++)
    {
        for (size_t RIndex = 0; RIndex < RSites.size(); RIndex++)
        {
            for (size_t KIndex = 0; KIndex < KWaveVectors.size(); KIndex++)
            {
                const double Rx = RSites.at(RIndex)(0);
                const double Ry = RSites.at(RIndex)(1);
                greenR(Rx, Ry, nn) += std::exp(im * dot(KWaveVectors.at(KIndex), RSites.at(RIndex))) * greenK(KIndex, nn);
            }
        }
    }
    return 1.0 / static_cast<double>(RSites.size()) * greenR;
}

std::valarray<cd_t> RtoK(const ClusterMatrixCD_t &greenR, const ClusterSites_t &RSites, const ClusterSites_t &KWaveVectors)
{
    assert(RSites.size() == KWaveVectors.size());
    std::valarray<cd_t> greenK(cd_t(0.0, 0.0), KWaveVectors.size());
    const cd_t im = cd_t(0.0, 1.0);
    for (size_t RIndex = 0; RIndex < RSites.size(); RIndex++)
    {
        for (size_t KIndex = 0; KIndex < KWaveVectors.size(); KIndex++)
        {
            const double Rx = RSites.at(RIndex)(0);
            const double Ry = RSites.at(RIndex)(1);
            greenK[KIndex] += std::exp(-im * dot(KWaveVectors.at(KIndex), RSites.at(RIndex))) * greenR(Rx, Ry);
        }
    }
    return greenK;
}

} // namespace FourierDCA