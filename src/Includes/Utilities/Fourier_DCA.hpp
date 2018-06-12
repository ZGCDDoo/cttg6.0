#pragma once

#include "Integrator.hpp"
#include "Utilities.hpp"
#include "MPIUtilities.hpp"
#include "GreenMat.hpp"
#include "ABC_SelfConsistency.hpp"

namespace FourierDCA
{
using DataK_t = ClusterMatrixCD_t; //first index = K, second index=iwn

DataK_t RtoK(const ClusterCubeCD_t &greenR, const ClusterSites_t &RSites, const ClusterSites &KWaveVectors)
{
    assert(RSites.size() == KWaveVectors.size());
    DataK_t greenK(KWaveVectors.size(), greenR.n_slices());

    //for each matsubara freq.
    const cd_t im = cd_t(0.0, 1.0);
    for (size_t nn = 0; nn < greenR.n_slices; n++)
    {
        for (size_t RIndex = 0; RIndex < RSites.size(); RIndex++)
        {
            for (size_t KIndex = 0; KIndex < KWaveVectors.size(); KIndex++)
            {
                const double Rx = RSites.at(RIndex)(0);
                const double Ry = RSites.at(RIndex)(1);
                greenK(KIndex, nn) += std::exp(-im * dot(K.at(KIndex), RSites_.at(rIndex))) * greenR(Rx, Ry, nn);
            }
        }
    }
    return greenK;
}

ClusterCubeCD_t KtoR(const DataK_t &greenK, const ClusterSites_t &RSites, const ClusterSites &KWaveVectors)
{
    assert(RSites.size() == KWaveVectors.size());
    ClusterCube_t greenR(KWaveVectors.size(), RSites.size(); greenK.at(0).size());
    greenR.zeros();

    //for each matsubara freq.
    const cd_t im = cd_t(0.0, 1.0);
    for (size_t nn = 0; nn < greenR.n_slices; n++)
    {
        for (size_t RIndex = 0; RIndex < RSites.size(); RIndex++)
        {
            for (size_t KIndex = 0; KIndex < KWaveVectors.size(); KIndex++)
            {
                const double Rx = RSites.at(RIndex)(0);
                const double Ry = RSites.at(RIndex)(1);
                greenR(Rx, Ry, nn) += std::exp(im * dot(K.at(KIndex), RSites_.at(rIndex))) * greenK(KIndex, nn);
            }
        }
    }
    return greenR;
}

} // namespace FourierDCA