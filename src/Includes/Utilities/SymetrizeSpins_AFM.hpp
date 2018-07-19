#pragma once

#include "../Models/ModelSquare2x2_AFM.hpp"
#include "../Models/ModelTriangle2x2.hpp"
#include "Utilities.hpp"
#include "IO.hpp"
#include "Conventions.hpp"

namespace SymmetrizeSpins
{

//This should simply be get downEquivalentSites

std::vector<size_t> GetDownEquivalentSites(const std::string &modelType)
{

    if (modelType == "Square2x2_AFM")
    {

        return IO::IOSquare2x2_AFM().downEquivalentSites();
    }
    else if (modelType == "Triangle2x2_AFM")
    {
        return IO::IOTriangle2x2().downEquivalentSites();
    }

    throw std::runtime_error("Miseria: modelType error in params file. Stupido !");
}

void SymmetrizeUpAndDown(const Json &jj)
{

    //Here we read the greens.dat files as specified in conventions.hpp and average them accordingly to the symmetries of up and down spins
    // remembert that the 0th column is the frequencies and that we have the real and imaginary parts in adjacent columns.

    const std::string modelType = jj["modelType"].get<std::string>();
    const std::vector<size_t> downEquivalentSites = GetDownEquivalentSites(modelType);

    if (downEquivalentSites.empty())
    {
        return;
    }

    const Conventions::NameVector_t nameVecGreens = Conventions::BuildGreensVectorNames();

    for (size_t ii = 0; ii < nameVecGreens.size() - 1; ii++)
    {
        ClusterMatrix_t tmpUp;
        assert(tmpUp.load(nameVecGreens.at(ii)));
        ClusterMatrix_t tmpDown;

        assert(tmpDown.load(nameVecGreens.at(ii + 1)));
        assert(tmpDown.n_cols == tmpUp.n_cols);
        assert(tmpUp.n_cols == downEquivalentSites.size());

        //Remember, the first column is the matsubara axis, so start at 1 for the column index
        for (size_t ii = 0; ii < downEquivalentSites.size(); ii++)
        {
            const size_t upIndex = 1 + 2 * downEquivalentSites.at(ii);
            const size_t downIndex = 1 + 2 * ii;
            tmpUp.col(upIndex) += tmpDown.col(downIndex);
            tmpUp.col(upIndex + 1) += tmpDown.col(downIndex + 1);

            tmpUp.col(upIndex) /= 2.0;
            tmpUp.col(upIndex + 1) /= 2.0;

            tmpDown.col(downIndex) = tmpUp.col(upIndex);
            tmpDown.col(downIndex + 1) = tmpUp.col(upIndex + 1);
        }

        assert(tmpUp.save(nameVecGreens.at(ii), arma::raw_ascii));
        assert(tmpDown.save(nameVecGreens.at(ii + 1), arma::raw_ascii));
    }

    mpiUt::Print("After symmetrize Spins.");
}

} // namespace SymmetrizeSpins
