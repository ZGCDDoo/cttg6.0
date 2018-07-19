#pragma once

#include "../Models/ModelSquare2x2_AFM.hpp"
#include "../Models/ModelTriangle2x2.hpp"
#include "Utilities.hpp"
#include "IO.hpp"
#include "Conventions.hpp"

namespace SymmetrizeSpins
{

std::unique_ptr<IO::Base_IOModel> GetBaseIOPtr(const std::string &modelType)
{

    const std::string modelType = jj["modelType"].get<std::string>();

    if (modelType == "Square2x2_AFM")
    {

        using IOModel_t = IO::IOSquare2x2_AFM;
        return std::make_unique<IOModel_t>()
    }
    else if (modelType == "Triangle2x2_AFM")
    {
        using IOModel_t = IO::IOTriangle2x2;
        return std::make_unique<IOModel_t>()
    }

    throw std::runtime_error("Miseria: modelType error in params file. Stupido !");
}

void SymmetrizeUpAndDown(const Json &jj)
{

    //Here we read the greens.dat files as specified in conventions.hpp and average them accordingly to the symmetries of up and down spins
    // remembert that the 0th column is the frequencies and that we have the real and imaginary parts in adjacent columns.

        const std::string modelType = jj["modelType"].get<std::string>();
    std::unique_ptr<IO::Base_IOModel> ioModelPtr = GetBaseIOPtr(modelType);

    if (ioModelPtr->downEquivalentSites().empty())
    {
        return;
    }

    const Conventions::NameVector_t nameVecGreens = Conventions::BuildGreensVectorNames();

    for (size_t ii = 0; ii < nameVecGreens.size() - 1; ii++)
    {
        const ClusterMatrixCD_t tmpUp = ioModelPtr->ReadGreenDat(nameVecGreens.at(ii));
        const ClusterMatrixCD_t tmpDown = ioModelPtr->ReadGreenDat(nameVecGreens.at(ii + 1));

        //Remember, the first column is the matsubara axis, so start at 1 for the column index
        for (size_t ii = 0; ii < ioModel_.downEquivalentSites().size(); ii++)
        {
            const size_t upIndex = ioModel_.downEquivalentSites().at(ii);
            tmpUp.col(upIndex) += tmpDown.col(ii);
            tmpUp.col(upIndex) /= 2.0;
            tmpDown.col(ii) = tmpUp.col(upIndex);
        }
    }
    mpiUt::Print("After symmetrize Spins.");
}

} // namespace SymmetrizeSpins
