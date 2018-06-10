#pragma once

#include "SelfConsistency.hpp"

namespace SelfCon
{

std::unique_ptr<ABC_SelfConsistency> SelfConsistencyBuilder(const Json &jj)
{
    const std::string modelType = jj["modelType"].get<std::string>();

    if (modelType == "SIAM_Square")
    {
        const size_t Nx = 1;
        using Model_t = Models::SIAM_Square;
        using IOModel_t = IO::IOSIAM;
        using H0_t = Models::H0Square<Nx, Nx>;

        Model_t model(jj);
        IOModel_t ioModel;
        const ClusterCubeCD_t greenImpurityUp = ioModel.ReadGreenDat("greenUp.dat");

        using SelfCon_t = SelfCon::SelfConsistency<IOModel_t, Model_t, H0_t>;
        return std::make_unique<SelfCon_t>(SelfCon_t(jj, model, greenImpurityUp, "Up"));
    }
    // else if (modelType == "Square2x2")
    // {
    //     using Model_t = Models::ModelSquare2x2;
    //     using IOModel_t = IO::IOSquare2x2;
    //     using MarkovInt_t = Markov::MarkovChain<IOModel_t, Model_t>;
    //     using MarkovAux_t = Markov::MarkovChainAux<IOModel_t, Model_t>;

    //     }
    //     else
    //     {
    //         return std::make_unique<MC::MonteCarlo<MarkovAux_t>>(std::make_shared<MarkovAux_t>(jj, seed), jj);
    //     }
    // }
    // else if (modelType == "Triangle2x2")
    // {
    //     using Model_t = Models::ModelTriangle2x2;
    //     using IOModel_t = IO::IOTriangle2x2;
    //     using MarkovInt_t = Markov::MarkovChain<IOModel_t, Model_t>;
    //     using MarkovAux_t = Markov::MarkovChainAux<IOModel_t, Model_t>;

    //     if (solverType == Int)
    //     {
    //         return std::make_unique<MC::MonteCarlo<MarkovInt_t>>(std::make_shared<MarkovInt_t>(jj, seed), jj);
    //     }
    //     else
    //     {
    //         return std::make_unique<MC::MonteCarlo<MarkovAux_t>>(std::make_shared<MarkovAux_t>(jj, seed), jj);
    //     }
    // }
    // else if (modelType == "Square4x4")
    // {
    //     using Model_t = Models::ModelSquare4x4;
    //     using IOModel_t = IO::IOSquare4x4;

    //     using MarkovInt_t = Markov::MarkovChain<IOModel_t, Model_t>;
    //     using MarkovAux_t = Markov::MarkovChainAux<IOModel_t, Model_t>;

    //     if (solverType == Int)
    //     {
    //         return std::make_unique<MC::MonteCarlo<MarkovInt_t>>(std::make_shared<MarkovInt_t>(jj, seed), jj);
    //     }
    //     else
    //     {
    //         return std::make_unique<MC::MonteCarlo<MarkovAux_t>>(std::make_shared<MarkovAux_t>(jj, seed), jj);
    //     }
    // }

    return NULL;
}

} // namespace SelfCon