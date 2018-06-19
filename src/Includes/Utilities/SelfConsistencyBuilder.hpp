#pragma once

#ifdef DCA
#include "SelfConsistency_DCA.hpp"
#else
#include "SelfConsistency.hpp"
#endif

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
    else if (modelType == "Square2x2")
    {
        const size_t Nx = 2;
        using Model_t = Models::ModelSquare2x2;
        using IOModel_t = IO::IOSquare2x2;
        using H0_t = Models::H0Square<Nx, Nx>;

        Model_t model(jj);
        IOModel_t ioModel;
        const ClusterCubeCD_t greenImpurityUp = ioModel.ReadGreenDat("greenUp.dat");

        using SelfCon_t = SelfCon::SelfConsistency<IOModel_t, Model_t, H0_t>;
        return std::make_unique<SelfCon_t>(SelfCon_t(jj, model, greenImpurityUp, "Up"));
    }
    else if (modelType == "Triangle2x2")
    {
        const size_t Nx = 2;
        using Model_t = Models::ModelTriangle2x2;
        using IOModel_t = IO::IOTriangle2x2;
        using H0_t = Models::H0Square<Nx, Nx>;

        Model_t model(jj);
        IOModel_t ioModel;
        const ClusterCubeCD_t greenImpurityUp = ioModel.ReadGreenDat("greenUp.dat");

        using SelfCon_t = SelfCon::SelfConsistency<IOModel_t, Model_t, H0_t>;
        return std::make_unique<SelfCon_t>(SelfCon_t(jj, model, greenImpurityUp, "Up"));
    }
    else if (modelType == "Square4x4")
    {
        const size_t Nx = 4;
        using Model_t = Models::ModelSquare4x4;
        using IOModel_t = IO::IOSquare4x4;
        using H0_t = Models::H0Square<Nx, Nx>;

        Model_t model(jj);
        IOModel_t ioModel;
        const ClusterCubeCD_t greenImpurityUp = ioModel.ReadGreenDat("greenUp.dat");

        using SelfCon_t = SelfCon::SelfConsistency<IOModel_t, Model_t, H0_t>;
        return std::make_unique<SelfCon_t>(SelfCon_t(jj, model, greenImpurityUp, "Up"));
    }
    else if (modelType == "Square6x6")
    {
        const size_t Nx = 6;
        using Model_t = Models::ModelSquare6x6;
        using IOModel_t = IO::IOSquare6x6;
        using H0_t = Models::H0Square<Nx, Nx>;

        Model_t model(jj);
        IOModel_t ioModel;
        const ClusterCubeCD_t greenImpurityUp = ioModel.ReadGreenDat("greenUp.dat");

        using SelfCon_t = SelfCon::SelfConsistency<IOModel_t, Model_t, H0_t>;
        return std::make_unique<SelfCon_t>(SelfCon_t(jj, model, greenImpurityUp, "Up"));
    }
    else if (modelType == "Square8x8")
    {
        const size_t Nx = 8;
        using Model_t = Models::ModelSquare8x8;
        using IOModel_t = IO::IOSquare8x8;
        using H0_t = Models::H0Square<Nx, Nx>;

        Model_t model(jj);
        IOModel_t ioModel;
        const ClusterCubeCD_t greenImpurityUp = ioModel.ReadGreenDat("greenUp.dat");

        using SelfCon_t = SelfCon::SelfConsistency<IOModel_t, Model_t, H0_t>;
        return std::make_unique<SelfCon_t>(SelfCon_t(jj, model, greenImpurityUp, "Up"));
    }

    return NULL;
} // namespace SelfCon

} // namespace SelfCon