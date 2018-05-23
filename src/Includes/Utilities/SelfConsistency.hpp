#pragma once

#include "Integrator.hpp"
#include "Utilities.hpp"
#include "GreenMat.hpp"

namespace SelfCon
{

template <typename TH0>
struct GreenLattice
{

  public:
    static const size_t Nc = TH0::Nc;
    static const ClusterMatrixCD_t II;
    static const size_t n_cols = TH0::n_cols;
    static const size_t n_rows = TH0::n_rows;
    static const size_t Nx = TH0::Nx;
    static const size_t Ny = TH0::Ny;

    GreenLattice(cd_t zz, ClusterMatrixCD_t selfEnergy, TH0 h0) : zz_(zz), selfEnergy_(selfEnergy), h0_(h0){};

    ClusterMatrixCD_t operator()(const double &kx, const double &ky)
    {
        return ((zz_ * II - h0_(kx, ky) - selfEnergy_).i());
    }

  private:
    const cd_t zz_;
    ClusterMatrixCD_t selfEnergy_;
    TH0 h0_;
};
template <typename TH0>
const ClusterMatrixCD_t GreenLattice<TH0>::II = ClusterMatrixCD_t(TH0::Nc, TH0::Nc).eye();

template <typename TIOModel, typename TModel, typename TH0>
class SelfConsistency
{

  public:
    static const size_t Nc = TModel::Nc;
    static const ClusterMatrixCD_t II;
    const double factNSelfCon = 2.0;

    SelfConsistency(const Json &jj, const TModel &model, const ClusterCubeCD_t &greenImpurity, const std::string &spinName) : model_(model),
                                                                                                                              ioModel_(TIOModel()),
                                                                                                                              greenImpurity_(greenImpurity),
                                                                                                                              hybridization_(model_.hybridizationMatUp()),
                                                                                                                              selfEnergy_(),
                                                                                                                              hybNext_(),
                                                                                                                              spinName_(spinName),
                                                                                                                              weights_(cd_t(jj["WEIGHTSR"].get<double>(), jj["WEIGHTSI"].get<double>()))
    {

        const size_t NGreen = greenImpurity_.n_slices;
        size_t NSelfConTmp = std::max<double>(0.5 * (jj["ESelfCon"].get<double>() * model_.beta() / M_PI - 1.0),
                                              0.5 * (200.0 * model_.beta() / M_PI - 1.0));
        if (NGreen >= NSelfConTmp)
        {
            NSelfConTmp = factNSelfCon * static_cast<double>(NGreen);
        }
        const size_t NSelfCon = NSelfConTmp;
        assert(NSelfCon > NGreen);
        //Patcher la hyb si necessaire
        hybridization_.PatchHF(NSelfCon, model_.beta());
        const size_t NHyb = hybridization_.n_slices();
        assert(NHyb >= NSelfCon);

        selfEnergy_.resize(Nc, Nc, NSelfCon);

        //0.) Extraire la self jusqu'a NGreen
        for (size_t nn = 0; nn < NGreen; nn++)
        {
            cd_t zz = cd_t(model_.mu(), (2.0 * nn + 1.0) * M_PI / model_.beta());
            selfEnergy_.slice(nn) = -greenImpurity_.slice(nn).i() + zz * ClusterMatrixCD_t(Nc, Nc).eye() - model_.tLoc() - hybridization_.slice(nn);
        }

        //1.) Patcher la self par HF de NGreen à NSelfCon
        ClusterMatrix_t nUpMatrix;
        assert(nUpMatrix.load("nUpMatrix.dat"));
        ClusterMatrix_t nDownMatrix;
        assert(nDownMatrix.load("nDownMatrix.dat"));
        ClusterMatrixCD_t nMatrix(nUpMatrix + nDownMatrix, ClusterMatrix_t(Nc, Nc).zeros());

        for (size_t nn = NGreen; nn < NSelfCon; nn++)
        {
            cd_t iwn = cd_t(0.0, (2.0 * nn + 1.0) * M_PI / model_.beta());
            selfEnergy_.slice(nn) = 0.5 * model_.U() * nMatrix + 1.0 / iwn * model_.U() * model_.U() * nMatrix / 2.0 * (II - nMatrix / 2.0);
        }

        std::cout << "In Selfonsistency constructor " << std::endl;
        Save("self" + spinName_, selfEnergy_);
        std::cout << "In Selfonsistency constructor, after save selfenery " << std::endl;
    }

    void DoSCGrid()
    {

        std::cout << "In Selfonsistency DOSC " << std::endl;
        const size_t NSelfCon = selfEnergy_.n_slices;
        ClusterCubeCD_t gImpUpNext(Nc, Nc, NSelfCon);
        gImpUpNext.zeros();
        hybNext_.resize(Nc, Nc, NSelfCon);
        hybNext_.zeros();
        ClusterCubeCD_t tKTildeGrid;
        assert(tKTildeGrid.load("tktilde.arma", arma::arma_ascii));
        size_t ktildepts = tKTildeGrid.n_slices;

        for (size_t nn = 0; nn < NSelfCon; nn++)
        {
            cd_t zz = cd_t(model_.mu(), (2.0 * nn + 1.0) * M_PI / model_.beta());
            for (size_t ktildeindex = 0; ktildeindex < ktildepts; ktildeindex++)
            {
                gImpUpNext.slice(nn) += 1.0 / (static_cast<double>(ktildepts)) * ((zz * ClusterMatrixCD_t(Nc, Nc).eye() - tKTildeGrid.slice(ktildeindex) - selfEnergy_.slice(nn)).i());
            }
            hybNext_.slice(nn) = -gImpUpNext.slice(nn).i() - selfEnergy_.slice(nn) + zz * ClusterMatrixCD_t(Nc, Nc).eye() - model_.tLoc();
        }

        hybNext_ *= (1.0 - weights_);
        hybNext_ += weights_ * hybridization_.data();
        Save("green" + spinName_, gImpUpNext);
        Save("hybNext" + spinName_, hybNext_);

        std::cout << "After Selfonsistency DOSC " << std::endl;
        return;
    }

    void DoSC()
    {
        //std::cout << "In Selfonsistency DOSC " << std::endl;
        const size_t NMat = greenImpurity_.n_slices;
        cd_t zz;
        ClusterCubeCD_t gImpUpPrime(Nc, Nc, NMat);
        ClusterCubeCD_t hybridizationNext(Nc, Nc, NMat);

        //std::cout << "In Selfonsistency DOSC before loop " << std::endl;
        for (size_t nn = 0; nn < NMat; nn++)
        {
            zz = cd_t(model_.mu(), (2.0 * nn + 1.0) * M_PI / model_.beta());
            GreenLattice<TH0> glattice(zz, selfEnergy_.slice(nn), model_.h0());
            //std::cout << "In Selfonsistency DOSC in  loop  after glattice" << std::endl;
            gImpUpPrime.slice(nn) = Integrator::GridKTilde(glattice);
            hybridizationNext.slice(nn) = -gImpUpPrime.slice(nn).i() - selfEnergy_.slice(nn) + zz * ClusterMatrixCD_t(Nc, Nc).eye() - model_.tLoc();
            //std::cout << "In Selfonsistency DOSC in end of loop " << std::endl;
        }

        Save("hybNext" + spinName_, hybridizationNext, true);
        hybNext_ = hybridizationNext;
        return;
    }

    void Save(std::string fname, ClusterCubeCD_t green, bool saveArma = false)
    {
        const size_t NMat = green.n_slices;
        ClusterMatrixCD_t greenOut(NMat, ioModel_.indepSites().size());

        std::ofstream fout;
        fout.open(fname + std::string(".dat"), std::ios::out);
        double iwn;
        for (size_t nn = 0; nn < green.n_slices; nn++)
        {
            iwn = (2.0 * nn + 1.0) * M_PI / model_.beta();
            fout << iwn << " ";

            for (Site_t ii = 0; ii < ioModel_.indepSites().size(); ii++)
            {
                Site_t s1 = ioModel_.indepSites().at(ii).first;
                Site_t s2 = ioModel_.indepSites().at(ii).second;

                greenOut(nn, ii) = green(s1, s2, nn);
                fout << green(s1, s2, nn).real()
                     << " "
                     << green(s1, s2, nn).imag()
                     << " ";
            }
            fout << "\n";
        }

        fout.close();
        if (saveArma)
        {
            greenOut.save(fname + std::string(".arma"), arma::arma_ascii);
        }
        return;
    }

    ClusterCubeCD_t
    hybNext() const
    {
        return hybNext_;
    };

  private:
    TModel model_;
    TIOModel ioModel_;

    const ClusterCubeCD_t greenImpurity_;
    GreenMat::HybridizationMat hybridization_;
    ClusterCubeCD_t selfEnergy_;
    ClusterCubeCD_t hybNext_;
    const std::string spinName_;
    const cd_t weights_;
};
template <typename TIOModel, typename TModel, typename TH0>
const ClusterMatrixCD_t SelfConsistency<TIOModel, TModel, TH0>::II = ClusterMatrixCD_t(TH0::Nc, TH0::Nc).eye();
}