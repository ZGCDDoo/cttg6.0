#pragma once

#include "Integrator.hpp"
#include "Utilities.hpp"
#include "MPIUtilities.hpp"
#include "GreenMat.hpp"
#include "ABC_SelfConsistency.hpp"
#include "Fourier_DCA.hpp"

namespace SelfCon
{

template <typename TH0>
struct GreenLattice
{

  public:
    static const size_t Nc;
    static const ClusterMatrixCD_t II;

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

template <typename TH0>
const size_t GreenLattice<TH0>::Nc = TH0::Nc;

template <typename TIOModel, typename TModel, typename TH0>
class SelfConsistency : public ABC_SelfConsistency
{

  public:
    static const size_t Nc;
    static const ClusterMatrixCD_t II;
    static const double factNSelfCon;

    SelfConsistency(const Json &jj, const TModel &model, const ClusterCubeCD_t &greenImpurity, const std::string &spinName) : model_(model),
                                                                                                                              ioModel_(TIOModel()),
                                                                                                                              h0_(jj["t"].get<double>(), jj["tPrime"].get<double>(), jj["tPrimePrime"].get<double>()),
                                                                                                                              greenImpurity_(greenImpurity),
                                                                                                                              hybridization_(model_.hybridizationMatUp()),
                                                                                                                              selfEnergy_(),
                                                                                                                              spinName_(spinName),
                                                                                                                              weights_(cd_t(jj["WEIGHTSR"].get<double>(), jj["WEIGHTSI"].get<double>()))
    {

        mpiUt::Print("Start of SC constructor");
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

        selfEnergyK_ = FourierDCA::RtoK(selfEnergy_, h0_.RSites(), h0_.KWaveVectors());

        if (mpiUt::Rank() == mpiUt::master)
        {
            Save("self" + spinName_, selfEnergy_);
            SaveK("self" + spinName + "_K", selfEnergyK_);
            std::cout << "In Selfonsistency constructor, after save selfenery " << std::endl;
        }

        mpiUt::Print("After SC constructor");
    }

    void DoSCGrid() override
    {

        DoSCGridSerial();
    }

    void DoSCGridSerial()
    {

        if (mpiUt::Rank() == mpiUt::master)
        {
            std::cout << "In Selfonsistency DOSC serial" << std::endl;
            const size_t NSelfCon = selfEnergyK_.n_cols;
            const size_t NKPTS = 200;
            ClusterMatrixCD_t gImpUpKNext(h0_.KWaveVectors().size(), NSelfCon);
            gImpUpKNext.zeros();
            hybKNext_ = gImpUpKNext;

            const double kxCenter = M_PI / static_cast<double>(h0_.Nx);
            const double kyCenter = M_PI / static_cast<double>(h0_.Ny);
            const std::valarray<cd_t> tK = FourierDCA::RtoK(model_.tLoc(), h0_.RSites(), h0_.KWaveVectors());

            for (size_t KIndex = 0; KIndex < h0_.KWaveVectors().size(); KIndex++)
            {
                // std::cout << "Here1 in selfonc" << std::endl;
                const double Kx = h0_.KWaveVectors().at(KIndex)(0);
                const double Ky = h0_.KWaveVectors().at(KIndex)(1);
                // std::cout << "Here2 in selfonc" << std::endl;
                for (size_t nn = 0; nn < NSelfCon; nn++)
                {
                    const cd_t zz = cd_t(model_.mu(), (2.0 * nn + 1.0) * M_PI / model_.beta());
                    for (size_t kxindex = 0; kxindex < NKPTS; kxindex++)
                    {
                        const double kx = (Kx - kxCenter) + static_cast<double>(kxindex) / static_cast<double>(NKPTS) * 2.0 * kxCenter;
                        // std::cout << "Here3 in selfonc" << std::endl;
                        for (size_t kyindex = 0; kyindex < NKPTS; kyindex++)
                        {
                            const double ky = (Ky - kyCenter) + static_cast<double>(kyindex) / static_cast<double>(NKPTS) * 2.0 * kyCenter;
                            gImpUpKNext(KIndex, nn) += 1.0 / (zz - h0_.Eps0k(kx, ky) - selfEnergyK_(KIndex, nn));
                            // std::cout << "Here4 in selfonc" << std::endl;
                        }
                        // std::cout << "Here5 in selfonc" << std::endl;
                    }
                    std::cout << "Here6 in selfonc" << std::endl;
                    gImpUpKNext(KIndex, nn) *= static_cast<double>(Nc) / static_cast<double>(NKPTS * NKPTS);
                    hybKNext_(KIndex, nn) = -1.0 / gImpUpKNext(KIndex, nn) - selfEnergyK_(KIndex, nn) + zz - tK[KIndex];
                    std::cout << "Here7 in selfonc" << std::endl;
                }
            }

            //     hybNextK_ *= (1.0 - weights_);
            //     hybNextK_ += weights_ * hybridization_.data();

            const ClusterCubeCD_t greenRNext = FourierDCA::KtoR(gImpUpKNext, h0_.RSites(), h0_.KWaveVectors());
            const ClusterCubeCD_t hybRNext = FourierDCA::KtoR(hybKNext_, h0_.RSites(), h0_.KWaveVectors());
            SaveK("green" + spinName_ + "_K", gImpUpKNext);
            Save("green" + spinName_, greenRNext);

            SaveK("hybNext" + spinName_ + "_K", hybKNext_);
            Save("hybNext" + spinName_, hybRNext);

            std::cout << "After Selfonsistency DOSC serial" << std::endl;
        }
    }

    void Save(std::string fname, ClusterCubeCD_t green, bool saveArma = false)
    {
        const size_t NMat = green.n_slices;
        ClusterMatrixCD_t greenOut(NMat, ioModel_.indepSites().size());

        std::ofstream fout;
        fout.open(fname + std::string(".dat"), std::ios::out);
        for (size_t nn = 0; nn < green.n_slices; nn++)
        {
            const double iwn = (2.0 * nn + 1.0) * M_PI / model_.beta();
            fout << iwn << " ";

            for (Site_t ii = 0; ii < ioModel_.indepSites().size(); ii++)
            {
                const Site_t s1 = ioModel_.indepSites().at(ii).first;
                const Site_t s2 = ioModel_.indepSites().at(ii).second;

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
    }

    void SaveK(std::string fname, FourierDCA::DataK_t greenK)
    {

        std::ofstream fout;
        fout.open(fname + std::string(".dat"), std::ios::out);

        for (size_t nn = 0; nn < greenK.n_cols; nn++)
        {
            const double iwn = (2.0 * nn + 1.0) * M_PI / model_.beta();
            fout << iwn << " ";

            for (size_t KIndex = 0; KIndex < h0_.KWaveVectors().size(); KIndex++)
            {

                fout << greenK(KIndex, nn).real()
                     << " "
                     << greenK(KIndex, nn).imag()
                     << " ";
            }
            fout << "\n";
        }

        fout.close();
    }

    // ClusterCubeCD_t
    //     // hybNext() const
    //     // {
    //     //     return hybNext_;
    //     // };

  private:
    TModel model_;
    TIOModel ioModel_;
    TH0 h0_;

    const ClusterCubeCD_t greenImpurity_;
    GreenMat::HybridizationMat hybridization_;

    ClusterCubeCD_t selfEnergy_;
    FourierDCA::DataK_t selfEnergyK_;
    FourierDCA::DataK_t hybKNext_;

    const std::string spinName_;
    const cd_t weights_;
};
template <typename TIOModel, typename TModel, typename TH0>
const ClusterMatrixCD_t SelfConsistency<TIOModel, TModel, TH0>::II = ClusterMatrixCD_t(TH0::Nc, TH0::Nc).eye();

template <typename TIOModel, typename TModel, typename TH0>
const size_t SelfConsistency<TIOModel, TModel, TH0>::Nc = TH0::Nc;

template <typename TIOModel, typename TModel, typename TH0>
const double SelfConsistency<TIOModel, TModel, TH0>::factNSelfCon = 2;

} // namespace SelfCon