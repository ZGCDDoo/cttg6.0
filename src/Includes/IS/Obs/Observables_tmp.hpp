#pragma once

#include "../../Utilities/Utilities.hpp"
#include "../../Utilities/LinAlg.hpp"
#include "../../Utilities/MPIUtilities.hpp"

#include "GreenBinning.hpp"
#include "FillingAndDocc.hpp"
#include "GreenTauMesure.hpp"
#include "../ISData.hpp"
#include "../ISResult.hpp"

namespace Markov
{
namespace Obs
{
using Matrix_t = LinAlg::Matrix<double>;

template <typename TIOModel, typename TModel>
class Observables
{

      public:
        static const size_t N_MC_SUSC = 1; //number of samples to take for Monte Carlo integration in susceptibilites calculations
        Observables(){};
        Observables(const std::shared_ptr<ISDataCT<TIOModel, TModel>> &dataCT,
                    const Json &jj) : modelPtr_(new TModel(jj)),
                                      ioModel_(TIOModel()),
                                      dataCT_(dataCT),
                                      rng_(jj["SEED"].get<size_t>() + mpiUt::Rank() * mpiUt::Rank()),
                                      urngPtr_(new Utilities::UniformRngFibonacci3217_t(rng_, Utilities::UniformDistribution_t(0.0, 1.0))),
                                      greenBinningUp_(modelPtr_, dataCT_, jj, FermionSpin_t::Up),
#ifdef AFM
                                      greenBinningDown_(modelPtr_, dataCT_, jj, FermionSpin_t::Down),
#endif
                                      fillingAndDocc_(dataCT_, urngPtr_, jj["N_T_INV"].get<size_t>()),
                                      greenTauMesure_(dataCT_, urngPtr_),
                                      isPrecise_(jj["PRECISE"].get<bool>()),
                                      signMeas_(0.0),
                                      expOrder_(0.0),
                                      knightshift_(0.0),
                                      NMeas_(0)
        {

                mpiUt::Print("In Obs constructor ");
                mpiUt::Print("After Obs  constructor ");
        }

        //Getters
        double signMeas() const { return signMeas_; };
        double expOrder() const { return expOrder; };

        void Measure()
        {

                // mpiUt::Print("start of Measure");

                NMeas_++;
                signMeas_ += static_cast<double>(dataCT_->sign_);
                expOrder_ += static_cast<double>(dataCT_->vertices_.size()) * static_cast<double>(dataCT_->sign_);

                fillingAndDocc_.MeasureFillingAndDocc();

                // if (isPrecise_)
                // {

                //         MeasureGreen();
                // }
                // else
                // {
#ifndef AFM
                Maveraged_ = 0.5 * (*(dataCT_->MupPtr_) + *(dataCT_->MdownPtr_));
                greenBinningUp_.MeasureGreenBinning(Maveraged_);
#else
                greenBinningUp_.MeasureGreenBinning(*dataCT_->MupPtr_);
                greenBinningDown_.MeasureGreenBinning(*dataCT_->MdownPtr_);
#endif
                // }

                MeasureKnightShift();

                // mpiUt::Print("End of Measure");
                // return;
        }

        void MeasureKnightShift()
        {
                // mpiUt::Print("start of measureKnightShift");

                // const double fillingTotalUp = fillingAndDocc_.fillingUpTotalCurrent();
                // const double fillingTotalDown = fillingAndDocc_.fillingDownTotalCurrent();

                //should i multiply by beta?
                // knightshift_ += dataCT_->beta() * (fillingTotalUp * fillingTotalUp + fillingTotalDown * fillingTotalDown - 2.0 * fillingTotalUp * fillingTotalDown);

                double sum_gf = 0.0;
                for (size_t ii = 0; ii < N_MC_SUSC; ii++)
                {
                        double tau = 3.9920000000; //(*urngPtr_)() * dataCT_->beta();
                        greenTauMesure_.MeasureGreen(tau);

                        // https://stackoverflow.com/questions/610245/where-and-why-do-i-have-to-put-the-template-and-typename-keywords
                        ClusterMatrix_t greenTauUp = ioModel_.template IndepToFull<SiteVector_t, ClusterMatrix_t>(greenTauMesure_.greenUpCurrent());
                        ClusterMatrix_t greenTauDown = ioModel_.template IndepToFull<SiteVector_t, ClusterMatrix_t>(greenTauMesure_.greenDownCurrent());

                        //greenTauMesure_.MeasureGreen(-tau);
                        ClusterMatrix_t greenMTauUp = ioModel_.template IndepToFull<SiteVector_t, ClusterMatrix_t>(greenTauMesure_.greenUpMCurrent());
                        ClusterMatrix_t greenMTauDown = ioModel_.template IndepToFull<SiteVector_t, ClusterMatrix_t>(greenTauMesure_.greenDownMCurrent());
                        //calculate sum_tau sum_ij\sigma g_ij\sigma(tau)g_ji\sigma(-tau)

                        for (size_t i = 0; i < ioModel_.Nc; i++)
                        {
                                for (size_t j = 0; j < ioModel_.Nc; j++)
                                {
                                        sum_gf += greenTauUp(i, j) * greenMTauUp(j, i) + greenTauDown(i, j) * greenMTauDown(j, i);
                                }
                        }
                }
                sum_gf *= dataCT_->beta() / static_cast<double>(N_MC_SUSC);
                knightshift_ -= sum_gf;

                // mpiUt::Print("End of measureKnightShift");
        }

        void Save()
        {
                mpiUt::Print("Start of Observables.Save()");
                signMeas_ /= NMeas_;

                fillingAndDocc_.Finalize(signMeas_, NMeas_);
                std::map<std::string, double> obsScal;

                obsScal = fillingAndDocc_.GetObs();

                obsScal["sign"] = signMeas_;
                obsScal["NMeas"] = NMeas_;

                //dont forget that the following obs have not been finalized (multiplied by following factor)
                const double fact = 1.0 / (NMeas_ * signMeas_);
                obsScal["k"] = fact * expOrder_;
                obsScal["knightshift"] = fact * knightshift_;

                //The greens
                // if (isPrecise_)
                // {
                //         FinalizeGreens();
                // }
                // else
                // {

                ClusterMatrixCD_t greenMatsubaraUp = ioModel_.FullCubeToIndep(greenBinningUp_.FinalizeGreenBinning(signMeas_, NMeas_));
#ifdef AFM
                ClusterMatrixCD_t greenMatsubaraDown = ioModel_.FullCubeToIndep(greenBinningDown_.FinalizeGreenBinning(signMeas_, NMeas_));
#endif
                // }

#ifndef AFM
                Result::ISResult isResult(obsScal, greenMatsubaraUp, greenMatsubaraUp, fillingAndDocc_.fillingUp(), fillingAndDocc_.fillingDown());
#else
                Result::ISResult isResult(obsScal, greenMatsubaraUp, greenMatsubaraDown, fillingAndDocc_.fillingUp(), fillingAndDocc_.fillingDown());
#endif
                std::vector<Result::ISResult> isResultVec;
#ifdef HAVEMPI

                mpi::communicator world;
                if (mpiUt::Rank() == mpiUt::master)
                {
                        mpi::gather(world, isResult, isResultVec, mpiUt::master);
                }
                else
                {
                        mpi::gather(world, isResult, mpiUt::master);
                }
                if (mpiUt::Rank() == mpiUt::master)
                {
                        mpiUt::IOResult<TIOModel>::SaveISResults(isResultVec, dataCT_->beta_);
                }

#else
                isResultVec.push_back(isResult);
                mpiUt::IOResult<TIOModel>::SaveISResults(isResultVec, dataCT_->beta_);
#endif

                mpiUt::Print("End of Observables.Save()");
                return;
        }

      private:
        std::shared_ptr<TModel> modelPtr_;
        TIOModel ioModel_;
        std::shared_ptr<ISDataCT<TIOModel, TModel>> dataCT_;
        Utilities::EngineTypeFibonacci3217_t rng_;
        std::shared_ptr<Utilities::UniformRngFibonacci3217_t> urngPtr_;

        GreenBinning<TIOModel, TModel> greenBinningUp_;
#ifdef AFM
        GreenBinning<TIOModel, TModel> greenBinningDown_;
#endif
        FillingAndDocc<TIOModel, TModel> fillingAndDocc_;
        GreenTauMesure<TIOModel, TModel> greenTauMesure_;

        Matrix_t Maveraged_;
        const bool isPrecise_;

        //=======Measured quantities
        double signMeas_;
        double expOrder_;
        double knightshift_;
        size_t NMeas_;
};

} // namespace Obs
} // namespace Markov
