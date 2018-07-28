#pragma once
#include "../Utilities/Utilities.hpp"
#include "../Utilities/Fourier_DCA.hpp"

namespace Models
{
template <size_t TNX, size_t TNY>
class ABC_H0
{

  public:
    static const size_t Nc = TNX * TNY;
    static const size_t n_cols = Nc;
    static const size_t n_rows = Nc;
    static const size_t Nx = TNX;
    static const size_t Ny = TNY;
    const size_t NKPTS = 1000;

    ABC_H0(const double &t, const double &tp, const double &tpp) : RSites_(Nc),
                                                                   KWaveVectors_(Nc),
                                                                   t_(t),
                                                                   tPrime_(tp),
                                                                   tPrimePrime_(tpp)

    {
        std::cout << "End of ABC_H0 Constructor " << std::endl;
        assert(TNX == TNY);

        for (size_t i = 0; i < TNX; i++)
        {
            for (size_t j = 0; j < TNY; j++)
            {
                const size_t index = i + TNY * j;
                RSites_.at(index) = {static_cast<double>(i), static_cast<double>(j)};
                KWaveVectors_.at(index) = {static_cast<double>(i) * 2.0 * M_PI / static_cast<double>(TNX), static_cast<double>(j) * 2.0 * M_PI / static_cast<double>(TNY)};
                KWaveVectors_.at(index).print();
                std::cout << "\n";
            }
        }

        assert(KWaveVectors_.size() == Nc);
        assert(KWaveVectors_.size() == RSites_.size());

        std::cout << "End of ABC_H0 Constructor " << std::endl;
    }

    virtual ~ABC_H0() = 0;

    double t() const { return t_; };
    double tPrime() const { return tPrime_; };
    double tPrimePrime() const { return tPrimePrime_; };
    ClusterSites_t RSites() const { return RSites_; };
    ClusterSites_t KWaveVectors() const { return KWaveVectors_; };

    virtual double Eps0k(const double &kx, const double &ky) const = 0;

    void SaveTKTildeAndHybFM()
    {
        //check if  file exists:
        using boost::filesystem::exists;
        if (exists("tloc.arma") && exists("hybFM.arma"))
        {
            ClusterMatrixCD_t tmp;
            tmp.load("tloc.arma");
            if (tmp.n_cols == Nc)
            {
                return;
            }
        }
        std::cout << "Calculating  tloc and hybFM. " << std::endl;

        ClusterMatrixCD_t hybFM(Nc, Nc);
        ClusterMatrixCD_t epsKBar(Nc, Nc);

        hybFM.zeros();
        epsKBar.zeros();

        assert(KWaveVectors_.size() == epsKBar.n_cols);

        for (size_t Kindex = 0; Kindex < KWaveVectors_.size(); Kindex++)
        {
            const double Kx = KWaveVectors_.at(Kindex)(0);
            const double Ky = KWaveVectors_.at(Kindex)(1);
            for (size_t kx = 1; kx < NKPTS; kx++)
            {
                const double kTildeX = -M_PI / static_cast<double>(Nx) + static_cast<double>(kx) / static_cast<double>(NKPTS) * 2.0 * M_PI / static_cast<double>(Nx);
                for (size_t ky = 1; ky < NKPTS; ky++)
                {
                    const double kTildeY = -M_PI / static_cast<double>(Ny) + static_cast<double>(ky) / static_cast<double>(NKPTS) * 2.0 * M_PI / static_cast<double>(Ny);
                    const double tmp = Eps0k(Kx + kTildeX, Ky + kTildeY);
                    epsKBar(Kindex, Kindex) += tmp;
                    hybFM(Kindex, Kindex) += tmp * tmp;
                }
            }
        }

        const size_t NKPTS_Squared = (NKPTS - 1) * (NKPTS - 1);
        epsKBar /= static_cast<double>(NKPTS_Squared);
        hybFM /= static_cast<double>(NKPTS_Squared);
        hybFM -= epsKBar * epsKBar;

        hybFM.save("hybFM.arma", arma::arma_ascii);
        epsKBar.save("tloc.arma", arma::arma_ascii);

        const ClusterMatrixCD_t hybFM_R = FourierDCA::KtoR(hybFM, RSites_, KWaveVectors_);
        const ClusterMatrixCD_t epsKBar_R = FourierDCA::KtoR(epsKBar, RSites_, KWaveVectors_);
        hybFM_R.save("hybFM_R.arma", arma::arma_ascii);
        epsKBar_R.save("tloc_R.arma", arma::arma_ascii);

        std::cout << "End of  SaveEpsKBarAndHybKFM" << std::endl;
    }

  protected:
    ClusterSites_t RSites_;
    ClusterSites_t KWaveVectors_;

    const double t_;
    const double tPrime_;
    const double tPrimePrime_;
};

template <size_t TNX, size_t TNY>
ABC_H0<TNX, TNY>::~ABC_H0() {} //destructors must exist
} // namespace Models