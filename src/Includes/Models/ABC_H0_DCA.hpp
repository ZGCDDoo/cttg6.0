#pragma once
#include "../Utilities/Utilities.hpp"

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
    const size_t NKPTS = 200;

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
                size_t index = i + TNY * j;
                RSites_.at(index) = {static_cast<double>(i), static_cast<double>(j)};
                KWaveVectors_.at(index) = {static_cast<double>(i) * 2.0 * M_PI / static_cast<double>(TNX), static_cast<double>(j) * 2.0 * M_PI / static_cast<double>(TNY)};
            }
        }

        std::cout << "End of ABC_H0 Constructor " << std::endl;
    }

    virtual ~ABC_H0() = 0;

    double t() const { return t_; };
    double tPrime() const { return tPrime_; };
    double tPrimePrime() const { return tPrimePrime_; };
    ClusterSites_t RSites() const { return RSites_; };
    ClusterSites_t KWaveVectors() const { return KWaveVectors_; };

    virtual double Eps0k(const double &kx, const double &ky) const = 0;

    ClusterMatrixCD_t operator()(const double &kTildeX, const double &kTildeY) //return t(ktilde)
    {
        const cd_t im = cd_t(0.0, 1.0);
        arma::vec ktilde = {kTildeX, kTildeY};
        ClusterMatrixCD_t HoppingKTilde(Nc, Nc);
        HoppingKTilde.zeros();

        for (size_t i = 0; i < Nc; i++)
        {
            for (size_t j = 0; j < Nc; j++)
            {
                for (const auto &K : this->KWaveVectors_)
                {
                    HoppingKTilde(i, j) += std::exp(im * dot(K, RSites_.at(i) - RSites_[j])) * Eps0k(K(0) + kTildeX, K(1) + kTildeY);
                }
            }
        }
        return (HoppingKTilde / static_cast<double>(Nc));
    }

    void SaveTKTildeAndHybFM(const size_t &kxpts)
    {
        std::cout << "Start of ABC_H0 SaveTKTildeAndHybFM " << std::endl;
        //0.) calculate tLoc
        CalculateEpsKBar(NKPTS);

        std::cout << "Here1 " << std::endl;
        ClusterMatrixCD_t tLoc(Nc, Nc);
        tLoc.zeros();
        std::cout << "Here2 " << std::endl;
        assert(KWaveVectors_.size() == epsKBar_.size());

        const cd_t im = cd_t(0.0, 1.0);
        for (size_t ii = 0; ii < Nc; ii++)
        {
            for (size_t jj = 0; jj < Nc; jj++)
            {
                for (size_t Kindex = 0; Kindex < KWaveVectors_.size(); Kindex++)
                {
                    tLoc(ii, jj) += std::exp(im * dot(KWaveVectors_.at(Kindex), RSites_.at(ii) - RSites_[jj])) * epsKBar_[Kindex];
                }
            }
        }

        std::cout << "Here3 " << std::endl;
        tLoc /= static_cast<double>(Nc);
        tLoc.save("tloc.arma", arma::arma_ascii);

        //1.) calculate tKTildeGrid
        ClusterCubeCD_t tKTildeGrid(Nc, Nc, kxpts * kxpts);
        tKTildeGrid.zeros();
        size_t sliceindex = 0;

        std::cout << "Here4 " << std::endl;
        for (size_t kx = 0; kx < kxpts; kx++)
        {
            const double kTildeX = static_cast<double>(kx) / static_cast<double>(kxpts) * 2.0 * M_PI / static_cast<double>(Nx);
            for (size_t ky = 0; ky < kxpts; ky++)
            {
                const double kTildeY = static_cast<double>(ky) / static_cast<double>(kxpts) * 2.0 * M_PI / static_cast<double>(Nx);
                tKTildeGrid.slice(sliceindex) = (*this)(kTildeX, kTildeY);
                sliceindex++;
            }
        }

        tKTildeGrid.save("tktilde.arma", arma::arma_ascii);

        //2.) calculate First moment of hyb
        ClusterMatrixCD_t hybFM(Nc, Nc);
        hybFM.zeros();

        const size_t Nkpts = tKTildeGrid.n_slices;
        for (size_t nn = 0; nn < Nkpts; nn++)
        {
            hybFM += tKTildeGrid.slice(nn) * tKTildeGrid.slice(nn);
        }
        hybFM /= Nkpts;
        hybFM -= tLoc * tLoc;
        hybFM.save("hybFM.arma", arma::arma_ascii);

        std::cout << "End of ABC_H0 SaveTKTildeAndHybFM " << std::endl;
    }

    void CalculateEpsKBar(const size_t &kxpts)
    {
        std::cout << "Start of ABC_H0 CalculateEpsKBar " << std::endl;
        std::valarray<double> epsKBar(0.0, KWaveVectors_.size());
        assert(KWaveVectors_.size() == epsKBar.size());
        for (size_t Kindex = 0; Kindex < KWaveVectors_.size(); Kindex++)
        {
            const double Kx = KWaveVectors_.at(Kindex)(0);
            const double Ky = KWaveVectors_.at(Kindex)(1);
            for (size_t kx = 0; kx < kxpts; kx++)
            {
                const double kTildeX = static_cast<double>(kx) / static_cast<double>(kxpts) * 2.0 * M_PI / static_cast<double>(Nx);
                for (size_t ky = 0; ky < kxpts; ky++)
                {
                    const double kTildeY = static_cast<double>(ky) / static_cast<double>(kxpts) * 2.0 * M_PI / static_cast<double>(Nx);
                    epsKBar[Kindex] += Eps0k(Kx + kTildeX, Ky + kTildeY);
                }
            }
        }

        epsKBar *= static_cast<double>(Nc) / static_cast<double>(kxpts * kxpts);
        epsKBar_ = epsKBar;
        std::cout << "End of ABC_H0 CalculateEpsKBar " << std::endl;
    }

  protected:
    ClusterSites_t RSites_;
    ClusterSites_t KWaveVectors_;
    std::valarray<double> epsKBar_;

    const double t_;
    const double tPrime_;
    const double tPrimePrime_;
};

template <size_t TNX, size_t TNY>
ABC_H0<TNX, TNY>::~ABC_H0() {} //destructors must exist
} // namespace Models