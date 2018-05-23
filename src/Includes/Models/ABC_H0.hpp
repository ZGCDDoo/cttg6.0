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

    // ABC_H0() : t_(0.0),
    //            tPrime_(0.0),
    //            tPrimePrime_(0.0),
    //            RSites_(),
    //            KWaveVectors_(){};
    ABC_H0(const double &t, const double &tp, const double &tpp) : RSites_(Nc),
                                                                   KWaveVectors_(Nc),
                                                                   t_(t),
                                                                   tPrime_(tp),
                                                                   tPrimePrime_(tpp)

    {
        assert(TNX == TNY);

        for (size_t i = 0; i < TNX; i++)
        {
            for (size_t j = 0; j < TNY; j++)
            {
                size_t index = i + TNY * j;
                RSites_[index] = {static_cast<double>(i), static_cast<double>(j)};
                KWaveVectors_[index] = {static_cast<double>(i) * 2.0 * M_PI / static_cast<double>(TNX), static_cast<double>(j) * 2.0 * M_PI / static_cast<double>(TNY)};
            }
        }
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
        cd_t im = cd_t(0.0, 1.0);
        arma::vec ktilde = {kTildeX, kTildeY};
        ClusterMatrixCD_t HoppingKTilde(Nc, Nc);
        HoppingKTilde.zeros();

        for (size_t i = 0; i < Nc; i++)
        {
            for (size_t j = 0; j < Nc; j++)
            {
                for (const auto &K : this->KWaveVectors_)
                {
                    HoppingKTilde(i, j) += std::exp(im * dot(K + ktilde, RSites_[i] - RSites_[j])) * Eps0k(K[0] + kTildeX, K[1] + kTildeY);
                }
            }
        }
        return (HoppingKTilde / static_cast<double>(Nc));
    }

    ClusterCubeCD_t SaveTKTilde(const size_t &kxpts)
    {
        // cd_t im = cd_t(0.0, 1.0);
        ClusterMatrixCD_t HoppingKTilde(Nc, Nc);
        HoppingKTilde.zeros();
        ClusterCubeCD_t tKTildeGrid(Nc, Nc, kxpts * kxpts);
        tKTildeGrid.zeros();

        size_t sliceindex = 0;
        for (size_t kx = 0; kx < kxpts; kx++)
        {
            double kTildeX = static_cast<double>(kx) / static_cast<double>(kxpts) * 2.0 * M_PI / static_cast<double>(Nx);
            for (size_t ky = 0; ky < kxpts; ky++)
            {
                double kTildeY = static_cast<double>(ky) / static_cast<double>(kxpts) * 2.0 * M_PI / static_cast<double>(Nx);
                tKTildeGrid.slice(sliceindex) = (*this)(kTildeX, kTildeY);
                sliceindex++;
            }
        }
        tKTildeGrid.save("tktilde.arma", arma::arma_ascii);
        return tKTildeGrid;
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
}