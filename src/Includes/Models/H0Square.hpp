#pragma once

#ifdef DCA
#include "ABC_H0_DCA.hpp"
#else
#include "ABC_H0.hpp"
#endif

namespace Models
{
template <size_t TNX, size_t TNY>
class H0Square : public ABC_H0<TNX, TNY>
{

  public:
    //H0Square(){};
    H0Square(const double &t, const double &tp, const double &tpp) : ABC_H0<TNX, TNY>(t, tp, tpp) {}

    double Eps0k(const double &kx, const double &ky) const
    {
        return (2.0 * (this->t_ * (cos(kx) + cos(ky))) + 2.0 * this->tPrime_ * (cos(kx + ky) + cos(kx - ky)) + 2.0 * this->tPrimePrime_ * (cos(2.0 * kx) + cos(2.0 * ky)));
    }
};
}