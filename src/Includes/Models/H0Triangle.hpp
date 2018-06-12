#pragma once

#ifdef DCA
#include "ABC_H0_DCA.hpp"
#else
#include "ABC_H0.hpp"
#endif

namespace Models
{
template <size_t TNX, size_t TNY>
class H0Triangle : public ABC_H0<TNX, TNY>
{

  public:
    H0Triangle(const double &t, const double &tp, const double &tpp) : ABC_H0<TNX, TNY>(t, tp, tpp) {}

    double Eps0k(const double &kx, const double &ky) const
    {
        return (
            2.0 * (this->t_ * (cos(kx) + cos(ky)) +
                   this->tPrime_ * (cos(kx + ky))));
    }
};
} // namespace Models