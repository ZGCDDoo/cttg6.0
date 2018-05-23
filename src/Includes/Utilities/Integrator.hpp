
#pragma once

#include "Utilities.hpp"

namespace Integrator
{
const size_t MAXEVALS = 100000;
const size_t KXPTS = 64;

template <typename TFct>
ClusterMatrixCD_t GridKTilde(TFct fct, size_t kxpts = KXPTS)
{
    double fact = 1.0 / (static_cast<double>(kxpts * kxpts));
    ClusterMatrixCD_t integral(fct.Nc, fct.Nc);
    integral.zeros();

    for (size_t ii = 0; ii < kxpts; ii++)
    {
        double kx = static_cast<double>(ii) / static_cast<double>(kxpts) * 2.0 * M_PI / static_cast<double>(fct.Nx);
        for (size_t jj = 0; jj < kxpts; jj++)
        {
            double ky = static_cast<double>(jj) / static_cast<double>(kxpts) * 2.0 * M_PI / static_cast<double>(fct.Ny);
            integral += fct(kx, ky);
        }
    }
    return (fact * integral);
}
}