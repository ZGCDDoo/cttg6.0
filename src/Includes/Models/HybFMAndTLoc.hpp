#pragma once

#include "../Utilities/Utilities.hpp"
#include "../Utilities/Integrator.hpp"

namespace Models
{

template <typename TH0>
class HybFMAndTLoc
{

  public:
    static void GetHybFMAndTLoc()
    {
        //Get TLoc = Int[t(ktilde)]
        TH0 h0();
        const ClusterMatrixCD_t tlocR = Integrator::CubatureKTilde<TH0>(h0);
        const ClusterMatrixCD_t tlocK = FourierDCA::RtoK(tlocR, h0.RSites(), h0.KWaveVectors());
        tlocR.save("tloc.arma", arma::arma_ascii);
        tlocK.save("tloc_K.arma", arma::arma_ascii);

        //Get Int[ t(ktilde)^2 ]
        TKTildeSquared fct;
        const ClusterMatrixCD_t tktildeSquaredIntegrated = Integrator::CubatureKTilde<TH0>(fct);
    }

    template <typename TH0>
    struct TKTildeSquared
    {

        TKTildeSquared() : h0_(){};
        const size_t n_rows = TH0::n_rows;
        const size_t n_cols = TH0::n_cols;

        ClusterMatrixCD_t operator()(const double &kTildeX, const double &kTildeY) //return t(ktilde)^2
        {
            const ClusterMatrixCD_t tmp = h0_(kTildeX, kTildeY);
            return (tmp * tmp);
        }
    }
}

} // namespace Models
