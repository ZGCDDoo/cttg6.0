#pragma once

#include "ABC_Model.hpp"
#include "H0Triangle.hpp"

namespace Models
{

class ModelTriangle4x4_DCA : public ABC_Model_2D<IO::IOTriangle4x4_DCA, H0Triangle<Nx4, Nx4>>
{

    using IOModel_t = IO::IOTriangle4x4_DCA;
    using H0_t = H0Triangle<Nx4, Nx4>;

  public:
    ModelTriangle4x4_DCA(const Json &jj) : ABC_Model_2D<IOModel_t, H0_t>(jj){};
};
} // namespace Models