#pragma once

#include "ABC_Model.hpp"
#include "H0Triangle.hpp"

namespace Models
{

class ModelTriangle3x2_DCA : public ABC_Model_2D<IO::IOTriangle3x2_DCA, H0Triangle<Nx3, Nx2>>
{

    using IOModel_t = IO::IOTriangle3x2_DCA;
    using H0_t = H0Triangle<Nx3, Nx2>;

  public:
    ModelTriangle3x2_DCA(const Json &jj) : ABC_Model_2D<IOModel_t, H0_t>(jj){};
};
} // namespace Models