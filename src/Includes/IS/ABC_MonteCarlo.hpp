#pragma once

#include "../Utilities/Utilities.hpp"

namespace MC
{

class ABC_MonteCarlo
{
  public:
    ABC_MonteCarlo(){};

    ~ABC_MonteCarlo(){};

    virtual void RunMonteCarlo() = 0;
}; // namespace MC
} // namespace MC
