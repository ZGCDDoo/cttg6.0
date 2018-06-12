#pragma once

#include "Utilities.hpp"

namespace SelfCon
{

class ABC_SelfConsistency
{
  public:
    ABC_SelfConsistency(){};
    ~ABC_SelfConsistency(){};
    virtual void DoSCGrid() = 0;

  private:
}; // class ABC_SelfConsistency

} // namespace SelfCon