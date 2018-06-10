#pragma once

#include "SIAM_Square.hpp"
#include "ModelSquare2x2.hpp"
#include "ModelTriangle2x2.hpp"
#include "ModelSquare4x4.hpp"

namespace Models
{

std::unique_ptr<ABC_Model_2D> ModelBuilder(const Json &jj)
{
    const std::string modelType = jj["modelType"].get<std::string>();

    if (modelType == "SIAM_Square")
    {
        using Model_t = Models::SIAM_Square;
        return std::make_unique<ABC_Model_2D>(Model_t(jj));
    }
    else if (modelType == "Square2x2")
    {
        using Model_t = Models::ModelSquare2x2;
        return std::make_unique<ABC_Model_2D>(Model_t(jj));
    }
    else if (modelType == "Triangle2x2")
    {
        using Model_t = Models::ModelTriangle2x2;
        return std::make_unique<ABC_Model_2D>(Model_t(jj));
    }
    else if (modelType == "Square4x4")
    {
        using Model_t = Models::ModelSquare4x4;
        return std::make_unique<ABC_Model_2D>(Model_t(jj));
    }

    return NULL;
}

} // namespace Models
