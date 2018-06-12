

#include <gtest/gtest.h>

#include "../src/Includes/Utilities/Fourier_DCA.hpp"
#include "../src/Includes/Models/H0Square.hpp"

using namespace FourierDCA;

const double DELTA = 1e-11;
TEST(FourierDCATests, KToR)
{
    const size_t Nx = 8;
    const size_t NMat = 100;
    const Models::H0Square<Nx, Nx> h0(-1.0, -0.32, 0.15);

    DataK_t greenK(Nx * Nx, NMat);
    greenK.randu();

    ClusterCubeCD_t greenRTest = KtoR(greenK, h0.RSites(), h0.KWaveVectors());
    DataK_t greenKTest = RtoK(greenRTest, h0.RSites(), h0.KWaveVectors());

    assert(greenKTest.n_cols == NMat);
    assert(greenKTest.n_rows == Nx * Nx);
    assert(greenK.n_rows == Nx * Nx);
    assert(greenKTest.n_cols == NMat);

    std::cout << "Here " << std::endl;
    for (size_t ii = 0; ii < greenKTest.n_rows; ii++)
    {
        for (size_t jj = 0; jj < greenKTest.n_cols; jj++)
        {
            // std::cout << "ii, jj = " << ii << ", " << jj << std::endl;
            ASSERT_NEAR(greenKTest(ii, jj).real(), greenK(ii, jj).real(), DELTA);
            ASSERT_NEAR(greenKTest(ii, jj).imag(), greenK(ii, jj).imag(), DELTA);
        }
    }
}

int main(int argc, char **argv)
{
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
