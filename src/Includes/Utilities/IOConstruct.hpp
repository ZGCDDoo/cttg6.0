#include "IO.hpp"

namespace IO
{

// Return the greensites, given a file name for the model at hand,
// in the file, the greensites are given as a complex matrix of ints.
GreenSites_t BuildGreenSites(const std::string &fname)
{
    using MatSites_t = arma::Matrix<std::complex<size_t>>;

    GreenSites_t greenSites;
    MatSites_t inGreen;
    inGreen.load(fname);

    for (size_t i = 0; i < inGreen.n_rows; ++i)
    {
        std::vector<std::pair<size_t, size_t>> tmpvec;
        for (size_t j = 0; j < inGreen.n_rows; ++j)
        {
            tmpvec.push_back(inGreen(i, j).real(), inGreen(i, j).imag());
        }

        greenSites.push_back(tmpvec);
    }
    return greenSites;
}

std::vector<std::pair<size_t, size_t>> BuildIndepSites(const GreenSites_t &greenSites)
{
    std::vector<std::pair<size_t, size_t>> indepSites;

    for (size_t i = 0; i < greenSites.size(); ++i)
    {
        for (size_t j = 0; j < greenSites.at(0.size(); ++j))
        {
            auto itsite = std::find(indepSites.begin(), indepSites.end(), greenSites.at(i).at(j));
            if (itsite == indepSites.end())
            {
                indepSites.emplace_back(greenSites.at(i).at(j));
            }
        }
    }

    return indepSites;
}

} // namespace IO