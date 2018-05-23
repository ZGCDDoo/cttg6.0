#pragma once

#include "Utilities.hpp"
#include "MPIUtilities.hpp"

namespace IO
{
const size_t Nx1 = 1;
const size_t Nx2 = 2;
const size_t Nx4 = 4;

template <size_t TNX, size_t TNY>
class Base_IOModel
{

  public:
    static const size_t Nc = TNX * TNY;
    static const size_t Nx = TNX;
    static const size_t Ny = TNY;
    Base_IOModel(){};

    void FinishConstructor()
    {
        ConstructfullSiteToIndepSite();
        SetnOfAssociatedSites();
        ConstructFillingSites();
        AssertSanity();
    }

    void ConstructFillingSites()
    {
        const size_t KK = indepSites_.size();
        for (size_t ii = 0; ii < KK; ii++)
        {
            Site_t s1 = indepSites_.at(ii).first;
            Site_t s2 = indepSites_.at(ii).second;
            if (s1 == s2)
            {
                fillingSites_.push_back(s1);
                fillingSitesIndex_.push_back(ii);
            }
        }
    }

    void ConstructfullSiteToIndepSite()
    {
        //construct equivalentSites_ also.
        equivalentSites_.clear();
        equivalentSites_.resize(indepSites_.size());

        for (size_t ii = 0; ii < indepSites_.size(); ii++)
        {
            std::pair<size_t, size_t> pairii = indepSites_.at(ii);
            equivalentSites_.at(ii).push_back(pairii);
        }

        for (Site_t s1 = 0; s1 < Nc; s1++)
        {
            for (Site_t s2 = 0; s2 < Nc; s2++)
            {

                std::pair<Site_t, Site_t> pairSites = GreenSites_.at(s1).at(s2);
                std::vector<std::pair<size_t, size_t>>::iterator llit = std::find(indepSites_.begin(), indepSites_.end(), pairSites);
                if (llit == indepSites_.end())
                {
                    //Try the pair with first and second exchanged:
                    llit = std::find(indepSites_.begin(), indepSites_.end(), std::make_pair(pairSites.second, pairSites.first));
                    if (llit == indepSites_.end())
                    {
                        throw std::runtime_error("Bad index in FindIndepSiteIndex!");
                    }
                }
                const size_t llDistance = std::distance(indepSites_.begin(), llit);
                fullSiteToIndepSite_.push_back(llDistance);
                equivalentSites_.at(llDistance).push_back(std::make_pair(s1, s2));

                if (s1 == s2)
                {
                    assert(indepSites_.at(llDistance).first == indepSites_.at(llDistance).second);
                }
            }
        }
    }

    size_t FindIndepSiteIndex(Site_t s1, Site_t s2)
    {
        return fullSiteToIndepSite_.at(s1 * Nc + s2);
    }

    void SetnOfAssociatedSites()
    {
        nOfAssociatedSites_.resize(indepSites_.size());
        // nOfFillingSites_ = 0;

        for (size_t ii = 0; ii < Nc; ii++)
        {
            for (size_t jj = 0; jj < Nc; jj++)
            {
                size_t ll = FindIndepSiteIndex(ii, jj);
                nOfAssociatedSites_.at(ll) = nOfAssociatedSites_.at(ll) + 1;

                // if (indepSites_.at(ll).first == indepSites_[ll].second)
                // {
                //     nOfFillingSites_ += 1;
                // }
            }
        }
    }

    //Save a green in .arma format.
    ClusterCubeCD_t
    ReadGreen(std::string filename)
    {
        mpiUt::Print("In IOModel READGREEN ");

        ClusterMatrixCD_t fileMat;
        ClusterMatrixCD_t tmp(Nc, Nc);
        fileMat.load(filename);
        assert(fileMat.n_cols == this->indepSites_.size());

        ClusterCubeCD_t cubetmp(Nc, Nc, fileMat.n_rows);

        for (size_t n = 0; n < cubetmp.n_slices; n++)
        {
            for (size_t ii = 0; ii < Nc; ii++)
            {
                for (size_t jj = 0; jj < Nc; jj++)
                {
                    size_t index = FindIndepSiteIndex(ii, jj);
                    tmp(ii, jj) = fileMat(n, index);
                }
            }

            cubetmp.slice(n) = tmp;
        }

        return cubetmp;
    }

    //Save a green in .dat format.
    ClusterCubeCD_t
    ReadGreenDat(std::string filename)
    {
        mpiUt::Print("In IOModel ReadGreenNDat ");

        ClusterMatrix_t fileMat;
        ClusterMatrixCD_t tmp(Nc, Nc);
        fileMat.load(filename);
        assert(fileMat.n_cols == 2 * this->indepSites_.size() + 1);
        fileMat.shed_col(0); // we dont want the matsubara frequencies

        ClusterCubeCD_t cubetmp(Nc, Nc, fileMat.n_rows);

        for (size_t n = 0; n < cubetmp.n_slices; n++)
        {
            for (size_t ii = 0; ii < Nc; ii++)
            {
                for (size_t jj = 0; jj < Nc; jj++)
                {
                    size_t index = FindIndepSiteIndex(ii, jj);
                    tmp(ii, jj) = cd_t(fileMat(n, 2 * index), fileMat(n, 2 * index + 1));
                }
            }

            cubetmp.slice(n) = tmp;
        }

        return cubetmp;
    }

    void SaveCube(std::string fname, ClusterCubeCD_t green, const double &beta, const bool &saveArma = false)
    {
        const size_t NMat = green.n_slices;
        ClusterMatrixCD_t greenOut(NMat, this->indepSites_.size());

        std::ofstream fout;
        fout.open(fname + std::string(".dat"), std::ios::out);
        double iwn;
        for (size_t nn = 0; nn < green.n_slices; nn++)
        {
            iwn = (2.0 * nn + 1.0) * M_PI / beta;
            fout << iwn << " ";

            for (Site_t ii = 0; ii < this->indepSites_.size(); ii++)
            {
                Site_t s1 = this->indepSites_.at(ii).first;
                Site_t s2 = this->indepSites_.at(ii).second;

                greenOut(nn, ii) = green(s1, s2, nn);
                fout << green(s1, s2, nn).real()
                     << " "
                     << green(s1, s2, nn).imag()
                     << " ";
            }
            fout << "\n";
        }

        fout.close();

        if (saveArma)
        {
            greenOut.save(fname + std::string(".arma"), arma::arma_ascii);
        }
        return;
    }

    //Save the independant values as tabular form
    // void SaveTabular(const std::string &fname, const ClusterMatrixCD_t &greenTab, const double &beta,
    //                  const bool &saveArma = true)
    // {

    //     assert(this->indepSites_.size() == greenTab.n_cols);
    //     mpiUt::IOResult::SaveTabular(fname, greenTab, beta, saveArma);
    // }

    //return the full matrix for only a vector with elements being the indep sites values.
    //i.e for one matsubara frequency, return the full matrix associated to the independant
    //elements of that matrix

    template <typename T1_t, typename T2_t = ClusterMatrixCD_t>
    T2_t IndepToFull(const T1_t &indepElements) //in practice will be a Sitevector_t or SitevectorCD_t
    {

        T2_t fullMatrix(Nc, Nc);

        for (size_t ii = 0; ii < Nc; ii++)
        {
            for (size_t jj = 0; jj < Nc; jj++)
            {
                size_t index = FindIndepSiteIndex(ii, jj);
                fullMatrix(ii, jj) = indepElements.at(index);
            }
        }

        return fullMatrix;
    }

    //from th full cube return the independant in tabular form
    template <typename T1_t>
    ClusterMatrixCD_t FullCubeToIndep(const T1_t &greenCube) //TT = {ClusterCube_t or clustercubeCD_t}
    {

        ClusterMatrixCD_t indepTabular(greenCube.n_slices, indepSites_.size());

        for (size_t i = 0; i < indepSites_.size(); i++)
        {
            Site_t s1 = indepSites_.at(i).first;
            Site_t s2 = indepSites_.at(i).second;
            for (size_t n = 0; n < greenCube.n_slices; n++)
            {
                indepTabular(n, i) = greenCube(s1, s2, n);
            }
        }

        return indepTabular;
    }

    void AssertSanity()
    {
        size_t sum = 0.0;
        assert(fillingSites_.size() == fillingSitesIndex_.size());

        for (size_t ii = 0; ii < fillingSitesIndex_.size(); ii++)
        {
            sum += nOfAssociatedSites_.at(fillingSitesIndex_.at(ii));
        }
        assert(sum == Nc);

        //make sure fullSitetoIndepSite is ok
        assert(fullSiteToIndepSite_.size() == Nc * Nc);

        //make sure equivalentSites_ is ok
        assert(equivalentSites_.size() == indepSites_.size());
        for (size_t ii = 0; ii < equivalentSites_.size(); ii++)
        {
            std::pair<size_t, size_t> pairii = equivalentSites_.at(ii).at(0);
            for (size_t jj = 0; jj < equivalentSites_.at(ii).size(); jj++)
            {
                const size_t s1 = equivalentSites_.at(ii).at(jj).first;
                const size_t s2 = equivalentSites_.at(ii).at(jj).second;
                assert(pairii == GreenSites_.at(s1).at(s2));
            }
        }
    }

    std::pair<size_t, size_t> FindSitesRng(const size_t &s1, const size_t &s2, const double &rngDouble)
    {

        const size_t indepSiteIndex = FindIndepSiteIndex(s1, s2);
        const size_t equivalentSize = equivalentSites_.at(indepSiteIndex).size();
        const size_t intRng = rngDouble * equivalentSize;
        return equivalentSites_.at(indepSiteIndex).at(intRng);
    }

    //Getters
    std::vector<std::pair<size_t, size_t>> const indepSites() { return indepSites_; };
    std::vector<std::vector<std::pair<size_t, size_t>>> const GreenSites() { return GreenSites_; };
    std::vector<std::vector<std::pair<size_t, size_t>>> const equivalentSites() { return equivalentSites_; };
    std::vector<size_t> const nOfAssociatedSites() { return nOfAssociatedSites_; };
    std::vector<size_t> const fillingSites() { return fillingSites_; };
    std::vector<size_t> const fillingSitesIndex() { return fillingSitesIndex_; };
    std::vector<size_t> const downEquivalentSites() { return downEquivalentSites_; };

  protected:
    std::vector<std::pair<size_t, size_t>> indepSites_;
    std::vector<std::vector<std::pair<size_t, size_t>>> GreenSites_;
    std::vector<std::vector<std::pair<size_t, size_t>>> equivalentSites_; // for example for square lattice equivalentsites.at(0) = {{0.0}, {1,1} , {2,2}, {3,3}}
    std::vector<size_t> nOfAssociatedSites_;
    std::vector<size_t> fullSiteToIndepSite_;
    std::vector<size_t> fillingSites_;
    std::vector<size_t> fillingSitesIndex_; //the indexes of the fillingsites in the indepSites_
    std::vector<size_t> downEquivalentSites_;
};

class IOTriangle2x2 : public Base_IOModel<Nx2, Nx2>
{
  public:
    IOTriangle2x2() : Base_IOModel<Nx2, Nx2>()
    {
        this->indepSites_ = {
            {0, 0}, {1, 1}, {0, 1}, {0, 3}, {1, 2}};

        this->GreenSites_ = {
            {{0, 0}, {0, 1}, {0, 1}, {0, 3}},
            {{0, 1}, {1, 1}, {1, 2}, {0, 1}},
            {{0, 1}, {1, 2}, {1, 1}, {0, 1}},
            {{0, 3}, {0, 1}, {0, 1}, {0, 0}}};

        FinishConstructor();
    }
};

class IOSquare2x2 : public Base_IOModel<Nx2, Nx2>
{
  public:
    IOSquare2x2() : Base_IOModel<Nx2, Nx2>()
    {
        this->indepSites_ = {
            {0, 0}, {0, 1}, {0, 3}};

        this->GreenSites_ = {
            {{0, 0}, {0, 1}, {0, 1}, {0, 3}},
            {{0, 1}, {0, 0}, {0, 3}, {0, 1}},
            {{0, 1}, {0, 3}, {0, 0}, {0, 1}},
            {{0, 3}, {0, 1}, {0, 1}, {0, 0}}};

        FinishConstructor();
    }
};

class IOSquare2x2AFM : public Base_IOModel<Nx2, Nx2>
{
  public:
    IOSquare2x2AFM() : Base_IOModel<Nx2, Nx2>()
    {
        this->indepSites_ = {
            {0, 0}, {1, 1}, {0, 1}, {0, 3}, {1, 2}};

        this->GreenSites_ = {
            {{0, 0}, {0, 1}, {0, 1}, {0, 3}},
            {{0, 1}, {1, 1}, {1, 2}, {0, 1}},
            {{0, 1}, {1, 2}, {1, 1}, {0, 1}},
            {{0, 3}, {0, 1}, {0, 1}, {0, 0}}};

        this->downEquivalentSites_ = {1, 0, 2, 4, 3};
        FinishConstructor();
    }
};

class IOSIAM : public Base_IOModel<Nx1, Nx1>
{
  public:
    IOSIAM() : Base_IOModel<Nx1, Nx1>()
    {
        this->indepSites_ = {
            {0, 0}};

        this->GreenSites_ = {
            {{0, 0}}};

        FinishConstructor();
    }
};

class IOSquare4x4 : public Base_IOModel<Nx4, Nx4>
{
  public:
    IOSquare4x4() : Base_IOModel<Nx4, Nx4>()
    {
        this->indepSites_ = {
            {0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 5}, {0, 6}, {0, 7}, {0, 10}, {0, 11}, {0, 15}, {1, 1}, {1, 2}, {1, 4}, {1, 5}, {1, 6}, {1, 7}, {1, 9}, {1, 10}, {1, 11}, {1, 13}, {1, 14}, {5, 5}, {5, 6}, {5, 10}};

        this->GreenSites_ = {
            {{0, 0}, {0, 1}, {0, 2}, {0, 3}, {0, 1}, {0, 5}, {0, 6}, {0, 7}, {0, 2}, {0, 6}, {0, 10}, {0, 11}, {0, 3}, {0, 7}, {0, 11}, {0, 15}},
            {{0, 1}, {1, 1}, {1, 2}, {0, 2}, {1, 4}, {1, 5}, {1, 6}, {1, 7}, {1, 7}, {1, 9}, {1, 10}, {1, 11}, {0, 7}, {1, 13}, {1, 14}, {0, 11}},
            {{0, 2}, {1, 2}, {1, 1}, {0, 1}, {1, 7}, {1, 6}, {1, 5}, {1, 4}, {1, 11}, {1, 10}, {1, 9}, {1, 7}, {0, 11}, {1, 14}, {1, 13}, {0, 7}},
            {{0, 3}, {0, 2}, {0, 1}, {0, 0}, {0, 7}, {0, 6}, {0, 5}, {0, 1}, {0, 11}, {0, 10}, {0, 6}, {0, 2}, {0, 15}, {0, 11}, {0, 7}, {0, 3}},
            {{0, 1}, {1, 4}, {1, 7}, {0, 7}, {1, 1}, {1, 5}, {1, 9}, {1, 13}, {1, 2}, {1, 6}, {1, 10}, {1, 14}, {0, 2}, {1, 7}, {1, 11}, {0, 11}},
            {{0, 5}, {1, 5}, {1, 6}, {0, 6}, {1, 5}, {5, 5}, {5, 6}, {1, 9}, {1, 6}, {5, 6}, {5, 10}, {1, 10}, {0, 6}, {1, 9}, {1, 10}, {0, 10}},
            {{0, 6}, {1, 6}, {1, 5}, {0, 5}, {1, 9}, {5, 6}, {5, 5}, {1, 5}, {1, 10}, {5, 10}, {5, 6}, {1, 6}, {0, 10}, {1, 10}, {1, 9}, {0, 6}},
            {{0, 7}, {1, 7}, {1, 4}, {0, 1}, {1, 13}, {1, 9}, {1, 5}, {1, 1}, {1, 14}, {1, 10}, {1, 6}, {1, 2}, {0, 11}, {1, 11}, {1, 7}, {0, 2}},
            {{0, 2}, {1, 7}, {1, 11}, {0, 11}, {1, 2}, {1, 6}, {1, 10}, {1, 14}, {1, 1}, {1, 5}, {1, 9}, {1, 13}, {0, 1}, {1, 4}, {1, 7}, {0, 7}},
            {{0, 6}, {1, 9}, {1, 10}, {0, 10}, {1, 6}, {5, 6}, {5, 10}, {1, 10}, {1, 5}, {5, 5}, {5, 6}, {1, 9}, {0, 5}, {1, 5}, {1, 6}, {0, 6}},
            {{0, 10}, {1, 10}, {1, 9}, {0, 6}, {1, 10}, {5, 10}, {5, 6}, {1, 6}, {1, 9}, {5, 6}, {5, 5}, {1, 5}, {0, 6}, {1, 6}, {1, 5}, {0, 5}},
            {{0, 11}, {1, 11}, {1, 7}, {0, 2}, {1, 14}, {1, 10}, {1, 6}, {1, 2}, {1, 13}, {1, 9}, {1, 5}, {1, 1}, {0, 7}, {1, 7}, {1, 4}, {0, 1}},
            {{0, 3}, {0, 7}, {0, 11}, {0, 15}, {0, 2}, {0, 6}, {0, 10}, {0, 11}, {0, 1}, {0, 5}, {0, 6}, {0, 7}, {0, 0}, {0, 1}, {0, 2}, {0, 3}},
            {{0, 7}, {1, 13}, {1, 14}, {0, 11}, {1, 7}, {1, 9}, {1, 10}, {1, 11}, {1, 4}, {1, 5}, {1, 6}, {1, 7}, {0, 1}, {1, 1}, {1, 2}, {0, 2}},
            {{0, 11}, {1, 14}, {1, 13}, {0, 7}, {1, 11}, {1, 10}, {1, 9}, {1, 7}, {1, 7}, {1, 6}, {1, 5}, {1, 4}, {0, 2}, {1, 2}, {1, 1}, {0, 1}},
            {{0, 15}, {0, 11}, {0, 7}, {0, 3}, {0, 11}, {0, 10}, {0, 6}, {0, 2}, {0, 7}, {0, 6}, {0, 5}, {0, 1}, {0, 3}, {0, 2}, {0, 1}, {0, 0}},
        };

        FinishConstructor();
    }
};
} // namespace IO