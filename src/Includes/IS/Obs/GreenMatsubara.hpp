#pragma once

namespace Markov
{
namespace Obs
{

class GreenMatsubara
{

    GreenMatsubara()
    {
        if (isPrecise_)
        {
            expsCached_ = Cache::ComputeExps(dataCT_->beta(), greenMatsubaraUp_.n_rows, NTauExp_);

#ifdef AFM
            MupAverages_.resize(NMat_);
            MdownAverages_.resize(NMat_);
#endif
#ifndef AFM
            Maverages_.resize(NMat_);
#endif
            const size_t Nc = ioModel_->Nc;
            for (size_t nn = 0; nn < NMat_; nn++)
            {
#ifdef AFM
                MupAverages_.at(nn).resize(Nc, Nc);
                MupAverages_.at(nn).zeros();
                MdownAverages_.at(nn).resize(Nc, Nc);
                MdownAverages_.at(nn).zeros();
#endif
#ifndef AFM
                Maverages_.at(nn).resize(Nc, Nc);
                Maverages_.at(nn).zeros();
#endif
            }
        }
    }

    void MeasureGreen() //Part not needed to be sampled are added to the greens in Save()
    {
#ifndef AFM
        Maveraged_ = 0.5 * (*MupPtr_ + *MdownPtr_);
#endif

        const double epsilon = 1e-12;
        double fact = dataCT_->sign_;

        // std::cout << "Maveraged_.print() " << std::endl;
        // Maveraged_.print();
        // std::cout << " " << std::endl;
#ifdef VERBOSE
        std::cout << "IN obs  Measuregreen " << std::endl;
#endif
        //size_t NMatMax = (NMeas_ % size_t(factHF_) == 0) ? NMat_ : indexHF_;
        //for (size_t nn = 0; nn < std::min(NMatMax, NMat_); nn++)
        for (size_t nn = 0; nn < NMat_; nn++)
        {
            for (size_t kk = 0; kk < dataCT_->vertices_.size(); kk++)
            {
                Site_t sitekk = dataCT_->vertices_[kk].site();
                for (size_t mm = 0; mm < dataCT_->vertices_.size(); mm++)
                {
                    Site_t sitemm = dataCT_->vertices_[mm].site();
                    Tau_t tau = dataCT_->vertices_[kk].tau() + epsilon - dataCT_->vertices_[mm].tau();

                    //epsilon added for tau < 0 in case tau == -beta and static_cast to size_t for negatif number is undefined.

                    //tau_double_index
                    double tdi = (tau < 0) ? ((tau + dataCT_->beta_ + epsilon) / dataCT_->beta_ * static_cast<double>(NTauExp_)) : (static_cast<size_t>(NTauExp_ + tau / dataCT_->beta_ * static_cast<double>(NTauExp_)));
                    //tau_size_index
                    size_t tsi = static_cast<size_t>(tdi);
                    cd_t expkkmm = ((1.0 - (tdi - tsi)) * expsCached_[nn].at(tsi) + (tdi - tsi) * expsCached_[nn][tsi]);
#ifndef AFM
                    Maverages_[nn](sitekk, sitemm) += 0.5 * fact * expkkmm * Maveraged_(kk, mm);
                    //Add the transpose of the configurations:
                    Maverages_[nn](sitekk, sitemm) += 0.5 * fact * std::conj(expkkmm) * Maveraged_(mm, kk);
#endif
#ifdef AFM
                    MupAverages_[nn](sitekk, sitemm) += 0.5 * fact * expkkmm * (*MupPtr_)(kk, mm);
                    MdownAverages_[nn](sitekk, sitemm) += 0.5 * fact * expkkmm * (*MdownPtr_)(kk, mm);

                    //Add the transpose of the configurations:
                    MupAverages_[nn](sitekk, sitemm) += 0.5 * fact * std::conj(expkkmm) * (*MupPtr_)(mm, kk);
                    MdownAverages_[nn](sitekk, sitemm) += 0.5 * fact * std::conj(expkkmm) * (*MdownPtr_)(mm, kk);
#endif
                }
            }
        }
#ifdef VERBOSE
        std::cout << "IN Obs  after Measuregreen " << std::endl;
#endif
        return;
    }

    void FinalizeGreens()
    {

        //First use the fact that some sites are equivalent for getting  smaller MsigmaAverages matrices. Then do the dot product.

        const size_t KK = ioModel_.indepSites().size();
        const size_t Nc = model_.Nc;

        //ClusterMatrix_t MupReduced(KK, KK); //MupAverageReduced
        //ClusterMatrix_t MdownReduced(KK, KK);

        SiteVectorCD_t vecLeftUp(Nc);
        SiteVectorCD_t vecRightUp(Nc);
#ifdef AFM
        SiteVectorCD_t vecLeftDown(Nc);
        SiteVectorCD_t vecRightDown(Nc);
#endif
        for (size_t nn = 0; nn < NMat_; nn++)
        {

            //The dot product
            for (Site_t ll = 0; ll < KK; ll++)
            {
                for (Site_t pp = 0; pp < Nc; pp++)
                {
                    Site_t sitell1 = ioModel_.indepSites().at(ll).first;
                    Site_t sitell2 = ioModel_.indepSites().at(ll).second;

                    vecLeftUp(pp) = model_.greenCluster0MatUp().data()(sitell1, pp, nn); //if time invariance the two lines are equal.
                    vecRightUp(pp) = model_.greenCluster0MatUp().data()(pp, sitell2, nn);
#ifdef AFM
                    vecLeftDown(pp) = model_.greenCluster0MatDown().data()(sitell1, pp, nn); //if time invariance the two lines are equal.
                    vecRightDown(pp) = model_.greenCluster0MatDown().data()(pp, sitell2, nn);
#endif
                }
#ifdef AFM
                greenMatsubaraUp_(nn, ll) -= LinAlg::Dot(vecLeftUp, MupAverages_.at(nn), vecRightUp);
                greenMatsubaraDown_(nn, ll) -= LinAlg::Dot(vecLeftDown, MdownAverages_.at(nn), vecRightDown);
#endif
#ifndef AFM
                greenMatsubaraUp_(nn, ll) -= LinAlg::Dot(vecLeftUp, Maverages_.at(nn), vecRightUp);
#endif
            }
        }

        greenMatsubaraUp_ /= (NMeas_ * dataCT_->beta_ * signMeas_);
#ifdef AFM
        greenMatsubaraDown_ /= (NMeas_ * dataCT_->beta_ * signMeas_);
#endif

        //Add the static, constant part to the greens
        //(the part not needing to be averaged be MonteCarlo Sampling)
        AddCstPartGreen();
    }

    void AddCstPartGreen()
    {

        cd_t tmp;
        for (Site_t ii = 0; ii < ioModel_.indepSites().size(); ii++)
        {
            Site_t s1 = ioModel_.indepSites().at(ii).first;
            Site_t s2 = ioModel_.indepSites().at(ii).second;
            for (size_t nn = 0; nn < NMat_; nn++)
            {
                greenMatsubaraUp_(nn, ii) += model_.greenCluster0MatUp().data()(s1, s2, nn);
#ifdef AFM
                greenMatsubaraDown_(nn, ii) += model_.greenCluster0MatDown().data()(s1, s2, nn);
#endif
            }
        }
    }

  private:
    std::vector<ClusterMatrixCD_t> Maverages_;
    std::vector<ClusterMatrixCD_t> MupAverages_; //For calculating the green functions in matsubara space
    std::vector<ClusterMatrixCD_t> MdownAverages_;
};

} // namespace Obs
} // namespace Markov