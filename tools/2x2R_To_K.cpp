#define DCA

#include "../src/Includes/Models/ModelSquare2x2.hpp"
#include "../src/Includes/Utilities/IO.hpp"

using Model_t = Models::ModelSquare2x2;
using H0_t = Models::H0Square<2, 2>;
using IOModel_t = IO::IOSquare2x2;

int main(int argc, char **argv)
{
    if (argc != 4)
    {
        throw std::runtime_error("Miseria: Wrong number of input parameters. Stupido !");
    }

    const std::string fname_hyb = argv[1];
    const double BETA = atof(argv[2]);
    std::cout << "BETA = " << BETA << std::endl;
    std::cout << "fname = " << fname_hyb << std::endl;

    const std::string fname_params = argv[3];
    Json jj;

    std::ifstream fin(fname_params);
    fin >> jj;
    fin.close();
    jj["beta"] = BETA;

    //0.) Read The real-space 2x2 Hyb
    const IOModel_t ioModel;
    // const Model_t model(jj);
    const H0_t h0(jj["t"], jj["tPrime"], jj["tPrimePrime"]);
    const ClusterCubeCD_t hyb_R = ioModel.ReadGreenDat(fname_hyb);
    const ClusterCubeCD_t hyb_K = FourierDCA::RtoK(hyb_R, h0.RSites(), h0.KWaveVectors());
    ioModel.SaveK("2x2_Converted_R_to_K.dat", hyb_K, BETA);

    return EXIT_SUCCESS;
}