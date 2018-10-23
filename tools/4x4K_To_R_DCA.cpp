#define DCA

#include "../src/Includes/Models/ModelSquare4x4_DCA.hpp"
#include "../src/Includes/Utilities/IO.hpp"

using Model_t = Models::ModelSquare4x4_DCA;
using H0_t = Models::H0Square<4, 4>;
using IOModel_t = IO::IOSquare4x4_DCA;

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

    //0.) Read The DCA 4x4 Hyb
    const IOModel_t ioModel;
    //const Model_t model(jj);
    const H0_t h0(jj["t"], jj["tPrime"], jj["tPrimePrime"]);
    const ClusterCubeCD_t hyb_K = ioModel.ReadGreenKDat(fname_hyb);
    const ClusterCubeCD_t hyb_R = FourierDCA::KtoR(hyb_K, h0.RSites(), h0.KWaveVectors());
    ioModel.SaveCube("4x4_Converted_K_to_R.dat", hyb_R, BETA);

    return EXIT_SUCCESS;
}