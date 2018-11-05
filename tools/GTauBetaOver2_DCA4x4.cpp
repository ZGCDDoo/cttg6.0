#define DCA

#include "../src/Includes/Models/ModelSquare4x4_DCA.hpp"
#include "../src/Includes/Utilities/IO.hpp"
#include "../src/Includes/IS/Obs/KineticEnergy.hpp"

using Model_t = Models::ModelSquare4x4_DCA;
using H0_t = Models::H0Square<4, 4>;
using IOModel_t = IO::IOSquare4x4_DCA;
using Kinetic_t = Markov::Obs::KineticEnergy<Model_t, IOModel_t>;

int main(int argc, char **argv)
{
    if (argc != 3)
    {
        throw std::runtime_error("Miseria: Wrong number of input parameters. Stupido !");
    }

    const std::string fname_green = argv[1];

    const std::string fname_params = argv[2];
    Json jj;

    std::ifstream fin(fname_params);
    fin >> jj;
    fin.close();

    //0.) Read The DCA 4x4 green
    const IOModel_t ioModel;
    const ClusterCubeCD_t green = ioModel.ReadGreenKDat(fname_green);
    const std::shared_ptr<Model_t> modelPtr = std::make_shared<Model_t>(Model_t(jj));

    const Kinetic_t kEnergy(modelPtr, green);
    kEnergy.GetKineticEnergy();

    return EXIT_SUCCESS;
}