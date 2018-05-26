
#include "Includes/IS/MonteCarlo.hpp"
#include "Includes/IS/MarkovChain.hpp"

#include "Includes/Models/ModelSquare4x4.hpp"

#include "Includes/Utilities/SelfConsistency.hpp"
#include "Includes/Utilities/FS.hpp"

int main(int argc, char **argv)
{

#ifdef HAVEMPI
    mpi::environment env;
    mpi::communicator world;
#endif

    if (argc != 3)
    {
        throw std::runtime_error("Miseria: Wrong number of input parameters. Stupido !");
    }

    const std::string paramsName = argv[1];
    const int ITER = atoi(argv[2]);
    const std::string fname_params = paramsName + std::to_string(ITER) + std::string(".json");
    Json jj;

    const size_t Nx = 4;

    using Model_t = Models::ModelSquare4x4;
    using IOModel_t = IO::IOSquare4x4;
    using H0_t = Models::H0Square<Nx, Nx>;
    using Markov_t = Markov::MarkovChain<IOModel_t, Model_t>;

    //init a model, to make sure all the files are present and that not all proc write to the same files

#ifndef HAVEMPI
    std::ifstream fin(fname_params);
    fin >> jj;
    fin.close();
    std::cout << "Iter = " << ITER << std::endl;
    const size_t seed = jj["SEED"].get<size_t>();
    MC::MonteCarlo<Markov_t> monteCarloMachine(std::make_shared<Markov_t>(jj, seed), jj);
    monteCarloMachine.RunMonteCarlo();
    Model_t model(jj);
    IOModel_t ioModel;
    const ClusterCubeCD_t greenImpurityUp = ioModel.ReadGreenDat("greenUp.dat");
    SelfCon::SelfConsistency<IOModel_t, Model_t, H0_t> selfconUp(jj, model, greenImpurityUp,
                                                                 "Up");
    selfconUp.DoSCGrid();
    IO::FS::PrepareNextIter(paramsName, ITER);
#endif

#ifdef HAVEMPI

    std::string jjStr;

    if (mpiUt::Rank() == mpiUt::master)
    {
        std::ifstream fin(fname_params);
        fin >> jj;
        jjStr = jj.dump();
        fin.close();
    }

    mpi::broadcast(world, jjStr, mpiUt::master);
    jj = Json::parse(jjStr);

    if (mpiUt::Rank() == mpiUt::master)
    {
        std::cout << "Iter = " << ITER << std::endl;
        const Model_t modelDummy(jj);
    }

    world.barrier();
    //wait_all

    const size_t seed = jj["SEED"].get<size_t>() + world.rank();
    {
        MC::MonteCarlo<Markov_t> monteCarloMachine(std::make_shared<Markov_t>(jj, seed), jj);
        monteCarloMachine.RunMonteCarlo();
    }
    world.barrier();

    Model_t model(jj);
    IOModel_t ioModel;
    const ClusterCubeCD_t greenImpurityUp = ioModel.ReadGreenDat("greenUp.dat");
    SelfCon::SelfConsistency<IOModel_t, Model_t, H0_t> selfconUp(jj, model, greenImpurityUp,
                                                                 "Up");
    selfconUp.DoSCGrid();

    if (mpiUt::Rank() == mpiUt::master)
    {
        IO::FS::PrepareNextIter(paramsName, ITER);
    }

#endif
    return EXIT_SUCCESS;
}
