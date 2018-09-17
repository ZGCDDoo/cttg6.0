#define SUBMATRIX
#define AFM

#include "Includes/IS/MonteCarloSubMatrixBuilder_AFM.hpp"
#include "Includes/Utilities/SelfConsistencyBuilder_AFM.hpp"
#include "Includes/Utilities/FS.hpp"
#include "Includes/Utilities/SymetrizeSpins_AFM.hpp"
#include "Includes/PrintVersion.hpp"

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

#ifndef HAVEMPI
    PrintVersion::PrintVersion();
    std::ifstream fin(fname_params);
    fin >> jj;
    fin.close();
    std::cout << "Iter = " << ITER << std::endl;
    const size_t seed = jj["SEED"].get<size_t>();

    //init a model, to make sure all the files are present and that not all proc write to the same files

    const std::unique_ptr<MC::ABC_MonteCarlo> monteCarloMachinePtr = MC::MonteCarloBuilder_AFM(jj, seed);

    monteCarloMachinePtr->RunMonteCarlo();

    const std::unique_ptr<SelfCon::ABC_SelfConsistency> selfconUpPtr = SelfCon::SelfConsistencyBuilder_AFM(jj, FermionSpin_t::Up);
    selfconUpPtr->DoSCGrid();

    const std::unique_ptr<SelfCon::ABC_SelfConsistency> selfconDownPtr = SelfCon::SelfConsistencyBuilder_AFM(jj, FermionSpin_t::Down);
    selfconDownPtr->DoSCGrid();

    SymmetrizeSpins::SymmetrizeUpAndDown(jj);
    IO::FS::PrepareNextIter(paramsName, ITER);

#endif

#ifdef HAVEMPI

    std::string jjStr;

    if (mpiUt::Rank() == mpiUt::master)
    {
        PrintVersion::PrintVersion();
        Logging::Info("ITER = " + std::to_string(ITER));
        std::ifstream fin(fname_params);
        fin >> jj;
        jjStr = jj.dump();
        fin.close();
    }

    mpi::broadcast(world, jjStr, mpiUt::master);
    jj = Json::parse(jjStr);
    world.barrier();
    //wait_all
    const size_t rank = world.rank();
    const size_t seed = jj["SEED"].get<size_t>() + 2797 * rank;

    {
        std::unique_ptr<MC::ABC_MonteCarlo> monteCarloMachinePtr = MC::MonteCarloBuilder_AFM(jj, seed);

        monteCarloMachinePtr->RunMonteCarlo();
    }

    world.barrier();

    const std::unique_ptr<SelfCon::ABC_SelfConsistency> selfconUpPtr = SelfCon::SelfConsistencyBuilder_AFM(jj, FermionSpin_t::Up);
    selfconUpPtr->DoSCGrid();

    const std::unique_ptr<SelfCon::ABC_SelfConsistency> selfconDownPtr = SelfCon::SelfConsistencyBuilder_AFM(jj, FermionSpin_t::Down);
    selfconDownPtr->DoSCGrid();

    if (mpiUt::Rank() == mpiUt::master)
    {
        SymmetrizeSpins::SymmetrizeUpAndDown(jj);
        IO::FS::PrepareNextIter(paramsName, ITER);
    }
#endif

    return EXIT_SUCCESS;
}
