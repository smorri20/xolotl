/**
 * Main.c, currently only able to load clusters
 */
#include <cstdlib>
#include <iostream>
#include <fstream>
#include <cassert>
#include <Reactant.h>
#include <PetscSolver.h>
#include <mpi.h>
#include <MPIUtils.h>
#include <Options.h>
#include <xolotlPerf.h>
#include <IMaterialFactory.h>
#include <TemperatureHandlerFactory.h>
#include <VizHandlerRegistryFactory.h>
#include <INetworkLoader.h>
#include <IReactionNetwork.h>
#include <SolverHandlerFactory.h>
#include <ISolverHandler.h>
#include <IReactionHandlerFactory.h>
#include <xolotlMemUsage.h>
#include <ctime>

using namespace std;
using std::shared_ptr;

namespace xperf = xolotlPerf;
namespace xmem = xolotlMemUsage;

//! This operation prints the start message
void printStartMessage() {
	std::cout << "Starting Xolotl Plasma-Surface Interactions Simulator"
			<< std::endl;
	// TODO! Print copyright message
	// Print date and time
	std::time_t currentTime = std::time(NULL);
	std::cout << std::asctime(std::localtime(&currentTime)); // << std::endl;
}

std::shared_ptr<xolotlFactory::IMaterialFactory> initMaterial(
		Options &options) {
	// Create the material factory
	auto materialFactory =
			xolotlFactory::IMaterialFactory::createMaterialFactory(
					options.getMaterial(), options.getDimensionNumber());

	// Initialize it with the options
	materialFactory->initializeMaterial(options);

	return materialFactory;
}

bool initTemp(Options &options) {

	bool tempInitOK = xolotlFactory::initializeTempHandler(options);
	if (!tempInitOK) {
		std::cerr << "Unable to initialize requested temperature.  Aborting"
				<< std::endl;
		return EXIT_FAILURE;
	} else
		return tempInitOK;
}

bool initViz(bool opts) {

	bool vizInitOK = xolotlFactory::initializeVizHandler(opts);
	if (!vizInitOK) {
		std::cerr
				<< "Unable to initialize requested visualization infrastructure. "
				<< "Aborting" << std::endl;
		return EXIT_FAILURE;
	} else
		return vizInitOK;
}

std::shared_ptr<xolotlSolver::PetscSolver> setUpSolver(
		std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry,
		std::shared_ptr<xolotlFactory::IMaterialFactory> material,
		std::shared_ptr<xolotlCore::ITemperatureHandler> tempHandler,
		std::shared_ptr<xolotlCore::IReactionNetwork> networkHandler,
		std::shared_ptr<xolotlSolver::ISolverHandler> solvHandler,
		Options &options) {
	// Initialize the solver handler
	solvHandler->initializeHandlers(material, tempHandler, networkHandler, options);

	// Setup the solver
	auto solverInitTimer = handlerRegistry->getTimer("initSolver");
	solverInitTimer->start();
	std::shared_ptr<xolotlSolver::PetscSolver> solver = std::make_shared<
			xolotlSolver::PetscSolver>(handlerRegistry);
	solver->setCommandLineOptions(options.getPetscArgc(),
			options.getPetscArgv());
	solver->initialize(solvHandler);
	solverInitTimer->stop();

	return solver;
}

void launchPetscSolver(std::shared_ptr<xolotlSolver::PetscSolver> solver,
		std::shared_ptr<xolotlPerf::IHandlerRegistry> handlerRegistry) {

	xperf::IHardwareCounter::SpecType hwctrSpec;
	hwctrSpec.push_back(xperf::IHardwareCounter::FPOps);
	hwctrSpec.push_back(xperf::IHardwareCounter::Cycles);
	hwctrSpec.push_back(xperf::IHardwareCounter::L3CacheMisses);

	// Launch the PetscSolver
	auto solverTimer = handlerRegistry->getTimer("solve");
	auto solverHwctr = handlerRegistry->getHardwareCounter("solve", hwctrSpec);
	solverTimer->start();
	solverHwctr->start();
	solver->solve();
	solverHwctr->stop();
	solverTimer->stop();
}

//! Main program
int main(int argc, char **argv) {

	// Local Declarations
	int rank;
    int ret = EXIT_SUCCESS;

    // Initialize MPI. We do this instead of leaving it to some
    // other package (e.g., PETSc), because we want to avoid problems
    // with overlapping Timer scopes.
    MPI_Init(&argc, &argv);

	// Check the command line arguments.
	// Skip the executable name before parsing
	argc -= 1; // one for the executable name
	argv += 1; // one for the executable name
	Options opts;
	opts.readParams(argv);
	if (!opts.shouldRun()) {
		return opts.getExitCode();
	}

	// Skip the name of the parameter file that was just used.
	// The arguments should be empty now.
	argc -= 1;
	argv += 1;

	// Extract the argument for the file name
	std::string networkFilename = opts.getNetworkFilename();
	assert(!networkFilename.empty());

	try {
		// Set up our performance data infrastructure and memory
        // usage monitoring infrastructure.
		xperf::initialize(opts.getPerfHandlerType());
		auto handlerRegistry = xolotlPerf::getHandlerRegistry();
        xmem::initialize(opts.getMemUsageHandlerType(),
                            opts.getMemUsageSamplingInterval(),
                            opts.getMemUsageProfileFilename());
        auto memUsageHandlerRegistry = xmem::getHandlerRegistry();

        // Start of "real" work.
        // We start a new scope so that any memory allocated
        // on the stack during the "real" work can be recovered
        // before aggregating and outputting performance data 
        // and memory usage data.
        {

		// Get the MPI rank
		MPI_Comm_rank(MPI_COMM_WORLD, &rank);

		if (rank == 0) {
			// Print the start message
			printStartMessage();
		}

		// Set up the material infrastructure that is used to calculate flux
		auto material = initMaterial(opts);
		// Set up the temperature infrastructure
		bool tempInitOK = initTemp(opts);
		if (!tempInitOK) {
            throw std::runtime_error("Unable to initialize temperature.");
        }
		// Set up the visualization infrastructure.
		bool vizInitOK = initViz(opts.useVizStandardHandlers());
		if (!vizInitOK) {
            throw std::runtime_error("Unable to initialize visualization infrastructure.");
        }

		// Access the temperature handler registry to get the temperature
		auto tempHandler = xolotlFactory::getTemperatureHandler();

		// Access our performance handler registry to obtain a Timer
		// measuring the runtime of the entire program.
		auto totalTimer = handlerRegistry->getTimer("total");
		totalTimer->start();

        auto totalMemUsageSampler = memUsageHandlerRegistry->getMemUsageSampler("total");
        totalMemUsageSampler->start();

		// Initialize and get the solver handler
		bool dimOK = xolotlFactory::initializeDimension(opts);
		if (!dimOK) {
            throw std::runtime_error("Unable to initialize dimension from inputs.");
        }
		auto solvHandler = xolotlFactory::getSolverHandler();

		// Create the network handler factory
		auto networkFactory =
					xolotlFactory::IReactionHandlerFactory::createNetworkFactory(opts.getMaterial());

		// Setup and load the network
		auto networkLoadTimer = handlerRegistry->getTimer("loadNetwork");
		networkLoadTimer->start();
		networkFactory->initializeReactionNetwork(opts, handlerRegistry);
		networkLoadTimer->stop();

		// Get the network handler
		auto networkHandler = networkFactory->getNetworkHandler();

		// Setup the solver
		auto solver = setUpSolver(handlerRegistry, material, tempHandler, networkHandler,
				solvHandler, opts);

		// Launch the PetscSolver
        auto solverMemUsageSampler = memUsageHandlerRegistry->getMemUsageSampler("solver");
        solverMemUsageSampler->start();
		launchPetscSolver(solver, handlerRegistry);

		// Finalize our use of the solver.
		auto solverFinalizeTimer = handlerRegistry->getTimer("solverFinalize");
		solverFinalizeTimer->start();
		solver->finalize();
		solverFinalizeTimer->stop();
        solverMemUsageSampler->stop();

        totalMemUsageSampler->stop();
		totalTimer->stop();

        } // end of "real" work.

		// Report statistics about the performance data collected during
		// the run we just completed.
        auto perfStats = handlerRegistry->collectStatistics();
        if (rank == 0) {
            handlerRegistry->reportStatistics(std::cout, perfStats);
        }

        // Report on memory usage data collected during the run we
        // just completed.
        auto memUsageData = memUsageHandlerRegistry->collectData();
        if (rank == 0) {
            memUsageHandlerRegistry->reportData(std::cout, memUsageData);
        }

	} catch (const std::exception& e) {
		std::cerr << e.what() << std::endl;
		std::cerr << "Aborting." << std::endl;
		ret = EXIT_FAILURE;
	} catch (const std::string& error) {
		std::cout << error << std::endl;
		std::cout << "Aborting." << std::endl;
		ret = EXIT_FAILURE;
	}

	// finalize our use of MPI
	MPI_Finalize();

	return ret;
}
