#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <PSIClusterNetworkLoader.h>
#include <memory>
#include <typeinfo>
#include <limits>
#include <string.h>
#include <PSIClusterNetworkLoader.h>
#include <PSIClusterReactionNetwork.h>
#include <PetscSolver.h>
#include <XolotlConfig.h>
#include <xolotlPerf.h>
#include <DummyHandlerRegistry.h>
#include <HDF5NetworkLoader.h>
#include <Options.h>
#include <PetscSolver1DHandler.h>
#include <IMaterialFactory.h>
#include <TemperatureHandlerFactory.h>
#include <VizHandlerRegistryFactory.h>

using namespace std;
using namespace xolotlCore;

/**
 * The test suite configuration
 */BOOST_AUTO_TEST_SUITE (PetscSolverTester_testSuite)

/**
 * This operation checks the concentration of clusters after solving a test case
 * in 1D.
 */
BOOST_AUTO_TEST_CASE(checkPetscSolver1DHandler) {

	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Local Declarations
	string sourceDir(XolotlSourceDirectory);

	// Create a fake command line to read the options
	argc = 1;
	argv = new char*[2];
	std::string parameterFile = sourceDir + "/tests/testfiles/param_good.txt";
	argv[0] = new char[parameterFile.length() + 1];
	strcpy(argv[0], parameterFile.c_str());
	argv[1] = 0; // null-terminate the array

	// Read the options
	Options opts;
	opts.readParams(argc, argv);

	// Create the network loader
	std::shared_ptr<HDF5NetworkLoader> loader = std::make_shared<HDF5NetworkLoader>(
			make_shared<xolotlPerf::DummyHandlerRegistry>());

	// Create the path to the network file
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string networkFilename = sourceDir + pathToFile;

	BOOST_TEST_MESSAGE(
			"PetscSolverTester Message: Network filename is: " << networkFilename);

	// Give the filename to the network loader
	loader->setFilename(networkFilename);

	// Create the solver
	std::shared_ptr<xolotlSolver::PetscSolver> solver = std::make_shared<
			xolotlSolver::PetscSolver>(make_shared<xolotlPerf::DummyHandlerRegistry>());

	// Create the material factory
	auto materialFactory = xolotlFactory::IMaterialFactory::createMaterialFactory(opts.getMaterial());

	// Initialize and get the temperature handler
	bool tempInitOK = xolotlFactory::initializeTempHandler(opts);
	auto tempHandler = xolotlFactory::getTemperatureHandler();

	// Set up our dummy performance and visualization infrastructures
    xolotlPerf::initialize(xolotlPerf::toPerfRegistryType("dummy"));
    xolotlFactory::initializeVizHandler(false);

    // Create a solver handler and initialize it
	auto solvHandler = std::make_shared<xolotlSolver::PetscSolver1DHandler> ();
	solvHandler->initializeHandlers(materialFactory, tempHandler, opts);

    // Set the solver command line to give the PETSc options and initialize it
    solver->setCommandLineOptions(opts.getPetscArgc(), opts.getPetscArgv());
	solver->initialize(solvHandler);

	// Give it the network loader
	solver->setNetworkLoader(loader);

	// Solve and finalize
	solver->solve();
	solver->finalize();

	// Check the concentrations left in the network
	auto network = solvHandler->getNetwork();
	double concs[network->getAll()->size()];
	network->fillConcentrationsArray(concs);

    BOOST_REQUIRE_CLOSE(concs[0], 0.0, 0.01);
    BOOST_REQUIRE_CLOSE(concs[1], 1.281e-11, 0.01);
    BOOST_REQUIRE_CLOSE(concs[2], 5.338e-22, 0.01);
    BOOST_REQUIRE_CLOSE(concs[3], 0.025, 0.01);
    BOOST_REQUIRE_CLOSE(concs[4], 2.2744e-24, 0.01);
}

BOOST_AUTO_TEST_SUITE_END()
