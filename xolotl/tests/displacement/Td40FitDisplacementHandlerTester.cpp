#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <Td40FitDisplacementHandler.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <DummyHandlerRegistry.h>
#include <mpi.h>

using namespace std;
using namespace xolotlCore;

/**
 * The test suite is responsible for testing the WFitFluxHandler.
 */
BOOST_AUTO_TEST_SUITE (Td40FitDisplacementHandlerTester_testSuite)

BOOST_AUTO_TEST_CASE(checkGetInitialDisplacement) {
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the network loader
	HDF5NetworkLoader loader =
			HDF5NetworkLoader(make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Load the network
	auto network = (PSIClusterReactionNetwork *) loader.load().get();

	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 2.50025;

	// Create the flux handler
    auto testFitDisplacement = make_shared<Td40FitDisplacementHandler>();
    // Set the factor to change the helium flux
    testFitDisplacement->setKrFluenceAmplitude(1.0);
    // Initialize the flux handler
    testFitDisplacement->initializeDisplacementHandler(network, nGridpts, step);

	// Get the flux vector
	auto testDisplacementVec = testFitDisplacement->getInitialDisplacementVec();

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(testDisplacementVec[1], 0.23219583, 0.01);
	BOOST_REQUIRE_CLOSE(testDisplacementVec[2], 0.12701147, 0.01);
	BOOST_REQUIRE_CLOSE(testDisplacementVec[3], 0.0407527, 0.01);

	return;
}

BOOST_AUTO_TEST_CASE(checkDisplacementIndex) {
	// Create the network loader
	HDF5NetworkLoader loader =
			HDF5NetworkLoader(make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Load the network
	auto network = (PSIClusterReactionNetwork *) loader.load().get();

	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 2.50025;

	// Create the flux handler
    auto testFitDisplacement = make_shared<Td40FitDisplacementHandler>();
    // Set the factor to change the helium flux
    testFitDisplacement->setKrFluenceAmplitude(1.0);
    // Initialize the flux handler
    testFitDisplacement->initializeDisplacementHandler(network, nGridpts, step);

    // Check the value of the index of the cluster for the flux
    BOOST_REQUIRE_EQUAL(testFitDisplacement->getInitialDisplacementClusterIndex(), 8);

	return;
}

BOOST_AUTO_TEST_CASE(checkKrFluenceAmplitude) {
	// Create the network loader
	HDF5NetworkLoader loader =
			HDF5NetworkLoader(make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten_diminutive.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Load the network
	auto network = (PSIClusterReactionNetwork *) loader.load().get();

	// Specify the number of grid points that will be used
	int nGridpts = 5;
	// Specify the step size between grid points
	double step = 2.50025;

	// Create the flux handler
    auto testFitDisplacement = make_shared<Td40FitDisplacementHandler>();
    // Set the factor to change the helium flux
    testFitDisplacement->setKrFluenceAmplitude(1.0);
    // Set the factor to change the helium flux
    testFitDisplacement->setKrFluenceAmplitude(2.5);
    // Initialize the flux handler
    testFitDisplacement->initializeDisplacementHandler(network, nGridpts, step);

    // Check the value of the helium flux
    BOOST_REQUIRE_EQUAL(testFitDisplacement->getKrFluenceAmplitude(), 2.5);

	// Get the flux vector
	auto testDisplacementVec = testFitDisplacement->getInitialDisplacementVec();

	// Check the value at some grid points
	BOOST_REQUIRE_CLOSE(testDisplacementVec[1], 0.58048957, 0.01);
	BOOST_REQUIRE_CLOSE(testDisplacementVec[2], 0.31752868, 0.01);
	BOOST_REQUIRE_CLOSE(testDisplacementVec[3], 0.10188176, 0.01);

	return;
}

BOOST_AUTO_TEST_SUITE_END()
