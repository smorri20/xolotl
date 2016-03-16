#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <DummyBubbleBurstingHandler.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <DummyHandlerRegistry.h>
#include <mpi.h>

using namespace xolotlCore;

/**
 * This suite is responsible for testing the DummyBubbleBurstingHandler.
 */
BOOST_AUTO_TEST_SUITE(DummyBubbleBurstingHandler_testSuite)

/**
 * Method checking the initialization and the compute bubble bursting methods.
 */
BOOST_AUTO_TEST_CASE(checkBubbleBursting) {
	// Initialize MPI for HDF5
	int argc = 0;
	char **argv;
	MPI_Init(&argc, &argv);

	// Create the network loader
	HDF5NetworkLoader loader =
			HDF5NetworkLoader(make_shared<xolotlPerf::DummyHandlerRegistry>());
	// Define the filename to load the network from
	string sourceDir(XolotlSourceDirectory);
	string pathToFile("/tests/testfiles/tungsten.h5");
	string filename = sourceDir + pathToFile;
	// Give the filename to the network loader
	loader.setFilename(filename);

	// Load the network
	auto network = (PSIClusterReactionNetwork *) loader.load().get();
	// Get its size
	const int size = network->getAll()->size();
	// Set the temperature to 1000.0 K
	network->setTemperature(1000.0);

	// Suppose we have a grid with 13 grip points and distance of
	// 0.1 nm between grid points
	std::vector<double> grid;
	for (int l = 0; l < 13; l++) {
		grid.push_back((double) l * 0.1);
	}
	// Set the surface position
	int surfacePos = 0;

	// Create the bubble bursting handler
	DummyBubbleBurstingHandler bubbleBurstingHandler;

	// Initialize it
	bubbleBurstingHandler.initialize(surfacePos, network, grid);

	// The arrays of concentration
	double concentration[13*size];
	double newConcentration[13*size];

	// Initialize their values
	for (int i = 0; i < 13*size; i++) {
		concentration[i] = (double) i * i;
		newConcentration[i] = 0.0;
	}

	// Get pointers
	double *conc = &concentration[0];
	double *updatedConc = &newConcentration[0];

	// Get the offset for the grid point in the middle
	double *concOffset = conc + size;
	double *updatedConcOffset = updatedConc + size;

	// Compute the bubble bursting at this grid point
	bubbleBurstingHandler.computeBursting(network, 1,
			updatedConcOffset);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[34], 0.0, 0.01); // Does not burst (He10V2)
	BOOST_REQUIRE_CLOSE(updatedConcOffset[94], 0.0, 0.01); // Does not burst (He15V5)
	BOOST_REQUIRE_CLOSE(updatedConcOffset[97], 0.0, 0.01); // Does not burst (He18V5)
	BOOST_REQUIRE_CLOSE(updatedConcOffset[361], 0.0, 0.01); // Does not burst (V10)
	BOOST_REQUIRE_CLOSE(updatedConcOffset[306], 0.0, 0.01); // Does not burst (He45V10)

	// Initialize the indices and values to set in the Jacobian
	int nBubble = network->getAll(heVType).size();
	int indices[nBubble];
	double val[2*nBubble];
	// Get the pointer on them for the compute bursting method
	int *indicesPointer = &indices[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the bursting a the grid point 1
	int nBursting = bubbleBurstingHandler.computePartialsForBursting(network, valPointer,
			indicesPointer, 1);

	// Check that the number of bursting cluster is 0
	BOOST_REQUIRE_EQUAL(nBursting, 0);

	// Finalize MPI
	MPI_Finalize();

	return;
}

BOOST_AUTO_TEST_SUITE_END()
