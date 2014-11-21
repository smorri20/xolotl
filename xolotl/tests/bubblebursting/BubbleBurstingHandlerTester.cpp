#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Regression

#include <boost/test/included/unit_test.hpp>
#include <BubbleBurstingHandler.h>
#include <HDF5NetworkLoader.h>
#include <XolotlConfig.h>
#include <DummyHandlerRegistry.h>
#include <mpi.h>

using namespace std;
using namespace xolotlCore;

/**
 * This suite is responsible for testing the BubbleBurstingHandler.
 */
BOOST_AUTO_TEST_SUITE(BubbleBurstingHandler_testSuite)

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
	auto network = loader.load();
	// Get its size
	const int size = network->getAll()->size();
	// Set the temperature to 1000.0 K
	network->setTemperature(1000.0);

	// Suppose we have a grid with 3 grip points and distance of 0.5 nm between grid points
	double hx = 0.5;
	int nGrid = 3;
	// And the surface is on the left side
	int surfacePos = 0;

	// Create the bubble bursting handler
	BubbleBurstingHandler bubbleBurstingHandler;

	// Initialize it
	bubbleBurstingHandler.initialize(network, hx, nGrid, surfacePos);

	// The arrays of concentration
	double concentration[3*size];
	double newConcentration[3*size];

	// Initialize their values
	for (int i = 0; i < 3*size; i++) {
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
	bubbleBurstingHandler.computeBursting(network, 1, surfacePos,
			concOffset, updatedConcOffset);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[2900], 0.0, 0.01); // Does not burst (He47V36)
	BOOST_REQUIRE_CLOSE(updatedConcOffset[3146], -5.1929e21, 0.01); // (He148V37)
	BOOST_REQUIRE_CLOSE(updatedConcOffset[3454], -5.5718e21, 0.01); // (He154V39)
	BOOST_REQUIRE_CLOSE(updatedConcOffset[4302], 1.2256e24, 0.01); // (V45)
	BOOST_REQUIRE_CLOSE(updatedConcOffset[4481], -6.9316e21, 0.01); // (He179V45)

	// Initialize the rows, columns, and values to set in the Jacobian
	int nBubble = network->getAll(heVType).size();
	int row[2*nBubble], col[nBubble];
	double val[2*nBubble];
	// Get the pointer on them for the compute bursting method
	int *rowPointer = &row[0];
	int *colPointer = &col[0];
	double *valPointer = &val[0];

	// Compute the partial derivatives for the advection a the grid point 1
	int nBursting = bubbleBurstingHandler.computePartialsForBursting(network, valPointer,
			rowPointer, colPointer, 1, 0, surfacePos);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(row[0], 13895);
	BOOST_REQUIRE_EQUAL(row[1], 13894);
	BOOST_REQUIRE_EQUAL(row[2], 13896);
	BOOST_REQUIRE_EQUAL(row[3], 13894);

	BOOST_REQUIRE_EQUAL(col[0], 13895);
	BOOST_REQUIRE_EQUAL(col[1], 13896);

	// Check values
	BOOST_REQUIRE_CLOSE(val[0], -7.0311e13, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 7.0311e13, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], -7.0311e13, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], 7.0311e13, 0.01);

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
