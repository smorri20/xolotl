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

	// Suppose we have a grid with 3 grip points and distance of 0.5 nm between grid points
	double hx = 0.5;
	int nGrid = 3;

	// Create the bubble bursting handler
	BubbleBurstingHandler bubbleBurstingHandler;

	// Initialize it
	bubbleBurstingHandler.initialize(network, hx, nGrid);

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

	// Compute the advection at this grid point
	bubbleBurstingHandler.computeBursting(network, 1,
			concOffset, updatedConcOffset);

	// Check the new values of updatedConcOffset
	BOOST_REQUIRE_CLOSE(updatedConcOffset[2900], 0.0, 0.01); // Does not burst (He47V36)
	BOOST_REQUIRE_CLOSE(updatedConcOffset[3146], -7.3857e21, 0.01); // (He148V37)
	BOOST_REQUIRE_CLOSE(updatedConcOffset[3454], -7.9246e21, 0.01); // (He154V39)
	BOOST_REQUIRE_CLOSE(updatedConcOffset[4302], 1.7431e24, 0.01); // (V45)
	BOOST_REQUIRE_CLOSE(updatedConcOffset[4481], -9.8585e21, 0.01); // (He179V45)

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
			rowPointer, colPointer, 1, 0);

	// Check the values for the indices
	BOOST_REQUIRE_EQUAL(row[0], 13895);
	BOOST_REQUIRE_EQUAL(row[1], 13894);
	BOOST_REQUIRE_EQUAL(row[2], 13896);
	BOOST_REQUIRE_EQUAL(row[3], 13894);

	BOOST_REQUIRE_EQUAL(col[0], 13895);
	BOOST_REQUIRE_EQUAL(col[1], 13896);

	// Check values
	BOOST_REQUIRE_CLOSE(val[0], -1.0e14, 0.01);
	BOOST_REQUIRE_CLOSE(val[1], 1.0e14, 0.01);
	BOOST_REQUIRE_CLOSE(val[2], -1.0e14, 0.01);
	BOOST_REQUIRE_CLOSE(val[3], 1.0e14, 0.01);

	// Finalize MPI
	MPI_Finalize();
}

BOOST_AUTO_TEST_SUITE_END()
